"""Data loaders shared by the Python notebooks (modeling, scRNA dot plots).

Usage (utils/ must be on sys.path, as for config.py):

    from config import in_path          # noqa: F401  (imported here too)
    from data_loading import load_modeling_data, build_clinical_matrix
    X, y, metadata = load_modeling_data("all")          # all 91 samples (CV and Cohort 3 holdout)
    X, y, metadata = load_modeling_data("12mo")         # >=12-month sensitivity cohort
    X, y, metadata = load_modeling_data("dlbcl")        # de novo DLBCL only
    X_clin = build_clinical_matrix(metadata, "IPI")

Functions
    load_sample_metadata(order="paper")     91 cfRNA samples indexed by sample_id
    load_multisource_metadata()             128 cfRNA / TME / PBMC samples
    load_counts(type="norm", samples_as_rows=False)
    load_gene_key()
    ensembl_to_symbol(mat, gene_key, strip_version=False, duplicates="first")
    filter_pre_rnr(meta), filter_12mo(meta), filter_de_novo_dlbcl(meta)
    add_stratify_group(meta)
    load_modeling_data(cohort)              X (samples x genes), y, metadata for the ML analyses
    build_clinical_matrix(meta, var)        clinical feature matrix for one covariate
    clinical_feature_matrices(meta)         dict of all of them
"""
import numpy as np
import pandas as pd

try:
    from config import in_path
except ImportError:  # imported as a package module
    from .config import in_path

N_CFRNA = 91


# ------------------------------------------------------------------------------------------------
# Metadata
# ------------------------------------------------------------------------------------------------
def load_sample_metadata(order="paper"):
    """One row per cfRNA sample (n = 91), indexed by ``sample_id``.

    order="paper"     default sample order (all-sample CV and the Cohort 3 holdout).
    order="clinical"  order given by ``row_order_clinical`` (12-month and DLBCL analyses).

    Row order matters: train_test_split is positional, so the 101 Monte-Carlo splits depend on it.

    Clinical text columns (ipi_prelymphodepletion, ldh_day_minus30, stage_prior_cart, ...) are kept
    verbatim; ``km_include`` and ``in_12mo_cohort`` are bool. Column definitions:
    <input_dir>/metadata/column_dictionary.tsv.
    """
    if order not in ("paper", "clinical"):
        raise ValueError("order must be 'paper' or 'clinical'")
    m = pd.read_csv(in_path("metadata", "sample_metadata.tsv"), sep="\t",
                    dtype={"ipi_prelymphodepletion": object, "ldh_day_minus30": object})
    assert len(m) == N_CFRNA and m.sample_id.is_unique
    assert m.km_include.dtype == bool and m.in_12mo_cohort.dtype == bool
    if order == "clinical":
        m = m.sort_values("row_order_clinical", kind="mergesort")
    return m.set_index("sample_id")


def load_multisource_metadata():
    """cfRNA (91) + TME (Li et al. 2025, 20) + PBMC (Maurer et al. 2024, 17) samples; 128 rows."""
    m = pd.read_csv(in_path("metadata", "multisource_metadata.tsv"), sep="\t")
    assert len(m) == 128
    return m


# ------------------------------------------------------------------------------------------------
# Counts
# ------------------------------------------------------------------------------------------------
_COUNT_FILES = {
    "raw": "counts_raw_cart_paper.tsv.gz",               # raw counts, versioned Ensembl IDs x 91
    "norm": "counts_norm_cart_paper.tsv.gz",             # normalized, versioned Ensembl IDs x 91
    "tmmcpm_all": "counts_TMMCPM_cart_paper_ALLDATA.tsv.gz",  # TMM-CPM, gene symbols x 130 (cfRNA+PBMC+TME)
}


def load_counts(type="norm", samples_as_rows=False):
    """Gene x sample DataFrame (index = gene IDs), read as ``pd.read_csv(index_col=0, sep='\\t')``.

    samples_as_rows=True returns the transpose (samples x genes), as used by the modeling code.
    """
    if type not in _COUNT_FILES:
        raise ValueError(f"type must be one of {list(_COUNT_FILES)}")
    x = pd.read_csv(in_path("counts", _COUNT_FILES[type]), sep="\t", index_col=0)
    return x.T if samples_as_rows else x


def load_gene_key():
    """GENCODE key with columns gene_id (versioned), gene_name, gene_type."""
    return pd.read_csv(in_path("reference", "gencode.biotype.name.key.tsv.gz"), sep="\t")


def ensembl_to_symbol(mat, gene_key, strip_version=False, duplicates="first"):
    """Python port of ensembl_to_symbol() in data_loading.R (gene x sample input).

    duplicates="first": keep the input row order, drop unmapped rows and later rows with an
        already-seen symbol.
    duplicates="mean": inner join sorted by ID, average rows sharing a symbol, symbols in order of
        first appearance.
    """
    strip = (lambda s: s.str.replace(r"\..*", "", regex=True))
    ids = pd.Series(mat.index, index=mat.index)
    key = gene_key[["gene_id", "gene_name"]].copy()
    if strip_version:
        ids = strip(ids)
        key["gene_id"] = strip(key.gene_id)
    if duplicates == "first":
        sym_map = key.drop_duplicates("gene_id").set_index("gene_id").gene_name
        new = ids.map(sym_map)
        valid = new.notna() & ~new.duplicated()
        out = mat.loc[valid.values].copy()
        out.index = new[valid].values
        return out
    if duplicates == "mean":
        key = key.drop_duplicates("gene_id")
        df = mat.copy()
        df.insert(0, "_id", ids.values)
        df = df.merge(key, left_on="_id", right_on="gene_id", how="inner").sort_values("_id", kind="mergesort")
        out = df.groupby("gene_name", sort=False)[list(mat.columns)].mean()
        out.index.name = None
        return out
    raise ValueError("duplicates must be 'first' or 'mean'")


# ------------------------------------------------------------------------------------------------
# Cohort filters
# ------------------------------------------------------------------------------------------------
def filter_pre_rnr(meta):
    """timepoint in {leuk, PRE} and response in {R, NR} (all 91 samples pass)."""
    return meta[meta["timepoint"].isin(["leuk", "PRE"]) & meta["response"].isin(["R", "NR"])]


def filter_12mo(meta):
    """>=12-month cohort (n = 85): all NR plus R alive or surviving > 365 days.

    Precomputed as ``in_12mo_cohort``.
    """
    return meta[meta["in_12mo_cohort"]]


def filter_de_novo_dlbcl(meta):
    """De novo DLBCL only (lymphoma_malignancy == "De novo DLBCL"), n = 67."""
    return meta[meta["lymphoma_malignancy"].isin(["De novo DLBCL"])]


def add_stratify_group(meta):
    """stratify_group = response + "_" + cohort (e.g. "R_Cohort 1").

    StratifiedShuffleSplit only uses the sorted order of these groups, so any coding of the
    cohorts with the same sort order gives the same splits.
    """
    meta = meta.copy()
    meta["stratify_group"] = meta["response"].astype(str) + "_" + meta["cohort"].astype(str)
    return meta


def load_modeling_data(cohort="all"):
    """X (samples x genes, counts_norm), y (response) and metadata for the ML analyses.

    cohort="all"    all samples (paper order, n = 91)
    cohort="12mo"   >=12-month sensitivity cohort (clinical order, n = 85)
    cohort="dlbcl"  de novo DLBCL only (clinical order, n = 67)
    """
    X = load_counts("norm", samples_as_rows=True)
    if cohort == "all":
        metadata = filter_pre_rnr(load_sample_metadata(order="paper"))
    elif cohort == "12mo":
        m = load_sample_metadata(order="clinical")
        metadata = filter_12mo(m[m["timepoint"].isin(["leuk", "PRE"])])
    elif cohort == "dlbcl":
        metadata = filter_de_novo_dlbcl(filter_pre_rnr(load_sample_metadata(order="clinical")))
    else:
        raise ValueError("cohort must be 'all', '12mo' or 'dlbcl'")
    metadata = add_stratify_group(metadata)
    X = X.loc[metadata.index]
    y = metadata["response"]
    return X, y, metadata


# ------------------------------------------------------------------------------------------------
# Clinical feature matrices
# ------------------------------------------------------------------------------------------------
CLINICAL_VARS = ["bridging", "IPI", "LDH", "subtype", "histology", "stage"]


def build_clinical_matrix(meta, var):
    """One clinical feature matrix (rows = meta.index).

    var            metadata column          encoding -> output column(s)
      bridging     bridging_therapy         Yes 1.0 / No 0.0 -> "bridging_therapy"
      IPI          ipi_prelymphodepletion   "Unknown" -> NaN -> "IPI"
      LDH          ldh_prelymphodepletion   numeric -> "LDH"
      subtype      dlbcl_subtype            "Unknown" -> NaN, one-hot "subtype_*"
      histology    lymphoma_malignancy      one-hot "histology_*"
      stage        stage_prior_cart         I..IV -> 1..4, "Unstaged or Unknown" -> NaN -> "Stage"
    var="clinical_all" concatenates them in the order bridging, IPI, LDH, subtype, histology, stage.
    One-hot columns are built from the levels present in ``meta``.
    """
    if var == "bridging":
        return meta["bridging_therapy"].map({"Yes": 1.0, "No": 0.0}).to_frame(name="bridging_therapy")
    if var == "IPI":
        return pd.to_numeric(meta["ipi_prelymphodepletion"].replace("Unknown", np.nan),
                             errors="coerce").to_frame(name="IPI")
    if var == "LDH":
        return pd.to_numeric(meta["ldh_prelymphodepletion"], errors="coerce").to_frame(name="LDH")
    if var == "subtype":
        raw = meta["dlbcl_subtype"].replace("Unknown", np.nan)
        X = pd.get_dummies(raw, prefix="subtype").astype(float)
        X.loc[raw.isna()] = np.nan
        return X
    if var == "histology":
        raw = meta["lymphoma_malignancy"]
        X = pd.get_dummies(raw, prefix="histology").astype(float)
        X.loc[raw.isna()] = np.nan
        return X
    if var == "stage":
        return (meta["stage_prior_cart"]
                .map({"I": 1.0, "II": 2.0, "III": 3.0, "IV": 4.0, "Unstaged or Unknown": np.nan})
                .to_frame(name="Stage"))
    if var == "clinical_all":
        return pd.concat([build_clinical_matrix(meta, v) for v in CLINICAL_VARS], axis=1)
    raise ValueError(f"var must be one of {CLINICAL_VARS + ['clinical_all']}")


def clinical_feature_matrices(meta):
    """{"bridging": X_bridging, "IPI": X_IPI, ..., "stage": X_stage}."""
    return {v: build_clinical_matrix(meta, v) for v in CLINICAL_VARS}

