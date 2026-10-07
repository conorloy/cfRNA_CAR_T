# Data loaders shared by all R notebooks.
#
# Requires utils/config.R to be sourced first (provides in_path()).
# Inputs live in <input_dir> (config.yaml); see README.md > Input data.
#
#   load_sample_metadata(order, factors)   91 cfRNA samples, one row each (metadata/sample_metadata.tsv)
#   load_multisource_metadata(factors)     128 cfRNA / TME / PBMC samples (metadata/multisource_metadata.tsv)
#   load_counts(type)                      gene x sample data.frame, gene IDs as rownames
#   load_gene_key()                        GENCODE gene_id / gene_name / gene_type
#   ensembl_to_symbol(mat, gene_key, ...)  Ensembl-ID rows -> gene-symbol rows
#   symbol_to_ensembl(symbols, gene_key, available_ids)
#   recode_clinical(meta, variant)         clinical covariates as coded for each analysis
#   filter_12mo(meta), filter_de_novo_dlbcl(meta)
#
# Base R only; data.table is needed only for ensembl_to_symbol(duplicates = "mean").

RESPONSE_LEVELS <- c("NR", "R")  # the alphabetical default, as used by DESeq2 and pROC
COHORT_LEVELS <- c("Cohort 1", "Cohort 2", "Cohort 3")
MULTISOURCE_COHORT_LEVELS <- c(COHORT_LEVELS, "Maurer et al. 2024", "Li et al. 2025")
DATASET_LEVELS <- c("cfRNA", "PBMC", "TME")

# ---------------------------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------------------------

# One row per cfRNA sample (n = 91). Column definitions: <input_dir>/metadata/column_dictionary.tsv.
#   order = "paper"    default sample order (used by most notebooks)
#   order = "clinical" order given by row_order_clinical (the 12-month and DLBCL analyses);
#                      row order sets the random splits and point jitter of some panels
#   factors = TRUE     response -> factor(c("NR","R")), cohort -> factor(Cohort 1..3)
load_sample_metadata <- function(order = c("paper", "clinical"), factors = TRUE) {
  order <- match.arg(order)
  m <- read.delim(in_path("metadata", "sample_metadata.tsv"), na.strings = "NA",
                  stringsAsFactors = FALSE,
                  colClasses = c(ipi_prelymphodepletion = "character",
                                 ldh_day_minus30 = "character"))
  stopifnot(nrow(m) == 91, !anyDuplicated(m$sample_id),
            is.logical(m$km_include), is.logical(m$in_12mo_cohort))
  if (order == "clinical") m <- m[order(m$row_order_clinical), ]
  rownames(m) <- NULL
  if (factors) {
    m$response <- factor(m$response, levels = RESPONSE_LEVELS)
    m$cohort <- factor(m$cohort, levels = COHORT_LEVELS)
  }
  m
}

# cfRNA (91) + TME (Li et al. 2025, 20) + PBMC (Maurer et al. 2024, 17) = 128 samples.
# Columns: sample_id, dataset, response, cohort, subtype (GCB / Non-GCB / Unknown / NA; NA for all PBMC).
load_multisource_metadata <- function(factors = TRUE) {
  m <- read.delim(in_path("metadata", "multisource_metadata.tsv"), na.strings = "NA",
                  stringsAsFactors = FALSE)
  stopifnot(nrow(m) == 128)
  if (factors) {
    m$response <- factor(m$response, levels = RESPONSE_LEVELS)
    m$cohort <- factor(m$cohort, levels = MULTISOURCE_COHORT_LEVELS)
    m$dataset <- factor(m$dataset, levels = DATASET_LEVELS)
  }
  m
}

# ---------------------------------------------------------------------------------------------
# Counts
# ---------------------------------------------------------------------------------------------

# Gene x sample data.frame (read.delim(..., row.names = 1)).
#   "raw"        counts_raw_cart_paper.tsv.gz       raw counts, versioned Ensembl IDs x 91 samples (DESeq2)
#   "norm"       counts_norm_cart_paper.tsv.gz       normalized, versioned Ensembl IDs x 91 samples
#   "tmmcpm_all" counts_TMMCPM_cart_paper_ALLDATA   TMM-CPM, gene SYMBOLS (R make.names style, e.g.
#                                                   "RP11.34P13.7") x 130 samples (91 cfRNA + 17 PBMC +
#                                                   22 TME; 2 TME samples have no metadata
#                                                   and are dropped by subsetting on sample_id)
# All sample IDs are syntactically valid, so check.names has no effect on column names.
load_counts <- function(type = c("raw", "norm", "tmmcpm_all"), as_matrix = FALSE) {
  type <- match.arg(type)
  f <- switch(type,
              raw = "counts_raw_cart_paper.tsv.gz",
              norm = "counts_norm_cart_paper.tsv.gz",
              tmmcpm_all = "counts_TMMCPM_cart_paper_ALLDATA.tsv.gz")
  x <- read.delim(in_path("counts", f), row.names = 1, check.names = FALSE)
  if (as_matrix) x <- as.matrix(x)
  x
}

load_gene_key <- function() {
  read.delim(in_path("reference", "gencode.biotype.name.key.tsv.gz"), stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------------------------
# Ensembl -> symbol
# ---------------------------------------------------------------------------------------------

# Convert a gene x sample matrix/data.frame with Ensembl-ID rownames to gene-symbol rownames.
#
#   duplicates = "first" (default): look up each row's ID in gene_key, drop rows with no symbol, and
#       drop every later row whose symbol was already seen. The input row order is kept:
#           sym_map <- setNames(gene_key$gene_name, gene_key$gene_id)
#           new_rownames <- sym_map[rownames(counts)]
#           valid <- !is.na(new_rownames) & !duplicated(new_rownames)
#   duplicates = "mean": inner-join on the ID. Rows come out sorted by ID, as data.table::merge sorts them.
#       Rows that share a symbol are averaged, and symbols are kept in order of first appearance.
#       With strip_version = TRUE (used for the decision-curve analysis): version suffixes are removed from both the matrix IDs and gene_key$gene_id,
#       gene_key is de-duplicated by stripped ID (first entry kept, e.g. PAR genes), then
#       merge + counts_dt[, lapply(.SD, mean), by = gene_name].
#   strip_version: remove the ".N" Ensembl version suffix before matching.
#
# counts_TMMCPM_cart_paper_ALLDATA already has symbol rownames; no conversion is needed.
ensembl_to_symbol <- function(mat, gene_key, strip_version = FALSE,
                              duplicates = c("first", "mean")) {
  duplicates <- match.arg(duplicates)
  strip <- function(x) sub("\\..*", "", x)
  if (duplicates == "first") {
    key_id <- if (strip_version) strip(gene_key$gene_id) else gene_key$gene_id
    row_id <- if (strip_version) strip(rownames(mat)) else rownames(mat)
    sym_map <- setNames(gene_key$gene_name, key_id)
    new_rownames <- sym_map[row_id]
    valid <- !is.na(new_rownames) & !duplicated(new_rownames)
    out <- mat[valid, , drop = FALSE]
    rownames(out) <- unname(new_rownames[valid])
    return(out)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) stop("duplicates = 'mean' needs data.table")
  `:=` <- data.table::`:=`
  sample_cols <- colnames(mat)
  dt <- data.table::as.data.table(as.data.frame(mat, check.names = FALSE),
                                  keep.rownames = "ensembl_id_versioned")
  gk <- data.table::as.data.table(gene_key)
  if (strip_version) {
    dt[, ensembl_id := sub("\\..*", "", ensembl_id_versioned)]
    gk[, gene_id := sub("\\..*", "", gene_id)]
  } else {
    dt[, ensembl_id := ensembl_id_versioned]
  }
  gk <- unique(gk, by = "gene_id")
  dt <- merge(dt, gk[, c("gene_id", "gene_name"), with = FALSE], by.x = "ensembl_id", by.y = "gene_id")
  collapsed <- dt[, lapply(.SD, mean), by = gene_name, .SDcols = sample_cols]
  out <- as.data.frame(collapsed[, sample_cols, with = FALSE], check.names = FALSE)
  rownames(out) <- collapsed$gene_name
  out
}

# Reverse lookup for a set of symbols: keep gene_key rows whose symbol is in
# `symbols` and whose ID is in `available_ids` (e.g. rownames of the norm counts), one row per
# gene_id, then one row per symbol (first kept, with a warning).
symbol_to_ensembl <- function(symbols, gene_key, available_ids) {
  id_map <- gene_key[gene_key$gene_name %in% symbols & gene_key$gene_id %in% available_ids,
                     c("gene_id", "gene_name")]
  id_map <- id_map[!duplicated(id_map), ]
  dup <- unique(id_map$gene_name[duplicated(id_map$gene_name)])
  if (length(dup)) {
    warning("Genes mapping to multiple gene_ids (keeping first): ", paste(dup, collapse = ", "))
    id_map <- id_map[!duplicated(id_map$gene_name), ]
  }
  miss <- setdiff(symbols, id_map$gene_name)
  if (length(miss)) warning("Dropping genes not found via gene_key: ", paste(miss, collapse = ", "))
  rownames(id_map) <- NULL
  id_map
}

# ---------------------------------------------------------------------------------------------
# Clinical covariates
# ---------------------------------------------------------------------------------------------

.num <- function(x) suppressWarnings(as.numeric(x))  # "Unknown" / "Not available" -> NA
.NON_DE_NOVO <- c("tFL", "HGBL", "tCLL", "tMZL", "FL", "PMBCL")

# Adds recoded clinical covariates to `meta` (a load_sample_metadata() data.frame). Rows and
# existing columns are unchanged.
#
# variant = "baseline"  MICE logistic/Cox models (SFig4, SFig9, Tables S3/S4):
#     IPI      numeric IPI ("Unknown" -> NA)
#     LDH      log2(LDH pre-lymphodepletion)
#     LDH_D30  log2(Day -30 LDH)  ("Not available" -> NA)
#     Bridging 1L Yes / 0L No / NA
#     Stage    1L III-IV / 0L I-II / NA ("Unstaged or Unknown", missing)
#     Lymphoma_malignancy 1L De novo DLBCL / 0L tFL,HGBL,tCLL,tMZL,FL,PMBCL / NA
# variant = "cox_mice"  univariate Cox forest (SFig2, MICE):
#     lines_of_therapy, IPI_prelymphodepletion, LDH_prelymphodepletion, LDH_Day30 (numeric, not logged),
#     stage_prior_cart_ord (ordered factor I<II<III<IV; a new name so that the verbatim
#     stage_prior_cart column is not overwritten), Malignancy (factor), Timepoint (factor),
#     Bridging (factor Yes/No), De_novo_DLBCL (factor Yes/No), stage_advanced
#     (factor "Advanced (III-IV)" / "Early (I-II)")
# variant = "dca"       decision curve (SFig10C):
#     LDH_Day30, LDH_prelymph, IPI (numeric), bridging_therapy (factor Yes/No),
#     stage (character "Advanced (III-IV)" / "Early (I-II)"),
#     histology (character "De novo DLBCL" / "Non De novo DLBCL")
# The encoding used by the Python models is clinical_feature_matrices() in data_loading.py.
recode_clinical <- function(meta, variant = c("baseline", "cox_mice", "dca")) {
  variant <- match.arg(variant)
  lm <- meta$lymphoma_malignancy
  st <- meta$stage_prior_cart
  br <- meta$bridging_therapy
  if (variant == "baseline") {
    meta$IPI <- .num(meta$ipi_prelymphodepletion)
    meta$LDH <- log2(.num(meta$ldh_prelymphodepletion))
    meta$LDH_D30 <- log2(.num(meta$ldh_day_minus30))
    meta$Bridging <- ifelse(br %in% "Yes", 1L, ifelse(br %in% "No", 0L, NA_integer_))
    meta$Stage <- ifelse(st %in% c("III", "IV"), 1L, ifelse(st %in% c("I", "II"), 0L, NA_integer_))
    meta$Lymphoma_malignancy <- ifelse(lm %in% "De novo DLBCL", 1L,
                                       ifelse(lm %in% .NON_DE_NOVO, 0L, NA_integer_))
  } else if (variant == "cox_mice") {
    meta$lines_of_therapy <- as.numeric(meta$lines_of_therapy)
    meta$IPI_prelymphodepletion <- .num(meta$ipi_prelymphodepletion)
    meta$LDH_prelymphodepletion <- as.numeric(meta$ldh_prelymphodepletion)
    meta$stage_prior_cart_ord <- factor(st, levels = c("I", "II", "III", "IV"), ordered = TRUE)
    meta$LDH_Day30 <- .num(meta$ldh_day_minus30)
    meta$Malignancy <- factor(lm)
    meta$Timepoint <- factor(meta$timepoint)
    meta$Bridging <- factor(br)
    meta$De_novo_DLBCL <- factor(ifelse(lm %in% "De novo DLBCL", "Yes",
                                        ifelse(lm %in% .NON_DE_NOVO, "No", NA_character_)))
    so <- meta$stage_prior_cart_ord
    meta$stage_advanced <- factor(ifelse(so %in% c("III", "IV"), "Advanced (III-IV)",
                                         ifelse(so %in% c("I", "II"), "Early (I-II)", NA_character_)))
  } else {
    meta$LDH_Day30 <- .num(meta$ldh_day_minus30)
    meta$LDH_prelymph <- .num(meta$ldh_prelymphodepletion)
    meta$bridging_therapy <- factor(br)
    meta$IPI <- .num(meta$ipi_prelymphodepletion)
    meta$stage <- ifelse(st %in% c("III", "IV"), "Advanced (III-IV)",
                         ifelse(st %in% c("I", "II"), "Early (I-II)", NA_character_))
    meta$histology <- ifelse(lm %in% "De novo DLBCL", "De novo DLBCL",
                             ifelse(lm %in% .NON_DE_NOVO, "Non De novo DLBCL", NA_character_))
  }
  meta
}

# ---------------------------------------------------------------------------------------------
# Cohort filters
# ---------------------------------------------------------------------------------------------

# >=12-month sensitivity cohort (n = 85): all NR plus R that are alive or survived > 365 days.
# use_flag = TRUE uses the precomputed in_12mo_cohort column; FALSE recomputes it as
#   (response_pd_new == "NR" | days_to_death > 365 | is.na(days_to_death)) &
#   response_pd_new %in% c("R","NR")
# Rows keep the order of `meta`; the paper's 12-month analyses use load_sample_metadata(order = "clinical").
filter_12mo <- function(meta, use_flag = TRUE) {
  keep <- if (use_flag) {
    meta$in_12mo_cohort
  } else {
    (meta$response_pd_new %in% "NR" | meta$days_to_death > 365 | is.na(meta$days_to_death)) &
      meta$response_pd_new %in% c("R", "NR")
  }
  keep[is.na(keep)] <- FALSE
  meta[keep, , drop = FALSE]
}

# De novo DLBCL only (SFig3B-C), n = 67.
filter_de_novo_dlbcl <- function(meta) {
  meta[meta$lymphoma_malignancy %in% "De novo DLBCL", , drop = FALSE]
}
