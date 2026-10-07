"""LymphoMAP archetype gene sets (LN, FMAC, TEX) for the Python notebooks.

The gene sets are defined once, in ``genesets.yaml`` next to this module, which
genesets.R also reads.

    from utils.genesets import GENESETS, LN_ENSEMBL, TEX_SYMBOLS

``GENESETS[name]`` is ``{"symbols": [...], "ensembl": [...]}`` for
name in ("LN", "FMAC", "TEX").

Notes
-----
* LN has 21 Ensembl IDs for 20 symbols, because SIGLEC5 has two GENCODE v38 IDs.
"""
from pathlib import Path

import yaml

GENESETS_YAML = Path(__file__).parent / "genesets.yaml"

with open(GENESETS_YAML) as _fh:
    GENESETS = yaml.safe_load(_fh)

LN_SYMBOLS = list(GENESETS["LN"]["symbols"])
LN_ENSEMBL = list(GENESETS["LN"]["ensembl"])
FMAC_SYMBOLS = list(GENESETS["FMAC"]["symbols"])
FMAC_ENSEMBL = list(GENESETS["FMAC"]["ensembl"])
TEX_SYMBOLS = list(GENESETS["TEX"]["symbols"])
TEX_ENSEMBL = list(GENESETS["TEX"]["ensembl"])

# Sanity checks on the yaml contents.
assert len(LN_ENSEMBL) == 21 and len(set(LN_ENSEMBL)) == 21, "LN should have 21 Ensembl IDs"
assert len(FMAC_ENSEMBL) == 20 and len(set(FMAC_ENSEMBL)) == 20, "FMAC should have 20 Ensembl IDs"
assert len(TEX_ENSEMBL) == 20 and len(set(TEX_ENSEMBL)) == 20, "TEX should have 20 Ensembl IDs"


def get_geneset(name, id_type="ensembl"):
    """Return a copy of gene set ``name`` ('LN', 'FMAC' or 'TEX') as 'ensembl' IDs or 'symbols'."""
    return list(GENESETS[name][id_type])
