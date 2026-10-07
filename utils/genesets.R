# LymphoMAP archetype gene sets (LN, FMAC, TEX), read from utils/genesets.yaml,
# the single source shared with genesets.py. Requires REPO_ROOT (config.R).

#' Named list of gene sets; each element has $symbols and $ensembl (GENCODE v38).
GENESETS <- lapply(yaml::read_yaml(file.path(REPO_ROOT, "utils", "genesets.yaml")),
                   function(gs) lapply(gs, as.character))

LN_SYMBOLS   <- GENESETS$LN$symbols
FMAC_SYMBOLS <- GENESETS$FMAC$symbols
TEX_SYMBOLS  <- GENESETS$TEX$symbols

LN_ENSEMBL   <- GENESETS$LN$ensembl
FMAC_ENSEMBL <- GENESETS$FMAC$ensembl
TEX_ENSEMBL  <- GENESETS$TEX$ensembl
