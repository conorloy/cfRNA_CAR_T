# Entry point for every R notebook:
#   source("../../utils/setup.R")
# Loads paths (config.R, which reads config.yaml and config.local.yaml), themes/palettes, gene sets and the shared helpers, and
# attaches the common packages. Packages used by one or two notebooks (DESeq2,
# Seurat, mice, pROC, survminer, pheatmap, ...) are loaded where they are used.

# Directory of this file: taken from source() when possible, otherwise found by
# walking up from the working directory to the folder holding config.yaml.
.find_utils_dir <- function() {
  for (i in rev(seq_len(sys.nframe()))) {
    f <- sys.frame(i)$ofile
    if (!is.null(f) && basename(f) == "setup.R") return(dirname(normalizePath(f)))
  }
  dir <- normalizePath(getwd())
  repeat {
    if (file.exists(file.path(dir, "config.yaml"))) return(file.path(dir, "utils"))
    parent <- dirname(dir)
    if (parent == dir) stop("Cannot locate utils/: config.yaml not found above ", getwd())
    dir <- parent
  }
}

# Startup messages and "package was built under R version x" warnings are muffled.
withCallingHandlers(suppressPackageStartupMessages({
  library(data.table)  # attached first so dplyr's first()/last()/between() take precedence
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(ggpubr)
}), warning = function(w) {
  if (grepl("was built under R version", conditionMessage(w))) invokeRestart("muffleWarning")
})

UTILS_DIR <- .find_utils_dir()
rm(.find_utils_dir)

# chdir = TRUE: config.R finds the repo root by walking up from the working directory.
source(file.path(UTILS_DIR, "config.R"), chdir = TRUE)
source(file.path(UTILS_DIR, "theme.R"))
source(file.path(UTILS_DIR, "genesets.R"))
source(file.path(UTILS_DIR, "geneset_scoring.R"))
source(file.path(UTILS_DIR, "plotting.R"))
source(file.path(UTILS_DIR, "survival.R"))
source(file.path(UTILS_DIR, "ml_eval.R"))
source(file.path(UTILS_DIR, "data_loading.R"))

options(repr.plot.width = 6, repr.plot.height = 6)

# Seeded steps (gene-set scoring, jitter, MICE) set their own seeds; this only
# fixes anything left unseeded, such as plain geom_jitter().
set.seed(42)
