# Gene-set scoring (adapted from Seurat::AddModuleScore). Used for both the
# cfRNA / bulk matrices and the single-cell pseudobulk scores.

#' Score one gene set per sample against expression-matched control genes.
#'
#' All genes in `pool` are ranked by their mean expression across samples and
#' cut into `nbin` equal-sized bins (`ggplot2::cut_number`). The 1e-30 random
#' jitter is inherited from Seurat. It is below double precision for genes with
#' a non-zero mean, but it does break ties among genes with a mean of exactly
#' zero, and it consumes random numbers, so it is kept for exact
#' reproducibility. For every gene in the set, `ctrl` control genes are
#' sampled without replacement from the same expression bin. The raw score of
#' a sample is mean(set genes) - mean(control genes); it is then min-max
#' rescaled to [0, 1] across the samples (columns) of `object`.
#'
#' `nbin` sets how finely control genes are matched on expression: more bins
#' give tighter matching but smaller bins to sample from (each bin must hold
#' at least `ctrl` genes). On the normalized cfRNA matrix (60,708 genes, many tied
#' low-expression means) nbin = 25 fails with "Insufficient data values to
#' produce 25 bins", which is why the notebooks scored that matrix with nbin = 6.
#' Because both the bins and the final rescaling
#' depend on which samples are in `object`, scores are only comparable
#' within one call.
#'
#' @param object  genes x samples numeric matrix or data.frame, gene symbols as rownames.
#' @param features character vector of genes; genes missing from `object` are
#'   dropped with a warning.
#' @param pool    genes used to build the expression bins (default: all rows).
#' @param nbin    number of expression bins (default 25).
#' @param ctrl    control genes sampled per set gene (default 100).
#' @param name    column name prefix; the output column is paste0(name, "1").
#' @param seed    passed to set.seed() before binning and sampling. This resets the
#'   global RNG state (it is not restored afterwards); NULL leaves the RNG alone.
#' @return data.frame with one column (e.g. "LN_set1") and samples as rownames.
AddGeneSetScore <- function(object, features, pool = NULL, nbin = 25, ctrl = 100,
                            name = "Set", seed = 123) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(features)) stop("Missing input feature list")

  missing_features <- setdiff(features, rownames(object))
  if (length(missing_features) > 0) {
    warning("The following features are not present in the object: ",
            paste(missing_features, collapse = ", "))
    warning(paste0("\n ", paste(missing_features, collapse = ", "),
                   " dropped for calculating the geneset score."))
  }
  features <- intersect(features, rownames(object))

  if (is.null(pool)) pool <- rownames(object)
  data_avg <- Matrix::rowMeans(object[pool, ])
  data_avg <- data_avg[order(data_avg)]
  data_cut <- ggplot2::cut_number(data_avg + rnorm(n = length(data_avg)) / 1e30,
                                  n = nbin, labels = FALSE, right = FALSE)
  names(data_cut) <- names(data_avg)

  ctrl_use <- NULL
  for (gene in features) {
    same_bin <- data_cut[which(data_cut == data_cut[gene])]
    ctrl_use <- c(ctrl_use, names(sample(x = same_bin, size = ctrl, replace = FALSE)))
  }
  ctrl_use <- unique(ctrl_use)

  ctrl_score    <- Matrix::colMeans(object[ctrl_use, ])
  feature_score <- Matrix::colMeans(object[features, , drop = FALSE])
  score <- feature_score - ctrl_score

  out <- data.frame((score - min(score)) / (max(score) - min(score)))
  colnames(out) <- paste0(name, "1")
  rownames(out) <- colnames(object)
  out
}
