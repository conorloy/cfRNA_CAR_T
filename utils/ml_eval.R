# Evaluation and plotting of the Monte Carlo cross-validated classifiers
# (outputs of code/01_compute: <method>/classifier_scores.csv and
# <method>/selected_features.csv). pROC is called as pROC::.

# ---- Display names ------------------------------------------------------------

#' Directory name of a model -> display label, for single gene sets and
#' clinical variables (SFig10A) and LN + clinical variable (SFig10B). Note the
#' different capitalisation of "Bridging therapy" and "LN + Bridging Therapy".
NAME_MAP <- c(
  FMAC            = "FMAC",
  TEX             = "TEX",
  LN              = "LN",
  LASSO           = "All Transcripts",
  RAND_30genes    = "Random Genes",
  LASSO_SHUFFLE   = "Shuffled Labels",
  bridging        = "Bridging therapy",
  IPI             = "IPI",
  LDH             = "LDH",
  subtype         = "Subtype (DLBCL)",
  histology       = "Histology",
  stage           = "Stage",
  LN_bridging     = "LN + Bridging Therapy",
  LN_IPI          = "LN + IPI",
  LN_histology    = "LN + Histology",
  LN_LDH          = "LN + LDH",
  LN_stage        = "LN + Stage",
  LN_subtype      = "LN + Subtype (DLBCL)",
  LN_clinical_all = "LN + Clinical (All)"
)

#' Display label -> model group.
GROUP_MAP <- c(
  "All Transcripts"       = "Gene set",
  "LN"                    = "Gene set",
  "FMAC"                  = "Gene set",
  "TEX"                   = "Gene set",
  "Bridging therapy"      = "Clinical",
  "IPI"                   = "Clinical",
  "LDH"                   = "Clinical",
  "Subtype (DLBCL)"       = "Clinical",
  "Histology"             = "Clinical",
  "Stage"                 = "Clinical",
  "Random Genes"          = "Control",
  "Shuffled Labels"       = "Control",
  "LN + Bridging Therapy" = "Gene Set + Clinical",
  "LN + IPI"              = "Gene Set + Clinical",
  "LN + LDH"              = "Gene Set + Clinical",
  "LN + Subtype (DLBCL)"  = "Gene Set + Clinical",
  "LN + Histology"        = "Gene Set + Clinical",
  "LN + Stage"            = "Gene Set + Clinical",
  "LN + Clinical (All)"   = "Gene Set + Clinical"
)

#' Forest-plot colour per display label: GENESET_COLOR for gene-set and control
#' models, CLINICAL_COLOR for models using clinical variables (SFig8, SFig10A/B).
FOREST_COLORS <- setNames(
  ifelse(GROUP_MAP %in% c("Clinical", "Gene Set + Clinical"), CLINICAL_COLOR, GENESET_COLOR),
  names(GROUP_MAP)
)

# ---- AUCs ---------------------------------------------------------------------

#' Train and test ROC AUC (percent, 0-100) for every Monte Carlo iteration of one model.
#'
#' With sens_spec = TRUE it also returns the per-split sensitivity/specificity
#' stored in the CSV (first row of each split after merging on sample_id; NA if
#' the CSV lacks those columns).
#' @param dir_path    folder containing <method_name>/classifier_scores.csv.
#' @param method_name model folder name, e.g. "LN", "LASSO".
#' @param mdf         metadata with sample_id and response (R/NR).
#' @return data.frame: iteration, train, test, [4 sens/spec columns], method.
compute_auc_df <- function(dir_path, method_name, mdf, sens_spec = FALSE) {
  scores <- read.csv(file.path(dir_path, method_name, "classifier_scores.csv"))
  scores_df <- merge(scores, mdf, by = "sample_id")
  has_sens_spec <- all(c("sensitivity", "specificity") %in% names(scores))
  if (sens_spec && !has_sens_spec) {
    warning(method_name, ": no sensitivity/specificity columns in classifier_scores.csv, filling with NA")
  }

  res <- lapply(0:max(scores$iteration), function(i) {
    train_mdf <- scores_df %>% dplyr::filter(iteration == i, train_or_test == "train")
    test_mdf  <- scores_df %>% dplyr::filter(iteration == i, train_or_test == "test")

    roc_train <- suppressMessages(pROC::roc(response = train_mdf$response,
                                            predictor = train_mdf$classifier_score,
                                            percent = TRUE, quiet = TRUE))
    roc_test  <- suppressMessages(pROC::roc(response = test_mdf$response,
                                            predictor = test_mdf$classifier_score,
                                            percent = TRUE, quiet = TRUE))

    row <- data.frame(iteration = i,
                      train     = as.numeric(roc_train$auc),
                      test      = as.numeric(roc_test$auc))
    if (sens_spec) {
      pick <- function(d, col) if (has_sens_spec) d[[col]][1] else NA
      row$train_sensitivity <- pick(train_mdf, "sensitivity")
      row$train_specificity <- pick(train_mdf, "specificity")
      row$test_sensitivity  <- pick(test_mdf, "sensitivity")
      row$test_specificity  <- pick(test_mdf, "specificity")
    }
    row$method <- method_name
    row
  })
  dplyr::bind_rows(res)
}

#' compute_auc_df() for several models, stacked.
#'
#' Used for the Fig3A violins, SFig7A and SFig3B/SFig17C.
roc_by_iteration <- function(dir_path, methods, mdf, sens_spec = FALSE) {
  dplyr::bind_rows(lapply(methods, compute_auc_df, dir_path = dir_path,
                          mdf = mdf, sens_spec = sens_spec))
}

#' Test-set AUCs (0-1) of all iterations in one classifier_scores.csv.
#'
#' Used for SFig8 and SFig10A/B. Iterations whose test set lacks one of the
#' classes are dropped.
#' @return numeric vector of AUCs, or NULL (with a message) if the file is missing.
compute_aucs <- function(csv_path, mdf) {
  if (!file.exists(csv_path)) {
    message("Not found: ", csv_path)
    return(NULL)
  }
  scores_df <- read.csv(csv_path) %>%
    dplyr::inner_join(mdf[, c("sample_id", "response")], by = "sample_id") %>%
    dplyr::filter(train_or_test == "test")

  iters <- sort(unique(scores_df$iteration))
  aucs  <- sapply(iters, function(i) {
    d <- dplyr::filter(scores_df, iteration == i)
    if (nrow(d) < 2 || dplyr::n_distinct(d$response) < 2) return(NA_real_)
    suppressMessages(as.numeric(pROC::roc(d$response, d$classifier_score, quiet = TRUE)$auc))
  })
  aucs[!is.na(aucs)]
}

#' Median and empirical 95% interval (plus quartiles and n) of AUCs per group.
auc_ci_summary <- function(df, value = "test", by = c("method", "group")) {
  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) %>%
    dplyr::summarise(
      median_auc = median(.data[[value]]),
      ci_low     = quantile(.data[[value]], 0.025),
      ci_high    = quantile(.data[[value]], 0.975),
      q25        = quantile(.data[[value]], 0.25),
      q75        = quantile(.data[[value]], 0.75),
      n          = dplyr::n(),
      .groups    = "drop"
    )
}

#' Horizontal AUC forest plot: point estimate with interval, dashed line at 0.5.
#'
#' Defaults give the SFig8 / SFig10A-B style (point size 0.5, x breaks every
#' 25%, 6 pt axis titles). SFig7A is forest_plot(df, colors = METHOD_COLORS,
#' x_lab = "Test AUC (95% empirical CI)", point_size = 0.4, x_breaks = NULL,
#' small_axis_title = FALSE).
#' @param df     one row per model; the factor levels of `y` set the row order.
#' @param colors named colours for the levels of `y`.
forest_plot <- function(df, x = "median_auc", y = "method", xmin = "ci_low", xmax = "ci_high",
                        colors = FOREST_COLORS, x_lab = "Test AUC", y_lab = NULL,
                        point_size = 0.5, x_breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                        small_axis_title = TRUE) {
  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]], xmin = .data[[xmin]],
                      xmax = .data[[xmax]], color = .data[[y]])) +
    geom_pointrange(linewidth = 0.5, size = point_size) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    scale_color_manual(values = colors) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                       breaks = if (is.null(x_breaks)) waiver() else x_breaks) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(x = x_lab, y = y_lab) +
    theme_cfrna_print() +
    theme_axis_lines()
  if (small_axis_title) {
    p <- p + theme(axis.title = element_text(size = 6), axis.text = element_text(size = 6))
  }
  p
}

# ---- LASSO feature selection --------------------------------------------------

#' Read selected_features.csv of a LASSO run and add gene_name/gene_type.
#' @param gene_key data.frame with gene_id, gene_name, gene_type (GENCODE key).
read_lasso_features <- function(path, gene_key) {
  read.csv(path) %>%
    dplyr::rename(gene_id = feature) %>%
    dplyr::left_join(gene_key[, c("gene_id", "gene_name", "gene_type")], by = "gene_id")
}

#' Selection count (n) and percentage of iterations (pct) per gene, most frequent first.
lasso_selection_table <- function(features) {
  num_iter <- max(features$iteration) + 1
  features %>%
    dplyr::group_by(gene_name) %>%
    dplyr::summarise(n = dplyr::n(), pct = (dplyr::n() / num_iter) * 100) %>%
    dplyr::arrange(dplyr::desc(n))
}

#' Bar plot of the top_n most frequently selected transcripts, labelled with
#' the percentage of models (Fig3B, SFig3C, SFig12A).
#' @return ggplot; attr(p, "gene_order") holds the top genes, most frequent first.
lasso_selection_barplot <- function(features, top_n = 10,
                                    title = "LASSO Transcript Selection") {
  lfts <- lasso_selection_table(features)
  gene_order <- lfts %>% head(top_n) %>% dplyr::pull(gene_name)

  p <- lfts %>%
    dplyr::filter(gene_name %in% gene_order) %>%
    dplyr::mutate(gene_name = factor(gene_name, levels = rev(gene_order))) %>%
    ggplot(aes(y = reorder(gene_name, n), x = n)) +
    geom_bar(stat = "identity", fill = GENESET_COLOR) +
    geom_text(aes(label = sprintf("%.1f%%", pct)), color = "black", hjust = -0.2, size = 2) +
    theme_cfrna_print() +
    labs(x = "Number of Models", y = "", title = title) +
    theme(axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 6),
          plot.title  = element_text(size = 8, hjust = 0.5)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.35)))
  attr(p, "gene_order") <- gene_order
  p
}

#' Box plot of LASSO coefficients across iterations for `gene_order`
#' (right half of Fig3B and SFig3C).
lasso_coef_boxplot <- function(features, gene_order) {
  features %>%
    dplyr::filter(gene_name %in% gene_order) %>%
    dplyr::mutate(gene_name = factor(gene_name, levels = rev(gene_order))) %>%
    ggplot(aes(y = gene_name, x = coefficient)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    geom_boxplot(width = 0.6, outlier.size = 0.5) +
    labs(title = "Coefficients", x = "Coefficient", y = NULL) +
    theme_cfrna_print() +
    theme(panel.grid.major = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 8))
}

#' Fraction of iterations in which each pair of the top_n transcripts was
#' selected together (diagonal = selection frequency).
#' @return top_n x top_n matrix, genes ordered by selection frequency.
coselection_matrix <- function(features, top_n = 10) {
  num_iter <- max(features$iteration) + 1
  top_genes <- features %>%
    dplyr::group_by(gene_name) %>%
    dplyr::summarise(freq = dplyr::n() / num_iter, .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(freq)) %>%
    head(top_n) %>%
    dplyr::pull(gene_name)

  full_mat <- features %>%
    dplyr::filter(gene_name %in% top_genes) %>%
    dplyr::distinct(iteration, gene_name) %>%
    dplyr::mutate(iteration = factor(iteration, levels = 0:(num_iter - 1)),
                  gene_name = factor(gene_name, levels = top_genes)) %>%
    with(table(iteration, gene_name)) %>%
    unclass()
  crossprod(full_mat) / num_iter
}

#' Tile heatmap of a co-selection matrix, rows and columns ordered by
#' hierarchical clustering on 1 - frequency (SFig11).
coselection_heatmap <- function(cosel_mat, n_iter = 101) {
  ord <- with(hclust(as.dist(1 - cosel_mat)), labels[order])
  cosel_long <- reshape2::melt(cosel_mat[ord, ord], varnames = c("gene1", "gene2"),
                               value.name = "freq") %>%
    dplyr::filter(!is.na(freq)) %>%
    dplyr::mutate(gene1 = factor(gene1, levels = ord),
                  gene2 = factor(gene2, levels = ord))

  ggplot(cosel_long, aes(x = gene1, y = gene2, fill = freq)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", freq)), size = 2, color = "grey20") +
    scale_fill_gradientn(colors = c("white", "#DEEBF7", "#9ECAE1", "#3182BD", "#08519C"),
                         limits = c(0, 1), name = "Co-selection\nfrequency") +
    coord_fixed() +
    labs(x = NULL, y = NULL,
         title = sprintf("Feature co-selection frequency (%d MC iterations)", n_iter)) +
    theme_minimal(base_size = 8) +
    theme(axis.text.x      = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y      = element_text(size = 7),
          panel.grid       = element_blank(),
          plot.title       = element_text(size = 8, hjust = 0.5),
          legend.key.width = unit(0.3, "cm"))
}
