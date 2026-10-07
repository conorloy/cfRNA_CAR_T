# Shared plot builders and the figure writer. Requires ggplot2 and ggpubr.

#' Arguments for every Wilcoxon rank-sum test of the paper: the normal
#' approximation (exact = FALSE). Pass as
#' `method.args = WILCOX_ARGS` to ggpubr, or `exact = FALSE` to wilcox.test().
WILCOX_ARGS <- list(exact = FALSE)

#' Box plot + dodged jittered points of a value by x, coloured by group
#' (default: R vs NR), optionally faceted and with a Wilcoxon p-value.
#'
#' Arguments that give each paper panel:
#'   Fig2A LN box: x = "cohort", facet = "dataset", strip_size = 8
#'     (the p-value brackets are added in the notebook with geom_bracket)
#'   SFig6A FMAC/TEX: facet = "dataset", strip_size = 6,
#'     box_linewidth = 0.25, jitter_width = 0.2
#'   SFig3A: x = facet = "Lymphoma", strip_size = 6,
#'     strip_margin = TRUE, show_title = FALSE, stat_method = "wilcox.test"
#'   SFig5 (x3): x = "celltype", jitter_width = 0.2,
#'     title_size = 8, axis_lines = TRUE
#'   Fig4B, SFig16: x = "cohort", title = gene, y_lab = "TMM CPM",
#'     jitter_width = 0.2, title_size = 8, axis_lines = TRUE, shade_x = c(2.5, 3.5)
#'     (Fig4B adds stat_method = "wilcox.test", stat_label_y = max(y) + 5% of range)
#'   Fig4C: facet = "dataset", strip_size = 6, axis_lines = TRUE,
#'     stat_method = "wilcox.test"
#'   SFig14: group = "group", box_dodge = NULL,
#'     stat_method = "wilcox.test", stat_size = 2.5, stat_family = "Helvetica",
#'     stat_label_y = 1.05 * max(y); then add a theme_classic() theme
#'     (a complete theme replaces the one set here)
#'   SFig19: x = group = "subtype", jitter_dodge = FALSE, jitter_width = 0.2,
#'     box_linewidth = 0.25, facet = "dataset", strip_size = 6, strip_margin = TRUE,
#'     stat_method = "wilcox.test", stat_size = 2.5, stat_label_x_npc = "center",
#'     palette = c(GCB = "#5B7F62", "Non-GCB" = "#C17B4A")
#'
#' @param df            data.frame in long format.
#' @param x,y,group     column names (x axis, value, colour/fill).
#' @param palette       named colours for the levels of `group`.
#' @param y_lab,title   axis label and plot title.
#' @param facet         column for facet_wrap(scales = "free", nrow = 1), or NULL.
#' @param box_linewidth outline width of the boxes.
#' @param box_dodge     dodge width of the boxes; NULL keeps ggplot's default dodging.
#' @param jitter_width  horizontal jitter of the points.
#' @param jitter_dodge  TRUE: position_jitterdodge (seed 42, no vertical jitter);
#'                      FALSE: plain geom_jitter (one group per x).
#' @param point_size    point size.
#' @param strip_size    facet label size; NULL leaves strips blank.
#' @param strip_margin  add a 6 pt bottom margin under strips, transparent strip background.
#' @param show_title    FALSE blanks the plot title.
#' @param title_size    title size; NULL keeps the theme default (8 pt).
#' @param axis_lines    apply theme_axis_lines().
#' @param shade_x       c(xmin, xmax) of a grey background band (e.g. the validation cohort).
#' @param stat_method   ggpubr::stat_compare_means test between groups (e.g.
#'                      "wilcox.test"); NULL for none. The label is the p-value.
#'                      Wilcoxon p-values use the normal approximation
#'                      (exact = FALSE).
#' @param stat_size,stat_label_y,stat_family size, y position and font of the p-value.
#' @param stat_label_x_npc horizontal position of the p-value in npc units
#'                      (e.g. "center"); NULL keeps ggpubr's default.
response_boxplot <- function(df, x, y, group = "response", palette = RESPONSE_COLORS,
                             y_lab = NULL, title = NULL, facet = NULL,
                             box_linewidth = 0.4, box_dodge = 0.8,
                             jitter_width = 0.3, jitter_dodge = TRUE, point_size = 1.25,
                             strip_size = NULL, strip_margin = FALSE,
                             show_title = TRUE, title_size = NULL, axis_lines = FALSE,
                             shade_x = NULL, stat_method = NULL, stat_size = 2,
                             stat_label_y = NULL, stat_family = NULL,
                             stat_label_x_npc = NULL) {
  fill_col <- aes(fill = .data[[group]], color = .data[[group]])
  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]]))

  # geom_blank() trains a discrete x scale first; with geom_rect as the first
  # layer ggplot2 >= 3.4 errors. One rect is drawn per row, so the stacked alpha gives a near-opaque grey band.
  if (!is.null(shade_x)) {
    p <- p + geom_blank() +
      geom_rect(aes(xmin = !!shade_x[1], xmax = !!shade_x[2], ymin = -Inf, ymax = Inf),
                fill = "grey90", alpha = 0.5)
  }

  box_args <- list(fill_col, outlier.shape = NA, width = 0.7, alpha = 0.25,
                   linewidth = box_linewidth)
  if (!is.null(box_dodge)) box_args$position <- position_dodge(width = box_dodge)
  p <- p + do.call(geom_boxplot, box_args)

  if (jitter_dodge) {
    p <- p + geom_point(fill_col, size = point_size, alpha = 0.75, stroke = NA,
                        position = position_jitterdodge(jitter.height = 0,
                                                        jitter.width = jitter_width,
                                                        seed = 42))
  } else {
    p <- p + geom_jitter(fill_col, size = point_size, alpha = 0.75, stroke = NA,
                         width = jitter_width)
  }

  if (!is.null(stat_method)) {
    stat_args <- list(aes(group = .data[[group]], label = after_stat(p.format)),
                      method = stat_method, size = stat_size)
    if (stat_method == "wilcox.test") stat_args$method.args <- WILCOX_ARGS
    if (!is.null(stat_label_y)) stat_args$label.y <- stat_label_y
    if (!is.null(stat_family)) stat_args$family <- stat_family
    if (!is.null(stat_label_x_npc)) stat_args$label.x.npc <- stat_label_x_npc
    p <- p + do.call(ggpubr::stat_compare_means, stat_args)
  }

  title_el <- if (!show_title) element_blank()
              else if (is.null(title_size)) element_text(hjust = 0.5)
              else element_text(hjust = 0.5, size = title_size)
  p <- p +
    labs(y = y_lab, title = title) +
    theme_cfrna_print() +
    theme(plot.title = title_el,
          panel.grid.major.x = element_blank(),
          axis.title.x = element_blank())

  if (!is.null(strip_size)) {
    p <- p + if (strip_margin) {
      theme(strip.text = element_text(size = strip_size, margin = margin(b = 6)),
            strip.background = element_rect(fill = NA, color = NA))
    } else {
      theme(strip.text = element_text(size = strip_size))
    }
  }
  if (axis_lines) p <- p + theme_axis_lines()

  p <- p + scale_fill_manual(values = palette) + scale_color_manual(values = palette)
  if (!is.null(facet)) p <- p + facet_wrap(facet, scales = "free", nrow = 1)
  p
}

#' Write a figure panel to plot_path(name) and return the path invisibly.
#'
#' `name` is the file name, e.g. "Fig2_LN_module.pdf" (".pdf" is added when no
#' extension is given); spaces are not allowed. The plot must be passed
#' explicitly: ggplot/patchwork objects, grobs/gtables, ggsurvplot objects
#' (their $plot is saved) and pheatmap objects (their $gtable is saved).
save_plot <- function(p, name, width, height, ...) {
  if (grepl("\\s", name)) stop("Plot file name must not contain spaces: '", name, "'")
  if (!nzchar(tools::file_ext(name))) name <- paste0(name, ".pdf")
  if (inherits(p, "ggsurvplot")) p <- p$plot
  if (inherits(p, "pheatmap")) p <- p$gtable
  path <- plot_path(name)
  ggsave(path, plot = p, width = width, height = height, ...)
  invisible(path)
}
