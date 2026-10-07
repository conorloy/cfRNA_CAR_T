# ggplot2 themes and colour palettes shared by all figure notebooks.
# Requires ggplot2 (attached by setup.R).

# ---- Themes -----------------------------------------------------------------

#' Print theme used for every paper panel: Helvetica, 8 pt titles, 6 pt axis
#' text, black border, no grid, no strips, no legend.
theme_cfrna_print <- function() {
  theme_bw() %+replace%
    theme(
      panel.border     = element_rect(color = "black", fill = NA),
      title            = element_text(family = "Helvetica", size = 8, color = "black"),
      axis.title       = element_text(family = "Helvetica", size = 8, color = "black"),
      axis.ticks       = element_line(color = "black"),
      axis.text        = element_text(family = "Helvetica", size = 6, color = "black"),
      legend.position  = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      strip.text       = element_blank(),
      strip.background = element_blank()
    )
}

#' theme_cfrna_print() for Kaplan-Meier and PCA panels.
#'
#'   major_grid = TRUE  -> Fig1 panels (PFS by cohort, PCA): adds light grey95
#'                         major gridlines.
#'   major_grid = FALSE -> Fig2 LN-split KM and SFig6B: identical to
#'                         theme_cfrna_print().
theme_cfrna_print_no_minor <- function(major_grid = TRUE) {
  th <- theme_cfrna_print()
  if (major_grid) {
    th <- th + theme(panel.grid.major = element_line(color = "grey95", linewidth = 0.5))
  }
  th
}

#' Remove the panel border and draw thin axis lines instead.
#'
#' Equivalent to
#'   theme(panel.border = element_blank(), axis.line = element_line(size = .25))
#' (`size` and `linewidth` are equivalent for lines in ggplot2 >= 3.4).
#' axis = "x" draws only the x axis line (the LASSO coefficient panels).
theme_axis_lines <- function(linewidth = 0.25, axis = c("xy", "x")) {
  axis <- match.arg(axis)
  line <- element_line(linewidth = linewidth)
  if (axis == "x") {
    theme(panel.border = element_blank(), axis.line.x = line)
  } else {
    theme(panel.border = element_blank(), axis.line = line)
  }
}

# ---- Palettes ---------------------------------------------------------------

#' Responders vs non-responders. Used in every R-vs-NR box plot (Fig2, Fig4,
#' SFig3A, SFig5, SFig6A, SFig14, SFig16), the volcano plot (Fig4, up = R,
#' down = NR), the GCB/non-GCB boxes (SFig19) and the SFig18A Venn.
RESPONSE_COLORS <- c(R = "#5B7F62", NR = "#C17B4A")

#' Cohorts 1-3. Fig1 PFS curves, read/feature-count violins and PCA, and the
#' Fig2 LN-split KM curves (each colour repeated for Low/High). The same three
#' colours are reused for sequencing batches 1-3 in SFig1A.
COHORT_COLORS <- c("Cohort 1" = "#8B6F4E", "Cohort 2" = "#6B8E8C", "Cohort 3" = "#D4A574")

#' Distinct colour per classifier / gene set. SFig7A CI forest and SFig7B
#' calibration; SFig9A adjustment forest (keyed by signature column; see
#' SIGNATURE_COLORS).
METHOD_COLORS <- c("All Transcripts" = "#555555", "LN" = "#5D7A8C",
                   "FMAC" = "#8B7355", "TEX" = "#D4A574")

#' Fig3A test-AUC violins. Differs from METHOD_COLORS:
#' All Transcripts is "grey" (not #555555) and Random Genes is added.
METHOD_COLORS_VIOLIN <- c("All Transcripts" = "grey", "LN" = "#5D7A8C", "FMAC" = "#8B7355",
                          "TEX" = "#D4A574", "Random Genes" = "#8d5519")

#' METHOD_COLORS keyed by score column name (SFig9A).
SIGNATURE_COLORS <- c(LN_set1 = "#5D7A8C", FMAC_set1 = "#8B7355", TEX_set1 = "#D4A574")

#' Two-colour scheme of the AUC forest plots: gene-set and control models vs
#' models that use clinical variables (SFig8, SFig10A/B). GENESET_COLOR is also the fill of the LASSO selection bars
#' (Fig3, SFig3C, SFig12A) and the LN-vs-LDH regression line (SFig4).
GENESET_COLOR  <- "#5D7A8C"
CLINICAL_COLOR <- "#CE771B"

#' Univariate forest plot, signature scores vs baseline covariates (SFig9B).
UNIVARIATE_GROUP_COLORS <- c("Signature score" = "#0e5785ff", "Baseline covariate" = "#555555")

#' "All patients" vs sensitivity subset violins: SFig3B ("De novo DLBCL only")
#' and SFig17C (">12mo follow-up"). Name it with
#' setNames(SENSITIVITY_COLORS, c("All patients", "<subset label>")).
SENSITIVITY_COLORS <- c("#555555", "#2c7bb6")

#' LymphoMAP archetype groups, SFig6B KM curves. The order
#' matches the alphabetical factor levels (FMAC orange, LN blue, TEX green).
LYMPHOMAP_COLORS <- c(FMAC = "#E69F00", LN = "#56B4E9", TEX = "#009E73")

#' Blue of the SFig2 Cox forest plot boxes and the SFig18B scatter points.
FOREST_BLUE <- "#377EB8"
