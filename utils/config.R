# Path configuration shared by all R notebooks.
# Reads <repo>/config.yaml, then the optional, untracked <repo>/config.local.yaml, whose keys
# override those of config.yaml. Exposes in_path(), res_path(), res_in() and plot_path().

find_repo_root <- function(start = getwd()) {
  dir <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(dir, "config.yaml"))) return(dir)
    parent <- dirname(dir)
    if (parent == dir) stop("config.yaml not found above ", start)
    dir <- parent
  }
}

REPO_ROOT <- find_repo_root()

.resolve <- function(p) {
  if (grepl("^(/|~)", p)) path.expand(p) else file.path(REPO_ROOT, p)
}

.cfg <- yaml::read_yaml(file.path(REPO_ROOT, "config.yaml"))
.local_cfg <- file.path(REPO_ROOT, "config.local.yaml")
if (file.exists(.local_cfg)) {
  .local <- yaml::read_yaml(.local_cfg)
  if (is.list(.local)) .cfg <- utils::modifyList(.cfg, .local)
  rm(.local)
}
rm(.local_cfg)
INPUT_DIR   <- normalizePath(.resolve(.cfg$input_dir), mustWork = FALSE)
RESULTS_DIR <- normalizePath(.resolve(.cfg$results_dir), mustWork = FALSE)
PLOTS_DIR   <- normalizePath(.resolve(.cfg$plots_dir), mustWork = FALSE)

if (!dir.exists(INPUT_DIR)) warning("input_dir does not exist: ", INPUT_DIR,
                                    " (edit config.yaml or config.local.yaml)")

# File in the input data directory; errors early if it is missing.
in_path <- function(...) {
  p <- file.path(INPUT_DIR, ...)
  if (!file.exists(p)) stop("Input file not found: ", p)
  p
}

# File in the results directory, for writing; parent folders are created on demand.
res_path <- function(...) {
  p <- file.path(RESULTS_DIR, ...)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  p
}

# Existing file or folder in the results directory, for reading; errors if it is missing
# and creates nothing.
res_in <- function(...) {
  p <- file.path(RESULTS_DIR, ...)
  if (!file.exists(p)) stop("Result not found: ", p,
                            " (run the code/01_compute or code/00_scrna_prep notebook that writes it)")
  p
}

# File in the plots directory.
plot_path <- function(...) {
  p <- file.path(PLOTS_DIR, ...)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  p
}
