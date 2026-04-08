
# List of packages
packages <- c("Rcpp", "mgcv", "sandwich", "pcse", "foreach",
              "doParallel", "lfe", "lmtest", "Lmoments","ggplot2",
              "plotrix", "grid", "gridExtra", "ggplotify", "ggpubr", "labelled",
              "RColorBrewer", "grDevices", "gtable", "MASS", "mvtnorm",
              "pROC", "ModelMetrics", "foreign","patchwork", "rmarkdown", "DT",
              "xgboost", "interflex")

# Install and load each package
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

# Reproducibility under parallel computing.
#
# interflex's bootstrap / CV / DML paths dispatch via foreach + future.apply
# with `.options.future = list(seed = TRUE)`. future.apply derives a
# deterministic per-iteration L'Ecuyer-CMRG stream from the caller's RNG state,
# so a single `set.seed()` before each stochastic interflex() call is
# sufficient to make parallel results reproducible. We enforce this globally
# by switching to L'Ecuyer-CMRG and installing a knitr chunk hook that seeds
# every chunk from a label-derived seed, so results are stable regardless of
# execution order, cache state, or cores count.
RNGkind("L'Ecuyer-CMRG", "default", "default")
set.seed(1234)

if (requireNamespace("knitr", quietly = TRUE)) {
  .interflex_chunk_seed <- function(before, options, envir) {
    if (isTRUE(before)) {
      label <- if (is.null(options$label)) "" else as.character(options$label)
      # Deterministic int from label so each chunk gets a distinct seed.
      offset <- if (nzchar(label)) sum(utf8ToInt(label)) else 0L
      set.seed(1234L + offset)
    }
  }
  knitr::knit_hooks$set(chunk_seed = .interflex_chunk_seed)
  knitr::opts_chunk$set(chunk_seed = TRUE)
}

