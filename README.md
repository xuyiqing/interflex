# interflex

<!-- badges: start -->
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/interflex)](https://www.datasciencemeta.com/rpackages)
<!-- badges: end -->

**interflex** performs estimation, diagnostics, and visualization of
**conditional marginal effects (CME)** and **group average treatment effects
(GATE)** of a treatment on an outcome across different values of a
moderator. It supports linear, kernel, binning, lasso, and double/debiased
machine-learning (DML) estimators, with optional fixed effects, multi-arm
treatments, and bootstrap or simulation-based inference.

**Maintainer:** [Yiqing Xu](https://yiqingxu.org/) (yiqingxu@stanford.edu)

**Authors:** Jens Hainmueller, Jonathan Mummolo, Yiqing Xu, Jiehan Liu,
Ziyi Liu, Licheng Liu, Tianzhu Qin

**Repos:** [GitHub](https://github.com/xuyiqing/interflex) ·
[CRAN](https://CRAN.R-project.org/package=interflex)

**Documentation:** [User manual (Quarto book)](https://yiqingxu.org/packages/interflex/)
covering installation, classic estimators, extensions, lasso, DML, discrete
moderators, and a full plot-options cyclopedia.

**Reference:** Hainmueller, J., Mummolo, J., & Xu, Y. (2019). [How Much
Should We Trust Estimates from Multiplicative Interaction Models? Simple
Tools to Improve Empirical Practice](https://www.cambridge.org/core/journals/political-analysis/article/how-much-should-we-trust-estimates-from-multiplicative-interaction-models-simple-tools-to-improve-empirical-practice/D8CAACB473F9B1EE256F43B38E458706).
*Political Analysis*, 27(2), 163–192.

---

## Installation

Install the development version from GitHub:

```r
# install.packages("devtools")  # if needed
devtools::install_github("xuyiqing/interflex")
```

Or install the released version from CRAN:

```r
install.packages("interflex")
```

`interflex` depends on a number of CRAN packages that are installed
automatically. If you need to install them manually:

```r
pkgs <- c(
  "Rcpp", "mgcv", "sandwich", "pcse", "foreach", "doParallel", "doFuture",
  "doRNG", "future.apply", "lfe", "lmtest", "Lmoments", "ggplot2", "ggpubr",
  "ggplotify", "gridExtra", "RColorBrewer", "MASS", "mvtnorm", "pROC",
  "ModelMetrics", "patchwork", "xgboost"
)
install.packages(setdiff(pkgs, rownames(installed.packages())))
```

### macOS toolchain notes

If you hit `-lgfortran` / `-lquadmath` errors during compilation, install
the official R macOS toolchain from
[CRAN tools](https://cran.r-project.org/bin/macosx/tools/).

If you hit `clang: error: unsupported option '-fopenmp'`, update R and the
Xcode command-line tools, or install the
[r-macos-rtools](https://github.com/coatless-mac/r-macos-rtools/releases)
bundle.

---

## Quick start

```r
library(interflex)
data(interflex)

# Linear estimator on simulated data with two covariates
fit <- interflex(
  estimator = "linear",
  data      = s5,
  Y         = "Y",
  D         = "D",
  X         = "X",
  Z         = c("Z1", "Z2"),
  vartype   = "bootstrap"
)

# Visualize the conditional marginal effect curve
plot(fit)
```

For an empirical example using `app_hma2015`:

```r
fit <- interflex(
  estimator = "linear",
  data      = app_hma2015,
  Y         = "totangry",
  D         = "threat",
  X         = "pidentity",
  Z         = c("issuestr2", "knowledge", "educ", "male", "age10"),
  Ylabel    = "Anger",
  Dlabel    = "Threat",
  Xlabel    = "Partisanship",
  vartype   = "bootstrap"
)

plot(fit, xlim = c(0.25, 1), ylim = c(-0.2, 0.6))
```

The companion [user manual](https://yiqingxu.org/packages/interflex/) walks
through every estimator, the diagnostic plots, multi-arm treatments,
fixed-effect models, lasso and DML, discrete moderators, and the full set
of plot options.

---

## Reporting bugs

Please open an issue on
[GitHub](https://github.com/xuyiqing/interflex/issues) with a minimal
reproducible example, or email **yiqingxu \[at\] stanford.edu**.
