
<!-- README.md is generated from README.Rmd. Please edit that file -->

# interflex

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://www.tidyverse.org/lifecycle/#stablel)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![downloads:
CRAN](https://cranlogs.r-pkg.org/badges/grand-total/interflex)](https://www.datasciencemeta.com/rpackages)
<!-- badges: end -->

**interflex** estimates, visualizes, and interprets conditional marginal
effects and performs diagnostics.

**Authors:** Jens Hainmueller, Jiehan Liu, Licheng Liu, Ziyi Liu,
Jonathan Mummolo, Tianzhu Qin, and [Yiqing Xu](https://yiqingxu.org/)
(maintainer)

**Date:** Feburary 12, 2025

**Repos:** [Github](https://github.com/xuyiqing/interflex) (1.3.2)
[CRAN](https://cran.r-project.org/web/packages/interflex/index.html)
(1.2.6)

**Tutoralis:** See tutorials for cases with
[continuous](https://yiqingxu.org/packages/interflex/articles/continuous.html)
and
[discrete](https://yiqingxu.org/packages/interflex/articles/discrete.html)
outcomes, as well as a tutorial for double/debiased machine learning
([DML](https://yiqingxu.org/packages/interflex/articles/dml.html))
estimators.

**Examples:** R code used in the tutorials can be downloaded from
[here](examples.R).

**Reference:** [*How Much Should We Trust Estimates from Multiplicative
Interaction Models? Simple Tools to Improve Empirical
Practice*](http://bit.ly/HMX2019). Political Analysis, Vol. 27, Iss. 2,
April 2019, pp. 163–192.

------------------------------------------------------------------------

## Installation

You can install the **interflex** package from CRAN:

``` r
install.packages('interflex', type = "source", 
                 repos = 'http://cran.us.r-project.org') 
```

Or you can install the up-to-date development version from Github:

``` r
# if not already installed
install.packages('devtools', repos = 'http://cran.us.r-project.org') 
# install from github
devtools::install_github('xuyiqing/interflex')
```

**interflex** depends on the following packages, which will be installed
automatically when **interflex** is being installed; if not, please
install them manually:

``` r
# Function to install packages
install_all <- function(packages) {
  installed_pkgs <- installed.packages()[, "Package"]
  for (pkg in packages) {
    if (!pkg %in% installed_pkgs) {
      install.packages(pkg, repos = 'http://cran.us.r-project.org')
    }
  }
}

# Packages to be installed
packages <- c("Rcpp", "mgcv", "sandwich", "pcse", "fixest", "foreach", "doParallel", 
              "lfe", "lmtest", "Lmoments", "ggplot2", "plotrix", "grid", 
              "gridExtra", "ggplotify", "ggpubr", "RColorBrewer", "grDevices", 
              "gtable", "MASS", "mvtnorm", "pROC", "ModelMetrics", "foreign",
              "patchwork", "rmarkdown")

# Install the packages
install_all(packages)
```

Mac users who encounter “-lgfortran” or “-lquadmath” error during
installation, please check out the solution
[here](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/).
Typing the following two lines of code in your **Terminal** should solve
this problem.

``` r
curl -OL http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```

Mac users who encounter **clang: error: unsupported option ‘-fopenmp’**,
please consider (1) updating your R and/or (2) installing new R macro
tools from
[Github](https://github.com/rmacoslib/r-macos-rtools/releases/tag/v3.1.0).

## Report bugs

Please report bugs to **yiqingxu \[at\] stanford.edu** with your sample
code and data file. Much appreciated!
