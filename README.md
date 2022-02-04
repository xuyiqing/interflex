
<!-- README.md is generated from README.Rmd. Please edit that file -->

# interflex

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://www.tidyverse.org/lifecycle/#stablel)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**interflex** produces flexible marginal effect estimates with
multiplicative interaction models.

**Authors:** Jens Hainmueller, Jonathan Mummolo, [Yiqing
Xu](https://yiqingxu.org/), and Ziyi Liu

**Date:** April 1, 2021

**Repos:** [Github](https://github.com/xuyiqing/interflex) (1.2.6)
[CRAN](https://cran.r-project.org/web/packages/interflex/index.html)
(1.2.6)

**Tutoralis:** See tutorials for cases with
[continuous](articles/continuous.html) and
[discrete](articles/discrete.html) outcomes.

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
install.packages('interflex', type = "source", repos = 'http://cran.us.r-project.org') 
```

Or you can install the up-to-date development version from Github:

``` r
install.packages('devtools', repos = 'http://cran.us.r-project.org') # if not already installed
devtools::install_github('xuyiqing/interflex')
```

**interflex** depends on the following packages, which will be installed
automatically when **interflex** is being installed; if not, please
install them manually:

``` r
require(Rcpp) # for processing C++ code
require(mgcv) # for GAM
require(sandwich) # for calculating standard errors
require(pcse) # in case panel-corrected standard errors will be used
require(foreach)  # for parallel computing in kernel estimation
require(doParallel) # for parallel computing in kernel estimation
require(lmtest) # for wald test
require(lfe) # for fixed effects estimations
require(Lmoments) # for L-kurtosis measure
require(ggplot2)  # for plotting
require(plotrix) # for plotting
require(grid) # for plotting
require(gridExtra) # for plotting
require(ggplotify) # for plotting
require(RColorBrewer) # for plotting
require(grDevices) # for plotting
require(gtable) # for plotting
require(MASS) # for wald test
require(mvtnorm) # for simulation
require(pROC) # for auc
require(ModelMetrics) # for cross entropy
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
