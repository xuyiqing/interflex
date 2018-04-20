\name{interflex-package}
\alias{interflex}
\docType{package}
\title{Multiplicative Interaction Models Diagnostics and Visualization}
\description{Producing Flexible Marginal Effect Estimates with Multiplicative Interaction Models}

\details{This package performs diagnostics and visualizations of
  multiplicative interaction models. Besides conventional linear
  interaction models, it provides two additional estimation
  strategies--linear regression based on pre-specified bins and locally
  linear regressions based on Gaussian kernels--to flexibly estimate the
  conditional marginal effect of a treatment variable on an outcome
  variable across different values of a moderatoring variable. These
  approaches relax the linear interaction effect assumption and
  safegiard against excessive extrapolation.}

\author{Jens Hainmueller; Jonathan Mummolo; Yiqing Xu (Maintainer).\cr
  Please report bugs to <yiqingxu@ucsd.edu>.}

\references{Jens Hainmueller; Jonathan Mummolo; Yiqing
Xu. 2018. "How Much Should We Trust Estimates from Multiplicative Interaction Models? Simple Tools to Improve Empirical Practice." Political Analysis, forthcoming.
Available at SSRN: \url{https://papers.ssrn.com/abstract_id=2739221}.
}

\seealso{
  \code{\link{inter.raw}}, \code{\link{inter.binning}}, 
  \code{\link{inter.gam}}, and \code{\link{inter.kernel}}
}
