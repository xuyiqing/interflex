\name{interflex-package}
\alias{interflex-package}
\docType{package}
\title{Multiplicative Interaction Models Diagnostics and Visualization}
\description{Producing Flexible Marginal Effect Estimates with Multiplicative Interaction Models}

\details{This package performs diagnostics and visualizations of multiplicative interaction models. Besides conventional linear interaction models, it provides two additional estimation strategies--linear regression based on pre-specified bins and locally linear regressions based on Gaussian kernels--to flexibly estimate the conditional marginal effect of a treatment variable on an outcome variable across different values of a moderating variable. These approaches relax the linear interaction effect assumption and safeguard against excessive extrapolation.}

\author{
  Jens Hainmueller; Jonathan Mummolo; Yiqing Xu (Maintainer); Ziyi Liu
}

\references{
  Jens Hainmueller; Jonathan Mummolo; Yiqing Xu. 2019. "How Much Should We Trust Estimates from Multiplicative Interaction Models? Simple Tools to Improve Empirical Practice." Political Analysis, Vol. 27, Iss. 2, April 2019, pp. 163--192. Available at: \url{http://bit.ly/HMX2019}.
}

\seealso{
  \code{\link{interflex}}, \code{\link{inter.raw}}, \code{\link{inter.gam}}, \code{\link{plot.interflex}}, and \code{\link{predict.interflex}}
}
