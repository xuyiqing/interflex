\name{inter.kernel}
\alias{inter.kernel}
\title{The Kernel Estimator}
\description{Implementing the kernel estimator, based on a
  semi-parametric model that relaxes the
  linear interaction effect (LIE) assumption of conventional
  multiplicative interaction models and safeguards against excessive extrapolation
}
\usage{inter.kernel(Y, D, X, Z = NULL, weights = NULL, FE = NULL,
          data, na.rm = FALSE, CI = TRUE, conf.level = 0.95,
          cl = NULL, neval = 50, nboots = 200, 
          parallel = TRUE, cores = 4, seed = 02139,
          bw = NULL, grid = 20, metric = "MSPE",
          Ylabel = NULL, Dlabel = NULL, Xlabel = NULL,
          main = NULL, xlim = NULL, ylim = NULL,
          Xdistr = "histogram", file = NULL)
}
\arguments{
  \item{Y}{a string, name of the outcome variable.}
  \item{D}{a string, name of the treatment variable.}
  \item{X}{a string, name of the moderating variable.}
  \item{Z}{a vector of variable names of control variables.}
  \item{weights}{a string, name of the weighting variable.}
  \item{FE}{a vector of variable names of fixed effects indicators.}
  \item{data}{a data.frame that stores Y, D, and X.}
  \item{na.rm}{a logical flag indicating whether to list-wise delete
    missing data. The algorithm will report an error if missing data
    exist.}
  \item{CI}{a logical flag indicating whether the confidence intervals
    need to be shown. If \bold{CI} = FALSE, the bootstrap procedure
    will not be conducted.}
  \item{conf.level}{a number that specifies the confidence interval. For
    example, 0.9 means a confidence interval of [5\%, 95\%]. The default is 0.95,
    i.e., [2.5\%, 97.5\%].}
  \item{cl}{a string specifying the name of clustering variable. If
    specified, the bootstrap procedure will be blocking on this variable.
    Ignored when \bold{CI} = FALSE.}
   \item{neval}{an integer specifying the number of evaluation
     points. A large number would produce smoother marginal effects
     estimates but the algorithm will take more time to run. The default
     is 50.}
  \item{nboots}{an integer specifying the number of bootstrap runs. The
    default is 200.}
  \item{parallel}{a logical flag indicating whether parallel computing
    will be used in bootstrapping and/or cross-validation. Ignored if
    \bold{CI} = FALSE and \bold{bw} is specified.}
  \item{cores}{an integer indicating the number of cores to be used in
    parallel computing. If not specified, the algorithm will use the
    maximum number of logical cores of your computer (warning: this
    could prevent you from multi-tasking on your computer).}
  \item{seed}{an integer that sets the seed in random number
    generation. Ignored if  \bold{CI} = FALSE and \bold{bw} is specified.}
  \item{bw}{a number that sets the bandwidth in Gaussian kernel reweighting.}
  \item{grid}{either an integer indicating the number of candidates in a grid
    search to find the optimal bandwidth, or a vector of positive
    numbers of candidate bandwidths. Ignored if \bold{bw} is specified.}
  \item{metric}{Either "MSPE" (Mean Squared Prediction Error) or "MAPE"
    (Mean Absolute Prediction Error), criterion to be used in
    cross-validation to select bandwidth. The default is "MSPE".}
  \item{Ylabel}{a string that controls the label of Y in the plot.}
  \item{Dlabel}{a string that controls the label of D in the plot.}
  \item{Xlabel}{a string that controls the label of X in the plot.} 
  \item{main}{a string that control the title of the plot.}
  \item{xlim}{a two-element numeric vector that controls the range of
    the x-axis to be shown in the plot.}
  \item{ylim}{a two-element numeric vector that controls the range of
    the y-axis to be shown in the plot (with small adjustments to
    improve aesthetics).}
  \item{Xdistr}{a string indicating the way the distribution of the
    moderator will be plotted, either "histogram" (or "hist") or
    "density". The default is "histogram".}
  \item{file}{save the output graph to the file.}	
}
\details{
  \bold{inter.kernel} implements a kernel smoothing estimator of the
  marginal effect, which is an application of semi-parametric smooth
  varying-coefficient models. This approach is more
  complicated than the simple binning estimator, but provides a
  generalization that allows researchers to flexibly estimate the
  functional form of the marginal effect of D on Y across the values
  of X by estimating a series of local effects with a Gaussian kernel
  re-weighting scheme.

  The kernel regression nests the standard linear interaction model
  as a special case. However, in the kernel regression the conditional
  effect of D on Y does not have to fall on a linear line as
  required by the linear interaction effect (LIE) assumption,
  but can vary freely across the range
  of X. In addition, if covariates Z are included in the model, the
  coefficients of those covariates are also allowed to vary freely
  across the range of X resulting in a very flexible estimator that
  also helps to guard against misspecification bias with respect to the
  covariates.

  \bold{inter.kernel} produces a smooth marginal effect of D on Y
  with respect to X. It estimates a series of locally
  linear regression using kernel re-weighting based on the distance
  between each value of X (the moderator) and each evaluation point.  If the bandwidth
  is not supplied, it is selected using a standard 5-fold
  cross-validation procedure. Standard errors and confidence intervals
  can be computed using a bootstrap. It also adds a histogram (or
  density plot)  at the bottom of the figure to help users assess the
  common support assumption based on the distribution of the moderator.

  \bold{inter.kernel} can accommodate linear fixed-effect models
  with multiple fixed-effect indicators. A fixed-effect model is
  estimated using a fast algorithm written in C++ at each evaluation
  point of the semi-parametric kernel estimation. 
}
\value{
  \item{bw}{specified or cross-validated bandwith.}
  \item{est}{marginal effects estimates (and confidence intervals if
    \bold{CI = TRUE})}
  \item{CV.out}{results from cross-validation if \bold{bw} is not specified.}
  \item{graph}{stores the graphic output, a \bold{ggplot2} object.}
}
\author{
  Jens Hainmueller; Jonathan Mummolo; Yiqing Xu (Maintainer)
}
\references{Jens Hainmueller; Jonathan Mummolo; Yiqing
Xu. 2016. "How Much Should We Trust Estimates from Multiplicative Interaction Models? Simple Tools to Improve Empirical Practice."
Available at SSRN: \url{https://papers.ssrn.com/abstract_id=2739221}.
}
\seealso{
  \code{\link{inter.raw}}, \code{\link{inter.gam}}, and
  \code{\link{inter.binning}}
}
\examples{
library(interflex)
data(interflex)
inter.kernel(Y = "Y", D = "D", X = "X", Z = "Z1", data = s3,
             nboots = 50, bw = 1.3, parallel = FALSE)
}
\keyword{graphics}


