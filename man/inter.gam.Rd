\name{inter.gam}
\alias{inter.gam}
\title{Visualization of a Generalized Additive Model (GAM)}
\description{Estimating and plotting a GAM before estimating a multiplicative interaction model}
\usage{inter.gam(data, Y, D, X, Z = NULL, weights = NULL, 
  FE = NULL, SE = FALSE, k = 10, angle = c(30, 100,-30,-120), 
  Ylabel = NULL, Dlabel = NULL, Xlabel = NULL)
}
\arguments{
  \item{data}{a data.frame that stores Y, D, X, and other variables.}
  \item{Y}{a string, name of the outcome variable.}
  \item{D}{a string, name of the treatment variable.}
  \item{X}{a string, name of the moderating variable.}
  \item{Z}{a vector of variable names of control variables.}
  \item{weights}{a string, name of the weighting variable.}
  \item{FE}{a vector of variable names of fixed effects indicators.}
  \item{SE}{a logical variable indicating whether to show the surfaces
    one standard error above and below the predicted surface.}
  \item{k}{an integer indicating the dimension of the basis used to
    represent the smooth term. Smaller values of k produce smoother surfaces.}
  \item{angle}{a vector of numbers indicating viewing angles (allows 1 to
    4 values).}
  \item{Ylabel}{a string that controls the label of the outcome variable Y in the plot.}
  \item{Dlabel}{a string that controls the label of the treatment variable D in the plot.}
  \item{Xlabel}{a string that controls the label of the moderating variable X in the plot.} 
}
\details{ 
  \bold{inter.gam} provides a convenient way to visualize the relationships between the outcome, treatment, and moderator using GAM when both the treatment and moderator are continuous (based on \code{\link[mgcv]{vis.gam}} of the \bold{mgcv} package).
}
\author{
  Jens Hainmueller; Jonathan Mummolo; Yiqing Xu (Maintainer)
}
\references{
  Jens Hainmueller; Jonathan Mummolo; Yiqing Xu. 2019. "How Much Should We Trust Estimates from Multiplicative Interaction Models? Simple Tools to Improve Empirical Practice." Political Analysis, Vol. 27, Iss. 2, April 2019, pp. 163--192. Available at SSRN: \url{https://papers.ssrn.com/abstract_id=2739221}.
}
\seealso{
  \code{\link{inter.raw}}, \code{\link{inter.binning}}, and \code{\link{inter.kernel}}
}
\examples{
library(interflex)
data(interflex)
inter.gam(Y="Y", D="D", X="X", Z=c("Z1"), data=s2)
}
\keyword{graphics}


