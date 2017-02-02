\name{inter.raw}
\alias{inter.raw}
\title{Plotting Raw Data}
\description{Plotting raw data before estimating a multiplicative interaction model}
\usage{inter.rawgplot(Y, D, X, weights = NULL, data,
              nbins  = 3, cutoffs = NULL, span = NULL,
              Ylabel = NULL, Dlabel = NULL, Xlabel = NULL)
}
\arguments{
  \item{Y}{a string, name of the outcome variable.}
  \item{D}{a string, name of the treatment variable.}
  \item{X}{a string, name of the moderating variable.}
  \item{weights}{a string, name of the weighting variable. It will
    influence both the linear and loess fits.}
  \item{data}{a data.frame that stores Y, D, and X.}
  \item{nbins}{an integer that determines the number of bins of the
    moderator (and hence the nubmer of subplots) when the treatment is
    continuous. The sample will be divided based on integral multiples
    of 100/nbins percentile.  \bold{nbins} will be automatically
    subtracted by 1 if more than one multiples of the percentile have
    the same value (for example, if the moderator has over 70\% zeros,
    both the 33 and 66 percentiles are zero). The default is 3. Ignored
    when \bold{cutoffs} is supplied.}
  \item{cutoffs}{a vector of numbers that determines how subgroups are
    divided based on the moderator X. When this option is supplied, \bold{nbins}
    will be ignored. The smallest number of the first interval and the
    largest number of the last interval do not need to be
    specified. Numbers smaller than the minimum or larger than the
    maximum of X will be ignored. \bold{nbins} equals the length of this
    vector plus 1. Ignored if the treatment is dichotomous.}
  \item{span}{hyperparameter for the loess fit. Large numbers produce
    smoother lines (from \bold{ggplot2}).}
  \item{Ylabel}{a string that controls the label of Y in the plot.}
  \item{Dlabel}{a string that controls the label of D in the plot.}
  \item{Xlabel}{a string that controls the label of X in the plot.} 
}
\details{ 
  \bold{inter.raw} provides a simple visual diagnostic, a scatterplot,
  to help researchers to detect potential problems with the linear
  interaction effect assumption and the lack of common support.
  
  If the treatment D is binary, it plots outcome Y against the moderator
  X separately for the sample of treatment group observations (D=1)
  and the sample of control group observations (D=0).  In each sample,
  it superimposes a linear regression line as well as loess fit in each group.

  To check whether there is sufficient common support in the data, we
  can simply compare the distribution of X in both groups and examine
  the range of X values for which there are a sufficient number of
  data points in both groups. The box plots near the center of the
  figures display quantiles of the moderator at each level of the
  treatment. The dot in the center denotes the median, the end points of
  the thick bars denote the 25th and 75th percentiles, and the end
  points of the thin bars denote the 5th and 95th percentiles. 

  If the treatment and moderator are continuous, we adopt a simple
  binning approach. Specifically, \bold{inter.raw} will split the
  sample into several subgroups based on the moderator. By default, it
  will split the sample into three roughly equal sized subgroups: low X
  (first tercile), medium X (second tercile), and high X (third
  tercile). For each of the subgroups, it plots Y against D while again
  adding the estimates for both the linear fit and the loess fit.
}
\value{
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
  \code{\link{inter.gam}}, \code{\link{inter.binning}},
  and \code{\link{inter.kernel}}
}
\examples{
library(interflex)
data(interflex)
inter.raw(Y = "Y", D = "D", X = "X", data = s1, 
              Ylabel = "Outcome", Dlabel = "Treatment", Xlabel="Moderator")
}
\keyword{graphics}


