\name{inter.binning}
\alias{inter.binning}
\title{The Binning Estimator}
\description{Implementing the binning estimator, a generalization of the
  multiplicative interaction model and conducting various diagnostic tests}
\usage{inter.binning(Y, D, X, Z = NULL, FE = NULL, weights = NULL,
          data, na.rm = FALSE, nbins = 3, cutoffs = NULL,
          vartype = "homoscedastic", 
          cl = NULL, time = NULL, pairwise = TRUE, figure = TRUE, main = NULL,
          Ylabel = NULL, Dlabel = NULL, Xlabel = NULL,
          xlim = NULL, ylim = NULL, interval = NULL, Xdistr = "histogram", wald = TRUE, theme.bw = FALSE)
}
\arguments{
  \item{Y}{a string, name of the outcome variable.}
  \item{D}{a string, name of the treatment variable.}
  \item{X}{a string, name of the moderating variable.}
  \item{Z}{a vector of variable names of control variables.}
  \item{weights}{a string, name of the weighting variable.}
  \item{FE}{a vector of variable names of fixed effects indicators.}
  \item{data}{a data.frame that stores Y, D, X, and other variables.}
  \item{na.rm}{a logical flag indicating whether to list-wise delete data. The algorithm will report error if missing data
    exist.}
  \item{nbins}{an integer that determines the number of bins of the
    moderator. The sample will be divided based on integral multiples
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
    maximum of X will be ignored. \bold{nbins} equals to the length of this
    vector plus 1. Ignored if the treatment is dichotomous.}
  \item{vartype}{a string that controls the variance-covariance
    estimator. Can be one of the following: "homoscedastic";
    "robust", "cluster", and "pcse" (panel corrected standard errors).
    The default is "homoscedastic". }
  \item{cl}{a string specifying the name of clustering variable for clustered
    standard errors. \bold{vartype} must be either "cluster" or "pcse".
    \bold{vartype} will be automatically set to "cluster" if \bold{cl}
    is supplied and \bold{vartype} is not "pcse". "pcse" is not allowed when \bold{FE} is specified.}
  \item{time}{a string specifying the name of time variable for panel corrected
    standard errors. Ignored when \bold{vartype} is not  "pcse".}
  \item{pairwise}{a logical flag indicating whether to switch on the
    \bold{pairwise} option when estimating panel corrected standard
    errors. Ignored when \bold{vartype} is not "pcse". The default is TRUE.}
  \item{Ylabel}{a string that controls the label of Y in the plot.}
  \item{Dlabel}{a string that controls the label of D in the plot.}
  \item{Xlabel}{a string that controls the label of X in the plot.} 
  \item{main}{a string that control the title of the plot.}
  \item{xlim}{a two-element numeric vector that controls the range of
    the x-axis to be shown in the plot.}
  \item{ylim}{a two-element numeric vector that controls the range of
    the y-axis to be shown in the plot (with small adjustments to
    improve aesthetics).}
  \item{bin.labs}{a logical flag indicating whether to show a label for each bin; 
    only meaningful when \bold{nbins} is 3 or 4. The default is TRUE.}  
  \item{Xdistr}{a string indicating the way the distribution of the
    moderator will be plotted, either "histogram" (or "hist") or
    "density". The default is "histogram".}
  \item{figure}{a logical flag indicating whether to draw the graph. 
    The default is TRUE.}
  \item{interval}{manually set the interval for the graph.}
  \item{wald}{wald test.}
  \item{theme.bw}{a logical flag specifying whether to use a black-white theme.}    
}
\details{
  \bold{inter.binning} implements the binning estimator. There are
  three steps to implement the estimator. First, we discretize the
  moderator X into several bins and create a dummy variable for each
  bin. The default is 3 bins, respectively corresponding to the three
  terciles.  Second, we pick an evaluation point within each bin, where
  we want to estimate the conditional marginal effect
  of D on Y. Typically, we choose the median value of X in each bin.
  Third, we estimate a model that includes interactions between the bin
  dummies and the treatment indicator, the bin dummies and the
  moderator X minus the evaluation points, as well as the triple
  interactions. The last two terms are to capture the effect of
  D on Y within each bin.

  The binning estimator has several key advantages over the standard
  multiplicative interaction model. First, the binning estimator is much more
  flexible as it jointly fits the interaction components of the standard
  model to each bin separately. The model does not impose the linear
  interaction effect (LIE) assumption. Instead, the conditional marginal effects can vary freely
  across the three bins and therefore can take on any non-linear or
  non-monotonic pattern that might describe the heterogeneity in the
  effect of D on Y across low, medium, or high levels of
  X.

  Second, since the bins are constructed based on the support of X,
  the binning ensures that the conditional marginal effects are
  estimated at typical values of the moderator and do not rely on
  excessive extrapolation or interpolation.

  Third, the binning estimator is easy to implement using any regression
  software and the standard errors for the conditional marginal effects
  are directly estimated by the regression so there are no need to
  compute linear combinations of coefficients to recover the conditional
  marginal effects.

  Fourth, the binning estimator actually provides a generalization that
  nests the standard multiplicative interaction model as a special
  case. It can therefore serve as a formal test on the validity of a
  global LIE assumption imposed by the standard model. So in the special
  case when the standard multiplicative interaction model is correct and
  therefore the global LIE assumption holds, then---as the sample size
  grows---the marginal effect estimates from the binning estimator
  converge in probability on the unbiased marginal effect estimates from the standard
  multiplicative interaction model.
  
  In the case of a binary treatment, \bold{inter.binning} displays at the
  bottom of the figure a stacked histogram that shows the distribution
  of the moderator X. In this histogram the total height of the
  stacked bars refers to the distribution of the moderator in the pooled
  sample and the red and white shaded bars refer to the distribution of
  the moderator in the treatment and control groups,
  respectively. Adding such a histogram makes it easy to judge the
  degree to which there is common support in the data. In the case of a
  continuous treatment, it displays a histogram at the bottom that
  simply shows the distribution of X in the entire sample.
}
\value{
  \item{est.binning}{stores the binning estimates}
  \item{binary.treatment}{whether the treatment variable is binary}
  \item{bin.size}{the percentage of observations within each bin}
  \item{treat.variation}{the variance of the treatment variable within each bin}
  \item{X.LKurtosis}{the L-kurtosis measure of the moderator}
  \item{correctOrder}{whether the three binning estimates are in the
    correct order (i.e. monotonically increasing or decreasing)}
  \item{p.twosided}{p values of pairwise t-tests of the binning estimates  when there are 2 or 3 bins.}
  \item{p.wald}{p-value of a Wald test. The NULL hypothesis is that the
  linear interaction model and the binning model are statistically equivalent.}
  \item{graph}{stores the graphic output, a \bold{ggplot2} object.}

}
\author{
  Jens Hainmueller; Jonathan Mummolo; Yiqing Xu (Maintainer)
}
\references{Jens Hainmueller; Jonathan Mummolo; Yiqing
Xu. 2018. "How Much Should We Trust Estimates from Multiplicative Interaction Models? Simple Tools to Improve Empirical Practice." Political Analysis, forthcoming.
Available at SSRN: \url{https://papers.ssrn.com/abstract_id=2739221}.
}


\seealso{
  \code{\link{inter.raw}}, \code{\link{inter.gam}}, and 
 \code{\link{inter.kernel}}
}
\examples{
library(interflex)
data(interflex)
inter.binning(Y = "Y", D = "D", X = "X", Z = "Z1",
                  data = s1, nbins = 3, vartype = "homoscedastic",
                  Ylabel = "Y", Dlabel = "Tr", Xlabel="X")

## specifying cutoffs
inter.binning(Y = "Y", D = "D", X = "X", Z = "Z1",
                  data = s2, cutoffs = c(1,2,4,5))

## with fixed effects
s4[,"wgt"] <- 1
inter.binning(Y = "Y", D = "D", X = "X", Z = "Z1",
                  weights = "wgt", FE = c("group","year"),
                  data = s4, cl = "group", vartype = "cluster")
}
\keyword{graphics}


