\name{inter.binning}
\alias{inter.binning}
\title{The Binning Estimator}
\description{Implementing the binning estimator, a generalization of the
  multiplicative interaction model and conducting various diagnostic tests}
\usage{inter.binning(data, Y, D, X, treat.type = NULL, base = NULL, 
  Z = NULL, FE = NULL, weights = NULL, full.moderate = TRUE, 
  na.rm = FALSE, Xunif = FALSE, nbins = 3, cutoffs = NULL, CI = TRUE, 
  vartype = "robust", nboots= 200, parallel = TRUE, cores =4, 
  cl = NULL, time = NULL, pairwise = TRUE, wald = TRUE, predict = FALSE, 
  D.ref = NULL, figure = TRUE, order = NULL, subtitles = NULL, 
  show.subtitles = NULL, Xdistr = "histogram", main = NULL, 
  Ylabel = NULL, Dlabel = NULL, Xlabel = NULL, xlab = NULL, 
  ylab = NULL, xlim = NULL, ylim = NULL, theme.bw = FALSE, 
  show.grid = TRUE, cex.main = NULL, cex.sub = NULL, cex.lab = NULL, cex.axis = NULL, 
  bin.labs = TRUE, interval = NULL, file = NULL, ncols = NULL, pool = FALSE, 
  color = NULL, jitter = FALSE, legend.title = NULL)
}
\arguments{
  \item{data}{a data.frame that stores Y, D, X, and other variables.}
  \item{Y}{a string, name of the outcome variable.}
  \item{D}{a string, name of the treatment variable.}
  \item{X}{a string, name of the moderating variable.}
  \item{treat.type}{a string that specifies the type of the treatment variable. Can be one of the following: "discrete" and "continuous". If not specified, it will be automatically set based on the treatment variable. (If D is in numeric format and the number of multiple arms of D is larger than 5, treat.type will be set to "continuous". Otherwise, it will be set to "discrete".)}
  \item{base}{a string specifying the base group if \bold{treat.type} is "discrete". If not specified, the base group will be selected based on numeric/character order of values in the treatment variable. For example, when treatment D can take values 0, 1, and 2, those observations with treatment 0 will be selected as the base group.}
  \item{Z}{a vector of variable names of control variables.}
  \item{FE}{a vector of variable names of fixed effects indicators.}
  \item{weights}{a string, name of the weighting variable.}
  \item{full.moderate}{a logical flag indicating whether to use fully moderated model. If it is set to TRUE, all interaction terms between moderator X and covariates Z(s) will be treated as extra covariates. The default is TRUE. }
  \item{na.rm}{a logical flag indicating whether to list-wise delete data. The algorithm will report error if missing data exist.}
  \item{Xunif}{a logical flag controlling whether to transform values of the moderator into their percentiles. The default is FALSE.}
  \item{nbins}{an integer that determines the number of bins of the moderator. The sample will be divided based on integral multiples of 100/nbins percentile.  \bold{nbins} will be automatically
    subtracted by 1 if more than one multiples of the percentile have the same value (for example, if the moderator has over 70\% zeros, both the 33 and 66 percentiles are zero). The default is 3. Ignored when \bold{cutoffs} is supplied. }
  \item{cutoffs}{a vector of numbers that determines how subgroups are divided based on the moderator X. When this option is supplied, \bold{nbins} will be ignored. The smallest number of the first interval and the largest number of the last interval do not need to be specified. Numbers smaller than the minimum or larger than the maximum of X will be ignored. \bold{nbins} equals to the length of this vector plus 1. Ignored if the treatment is dichotomous.}
  \item{CI}{a logical flag indicating whether the confidence intervals need to be shown.}
  \item{vartype}{a string that controls the variance-covariance estimator. Can be one of the following: "homoscedastic"; "robust", "cluster", "pcse" (panel corrected standard errors) and "bootstrap".
    The default is "robust".}
  \item{nboots}{an integer specifying the number of bootstrap runs when \bold{vartype} is "bootstrap". The default is 200.}
  \item{parallel}{a logical flag indicating whether parallel computing will be used in bootstrapping.}
  \item{cores}{an integer indicating the number of cores to be used in parallel computing. If not specified, the algorithm will use the maximum number of logical cores of your computer (warning: this could prevent you from multi-tasking on your computer).}
  \item{cl}{a string specifying the name of clustering variable for clustered standard errors. \bold{vartype} must be either "cluster" or "pcse". \bold{vartype} will be automatically set to "cluster" if \bold{cl} is supplied and \bold{vartype} is not "pcse". "pcse" is not allowed when \bold{FE} is specified.}
  \item{time}{a string specifying the name of time variable for panel corrected standard errors. Ignored when \bold{vartype} is not  "pcse".}
  \item{pairwise}{a logical flag indicating whether to switch on the \bold{pairwise} option when estimating panel corrected standard errors. Ignored when \bold{vartype} is not "pcse". The default is TRUE.}
  \item{wald}{a logical flag controlling whether to conduct a Wald test.}
  \item{predict}{a logical flag that specifies whether to estimate the expected value of Y conditional on X, D and the mean of all covariates. The default is FALSE. If it is set to TRUE, it will take a bit more time to estimate and save expected values of Y. Users can later use the command predict(...) to visualize the plots.}
  \item{D.ref}{a vector of numbers that specifies the reference values of D which the expected values of Y will be conditional on when \bold{treat.type} is "continuous". If it is not specified, the 0.25, 0.5, and 0.75 quantiles of D will be chosen.}
  \item{figure}{a logical flag controlling whether to draw a marginal effect plot.}
  \item{order}{a vector of strings that determines the order of subplots when visualizing marginal effects. It should contain all kinds of treatment arms except for the base group.}
  \item{subtitles}{a vector of strings that determines the subtitles of subplots when \bold{pool} is FALSE, or determines the label in the legend when \bold{pool} is TRUE. It should have the same length as the number of categories of treatment arms minus one.}
  \item{show.subtitles}{a logical flag controlling whether to show subtitles.}
  \item{Xdistr}{a string indicating the way the distribution of the moderator will be plotted: "histogram" (or "hist"), "density", or "none". The default is "histogram".}
  \item{main}{a string that controls the title of the plot.}
  \item{Ylabel}{a string that controls the label of the outcome variable Y in the plot.}
  \item{Dlabel}{a string that controls the label of the treatment variable D in the plot.}
  \item{Xlabel}{a string that controls the label of the moderating variable X in the plot.} 
  \item{xlab}{a string that specifies the label of the x-axis.}
  \item{ylab}{a string that specifies the label of the y-axis.}
  \item{xlim}{a two-element numeric vector that controls the range of the x-axis to be shown in the plot.}
  \item{ylim}{a two-element numeric vector that controls the range of the y-axis to be shown in the plot (with small adjustments to improve aesthetics).}
  \item{theme.bw}{a logical flag specifying whether to use a black-white theme.}
  \item{show.grid}{a logical flag indicating whether to show grid in the plot.}
  \item{cex.main}{a numeric value that controls the font size of the plot title.}  
  \item{cex.sub}{a numeric value that controls the font size of subtitles.}  
  \item{cex.lab}{a numeric value that controls the font size of the axis labels.}   
  \item{cex.axis}{a numeric value that controls the font size of the axis numbers.}
  \item{bin.labs}{a logical flag indicating whether to show a label for each bin in a binning plot; only meaningful when \bold{nbins} is 3 or 4. The default is TRUE.}
  \item{interval}{draw two vertical lines to demonstrate the interval used in replicated papers.}
  \item{file}{a string that specifies the filename in which the plot is saved.}
  \item{ncols}{a integral value that controls the number of columns in the plot if \bold{pool} is FALSE.}
  \item{pool}{a logical flag specifying whether to integrate marginal effects or expected values of Y for each treatment arm in one plot.}
  \item{color}{a vector of colors that determines the color of lines when drawing the pool plot.}
  \item{jitter}{a logical flag specifying whether to jitter error-bars to avoid overlapping. The default is FALSE.}
  \item{legend.title}{a string that specifies the title of the legend when \bold{pool} is TRUE.}
}
\details{
  \bold{inter.binning} implements the binning estimator. There are three steps to implement the estimator. First, we discretize the moderator X into several bins and create a dummy variable for each bin. The default is 3 bins, respectively corresponding to the three terciles.  Second, we pick an evaluation point within each bin, where we want to estimate the conditional marginal effect of D on Y. Typically, we choose the median value of X in each bin. Third, we estimate a model that includes interactions between the bin dummies and the treatment indicator, the bin dummies and the moderator X minus the evaluation points, as well as the triple interactions. The last two terms are to capture the effect of D on Y within each bin.

  The binning estimator has several key advantages over the standard multiplicative interaction model. First, the binning estimator is much more flexible as it jointly fits the interaction components of the standard model to each bin separately. The model does not impose the linear
  interaction effect (LIE) assumption. Instead, the conditional marginal effects can vary freely
  across the three bins and therefore can take on any non-linear or non-monotonic pattern that might describe the heterogeneity in the effect of D on Y across low, medium, or high levels of
  X.

  Second, since the bins are constructed based on the support of X, the binning ensures that the conditional marginal effects are estimated at typical values of the moderator and do not rely on
  excessive extrapolation or interpolation.

  Third, the binning estimator is easy to implement using any regression software and the standard errors for the conditional marginal effects are directly estimated by the regression so there are no need to compute linear combinations of coefficients to recover the conditional
  marginal effects.

  Fourth, the binning estimator actually provides a generalization that nests the standard multiplicative interaction model as a special case. It can therefore serve as a formal test on the validity of a global LIE assumption imposed by the standard model. So in the special case when the standard multiplicative interaction model is correct and therefore the global LIE assumption holds, then---as the sample size grows---the marginal effect estimates from the binning estimator converge in probability on the unbiased marginal effect estimates from the standard multiplicative interaction model.
  
  In the case of a discrete treatment, \bold{inter.binning} displays at the bottom of the figure a stacked histogram that shows the distribution of the moderator X. In this histogram the total height of the stacked bars refers to the distribution of the moderator in the pooled sample and the red and gray shaded bars refer to the distribution of the moderator in the treatment and control groups, respectively. Adding such a histogram makes it easy to judge the degree to which there is common support in the data. In the case of a continuous treatment, it displays a histogram at the bottom that simply shows the distribution of X in the entire sample.
}
\value{
  \item{est.lin}{stores the linear marginal effect estimates}
  \item{est.bin}{stores the binning estimates}
  \item{tests}{stores the results of various statistical tests. See below}
  \item{tests$bin.size}{the percentage of observations within each bin}
  \item{tests$X.LKurtosis}{the L-kurtosis measure of the moderator}
  \item{tests$correctOrder}{whether the three binning estimates are in the correct order (i.e. monotonically increasing or decreasing) when there are 3 bins.}
  \item{tests$p.twosided}{p values of pairwise t-tests of the binning estimates when there are 2 or 3 bins.}
  \item{tests$p.wald}{p-value of a Wald test. The NULL hypothesis is that the linear interaction model and the binning model are statistically equivalent.}
  \item{graph}{stores the graphic output, if \bold{pool} is TRUE, it is a ggplot object, otherwise it is a Grid object.}
  \item{est.predict}{stores the expected value of Y conditional on X, D and the mean of all covariates if \bold{predict = TRUE}.}
}
\author{
  Jens Hainmueller; Jonathan Mummolo; Yiqing Xu (Maintainer); Ziyi Liu
}
\references{
  Jens Hainmueller; Jonathan Mummolo; Yiqing Xu. 2019. "How Much Should We Trust Estimates from Multiplicative Interaction Models? Simple Tools to Improve Empirical Practice." Political Analysis, Vol. 27, Iss. 2, April 2019, pp. 163--192. Available at: \url{http://bit.ly/HMX2019}.
}
\seealso{
  \code{\link{interflex}}
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
  weights = "wgt", FE = c("group","year"), data = s4, 
  cl = "group", vartype = "cluster")
  
}
\keyword{graphics}


