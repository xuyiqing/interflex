\name{interflex}
\alias{interflex}
\title{Flexible Interactive Models}
\description{Producing flexible marginal effect estimates with linear, binning, and kernel estimators}
\usage{interflex(data,Y, D, X, estimator, treat.type = NULL, 
      base = NULL, Z = NULL, weights = NULL, 
      FE = NULL,cl = NULL, full.moderate = TRUE, 
      na.rm = FALSE, Xunif = FALSE, 
      CI = TRUE, conf.level = 0.95, nbins = 3, cutoffs = NULL, 
      vartype = "robust", time = NULL, 
      pairwise = TRUE, wald = TRUE, bin.labs = TRUE, 
      CV.method = NULL, kfold = 10, grid = 30, 
      neval = 50, nboots = 200, parallel = TRUE,
      cores = 4, seed = 02139, bw = NULL, bw.adaptive = TRUE, 
      quantile.eval = FALSE, metric = "MSPE", figure = TRUE, 
      order = NULL, subtitles = NULL,
      show.subtitles = NULL, Xdistr = "histogram", 
      main = NULL, Ylabel = NULL, Dlabel = NULL, Xlabel = NULL, 
      xlab = NULL,ylab = NULL, xlim = NULL, ylim = NULL, 
      theme.bw = FALSE, interval = NULL, show.grid = TRUE,
      cex.main = NULL, cex.sub = NULL, cex.lab = NULL, cex.axis = NULL, 
      file = NULL, ncols = NULL, pool = FALSE, color = NULL, 
      legend.title = NULL, jitter = FALSE,
      predict = FALSE, D.ref = NULL, diff.values = NULL,
	  percentile = FALSE)
}
\arguments{
  \item{data}{a data.frame that stores Y, D, and X.}
  \item{Y}{a string, name of the outcome variable.}
  \item{D}{a string, name of the treatment variable.}
  \item{X}{a string, name of the moderating variable.}
  \item{estimator}{a string that specifies the estimator. Can be one of the following: "linear", "binning" and "kernel".}
  \item{treat.type}{a string that specifies the type of the treatment variable. Can be one of the following: "discrete" and "continuous". If not specified, it will be automatically set based on the treatment variable. (If D is in numeric format and the number of multiple arms of D is larger than 5, treat.type will be set to "continuous". Otherwise, it will be set to "discrete".)}
  \item{base}{a string specifying the base group if \bold{treat.type} is "discrete". If not specified, the base group will be selected based on numeric/character order of values in the treatment variable. For example, when treatment D can take values 0, 1, and 2, those observations with treatment 0 will be selected as the base group.}
  \item{Z}{a vector of variable names of control variables.}
  \item{weights}{a string, name of the weighting variable.}
  \item{FE}{a vector of variable names of fixed effects indicators.}
  \item{cl}{a string specifying the name of clustering variable. If specified, the bootstrap and cross-validation procedure will be blocking on this variable. Ignored when \bold{CI} = FALSE. 
    When \bold{vartype} is "pcse", it specifies the name of unit variable for panel corrected standard errors.}
  \item{full.moderate}{a logical flag indicating whether to use fully moderated model. If it is set to TRUE, all interaction terms between moderator X and covariates Z(s) will be treated as extra covariates. The default is TRUE. }
  \item{na.rm}{a logical flag indicating whether to list-wise delete missing data. The algorithm will report an error if missing data exist.}
  \item{Xunif}{a logical flag controlling whether to transform values of the moderator into their rank order.}  
  \item{CI}{a logical flag indicating whether the confidence intervals need to be shown. If \bold{CI} = FALSE, the bootstrap procedure will not be conducted.}
  \item{conf.level}{a number that specifies the confidence interval. For example, 0.9 means a confidence interval of [5\%, 95\%]. The default is 0.95,
    i.e., [2.5\%, 97.5\%].}
  \item{nbins}{an integer that determines the number of bins of the moderator. The sample will be divided based on integral multiples of 100/nbins percentile.  \bold{nbins} will be automatically
    subtracted by 1 if more than one multiples of the percentile have the same value (for example, if the moderator has over 70\% zeros, both the 33 and 66 percentiles are zero). The default is 3. Ignored when \bold{cutoffs} is supplied. When nbins are set to 1, it only conducts a linear fit.}
  \item{cutoffs}{a vector of numbers that determines how subgroups are divided based on the moderator X. When this option is supplied, \bold{nbins} will be ignored. The smallest number of the first interval and the largest number of the last interval do not need to be specified. Numbers smaller than the minimum or larger than the maximum of X will be ignored. \bold{nbins} equals to the length of this vector plus 1. Ignored if the treatment is dichotomous.}
  \item{vartype}{a string that controls the variance-covariance estimator. Can be one of the following: "homoscedastic"; "robust", "cluster", "pcse" (panel corrected standard errors) and "bootstrap". The default is "robust".}
  \item{time}{a string specifying the name of time variable for panel corrected standard errors. Ignored when \bold{vartype} is not "pcse".}
  \item{pairwise}{a logical flag indicating whether to switch on the \bold{pairwise} option when estimating panel corrected standard errors. Ignored when \bold{vartype} is not "pcse". The default is TRUE.}
  \item{wald}{a logical flag controlling whether to conduct a Wald test.}
  \item{bin.labs}{a logical flag indicating whether to show a label for each bin in a binning plot; only meaningful when \bold{nbins} is 3 or 4. The default is TRUE.}
  \item{CV.method}{a string specified the method of sampling in cross-validation. Can be one of the following: "simple", "cluster", and "stratify".The default is "simple".}
  \item{kfold}{an integer specified the number of folds in cross-validation, default to 10. If it is larger than the number of observations in the dataset, it will conduct a leave-one-out CV.}
  \item{grid}{either an integer indicating the number of candidates in a grid search to find the optimal bandwidth, or a vector of positive numbers of candidate bandwidths. Ignored if \bold{bw} is specified.}
  \item{neval}{an integer specifying the number of evaluation points. A large number would produce smoother marginal effects estimates but the algorithm will take more time to run. The default is 50.}
  \item{nboots}{an integer specifying the number of bootstrap runs. The default is 200.}
  \item{parallel}{a logical flag indicating whether parallel computing will be used in bootstrapping and/or cross-validation. Ignored if \bold{CI} = FALSE and \bold{bw} is specified.}
  \item{cores}{an integer indicating the number of cores to be used in parallel computing. If not specified, the algorithm will use the maximum number of logical cores of your computer (warning: this could prevent you from multi-tasking on your computer).}
  \item{seed}{an integer that sets the seed in random number generation. Ignored if  \bold{CI} = FALSE and \bold{bw} is specified.}
  \item{bw}{a number that sets the bandwidth in Gaussian kernel reweighting.}
  \item{bw.adaptive}{a logical flag indicating whether to use adaptive bandwidth. If it is set to TRUE, a small bandwidth will be assigned to the observations in the high-density
	region and vice versa, \bold{bw} will indicate the bandwidth applied to the region of X with highest density of observations.}
  \item{quantile.eval}{a logical flag indicating whether to choose evaluation points based on quantiles rather than uniformly. If it is set to TRUE, denser evaluation points will be assigned to the region with high-density observations and vice versa. The default is FALSE.}
  \item{metric}{either "MSPE" (Mean Squared Prediction Error) or "MAPE" (Mean Absolute Prediction Error), criterion to be used in cross-validation to select bandwidth. The default is "MSPE".}
  \item{figure}{a logical flag controlling whether to draw a marginal effect plot.}
  \item{order}{a vector of strings that determines the order of treatment arms in the plot when visualizing marginal effects. It should contain all kinds of treatment arms except for the base group.}
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
  \item{interval}{draw two vertical lines to demonstrate the interval used in replicated papers.}
  \item{show.grid}{a logical flag indicating whether to show grid in the plot.}
  \item{cex.main}{a numeric value that controls the font size of the plot title.}   
  \item{cex.sub}{a numeric value that controls the font size of subtitles.}
  \item{cex.lab}{a numeric value that controls the font size of the axis labels.}   
  \item{cex.axis}{a numeric value that controls the font size of the axis numbers.}
  \item{file}{a string that specifies the filename in which the plot is saved.}
  \item{ncols}{an integral that controls the number of columns in the plot if \bold{pool} is FALSE.}
  \item{pool}{a logical flag specifying whether to draw marginal effects or expected values of Y for each treatment arm in a single plot. The default is FALSE.}
  \item{color}{a vector of colors that determines the color of lines when drawing the pool plot.}
  \item{legend.title}{a string that specifies the title of the legend when \bold{pool} is TRUE.}
  \item{jitter}{a logical flag specifying whether to jitter error-bars to avoid overlapping. The default is FALSE.}
  \item{predict}{a logical flag that specifies whether to estimate the expected value of Y conditional on X, D and the mean of all covariates. The default is FALSE. If it is set to TRUE, it will take a bit more time to estimate and save expected values of Y. Users can later use the command plot(predict=TRUE,...) to visualize the plots.}
  \item{D.ref}{a vector of numbers that specifies the reference values of D which the expected values of Y will be conditional on when \bold{treat.type} is "continuous". If it is not specified, the 0.25, 0.5, and 0.75 quantiles of D will be chosen.}
  \item{diff.values}{a numeric vector contain 2 or 3 elements which within the range of the moderator. The treatment effects at corresponding values of the moderator will be compared.}
  \item{percentile}{a logical flag indicating whether to take values of the moderator on a percentile scale.}
}
\details{ 
  \bold{interflex} is an umbrella function that allow users to call inter.kernel and inter.binning more conveniently.
}

\author{
  Jens Hainmueller; Jonathan Mummolo; Yiqing Xu (Maintainer); Ziyi Liu
}
\references{
  Jens Hainmueller; Jonathan Mummolo; Yiqing Xu. 2019. "How Much Should We Trust Estimates from Multiplicative Interaction Models? Simple Tools to Improve Empirical Practice." Political Analysis, Vol. 27, Iss. 2, April 2019, pp. 163--192. Available at: \url{http://bit.ly/HMX2019}.
}
\seealso{
  \code{\link{inter.raw}}, \code{\link{inter.gam}}, \code{\link{plot.interflex}}, and \code{\link{predict.interflex}}
}
\examples{
library(interflex)
data(interflex)
s5.binning <-interflex(estimator = 'binning', Y = "Y", D = "D", X = "X", 
Z = c("Z1","Z2"), data = s5, CI = FALSE)
}



