\name{inter.kernel}
\alias{inter.kernel}
\title{The Kernel Estimator}
\description{Implementing the kernel estimator, based on a semi-parametric model that relaxes the linear interaction effect assumption of conventional multiplicative interaction models and safeguards against excessive extrapolation
}
\usage{inter.kernel(data, Y, D, X, treat.type = NULL, base = NULL, 
  Z = NULL, weights = NULL, FE = NULL, full.moderate = TRUE, 
  na.rm = FALSE, Xunif = FALSE, CI = TRUE, conf.level = 0.95, cl = NULL,
  CV.method = NULL, kfold = 10, grid = 30, neval = 50, nboots = 200, parallel = TRUE, 
  cores = 4, seed = 02139, bw = NULL, bw.adaptive = TRUE, quantile.eval = FALSE, 
  metric = "MSPE", predict = FALSE, D.ref = NULL, figure = TRUE, order = NULL, 
  subtitles = NULL, show.subtitles = NULL, Xdistr = "histogram", main = NULL, 
  Ylabel = NULL, Dlabel = NULL, Xlabel = NULL, xlab = NULL, ylab = NULL, 
  xlim = NULL, ylim = NULL, theme.bw = FALSE, interval = NULL, show.grid = TRUE, 
  cex.main = NULL, cex.sub = NULL, cex.lab = NULL, cex.axis = NULL, 
  file = NULL, ncols = NULL, pool = FALSE, color = NULL, legend.title = NULL, 
  diff.values = NULL)
}
\arguments{
  \item{data}{a data.frame that stores Y, D, and X.}
  \item{Y}{a string, name of the outcome variable.}
  \item{D}{a string, name of the treatment variable.}
  \item{X}{a string, name of the moderating variable.}
  \item{treat.type}{a string that specifies the type of the treatment variable. Can be one of the following: "discrete" and "continuous". If not specified, it will be automatically set based on the treatment variable. (If D is in numeric format and the number of multiple arms of D is larger than 5, treat.type will be set to "continuous". Otherwise, it will be set to "discrete".)}
  \item{base}{a string specifying the base group if \bold{treat.type} is "discrete". If not specified, the base group will be selected based on numeric/character order of values in the treatment variable. For example, when treatment D can take values 0, 1, and 2, those observations with treatment 0 will be selected as the base group.}
  \item{Z}{a vector of variable names of control variables.}
  \item{weights}{a string, name of the weighting variable.}
  \item{FE}{a vector of variable names of fixed effects indicators.}
  \item{full.moderate}{a logical flag indicating whether to use fully moderated model. If it is set to TRUE, all interaction terms between moderator X and covariates Z(s) will be treated as extra covariates. The default is TRUE. }
  \item{na.rm}{a logical flag indicating whether to list-wise delete missing data. The algorithm will report an error if missing data exist.}
  \item{Xunif}{a logical flag controlling whether to transform values of the moderator into their percentiles.}  
  \item{CI}{a logical flag indicating whether the confidence intervals need to be shown. If \bold{CI} = FALSE, the bootstrap procedure will not be conducted.}
  \item{conf.level}{a number that specifies the confidence interval. For example, 0.9 means a confidence interval of [5\%, 95\%]. The default is 0.95,
    i.e., [2.5\%, 97.5\%].}
  \item{cl}{a string specifying the name of clustering variable. If specified, the bootstrap procedure will be blocking on this variable. Ignored when \bold{CI} = FALSE.}
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
  \item{quantile.eval}{a logical flag indicating whether to choose evaluation points based on quantiles rather than uniformly. If it is set to TRUE, denser evaluation points will be assigned to the region with high-density observations
	 and vice versa. The default is FALSE.}
  \item{metric}{either "MSPE" (Mean Squared Prediction Error) or "MAPE" (Mean Absolute Prediction Error), criterion to be used in cross-validation to select bandwidth. The default is "MSPE".}
  \item{predict}{a logical flag that specifies whether to estimate the expected value of Y conditional on X, D and the mean of all covariates. The default is FALSE. If it is set to TRUE, it will take a bit more time to estimate and save expected values of Y. Users can later use the command predict(...) to visualize the plots.}
  \item{D.ref}{a vector of numbers that specifies the reference values of D which the expected values of Y will be conditional on when \bold{treat.type} is "continuous". If it is not specified, the 0.25, 0.5, and 0.75 quantiles of D will be chosen.}
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
  \item{ncols}{a integral value that controls the number of columns in the plot if \bold{pool} is FALSE.}
  \item{pool}{a logical flag specifying whether to integrate marginal effects or expected values of Y for each treatment arm in one plot.}
  \item{color}{a vector of colors that determines the color of lines when drawing the pool plot.}
  \item{legend.title}{a string that specifies the title of the legend when \bold{pool} is TRUE.}
  \item{diff.values}{a three-element numeric vector. The difference between treatment effects at corresponding values of the moderator will be estimated. If not specified, the 0.25,0.5,0.75 quantiles of the moderator will be chosen.}
}
\details{
  \bold{inter.kernel} implements a kernel smoothing estimator of the marginal effect, which is an application of semi-parametric smooth varying-coefficient models. This approach is more complicated than the simple binning estimator, but provides a generalization that allows researchers to flexibly estimate the functional form of the marginal effect of D on Y across the values of X by estimating a series of local effects with a Gaussian kernel
  re-weighting scheme.

  The kernel regression nests the standard linear interaction model as a special case. However, in the kernel regression the conditional effect of D on Y does not have to fall on a linear line as required by the linear interaction effect (LIE) assumption, but can vary freely across the range of X. In addition, if covariates Z are included in the model, the coefficients of those covariates are also allowed to vary freely across the range of X resulting in a very flexible estimator that also helps to guard against misspecification bias with respect to the
  covariates.

  \bold{inter.kernel} produces a smooth marginal effect of D on Y with respect to X. It estimates a series of locally linear regression using kernel re-weighting based on the distance between each value of X (the moderator) and each evaluation point.  If the bandwidth is not supplied, it is selected using a standard n-fold cross-validation procedure. Standard errors and confidence intervals can be computed using a bootstrap. It also adds a histogram (or density plot)  at the bottom of the figure to help users assess the common support assumption based on the distribution of the moderator.

  \bold{inter.kernel} can accommodate linear fixed-effect models with multiple fixed-effect indicators. A fixed-effect model is estimated using a fast algorithm written in C++ at each evaluation point of the semi-parametric kernel estimation. 
  
  \bold{inter.kernel} can estimate the expected value of Y conditional on the moderator X and treatment D when the option "predict" is set to TRUE.
}
\value{
  \item{bw}{specified or cross-validated bandwidth.}
  \item{est}{marginal effects estimates (and confidence intervals if \bold{CI = TRUE})}
  \item{est.matrix}{stores the variance-covariance matrix of the series of point estimates if \bold{CI = TRUE}.}
  \item{CV.out}{results from cross-validation if \bold{bw} is not specified.}
  \item{graph}{stores the graphic output, a ggplot object.}
  \item{est.predict}{stores the expected value of Y conditional on X, D and the mean of all covariates if \bold{predict = TRUE}.}
  \item{difference.est}{stores the difference between treatment effects at corresponding values specified in \bold{diff.values}.}

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
inter.kernel(Y = "Y", D = "D", X = "X", Z = "Z1", 
			 data = s3, CI = FALSE, bw = 1.3, parallel = FALSE)
}
\keyword{graphics}


