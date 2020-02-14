\name{inter.linear}
\alias{inter.linear}
\title{The Linear Estimator}
\description{Implementing the linear interaction estimator}
\usage{inter.linear(data, Y, D, X, treat.type = NULL, base = NULL, 
  Z = NULL, FE = NULL, weights = NULL, full.moderate = TRUE, 
  na.rm = FALSE, Xunif = FALSE, CI = TRUE, vartype = "robust", 
  nboots= 200, parallel = TRUE, cores =4, cl = NULL, time = NULL, 
  pairwise = TRUE, predict = FALSE, D.ref = NULL, figure = TRUE, 
  order = NULL, subtitles = NULL, show.subtitles = NULL, Xdistr = "histogram", 
  main = NULL, Ylabel = NULL, Dlabel = NULL, Xlabel = NULL, 
  xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, theme.bw = FALSE, 
  show.grid = TRUE, cex.main = NULL, cex.sub = NULL, cex.lab = NULL, 
  cex.axis = NULL, interval = NULL, file = NULL, ncols = NULL, 
  pool = FALSE, color = NULL, legend.title = NULL, diff.values = NULL,
  percentile = FALSE)
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
  \item{Xunif}{a logical flag controlling whether to transform values of the moderator into their rank order. The default is FALSE.}
  \item{CI}{a logical flag indicating whether the confidence intervals need to be shown.}
  \item{vartype}{a string that controls the variance-covariance estimator. Can be one of the following: "homoscedastic"; "robust", "cluster", "pcse" (panel corrected standard errors) and "bootstrap". The default is "robust".}
  \item{nboots}{an integer specifying the number of bootstrap runs when \bold{vartype} is "bootstrap". The default is 200.}
  \item{parallel}{a logical flag indicating whether parallel computing will be used in bootstrapping.}
  \item{cores}{an integer indicating the number of cores to be used in parallel computing. If not specified, the algorithm will use the maximum number of logical cores of your computer (warning: this could prevent you from multi-tasking on your computer).}
  \item{cl}{a string specifying the name of clustering variable for clustered standard errors. \bold{vartype} must be either "cluster" or "pcse". \bold{vartype} will be automatically set to "cluster" if \bold{cl} is supplied and \bold{vartype} is not "pcse". "pcse" is not allowed when \bold{FE} is specified.}
  \item{time}{a string specifying the name of time variable for panel corrected standard errors. Ignored when \bold{vartype} is not  "pcse".}
  \item{pairwise}{a logical flag indicating whether to switch on the \bold{pairwise} option when estimating panel corrected standard errors. Ignored when \bold{vartype} is not "pcse". The default is TRUE.}
  \item{predict}{a logical flag that specifies whether to estimate the expected value of Y conditional on X, D and the mean of all covariates. The default is FALSE. If it is set to TRUE, it will take a bit more time to estimate and save expected values of Y. Users can later use the command plot(predict=TRUE,...) to visualize the plots.}
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
  \item{interval}{draw two vertical lines to demonstrate the interval used in replicated papers.}
  \item{file}{a string that specifies the filename in which the plot is saved.}
  \item{ncols}{a integral value that controls the number of columns in the plot if \bold{pool} is FALSE.}
  \item{pool}{a logical flag specifying whether to integrate marginal effects or expected values of Y for each treatment arm in one plot.}
  \item{color}{a vector of colors that determines the color of lines when drawing the pool plot.}
  \item{legend.title}{a string that specifies the title of the legend when \bold{pool} is TRUE.}
  \item{diff.values}{a numeric vector contain 2 or 3 elements which are between the minimum and the maximum of the moderator. The treatment effects at corresponding values of the moderator will be compared.}
  \item{percentile}{a logical flag indicating whether to take values of the moderator on a percentile scale.}
}
\details{
  \bold{inter.linear} implements the classic linear estimator. 
}
\value{
  \item{est.lin}{stores the linear marginal effect estimates.}
  \item{vcov.matrix}{stores the variance-covariance matrix of the series of point estimates.}
  \item{tests}{stores the results of various statistical tests. See below}
  \item{tests$X.LKurtosis}{the L-kurtosis measure of the moderator}
  \item{graph}{stores the graphic output, if \bold{pool} is TRUE, it is a ggplot object, otherwise it is a Grid object.}
  \item{est.predict}{stores the expected value of Y conditional on X, D and the mean of all covariates if \bold{predict = TRUE}.}
  \item{ttest.diffs}{stores the difference between treatment effects at corresponding values specified in \bold{diff.values}.}
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
inter.linear(Y = "Y", D = "D", X = "X", Z = "Z1", 
  data = s1, vartype = "homoscedastic", 
  Ylabel = "Y", Dlabel = "Tr", Xlabel="X")
 
inter.linear(Y = "Y", D = "D", X = "X", Z = c("Z1","Z2"),
  data = s5, base = 'A')

}
\keyword{graphics}


