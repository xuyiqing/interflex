\name{predict.interflex}
\alias{predict.interflex}
\title{Plotting Marginal Effect Estimates}
\description{Plotting expected outcomes given fixed values of the treatment and moderator after either the linear, binning or the kernel estimator is applied.}
\usage{\method{predict}{interflex}(object, to.show = NULL, order = NULL, subtitles = NULL, 
  show.subtitles = NULL, Xdistr = "histogram", CI = NULL, pool = FALSE, main = NULL, 
  Ylabel = NULL, Xlabel = NULL, xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, 
  theme.bw = FALSE, show.grid = TRUE, cex.main = NULL, cex.sub = NULL, 
  cex.lab = NULL, cex.axis = NULL, color = NULL, file = NULL,interval = NULL, 
  legend.title = NULL, ncols = NULL, \dots)
}
\arguments{
  \item{object}{an \bold{interflex} object (\bold{predict==TRUE}).}
  \item{to.show}{a vector of strings that determines which treatment arms will be displayed when visualizing expected values. It should should be a subset of all treatment arms including the base group.}
  \item{order}{a vector of strings that determines the order of treatment arms in the plot when visualizing expected values. It should contain all treatment arms or all elements in "to.show".}
  \item{subtitles}{a vector of strings that determines the subtitles of subplots when \bold{pool} is FALSE, or determines the label in the legend when \bold{pool} is TRUE. }
  \item{show.subtitles}{a logical flag controlling whether to show subtitles. Default to TRUE when the number of treatment arms is larger than 2.}
  \item{Xdistr}{a string indicating the way the distribution of the moderator will be plotted: "histogram" (or "hist"), "density", or "none". The default is "histogram".}
  \item{CI}{a logical flag indicating whether the confidence intervals need to be shown.}
  \item{pool}{a logical flag specifying whether to integrate expected values of Y for each treatment arm in one plot.}
  \item{main}{a string that controls the title of the plot.}
  \item{Ylabel}{a string that controls the label of the outcome variable Y in the plot.}
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
  \item{color}{a vector of colors that determines the color of lines when \bold{pool} is TRUE.}
  \item{file}{a string that specifies the filename in which the plot is saved.}
  \item{interval}{draw two vertical lines to demonstrate the interval used in replicated papers.}
  \item{legend.title}{a string that specifies the title of the legend when \bold{pool} is TRUE.}
  \item{ncols}{na integral that controls the number of columns in the plot if \bold{pool} is FALSE.}
  \item{\dots}{other argv.}
}
\details{
  \bold{predict.interflex} visualize expected outcomes given fixed values of the treatment and moderator after either the linear, binning or the kernel estimator is applied.  It allows users to flexibly change the look of a plot without re-estimating the model, hence saving time. 
}
\value{  
  \item{graph}{stores the graphic output, a \bold{ggplot2} object.}
}
\author{
  Jens Hainmueller; Jonathan Mummolo; Yiqing Xu (Maintainer); Ziyi Liu
}
\references{
  Jens Hainmueller; Jonathan Mummolo; Yiqing Xu. 2019. "How Much Should We Trust Estimates from Multiplicative Interaction Models? Simple Tools to Improve Empirical Practice." Political Analysis, Vol. 27, Iss. 2, April 2019, pp. 163--192. Available at: \url{http://bit.ly/HMX2019}.
}
\seealso{
  \code{\link{interflex}} and \code{\link{plot.interflex}}
}
\examples{
library(interflex)
data(interflex)
out <- inter.binning(data = s1, Y = "Y", D = "D", X = "X", Z = "Z1",predict = TRUE)
predict(out)
}
\keyword{graphics}
