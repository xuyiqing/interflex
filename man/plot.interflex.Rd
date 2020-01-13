\name{plot.interflex}
\alias{plot.interflex}
\title{Plotting Marginal Effect Estimates}
\description{Plotting marginal effect estimates after either the binning estimator or the kernel estimator is applied}
\usage{\method{plot}{interflex}(x, order = NULL, subtitles = NULL, show.subtitles = NULL,
  CI = NULL, diff.values = NULL, Xdistr = "histogram", main = NULL, 
  Ylabel = NULL, Dlabel = NULL, Xlabel = NULL, xlab = NULL, ylab = NULL, 
  xlim = NULL, ylim = NULL, theme.bw = FALSE, show.grid = TRUE, 
  cex.main = NULL, cex.sub = NULL, cex.lab = NULL, cex.axis = NULL, 
  bin.labs = TRUE, interval = NULL, file = NULL, ncols = NULL, pool = FALSE, 
  legend.title = NULL, color = NULL, jitter = FALSE, \dots)
}
\arguments{
  \item{x}{an \bold{interflex} object.}
  \item{order}{a vector of strings that determines the order of treatment arms in the plot when visualizing marginal effects. It should contain all kinds of treatment arms except for the base group.}
  \item{subtitles}{a vector of strings that determines the subtitles of subplots. Its length should be equal to the number of categories of treatment arms minus one. }
  \item{show.subtitles}{a logical flag controlling whether to show subtitles. Default to TRUE when the number of treatment arms is larger than 2.}
  \item{CI}{a logical flag indicating whether the confidence intervals need to be shown.}
  \item{diff.values}{a two-element numeric vector. The treatment effects at corresponding values of the moderator will be marked on the plot.}
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
  \item{pool}{a logical flag specifying whether to integrate marginal effects for each treatment arm in one plot.}
  \item{legend.title}{a string that specifies the title of the legend when \bold{pool} is TRUE.}
  \item{color}{a vector of colors that determines the color of lines when \bold{pool} is TRUE.}
  \item{jitter}{a logical flag specifying whether to jitter error-bars to avoid overlapping. The default is FALSE.}
  \item{\dots}{other argv.}
}
\details{
  \bold{plot.interflex} draws a margining effect plot based on a \bold{interflex} object after either \bold{inter.binning} or \bold{inter.kernel} is applied. It allows users to flexibly change the look of a plot without re-estimating the model, hence saving time. 
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
  \code{\link{interflex}} and \code{\link{predict.interflex}}
}
\examples{
library(interflex)
data(interflex)
out <- inter.binning(data = s1, Y = "Y", D = "D", X = "X", Z = "Z1", CI = FALSE)
plot(out)
}
\keyword{graphics}
