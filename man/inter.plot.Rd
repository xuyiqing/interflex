\name{inter.plot}
\alias{inter.plot}
\title{Plotting Marginal Effect Estimates}
\description{Plotting marginal effect estimates after either the binning estimator or the kernel estimator is applied}
\usage{inter.plot(out, CI = TRUE, Xdistr = "histogram", 
  main = NULL, Ylabel = NULL, Dlabel = NULL, Xlabel = NULL, 
  xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, 
  theme.bw = FALSE, show.grid = TRUE, 
  cex.main = NULL, cex.lab = NULL, cex.axis = NULL, 
  bin.labs = TRUE, interval = NULL, file = NULL)
}
\arguments{
  \item{out}{an \bold{interflex} object.}
  \item{CI}{a logical flag indicating whether the confidence intervals need to be shown.}
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
  \item{cex.lab}{a numeric value that controls the font size of the axis labels.}   
  \item{cex.axis}{a numeric value that controls the font size of the axis numbers.}
  \item{bin.labs}{a logical flag indicating whether to show a label for each bin in a binning plot; only meaningful when \bold{nbins} is 3 or 4. The default is TRUE.}
  \item{interval}{draw two vertical lines to demonstrate the interval used in replicated papers.}
  \item{file}{a string that specifies the filename in which the plot is saved.}
}
\details{
  \bold{inter.plot} draws a margining effect plot based on a \bold{interflex} object after either \bold{inter.binning} or \bold{inter.kernel} is applied. It allows users to flexibly change the look of a plot without re-estimating the model, hence saving time. 
}
\value{  
  \item{graph}{stores the graphic output, a \bold{ggplot2} object.}
}
\author{
  Jens Hainmueller; Jonathan Mummolo; Yiqing Xu (Maintainer)
}
\references{
  Jens Hainmueller; Jonathan Mummolo; Yiqing Xu. 2019. "How Much Should We Trust Estimates from Multiplicative Interaction Models? Simple Tools to Improve Empirical Practice." Political Analysis, Vol. 27, Iss. 2, April 2019, pp. 163--192. Available at SSRN: \url{https://papers.ssrn.com/abstract_id=2739221}.
}


\seealso{\code{\link{inter.binning}} and \code{\link{inter.kernel}}}
\examples{
library(interflex)
data(interflex)
out <- inter.binning(data = s1, Y = "Y", D = "D", X = "X", Z = "Z1")
inter.plot(out)
}
\keyword{graphics}
