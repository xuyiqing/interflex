\name{inter.test}
\alias{inter.test}
\title{Conduct t-Test after Estimation}
\description{Compare treatment effects at specific values of the moderator}
\usage{inter.test(out, diff.values, percentile=FALSE, k=16)
}
\arguments{
  \item{out}{an \bold{interflex} object.}
  \item{diff.values}{a numeric vector contain 2 or 3 elements which are within the range of the moderator. The treatment effects at corresponding values of the moderator will be compared.}
  \item{percentile}{a logical flag indicating whether to take values of the moderator on a percentile scale.}
  \item{k}{an integer specifying the dimension(s) of the bases used to represent the smooth term, default to 16.}
}
\details{
  \bold{inter.test} compare treatment effects at specific values of the moderator using marginal effects and vcov matrix derived from linear/kernel estimation. Based on GAM model(relies on \bold{mgcv} package), users can approximate the treatment effect and its variance using smooth functions without re-estimating the model, hence saving time.
}

\author{
  Jens Hainmueller; Jonathan Mummolo; Yiqing Xu; Ziyi Liu (Maintainer)
}
\references{
  Jens Hainmueller; Jonathan Mummolo; Yiqing Xu. 2019. "How Much Should We Trust Estimates from Multiplicative Interaction Models? Simple Tools to Improve Empirical Practice." Political Analysis, Vol. 27, Iss. 2, April 2019, pp. 163--192. Available at: \url{https://www.cambridge.org/core/journals/political-analysis/article/how-much-should-we-trust-estimates-from-multiplicative-interaction-models-simple-tools-to-improve-empirical-practice/D8CAACB473F9B1EE256F43B38E458706}.
}
\seealso{
  \code{\link{interflex}}, \code{\link{plot.interflex}} and \code{\link{predict.interflex}}
}


