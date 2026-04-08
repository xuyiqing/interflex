\encoding{UTF-8}
\name{app-data}
\alias{app-data}
\alias{app_adiguzel2023}
\alias{app_bb2024}
\alias{app_et2023}
\alias{app_hma2015}
\alias{app_vernby2013}
\docType{data}
\title{
  Replication Datasets from Published Interaction-Effect Studies
}
\description{
  A collection of data frames from published political science studies
  that feature conditional marginal effects. These datasets are shipped
  with \pkg{interflex} so users can reproduce the applied examples in
  the package's user manual (\url{https://yiqingxu.org/packages/interflex/})
  and experiment with the provided estimators on real data.
}
\usage{
data(interflex)
}
\format{
  Each object is a \code{data.frame}. Column sets differ across studies;
  the treatment variable, moderator, and outcome used in the manual
  examples are listed below.

  \describe{
    \item{\code{app_adiguzel2023}}{
      1,295 rows, 10 columns. Treatment: \code{Dakp_3}. Moderator:
      \code{lograyic09}. Used in the user manual to illustrate linear
      and kernel estimators on a discrete treatment with continuous
      moderator.
    }
    \item{\code{app_bb2024}}{
      810 rows, 22 columns. Treatment: \code{treat_rev}. Outcome:
      \code{sup_tax}. Moderator: \code{clim_belief}. Multi-arm
      treatment example.
    }
    \item{\code{app_et2023}}{
      570 rows, 39 columns. Treatment: \code{board_service}. Outcome:
      \code{seat_next_5yr}. Illustrates binning and fixed-effect
      estimators on legislative career data.
    }
    \item{\code{app_hma2015}}{
      1,482 rows, 11 columns. Treatment: \code{threat} (binary).
      Outcome: \code{totangry}. Moderator: \code{pidentity}. Based on
      Huddy, Mason, and Aarøe (2015); used throughout the manual as the
      primary binary-treatment running example.
    }
    \item{\code{app_vernby2013}}{
      183 rows, 9 columns. Treatment: \code{noncitvotsh}. Outcome:
      \code{school_diff}. Moderator: \code{Taxbase2}. Small-sample
      regional dataset used in the extensions chapter.
    }
  }
}
\source{
  Each dataset is redistributed with permission from its original
  study's replication archive. See the user manual for full citations
  and links to the source archives:
  \url{https://yiqingxu.org/packages/interflex/}.
}
\references{
  Adiguzel, F. S. (2023).
  Bernauer, J. and Bodnar, M. (2024).
  Eggers, A. C. and Tanaka, S. (2023).
  Huddy, L., Mason, L. and Aarøe, L. (2015). Expressive partisanship:
  Campaign involvement, political emotion, and partisan identity.
  \emph{American Political Science Review} 109(1), 1-17.
  Vernby, K. (2013). Inclusion and public policy: Evidence from
  Sweden's introduction of noncitizen suffrage.
  \emph{American Journal of Political Science} 57(1), 15-29.
}
\keyword{datasets}
