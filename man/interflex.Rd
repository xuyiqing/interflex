\name{interflex}
\alias{interflex}
\title{Flexible Interactive Models}
\description{Producing flexible marginal effect estimates with linear, binning, and kernel estimators}
\usage{
interflex(estimator, data, Y, D, X, treat.type = NULL,
  base = NULL, Z = NULL, IV = NULL, FE = NULL, full.moderate = FALSE,
  weights = NULL, na.rm = FALSE, Xunif = FALSE, CI = TRUE, neval = 50,
  X.eval = NULL, method = "linear", vartype = "delta",
  vcov.type = "robust", time = NULL, pairwise = TRUE, nboots = 200,
  nsimu = 1000, parallel = FALSE, cores = 4, cl = NULL, Z.ref = NULL,
  D.ref = NULL, nbins = 3, cutoffs = NULL, wald = TRUE, bw = NULL,
  kfold = 10, grid = 30, metric = NULL, model.y = "rf",
  param.y = NULL, param.grid.y = NULL,
  scoring.y = "neg_mean_squared_error", model.t = "rf", param.t = NULL,
  param.grid.t = NULL, scoring.t = "neg_mean_squared_error", CV = FALSE,
  n.folds = 10, n.jobs = -1, cf.n.folds = 5, gate = FALSE,
  grf.num.trees = 4000, signal = c("aipw"), estimand = c("ATE"),
  basis.type = c("bspline"), include.interactions = TRUE, poly.degree = 2,
  spline.df = 5, spline.degree = 3, lambda.seq = NULL,
  reduce.dimension = c("bspline"), figure = TRUE, bin.labs = TRUE,
  order = NULL, subtitles = NULL, show.subtitles = NULL,
  Xdistr = "histogram", main = NULL, Ylabel = NULL, Dlabel = NULL,
  Xlabel = NULL, xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL,
  theme.bw = FALSE, show.grid = TRUE, show.uniform.CI = TRUE,
  cex.main = NULL, cex.sub = NULL, cex.lab = NULL, cex.axis = NULL,
  interval = NULL, file = NULL, ncols = NULL, pool = FALSE, color = NULL,
  show.all = FALSE, legend.title = NULL, diff.values = NULL,
  percentile = FALSE, gam.k = 10, angle = c(30, 100, -30, -120),
  span = NULL, scale = 1.1, height = 7, width = 10, box.pos = "down")
}
\arguments{
  \item{estimator}{a string that specifies the estimator. Can be one of the following: "raw", "linear", "binning", "kernel", "gam", "grf", "dml", or "lasso".}
  \item{data}{a data.frame that stores \bold{Y}, \bold{D}, \bold{X} (and optionally \bold{Z}, \bold{IV}, \bold{FE}, \bold{weights}).}
  \item{Y}{a string, name of the outcome variable.}
  \item{D}{a string, name of the treatment variable.}
  \item{X}{a string, name of the moderating variable.}
  \item{treat.type}{a string specifying the type of the treatment. Must be "discrete" or "continuous". If not specified, it is set automatically: if \bold{D} is numeric and has more than 5 unique values, then "continuous"; otherwise "discrete".}
  \item{base}{a string specifying the baseline treatment group when \bold{treat.type} is "discrete". If not specified, the baseline is chosen as the first level in sorted order of \bold{D}.}
  \item{Z}{a character vector of covariate names.}
  \item{IV}{a character vector of instrument variable names (only supported when \bold{method} = "linear").}
  \item{FE}{a character vector of fixed-effect indicator names (only supported when \bold{method} = "linear").}
  \item{full.moderate}{a logical flag indicating whether to use a fully moderated model, i.e., include additional interactions between \bold{X} and other covariates \bold{Z}.}
  \item{weights}{a string, name of the weighting variable.}
  \item{na.rm}{a logical flag indicating whether to list-wise delete missing data in required variables. If \bold{na.rm} = FALSE, the function errors when missing values exist.}
  \item{Xunif}{a logical flag controlling whether to transform \bold{X} to its percentile rank (scaled to 0--100).}
  \item{CI}{a logical flag indicating whether to compute confidence intervals (if \bold{CI} = FALSE, CI-related bootstrap/simulation routines are skipped when applicable).}
  \item{neval}{an integer specifying the number of evaluation points/grid size used by some estimators (larger values produce smoother curves but take longer).}
  \item{X.eval}{a numeric vector of manually selected evaluation points for \bold{X} (used by some estimators such as \bold{estimator}="linear" or "kernel").}
  \item{method}{a string specifying the outcome model: one of "linear", "probit", "logit", "poisson", or "nbinom".}
  \item{vartype}{a string controlling uncertainty estimation for \bold{estimator} = "linear" or "binning": one of "delta", "simu", or "bootstrap".}
  \item{vcov.type}{a string controlling the variance-covariance estimator for linear models: one of "homoscedastic", "robust", "cluster", or "pcse" (panel-corrected SEs; only for linear models).}
  \item{time}{a string specifying the time variable name for \bold{vcov.type} = "pcse".}
  \item{pairwise}{a logical flag for \bold{vcov.type} = "pcse" controlling the pairwise option.}
  \item{nboots}{an integer specifying the number of bootstrap runs (when bootstrap-based inference is used).}
  \item{nsimu}{an integer specifying the number of simulation draws (when simulation-based inference is used).}
  \item{parallel}{a logical flag indicating whether to use parallel computing in bootstrap and/or cross-validation routines (when supported by the chosen estimator).}
  \item{cores}{an integer indicating the number of cores to be used when \bold{parallel} = TRUE.}
  \item{cl}{a string specifying the clustering/unit variable name for \bold{vcov.type} = "cluster" or \bold{vcov.type} = "pcse". It is also used to block some resampling routines when applicable.}
  \item{Z.ref}{a vector specifying the reference values of covariates when estimating marginal effects/predicted values in GLM-type models. Must have the same length and order as the original \bold{Z}. For factor covariates, provide a level name; for numeric covariates, provide a numeric value. By default, numeric covariates are set to their sample mean (factor covariates are handled internally via contrast coding).}
  \item{D.ref}{a numeric vector specifying reference values of \bold{D} when \bold{treat.type} = "continuous" (used for labeling/plotting). If not specified, it defaults to the median of \bold{D}.}
  \item{nbins}{an integer determining the number of bins when discretizing the moderator for \bold{estimator} = "binning" (and also used in \bold{estimator}="raw" for some plots). When quantile cutpoints are not unique, the effective number of bins may be reduced automatically.}
  \item{cutoffs}{a numeric vector of custom cutpoints for \bold{X}. If supplied, \bold{nbins} is ignored and the number of bins equals \code{length(cutoffs) + 1}.}
  \item{wald}{a logical flag controlling whether to conduct a Wald test (used by some estimators such as \bold{estimator}="binning").}
  \item{bw}{a positive number setting the bandwidth in Gaussian kernel reweighting for \bold{estimator} = "kernel".}
  \item{kfold}{an integer specifying the number of folds used in cross-validation for bandwidth selection (kernel estimator).}
  \item{grid}{either an integer indicating the number of candidate bandwidths in a grid search, or a numeric vector of candidate bandwidths (kernel estimator; ignored if \bold{bw} is specified).}
  \item{metric}{a string specifying the criterion used in cross-validation to select bandwidth: one of "MSE", "MAE", "Cross Entropy", or "AUC". Default is "MSE" for non-binary outcomes and "Cross Entropy" for binary outcomes.}
  \item{model.y}{outcome model used by \bold{estimator}="dml" (and also used by \bold{estimator}="lasso"). In the DML estimator, this string is passed to the underlying Python implementation (default "rf").}
  \item{param.y}{named list of hyper-parameters for the DML outcome model (passed to Python).}
  \item{param.grid.y}{named list specifying a hyper-parameter search grid for the DML outcome model (passed to Python; used when \bold{CV} = TRUE).}
  \item{scoring.y}{a string specifying the scoring criterion used for tuning the DML outcome model (passed to Python; default "neg_mean_squared_error").}
  \item{model.t}{treatment/propensity model used by \bold{estimator}="dml" (and also used by \bold{estimator}="lasso"). In the DML estimator, this string is passed to the underlying Python implementation (default "rf").}
  \item{param.t}{named list of hyper-parameters for the DML treatment model (passed to Python).}
  \item{param.grid.t}{named list specifying a hyper-parameter search grid for the DML treatment model (passed to Python; used when \bold{CV} = TRUE).}
  \item{scoring.t}{a string specifying the scoring criterion used for tuning the DML treatment model (passed to Python; default "neg_mean_squared_error").}
  \item{CV}{a logical flag indicating whether to tune the DML nuisance models via cross-validation (DML estimator).}
  \item{n.folds}{an integer specifying the number of folds used in DML cross-fitting (DML estimator).}
  \item{n.jobs}{an integer specifying the number of parallel jobs used in the Python implementation for DML (default -1 uses all available cores).}
  \item{cf.n.folds}{an integer specifying the number of folds for internal cross-fitting (DML estimator).}
  \item{gate}{a logical flag indicating whether to compute group-specific ATEs (GATE) in the DML estimator (if supported by the underlying Python routine).}
  \item{grf.num.trees}{an integer specifying the number of trees for the generalized random forest when \bold{estimator}="grf".}
  \item{signal}{(LASSO estimator) which signal to use when \bold{treat.type}="discrete": one of "outcome", "ipw", or "aipw".}
  \item{estimand}{(LASSO estimator) the target estimand when \bold{treat.type}="discrete": one of "ATE" or "ATT".}
  \item{basis.type}{(LASSO estimator) basis expansion for the moderator: one of "polynomial", "bspline", or "none".}
  \item{include.interactions}{(LASSO estimator) whether to include interactions in the basis expansion.}
  \item{poly.degree}{(LASSO estimator) degree of the polynomial basis when \bold{basis.type}="polynomial".}
  \item{spline.df}{(LASSO estimator) degrees of freedom for B-splines when \bold{basis.type}="bspline".}
  \item{spline.degree}{(LASSO estimator) degree of B-splines when \bold{basis.type}="bspline".}
  \item{lambda.seq}{(LASSO estimator) an optional sequence of \code{lambda} values passed to \code{glmnet}.}
  \item{reduce.dimension}{(LASSO estimator) method used for dimension reduction in the continuous-treatment case (e.g., "bspline" or "kernel").}
  \item{figure}{a logical flag controlling whether to draw a marginal effect plot.}
  \item{bin.labs}{a logical flag indicating whether to show a label for each bin in a binning plot; only meaningful when \bold{nbins} is 2--4.}
  \item{order}{controls the order of panels/series in plots. For discrete treatments, it should contain all treatment arms except the baseline group. For continuous treatments, it can be a numeric vector of the reference \bold{D} values (matching \bold{D.ref}).}
  \item{subtitles}{a vector controlling the subtitles of subplots when \bold{pool} = FALSE, or legend labels when \bold{pool} = TRUE.}
  \item{show.subtitles}{a logical flag controlling whether to show subtitles.}
  \item{Xdistr}{a string indicating how the distribution of the moderator will be plotted: "histogram" (or "hist"), "density", or "none".}
  \item{main}{a string controlling the title of the plot.}
  \item{Ylabel}{a string controlling the label of the outcome variable \bold{Y} in the plot.}
  \item{Dlabel}{a string controlling the label of the treatment variable \bold{D} in the plot.}
  \item{Xlabel}{a string controlling the label of the moderator \bold{X} in the plot.}
  \item{xlab}{a string specifying the label of the x-axis.}
  \item{ylab}{a string specifying the label of the y-axis.}
  \item{xlim}{a two-element numeric vector controlling the x-axis range.}
  \item{ylim}{a two-element numeric vector controlling the y-axis range.}
  \item{theme.bw}{a logical flag specifying whether to use a black-white theme.}
  \item{show.grid}{a logical flag indicating whether to show grid lines in the plot.}
  \item{show.uniform.CI}{a logical flag indicating whether to show uniform/simultaneous confidence intervals when they are available in the output (plotted as dashed lines).}
  \item{cex.main}{a numeric value that controls the font size of the plot title.}
  \item{cex.sub}{a numeric value that controls the font size of subtitles.}
  \item{cex.lab}{a numeric value that controls the font size of the axis labels.}
  \item{cex.axis}{a numeric value that controls the font size of axis tick labels.}
  \item{interval}{draw two vertical lines to demonstrate an interval used in replicated papers.}
  \item{file}{a string specifying the filename to save the plot (if not NULL).}
  \item{ncols}{an integer controlling the number of columns in the plot when \bold{pool} is FALSE.}
  \item{pool}{a logical flag specifying whether to draw all series in a single pooled plot.}
  \item{color}{a vector of colors that determines the color of lines in pooled plots.}
  \item{show.all}{a logical flag specifying whether to return all subplot objects (if TRUE, returns a list of ggplot objects).}
  \item{legend.title}{a string specifying the title of the legend when \bold{pool} is TRUE.}
  \item{diff.values}{a numeric vector (length 2 or 3) specifying moderator values at which treatment effects are compared (only available for some estimators; \bold{estimator}="binning" does not support this feature).}
  \item{percentile}{a logical flag indicating whether \bold{diff.values} should be interpreted on a percentile scale (between 0 and 1).}
  \item{gam.k}{(GAM estimator) dimension of the basis used to represent smooth terms.}
  \item{angle}{(GAM estimator) a vector of numbers indicating viewing angles (allows 1 to 4 values).}
  \item{span}{(RAW estimator) hyper-parameter for the loess fit (larger values produce smoother lines).}
  \item{scale}{multiplicative scaling factor for saving the graph.}
  \item{height}{plot height (in inches) when saving.}
  \item{width}{plot width (in inches) when saving.}
  \item{box.pos}{(RAW estimator) where to place the reference box in discrete-treatment raw plots: "down" (near minimum of \bold{Y}) or "up" (near maximum of \bold{Y}).}
}
\details{
  \bold{interflex} can implement the linear, binning and kernel estimator to estimate the marginal effects of D on Y across the supports of X in linear or generalized linear models. It also provides a few diagnostic tests for classic interaction models and convenient tools for visulization.

  When \bold{estimator}="linear", the program will estimate and visualize the treatment effects/marginal effects of D on Y across the supports of X and the average treatment effects/average effects using classic linear or generalized linear models.

  When \bold{estimator}="binning", There are three steps to implement the estimator. First, the program discretizes the moderator X into several bins and create a dummy variable for each bin. The default is 3 bins, respectively corresponding to the three terciles.  Second, it picks an evaluation point within each bin, where users want to estimate the conditional marginal effect of D on Y. Third, it estimates a model that includes interactions between the bin dummies and the treatment indicator, the bin dummies and the moderator X minus the evaluation points, as well as the triple interactions. The last two terms are to capture the effect of D on Y within each bin.

  The binning estimator has several key advantages over the standard multiplicative interaction model.

  First, the binning estimator is much more flexible as it jointly fits the interaction components of the standard model to each bin separately. The model does not impose the linear interaction effect (LIE) assumption. Instead, the conditional marginal effects can vary freely across the three bins and therefore can take on any non-linear or non-monotonic pattern that might describe the heterogeneity in the effect of D on Y across low, medium, or high levels of X.

  Second, since the bins are constructed based on the support of X, the binning ensures that the conditional marginal effects are estimated at typical values of the moderator and do not rely on excessive extrapolation or interpolation.

  Third, the binning estimator is easy to implement using any regression software and the standard errors for the conditional marginal effects are directly estimated by the regression so there are no need to compute linear combinations of coefficients to recover the conditional marginal effects.

  Fourth, the binning estimator actually provides a generalization that nests the standard multiplicative interaction model as a special case. It can therefore serve as a formal test on the validity of a global LIE assumption imposed by the standard model. So in the special case when the standard multiplicative interaction model is correct and therefore the global LIE assumption holds, then---as the sample size grows---the marginal effect estimates from the binning estimator converge in probability on the unbiased marginal effect estimates from the standard multiplicative interaction model.

  In the case of a discrete treatment, \bold{inter.binning} displays at the bottom of the figure a stacked histogram that shows the distribution of the moderator X. In this histogram the total height of the stacked bars refers to the distribution of the moderator in the pooled sample and the red and gray shaded bars refer to the distribution of the moderator in the treatment and control groups, respectively. Adding such a histogram makes it easy to judge the degree to which there is common support in the data. In the case of a continuous treatment, it displays a histogram at the bottom that simply shows the distribution of X in the entire sample.

  When \bold{estimator}="kernel", the program will implement a kernel smoothing estimator of the marginal effect, which is an application of semi-parametric smooth varying-coefficient models. This approach is more complicated than the simple binning estimator, but provides a generalization that allows researchers to flexibly estimate the functional form of the marginal effect of D on Y across the values of X by estimating a series of local effects with a Gaussian kernel re-weighting scheme.

  The kernel regression nests the standard linear interaction model as a special case. However, in the kernel regression the conditional effect of D on Y does not have to fall on a linear line as required by the linear interaction effect (LIE) assumption, but can vary freely across the range of X. In addition, if covariates Z are included in the model, the coefficients of those covariates are also allowed to vary freely across the range of X resulting in a very flexible estimator that also helps to guard against misspecification bias with respect to the covariates.

}
\value{
  \item{est.lin}{marginal effects/treatment effects estimated using the linear estimator.}
  \item{est.bin}{marginal effects/treatment effects estimated using the binning estimator.}
  \item{est.kernel}{marginal effects/treatment effects estimated using the kernel estimator.}
  \item{pred.lin}{the expected value of Y conditional on X, D and Z estimated using the linear estimator.}
  \item{pred.bin}{the expected value of Y conditional on X, D and Z estimated using the binning estimator.}
  \item{pred.kernel}{the expected value of Y conditional on X, D and Z estimated using the kernel estimator.}
  \item{diff.estimate}{stores the difference between marginal effects/treatment effects at corresponding values specified in \bold{diff.values}.}
  \item{vcov.matrix}{stores the variance-covariance matrix of the series of point estimates.}
  \item{Avg.estimate}{stores the average treatment effects/marginal effects.}
  \item{CV.output}{results from cross-validation if \bold{bw} is not specified.}
  \item{bw}{specified or cross-validated bandwidth.}
  \item{figure}{stores the graphic output.}
}


\author{
  Jens Hainmueller; Jonathan Mummolo; Yiqing Xu (Maintainer);  Ziyi Liu
}
\references{
  Jens Hainmueller; Jonathan Mummolo; Yiqing Xu. 2019. "How Much Should We Trust Estimates from Multiplicative Interaction Models? Simple Tools to Improve Empirical Practice." Political Analysis, Vol. 27, Iss. 2, April 2019, pp. 163--192. Available at: \url{https://www.cambridge.org/core/journals/political-analysis/article/how-much-should-we-trust-estimates-from-multiplicative-interaction-models-simple-tools-to-improve-empirical-practice/D8CAACB473F9B1EE256F43B38E458706}.
}

\examples{
library(interflex)
data(interflex)
s1.binning <-interflex(estimator = 'binning', Y = "Y", D = "D", X = "X",
data = s1)

s1.linear <-interflex(estimator = 'linear', Y = "Y", D = "D", X = "X",
data = s1)


}



