\name{interflex}
\alias{interflex}
\title{Flexible Interactive Models}
\description{Producing flexible marginal effect estimates with linear, binning, and kernel estimators}
\usage{interflex(estimator, data,Y, D, X, treat.type = NULL, 
      base = NULL, Z = NULL, IV = NULL, FE = NULL, 
	  full.moderate = FALSE,
      weights = NULL, na.rm = FALSE, Xunif = FALSE, 
      CI = TRUE, neval = 50, X.eval = NULL, 
	  method = "linear", vartype = "delta", 
	  vcov.type = "robust",time=NULL,
      pairwise = TRUE,nboots = 200, nsimu = 1000, 
      parallel =TRUE,cores = 4, cl = NULL, 
      Z.ref = NULL, D.ref = NULL,
      nbins = 3, cutoffs = NULL, wald = TRUE,
      bw = NULL, kfold = 10, grid = 30,
      metric = NULL, figure = TRUE, bin.labs = TRUE,
      order = NULL, subtitles = NULL, show.subtitles = NULL,
      Xdistr = "histogram", main = NULL,
      Ylabel = NULL, Dlabel = NULL, Xlabel = NULL,
      xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL,
      theme.bw = FALSE, show.grid = TRUE, cex.main = NULL,
      cex.sub = NULL, cex.lab = NULL, cex.axis = NULL,
      interval = NULL, file = NULL, ncols = NULL,
      pool = FALSE, color = NULL, show.all = FALSE,
      legend.title = NULL, diff.values = NULL, 
      percentile = FALSE, gam.k = 10,
	    angle=c(30, 100,-30,-120),
	    span=NULL, scale = 1.1, height = 7, width = 10)
}
\arguments{
  \item{estimator}{a string that specifies the estimator. Can be one of the following: "linear", "binning", "kernel" or "gam".}
  \item{data}{a data.frame that stores Y, D, and X.}
  \item{Y}{a string, name of the outcome variable.}
  \item{D}{a string, name of the treatment variable.}
  \item{X}{a string, name of the moderating variable.}
  \item{treat.type}{a string that specifies the type of the treatment variable. Can be one of the following: "discrete" and "continuous". If not specified, it will be automatically set based on the treatment variable. (If D is in numeric format and the number of multiple arms of D is larger than 5, treat.type will be set to "continuous". Otherwise, it will be set to "discrete".)}
  \item{base}{a string specifying the base group if \bold{treat.type} is "discrete". If not specified, the base group will be selected based on numeric/character order of values in the treatment variable. For example, when treatment D can take values 0, 1, and 2, those observations with treatment 0 will be selected as the base group.}
  \item{Z}{a vector of variable names of control variables.}
  \item{IV}{a vector of variable names for the instrument variables of the treatment, doesn't support glm model in this version}
  \item{FE}{a vector of variable names of fixed effects indicators, doesn't support glm model in this version}
  \item{full.moderate}{a logical flag indicating whether to use fully moderated model that is default to FALSE. If it is set to TRUE, the program will include additional interactions between the moderating variable and other covariates. }
  \item{weights}{a string, name of the weighting variable.}
  \item{na.rm}{a logical flag indicating whether to list-wise delete missing data. The algorithm will report an error if missing data exist.}
  \item{Xunif}{a logical flag controlling whether to transform values of the moderator into their rank order.}  
  \item{CI}{a logical flag indicating whether the confidence intervals need to be shown. If \bold{CI} = FALSE, the bootstrap procedure will not be conducted.}
  \item{neval}{an integer specifying the number of evaluation points. A large number would produce smoother marginal effects estimates but the algorithm will take more time to run.}
  \item{X.eval}{a vector of manually selected evaluation points.}
  \item{method}{a string that specifies the link function. Can be one of the following: "linear", "probit", "logit", "poisson", "nbinom". The default is "linear".}
  \item{vartype}{a string that controls the method used for uncertainty estimation when \bold{estimator} = "linear" or "binning". Can be one of the following: "delta", "simu" and "bootstrap". The default is "delta".}
  \item{vcov.type}{a string that controls the variance-covariance estimator. Can be one of the following: "homoscedastic"; "robust"; "cluster" and "pcse"(only for linear models). The default is "robust".}
  \item{time}{a string specifying the name of time variable for panel corrected standard errors. Ignored when \bold{vcov.type} is not "pcse".}
  \item{pairwise}{a logical flag indicating whether to switch on the \bold{pairwise} option when estimating panel corrected standard errors. Ignored when \bold{vcov.type} is not "pcse". The default is TRUE.}
  \item{nboots}{an integer specifying the number of bootstrap runs. The default is 200.}
  \item{nsimu}{an integer specifying the rounds of simulations. The default is 1000.}
  \item{parallel}{a logical flag indicating whether parallel computing will be used in bootstrapping and/or cross-validation. Ignored if \bold{CI} = FALSE and \bold{bw} is specified.}
  \item{cores}{an integer indicating the number of cores to be used in parallel computing. If not specified, the algorithm will use the maximum number of logical cores of your computer (warning: this could prevent you from multi-tasking on your computer).}
  \item{cl}{a string specifying the name of clustering variable. If specified, the bootstrap and cross-validation procedure will be blocking on this variable. Ignored when \bold{CI} = FALSE. 
    When \bold{vartype} is "cluster", it specifies the name of unit variable for estimating clustered standard errors.}
  \item{Z.ref}{a numeric vector specifies the value of covariates when estimating marginal effects in glm models, must have the same length and order as \bold{Z}. In the default setting, each element of \bold{Z.ref} equals to the mean of the corresponding covariate in the data.}
  \item{D.ref}{a numeric vector specifies the value of \bold{D} when estimating marginal effects when \bold{D} is continuous. The default is the mean of \bold{D} in the data.}
  \item{nbins}{an integer that determines the number of bins of the moderator. The sample will be divided based on integral multiples of 100/nbins percentile.  \bold{nbins} will be automatically
    subtracted by 1 if more than one multiples of the percentile have the same value (for example, if the moderator has over 70\% zeros, both the 33 and 66 percentiles are zero). The default is 3. Ignored when \bold{cutoffs} is supplied. When nbins are set to 1, it only conducts a linear fit.}
  \item{cutoffs}{a vector of numbers that determines how subgroups are divided based on the moderator X. When this option is supplied, \bold{nbins} will be ignored. The smallest number of the first interval and the largest number of the last interval do not need to be specified. Numbers smaller than the minimum or larger than the maximum of X will be ignored. \bold{nbins} equals to the length of this vector plus 1. Ignored if the treatment is dichotomous.}
  \item{wald}{a logical flag controlling whether to conduct a Wald test.}
  \item{bw}{a number that sets the bandwidth in Gaussian kernel reweighting.}
  \item{kfold}{an integer specified the number of folds in cross-validation, default to 10.}
  \item{grid}{either an integer indicating the number of candidates in a grid search to find the optimal bandwidth, or a vector of positive numbers of candidate bandwidths, default to 30. Ignored if \bold{bw} is specified.}
  \item{metric}{a string specifies the criterion to be used in cross-validation to select bandwidth, can be one of the following:"MSE", "MAE", "Cross Entropy" and "AUC". The default will be "MSPE" if the outcome is not binary, otherwise the default will be "Cross Entropy"}
  \item{figure}{a logical flag controlling whether to draw a marginal effect plot.}
  \item{bin.labs}{a logical flag indicating whether to show a label for each bin in a binning plot; only meaningful when \bold{nbins} is 3 or 4. The default is TRUE.}
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
  \item{show.grid}{a logical flag indicating whether to show grid in the plot.}
  \item{cex.main}{a numeric value that controls the font size of the plot title.}   
  \item{cex.sub}{a numeric value that controls the font size of subtitles.}
  \item{cex.lab}{a numeric value that controls the font size of the axis labels.}   
  \item{cex.axis}{a numeric value that controls the font size of the axis numbers.}
  \item{interval}{draw two vertical lines to demonstrate the interval used in replicated papers.}
  \item{file}{a string that specifies the filename in which the plot is saved.}
  \item{ncols}{an integral that controls the number of columns in the plot if \bold{pool} is FALSE.}
  \item{pool}{a logical flag specifying whether to draw marginal effects or expected values of Y for each treatment arm in a single plot. The default is FALSE.}
  \item{color}{a vector of colors that determines the color of lines when drawing the pool plot.}
  \item{show.all}{a logical flag specifying whether to return the plots for each kinds of treatment level. If setting to TRUE, the program will return a list of ggplot2 objects.}
  \item{legend.title}{a string that specifies the title of the legend when \bold{pool} is TRUE.}
  \item{diff.values}{a numeric vector contain 2 or 3 elements which within the range of the moderator. The treatment effects at corresponding values of the moderator will be compared.}
  \item{percentile}{a logical flag indicating whether to take values of the moderator on a percentile scale.}
  \item{gam.k}{an integer indicating the dimension of the basis used to represent the smooth term. Smaller values of k produce smoother surfaces.}
  \item{angle}{a vector of numbers indicating viewing angles (allows 1 to 4 values).}
  \item{span}{hyper-parameter for the lowess fit. Large numbers produce smoother lines (from \bold{ggplot2}).}
  \item{scale}{Multiplicative scaling factor of printed graph.}
  \item{height}{Plot size (Height).}
  \item{width}{Plot size (Width).}
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
  Jens Hainmueller; Jonathan Mummolo; Yiqing Xu. 2019. "How Much Should We Trust Estimates from Multiplicative Interaction Models? Simple Tools to Improve Empirical Practice." Political Analysis, Vol. 27, Iss. 2, April 2019, pp. 163--192. Available at: \url{http://bit.ly/HMX2019}.
}

\examples{
library(interflex)
data(interflex)
s1.binning <-interflex(estimator = 'binning', Y = "Y", D = "D", X = "X", 
data = s1)

s1.linear <-interflex(estimator = 'linear', Y = "Y", D = "D", X = "X", 
data = s1)


}



