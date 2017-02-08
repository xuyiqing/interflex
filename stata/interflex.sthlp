{smcl}
{* *! version 1.1  February 3, 2017 @ 01:34:55}{...}
{cmd:help interflex}
{hline}

{title:Title}

{p2colset 5 20 22 2}{...}  {p2col :{hi:interflex} {hline 2}}
 Multiplicative interaction models diagnostics and
 visualization{p_end} {p2colreset}{...}


{title:Syntax}

{p 8 17 2}
{cmdab:interflex}
{it:{help varname:outcome}}
{it:{help varname:treat}}
{it:{help varname:moderator}}
[{it:{help varlist:covar}}]
{ifin}
{it:{weight}} 
[{cmd:,} {it:options}]

{synoptset 23 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt fe(varlist)}} specify fixed effects variables {p_end}
{synopt:{opt type(string)}} specify the estimation strategy, including {opt binning} (default), {opt linear}, and {opt kernal} {p_end}
{synopt:{opt vce(vcetype)}} specify the variance–covariance estimator; vcetype can be {opt homoscedastic} (default), {opt robust}, {opt cluster}, {opt boot:strap}, or {opt off} {p_end}
{synopt:{opt cl:uster(varname)}} specify the clustering variable for clustered standard errors in linear and binning estimtions and for the block bootstrap procedure in the kernel estimation {p_end}
{synopt:{opt n:bins(integer)}} set the number of bins; default is {opt 3}; {opt type(binning)} is required; ignored if {opt cutoffs} is supplied {p_end}
{synopt:{opt cut:offs(numlist)}} set the cutoff points (in ascending order) to determine bins of the moderator in {opt binning} estimations {p_end}
{synopt:{opt sav:ing(string)}} specify a filename with which the graph will be stored; if no suffix is provided, the graph will be saved as a ".pdf" {p_end}

{syntab:Advanced}
{synopt:{opt reps(integer)}} specify the number of bootstrap runs; the default is {opt 50} {p_end}
{synopt:{opt seed(integer)}} set the random seed for cross validation or bootstrapping {p_end}
{synopt:{opt bw(real)}} set the bandwidth for {opt kernel} estimations {p_end}
{synopt:{opt grid(integer)}} set the number of candidates in a grid search to determine the optimal bandwidth; default is {opt 20}; {opt type(kernel)} is required; ignored if {cmd: bw} is supplied {p_end}
{synopt:{opt neval(integer)}} set the number of evaluation points in the kernel estimation; default is {opt 50}; {opt type(kernel)} is required {p_end}
{synopt:{opt ti:tle(string)}} set the title of the graph {p_end}
{synopt:{opt xr:ange(numlist)}} set the range of values to be shown on the x-axis {p_end}
{synopt:{opt yr:ange(numlist)}} set the range of values to be shown on the y-axis {p_end}
{synopt:{opt ylab:el(string)}} set the label for the outcome variable {p_end}
{synopt:{opt dlab:el(string)}} set the label for the treatment variable {p_end}
{synopt:{opt xlab:el(string)}} set the label for the moderating variable {p_end}
{synopt:{opt xd:istr(string)}} specify how the distribution of the moderator should be visualized; either {opt hist:ogram} (default) or {opt density} {p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
{p_end}


{title:Description}

{pstd} {opt interflex} performs diagnostics and visualizations of
multiplicative interaction models. Besides conventional linear
interaction models, it provides two additional estimation
strategies--linear regression based on pre-specified bins and locally
linear regressions based on Gaussian kernel reweighting--to estimate the
conditional marginal effect of a treatment variable on an outcome
variable across different values of a moderating variable. These
approaches relax the linear interaction effect assumption and
safeguard against excessive extrapolation. {p_end}


{title:Options}

{dlgtab:Main}


{phang} {opt fe(varlist)} specifies fixed effects variables in linear
fixed effect models. We use a new STATA command that implements a fast
routine of high dimensional linear fixed effects models {help reghdfe: reghdfe} developed by Sergio Correia. Install via {stata "ssc install reghdfe":ssc install reghdfe}.

{phang} {opt type(string)} specifies the type of estimation
strategy. The default is {cmd:type(binning)}, which plots the linear
marginal effects superimposed by the binning estimates (at low,
medium, and high levels of the moderator if there are 3 bins, for
instance). {cmd:type(linear)} plots the conventional linear margianl
effects (Brambor, Clark, and Golder 2006). {cmd:type(kernel)} plots the
marginal effects based on a kernel smoothing estimator. It produces
the marginal effect estimates of the treatment on the outcome at a
series of values of the moderator using a kernel-weighted locally
linear regressions.

{phang} {opt vce(string)} specifies the variance–covariance
 estimator. vcetype can be {opt homoscedastic} (default), {opt robust}, {opt cluster}, {opt boot:strap}, or {opt off}. {opt vce(off)} only applies to the kernel estimator. When supplied, uncertainty estimates will not be produced and hence confidence
 intervals will not appear in the graph.

{phang} {opt cl:uster(varname)} specifies the clustering variable for
clustered standard errors when {opt vce(cluster)} is chosen and for the blocking variable in bootstrap procedure when {opt vce(bootstrap)} is chosen. The current algorithm only allows one
way clustering. 

{phang} {opt n:bins(integer)} sets the number of bins for
{opt type(binning)} estimaion; the default is {opt 3}. When
{opt cutoffs(numlist)} are not specified, the binning estimates are positioned
at the median value of the moderator within each bin.  {opt nbins} will be automatically
    subtracted by 1 if more than one multiples of the percentile have
    the same value (for example, if the moderator has over 70% zeros,
    both the 33 and 66 percentiles are zero). Ignored
    when {opt cutoffs} is supplied.

{phang}{opt cut:offs(numlist)} sets the cutoff points for moderator bins in {opt type(binning)}
estimation. When supplied, {opt nbins} will be ignored. The smallest number of the first interval and the
    largest number of the last interval do not need to be specified. Numbers smaller than the minimum or larger than the
    maximum of X will be ignored. The number of bins equals to the length of numlist plus 1. Ignored if the treatment is dichotomous.

{phang}{opt sav:ing(filename)} specifies a filename with which the graph will be stored; if no suffix is provided, the graph will be saved as a ".pdf".
{dlgtab:Advanced}

{phang}{opt reps(integer)} specifies the number of bootstrap runs; the default is {opt 50}; 

{phang}{opt seed(integer)} sets the random seed for cross-validation and/or boostrapping in {opt kernel} estimations, the default is {opt seed(12345678)}. 

{phang} {opt bw(real)} sets the bandwidth for the {opt kernel} estimation; If {opt bw} is not specified, the optimal bandwidth
will be selected via a least-squares cross-validation procedure.

{phang} {opt grid(integer)} sets the number of candidates in a grid
search to determine the optimal bandwidth using cross-validation. Ignored when {opt bw} is supplied. The default is {opt 20}.

{phang}{opt neval(integer)} sets the number of evaluation points in {opt kernel} estimations.

{phang}{opt ti:tle(string)} sets the title of the graph.

{phang}{opt xr:ange(numlist)} sets the range of values to be shown on the x-axis. Only the first two numbers (in ascending order) will be considered.

{phang}{opt yr:ange(numlist)} sets the range of values to be shown on the y-axis. Only the first two numbers (in ascending order) will be considered.

{phang}{opt ylab:el(string)} sets the label for the outcome variable.

{phang}{opt dlab:el(string)} sets the label for the treatment variable.

{phang}{opt xlab:el(string)} sets the label for the moderating variable.

{phang}{opt xd:istr(string)} specifies the way in which the distribution of the moderator is visualized; either {opt hist:ogram} (default) or {opt density}.

{title:Examples}

We provide four simulated datasets to illustrate how {cmd:interflex} works. 

{pstd}{ul:Dichotomous or countinuous treatment with linear marginal effects}{p_end}

{pstd}Load sample1 (dichotomous treatment, linear) {p_end}
{p 4 8 2}{stata "sysuse sample1.dta, clear":. sysuse sample1.dta, clear}{p_end}

{pstd}Plot the raw data{p_end}
{p 4 8 2}{stata "twoway (sc Y X) (lowess Y X), by(D)":. twoway (sc Y X) (lowess Y X), by(D)}{p_end}

{pstd}Estimate a linear interaction model{p_end}
{p 4 8 2}{stata "interflex Y D X Z1, type(linear)":. interflex Y D X Z1, type(linear)}{p_end}

{pstd}Linear interaction model with robust standard errors{p_end}
{p 4 8 2}{stata "interflex Y D X Z1, type(linear) vce(robust)":. interflex Y D X Z1, type(linear) vce(robust)}{p_end}

{pstd}Linear interaction model with bootstrap standard errors{p_end}
{p 4 8 2}{stata "interflex Y D X Z1, type(linear) vce(boot)":. interflex Y D X Z1, type(linear) vce(boot)}{p_end}

{pstd}Superimpose three binning estimates (by default) with robust standard errors {p_end}
{p 4 8 2}{stata "interflex Y D X Z1, vce(r)":. interflex Y D X Z1, vce(r)}{p_end}
{p 6 6 2}When the linear interaction model is correct, the binning estimates are consistent with the linear marginal effects.{p_end}

{pstd}Estimate potentially nonlinear marginel effects based on the kernel estimator with specified bandwidth {p_end}
{p 4 8 2}{stata "interflex Y D X Z1, type(kernel) bw(5.5)":. interflex Y D X Z1, type(kernel) bw(5.5)}{p_end}

{pstd}Optimal bandwidth selected via cross-validation (THIS's GONNA TAKE A WHILE)) {p_end}
{p 4 8 2}{stata "interflex Y D X Z1, type(kernel)":. interflex Y D X Z1, type(kernel)}{p_end}
{p 6 6 2}Results are similar to what we got from the lienar
interaction model, as they should be. The optimal bandwithd is
relatively large.{p_end}

{pstd}Similar for the case of continuous treatment; first, plot the raw data{p_end}
{p 4 8 2}{stata "sysuse sample2.dta, clear":. sysuse sample2.dta, clear}{p_end}
{p 4 8 2}{stata "egen Xbin = cut(X), group(3)":. egen Xbin = cut(X), group(3)}{p_end}
{p 4 8 2}{stata "twoway (sc Y D) (lowess Y D), by(Xbin)":. twoway (sc Y D) (lowess Y D), by(Xbin)}{p_end}

{pstd}Kernel estimates{p_end}
{p 4 8 2}{stata "interflex Y D X Z1, type(kernel) bw(5.0)":. interflex Y D X Z1, type(kernel) bw(5.0)}{p_end}

{pstd}{ul:Nonlinear marginal effects}{p_end}

{pstd}Plot the raw data{p_end}
{p 4 8 2}{stata "sysuse sample3.dta, clear":. sysuse sample3.dta, clear}{p_end}
{p 4 8 2}{stata "twoway (sc Y X) (lowess Y X), by(D)":. twoway (sc Y X) (lowess Y X), by(D)}{p_end}

{pstd}Apply the binning estimator{p_end}
{p 4 8 2}{stata "interflex Y D X Z1, vce(r)":. interflex Y D X Z1, vce(r)}{p_end}
{p 6 6 2}When the linear interaction model is incorrect, the binning estimates diagree with the marginal effect estimates from the conventional linear model.{p_end}

{pstd}Apply the kernel estimator (bandwidth selected via cross-validation) {p_end}
{p 4 8 2}{stata "interflex Y D X Z1, type(kernel) bw(0.345)":. interflex Y D X Z1, type(kernel) bw(0.345)}{p_end}
{p 6 6 2}Marginal effect estimates using the kernel estimator are consistent with the true DGP.{p_end}

{pstd}{ul:Nonlinear marginal effects with additive two-way fixed effects}{p_end}

{pstd}Load sample4{p_end}
{p 4 8 2}{stata "sysuse sample4.dta, clear":. sysuse sample4.dta, clear}{p_end}
{p 6 6 2}The two fixed effect indicators are "group" and "year".{p_end}

{pstd}Plot the raw data{p_end}
{p 4 8 2}{stata "twoway (sc Y X) (lowess Y X), by(D)":. twoway (sc Y X) (lowess Y X), by(D)}{p_end}
{p 6 6 2}With fixed effects affecting the outcome, we cannot observe a clear pattern of marginal effects in the raw plot as before.{p_end}

{p 4 8 2}{stata "interflex Y D X Z1, cl(group)":. interflex Y D X Z1, cl(group)}{p_end}
{p 6 6 2}The binning estimates have wide confidence intervals if fixed effects are not controlled for.{p_end}

{pstd} Controlling for fixed effects with the binning estimator {p_end}
{p 4 8 2}{stata "interflex Y D X Z1, fe(group year) cl(group)":. interflex Y D X Z1, fe(group year) cl(group)}{p_end}

{pstd} Controlling for fixed effects with the kernel estimator (bandwith is selected via cross-validation) {p_end}
{p 4 8 2}{stata "interflex Y D X Z1, type(kernel) fe(group year) cl(group) bw(0.40)":. interflex Y D X Z1, type(kernel) fe(group year) cl(group) bw(0.40)}{p_end}


{title:Saved results}
{pstd}{cmd:interflex} returns the estimated marginal effect of treatment variable on dependent variable (as well as its standard error and the 95% confidence interval) at each evaluation point of moderating variable. {p_end}

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:e(margeff)}} estimated marginal effects {p_end}
{synopt:{cmd:e(estBin)}} binning estimates {p_end}
{synopt:{cmd:e(bandiwth)}} specified or cross-validated bandwith if {opt type(kernel)}{p_end}
{synopt:{cmd:e(CVout)}} mean squared prediction errors in cross-validation {p_end}

{title:Reference}
{p 4 8 2}

{pstd}Brambor, Thomas, William Roberts
Clark and Matt Golder (2006). "Understanding Interaction Models:
Improving Empirical Analyses." Political Analysis 14:63–82.{p_end}

{pstd}Hainmueller, Jens, Mummolo, Jonathan, and Xu, Yiqing
(2016). "How Much Should We Trust Estimates from Multiplicative
Interaction Models? Simple Tools to Improve Empirical Practice."
Available at SSRN: https://papers.ssrn.com/abstract_id=2739221.{p_end}

{title:Authors}

    Yiqing Xu (Maintainer), yiqingxu@ucsd.edu
    Univeristy of California, San Diego

    Jens Hainmueller, jhain@stanford.edu
    Stanford University

    Jonathan Mummolo, jmummolo@stanford.edu
    Stanford University

    Licheng Liu, liulch.16@sem.tsinghua.edu.cn
    Tsinghua University



