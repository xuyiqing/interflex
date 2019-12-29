interflex <- function(data,
						 Y, D, X,  
						 estimator,
						 treat.type = NULL, #discrete or continuous
						 base = NULL,
						 Z = NULL,
						 weights = NULL, # weighting variable
						 FE = NULL,
						 cl = NULL, ## cluster
						 full.moderate = TRUE, ## fully moderated model
						 na.rm = FALSE,
						 Xunif = FALSE, ## transform moderate into a uniform distribution
						 CI = TRUE, ## confidence intervals
						 conf.level = 0.95,
						 
						 #binning
						 nbins = 3,  # No. of X bins
                         cutoffs = NULL,
                         vartype = "robust", # variance type
                         time = NULL, # time variable for pcse
                         pairwise = TRUE, # pcse option
                         wald = TRUE,	
						 bin.labs = TRUE, 
						 
						 #kernel
						 CV.method = 'simple',
                         kfold = 10,
						 grid = 30, # either a number of a sequence of numbers
                         neval = 50,
                         nboots = 200,
                         parallel = TRUE,
                         cores = 4,
                         seed = 02139,
						 bw = NULL,
						 bw.adaptive = TRUE, # adaptive bandwidth
						 quantile.eval = FALSE,
						 metric = "MSPE",
						 
						 # plot
                         figure = TRUE,
						 order = NULL,
						 subtitles = NULL,
						 show.subtitles = NULL,
                         Xdistr = "histogram",
                         main = NULL,
                         Ylabel = NULL,
                         Dlabel = NULL,
                         Xlabel = NULL,
                         xlab = NULL,
                         ylab = NULL,
                         xlim = NULL,
                         ylim = NULL,                         
                         theme.bw = FALSE,
                         interval = NULL,
                         show.grid = TRUE,
                         cex.main = NULL,
						 cex.sub = NULL,
                         cex.lab = NULL,
                         cex.axis = NULL,
                         file = NULL,
                         ncols = NULL,
						 
						 # pool plot
						 pool = FALSE,
						 color = NULL,
						 legend.title = NULL,
						 jitter = FALSE,
						 
						 
						 #predict
						 predict = FALSE,
						 D.ref = NULL,
						 
						 #difference
						 diff.values = NULL
){
	if (!estimator %in% c("linear","binning","kernel") ){
    stop("\"estimator\" must be one of the following: \"linear\",\"binning\",\"kernel\".")
	}
	
	if(estimator=='linear'){
	out <- inter.linear(data=data,
                            Y=Y, # outcome
                            D=D, # treatment indicator
                            X=X, # moderator
                            treat.type = treat.type, #discrete or continuous
                            base = base, # base group when treatments are discrete
                            Z = Z, # covariates
                            FE = FE, # fixed effects
                            weights = weights, # weighting variable
                            full.moderate = full.moderate, # whether use fully moderated model
                            na.rm = na.rm,
                            Xunif = Xunif,
                            vartype = vartype, # variance type
							nboots = nboots,
							parallel = parallel,
							cores = cores,
                            ##  "homoscedastic" (default); "robust"; "cluster", "pcse"
                            cl = cl, # variable to be clustered on
							CI = CI,
                            time = time, # time variable for pcse
                            pairwise = pairwise, # pcse option
							predict = predict,
							D.ref = D.ref,
                            figure = figure,
	                        order = order,
	                        subtitles = subtitles,
	                        show.subtitles = show.subtitles,
                            Xdistr = Xdistr, # c("density","histogram","none")
                            main = main,
                            Ylabel = Ylabel,
                            Dlabel = Dlabel,
                            Xlabel = Xlabel ,
                            xlab = xlab,
                            ylab = ylab,
                            xlim = xlim,
                            ylim = ylim,
                            theme.bw = theme.bw,
                            show.grid = show.grid,
                            cex.main = cex.main,
	                        cex.sub = cex.sub,
                            cex.lab = cex.lab,
                            cex.axis = cex.axis,                      
                            interval = interval,
                            file = file,
                            ncols = ncols,
							color = color,
                            pool = pool,
							legend.title = legend.title,
							diff.values = diff.values)
	
	}
	if(estimator=='binning'){
	out <- inter.binning(data=data,
                            Y=Y, # outcome
                            D=D, # treatment indicator
                            X=X, # moderator
                            treat.type = treat.type, #discrete or continuous
                            base = base, # base group when treatments are discrete
                            Z = Z, # covariates
                            FE = FE, # fixed effects
                            weights = weights, # weighting variable
                            full.moderate = full.moderate, # whether use fully moderated model
                            na.rm = na.rm,
                            Xunif = Xunif,
                            nbins = nbins,  # No. of X bins
                            cutoffs = cutoffs,
                            vartype = vartype, # variance type
                            ##  "homoscedastic" (default); "robust"; "cluster", "pcse"
							nboots = nboots,
							parallel = parallel,
							cores = cores,
                            cl = cl, # variable to be clustered on
							CI = CI,
                            time = time, # time variable for pcse
                            pairwise = pairwise, # pcse option
                            wald = wald,
							predict = predict,
							D.ref = D.ref,
                            figure = figure,
							order = order,
							subtitles = subtitles,
							show.subtitles = show.subtitles,
                            Xdistr = Xdistr, # c("density","histogram","none")
                            main = main,
                            Ylabel = Ylabel,
                            Dlabel = Dlabel,
                            Xlabel = Xlabel ,
                            xlab = xlab,
                            ylab = ylab,
                            xlim = xlim,
                            ylim = ylim,
                            theme.bw = theme.bw,
                            show.grid = show.grid,
                            cex.main = cex.main,
	                        cex.sub = cex.sub,
                            cex.lab = cex.lab,
                            cex.axis = cex.axis,
                            bin.labs = bin.labs,                        
                            interval = interval,
                            file = file,
                            ncols = ncols,
							color = color,
                            pool = pool,
                            jitter = jitter,
							legend.title = legend.title)
	
	}
	if(estimator=='kernel'){
	out <- inter.kernel(data=data,
                         Y=Y, 
						 D=D,
						 X=X, 
                         treat.type = treat.type, #discrete or continuous
                         base = base,
                         Z = Z,
                         weights = weights, # weighting variable
                         FE = FE,
                         full.moderate = full.moderate,
						 CV.method = CV.method,
                         kfold = kfold,
                         na.rm = na.rm,
                         Xunif = Xunif, ## transform moderate into a uniform distribution
                         CI = CI,
                         conf.level = conf.level,
                         cl = cl,
                         neval = neval,
                         nboots = nboots,
                         parallel = parallel,
                         cores = cores,
                         seed = seed,
						 bw = bw,
						 bw.adaptive = bw.adaptive,
						 quantile.eval = quantile.eval,
                         figure = figure,   
						 order = order,
						 subtitles = subtitles,
						 show.subtitles = show.subtitles,
                         grid = grid, # either a number of a sequence of numbers
                         metric = metric,
                         Xdistr = Xdistr,
                         main = main,
                         Ylabel = Ylabel,
                         Dlabel = Dlabel,
                         Xlabel = Xlabel,
                         xlab = xlab,
                         ylab = ylab,
                         xlim = xlim,
                         ylim = ylim,                         
                         theme.bw = theme.bw,
                         interval = interval,
                         show.grid = show.grid ,
                         cex.main = cex.main,
						 cex.sub = cex.sub,
                         cex.lab = cex.lab,
                         cex.axis = cex.axis,
                         file = file,
                         ncols = ncols,
                         pool = pool,
						 color = color ,
						 #predict
						 predict = predict,
						 D.ref = D.ref,
						 legend.title = legend.title,
						 diff.values = diff.values
						 )
	}
	return(out)
}