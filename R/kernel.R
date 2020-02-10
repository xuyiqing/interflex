inter.kernel <- function(data,
                         Y, D, X, 
                         treat.type = NULL, #discrete or continuous
                         base = NULL,
                         Z = NULL,
                         weights = NULL, # weighting variable
                         FE = NULL,
                         full.moderate = TRUE,
						 na.rm = FALSE,
                         Xunif = FALSE, ## transform moderate into a uniform distribution
                         CI = TRUE,
                         conf.level = 0.95,
                         cl = NULL,
						 CV.method = NULL,
                         kfold = 10,
						 grid = 30, # either a number of a sequence of numbers
                         neval = 50,
                         nboots = 200,
                         parallel = TRUE,
                         cores = 4,
                         seed = 02139,
						 bw = NULL,
						 bw.adaptive = TRUE,
						 quantile.eval = FALSE,
                         metric = "MSPE",
						 predict = FALSE,
						 D.ref = NULL,
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
                         pool = FALSE,
						 color = NULL,
						 legend.title = NULL,
						 diff.values = NULL
){
  
  x <- NULL
  y <- NULL
  xmin <- NULL
  xmax <- NULL
  ymin <- NULL
  ymax <- NULL
  ME <- NULL
  CI_lower <- NULL
  CI_upper <- NULL
  diff.table <- NULL
 
if(TRUE){ #INPUT CHECK
## in case data is in tibble format
  if (is.data.frame(data) == FALSE) {
    data <- as.data.frame(data)
  }
  n<-dim(data)[1]
  ## Y
  if (is.character(Y) == FALSE) {
    stop("\"Y\" is not a string.")
  } else {
    Y <- Y[1]
  }
  ## D
  if (is.character(D) == FALSE) {
    stop("\"D\" is not a string.")
  } else {
    D <- D[1]
  }
  ## X
  if (is.character(X) == FALSE) {
    stop("\"X\" is not a string.")
  } else {
    X <- X[1]
  }
  ## Z
  if (is.null(Z) == FALSE) {
    for (i in 1:length(Z)) {
      if (is.character(Z[i]) == FALSE) {
        stop("Some element in Z is not a string.")
      }
    }
  }
  
  ## weights
  if (is.null(weights) == FALSE) {
    if (is.character(weights) == FALSE) {
      stop("\"weights\" is not a string.")
    } else {
      dataweights <- data[,weights]
    }   
  }
  
  ## FE
  if (is.null(FE) == FALSE) {
    requireNamespace("lfe")
    for (i in 1:length(FE)) {
      if (is.character(FE[i]) == FALSE) {
        stop("Some element in FE is not a string.")
      }
    }
  }
  
  ## full moderate model
  if (is.logical(full.moderate) == FALSE & is.numeric(full.moderate)==FALSE) {
    stop("full.moderate is not a logical flag.")
  }else{
    full <- full.moderate
  }
  
  ## na.rm
  if (is.logical(na.rm) == FALSE & is.numeric(na.rm)==FALSE) {
    stop("\"na.rm\" is not a logical flag.")
  }

  
  ## Xunif
  if (is.logical(Xunif) == FALSE & is.numeric(Xunif)==FALSE) {
    stop("\"Xunif\" is not a logical flag.")
  }
  
  ## CI
  if (is.logical(CI) == FALSE & is.numeric(CI)==FALSE) {
    stop("\"CI\" is not a logical flag.")
  }
  
  ## conf.level
  if (is.null(conf.level)==FALSE) {
    if (is.numeric(conf.level)==FALSE) {
      stop("\"conf.level\" should be a number between 0.5 and 1.")
    } else {
      if (conf.level<=0.5 | conf.level>1) {
        stop("\"conf.level\" should be a number between 0.5 and 1.")
      }
    } 
  }
  
  ## cl
  if (is.null(cl)==FALSE) {
    if (is.character(cl) == FALSE) {
      stop("\"cl\" is not a string.")
    } else {
      cl <- cl[1] 
    }
  }
  
  ## check missing values
  vars <- c(Y, D, X, Z, FE, cl, weights)
  if (na.rm == TRUE) {        
    data <- na.omit(data[,vars])
  } else {
    if (sum(is.na(data[,vars]))>0) {
      stop("Missing values. Try option na.rm = TRUE\n")
    }
  }
  
  
  ## CV.method
  if(is.null(CV.method)==FALSE){
	if (!CV.method %in% c("simple","cluster","stratify") ){
		stop("\"treat.type\" must be one of the following: \"simple\",\"cluster\",\"stratify\".")
	}
	if(CV.method=='cluster' & is.null(cl)==TRUE){
		stop("\"cl\" is not specified.")
	}
  }else{
	if(is.null(cl)==FALSE){
		CV.method <- 'cluster'
	}else{CV.method <- 'simple'}
  }
  
  ## kfold
  if (is.numeric(kfold)==FALSE) {
    stop("\"kfold\" is not a positive integer.")
  } else {
    kfold <- kfold[1]
    if (kfold%%1!= 0 | kfold<=0) {
      stop("\"kfold\" is not a positive integer.")
    }
	if (kfold<3) {
      stop("\"kfold\" should be greater than 3.")
    }
  }
  
  ## grid
  if (is.numeric(grid) == FALSE) {
    stop("\"grid\" should be numeric.")
  } else {
    if (length(grid)==1) {
      if (grid%%1 != 0 | grid<1) {
        stop("\"grid\" is not a positive integer.")
      }
    } else {
      grid <- grid[which(grid>0)]
    }
  }
  
  ## neval
  if (is.null(neval)==FALSE) {
    if (is.numeric(neval)==FALSE) {
      stop("\"neval\" is not a positive integer.")
    } else {
      neval <- neval[1]
      if (neval%%1!=0 | neval<=0) {
        stop("\"neval\" is not a positive integer.")
      }  
    } 
  } 
  
  ## nboots
  if (is.null(nboots) == FALSE) {
    if (is.numeric(nboots)==FALSE) {
      stop("\"nboots\" is not a positive integer.")
    } else {
      nboots <- nboots[1]
      if (nboots%%1 != 0 | nboots < 1) {
        stop("\"nboots\" is not a positive number.")
      }
    } 
  }
  
  ## parallel
  if (is.logical(parallel) == FALSE & is.numeric(parallel)==FALSE) {
    stop("\"paralell\" is not a logical flag.")
  }
  
  ## cores
  if (is.numeric(cores)==FALSE) {
    stop("\"cores\" is not a positive integer.")
  } else {
    cores <- cores[1]
    if (cores%%1!= 0 | cores<=0) {
      stop("\"cores\" is not a positive integer.")
    }
  }
  
  #seed
  if (is.numeric(seed)==FALSE) {
    stop("\"seed\" should be a number.")
  }
  
  #bw
  if (is.null(bw)==FALSE) {
    if (is.numeric(bw)==FALSE) {
      stop("\"bw\" should be a positive number.")
    } else {
      bw <- bw[1]
    } 
    if (bw<=0) {
      stop("\"bw\" should be a positive number.")
    }
  }
  
  #bw.adaptive
  if (is.logical(bw.adaptive) == FALSE & is.numeric(bw.adaptive)==FALSE) {
    stop("\"bw.adaptive\" is not a logical flag.")
  }
  
  #quantile.eval
  if (is.logical(quantile.eval) == FALSE & is.numeric(quantile.eval)==FALSE) {
    stop("\"quantile.eval\" is not a logical flag.")
  }
  
  #metric
  if (!metric%in%c("MSPE","MAPE")) {
    stop("\"metric\" should be either \"MSPE\" or \"MAPE\".")
  }
  
  # predict
  if (is.logical(predict) == FALSE & is.numeric(predict)==FALSE) {
    stop("\"predict\" is not a logical flag.")
  }
  
  # D.ref
  if(is.null(D.ref)==FALSE){
	if(is.numeric(D.ref)==FALSE){
      stop("\"D.ref\" is not a numeric variable.")
    }
	if(length(D.ref)>=10){
	  stop("Too many elements in \"D.ref\".")
	}
  }
  
    # figure
  if (is.logical(figure) == FALSE & is.numeric(figure)==FALSE) {
    stop("\"figure\" is not a logical flag.")
  }
  
  # show.subtitles
  if(is.null(show.subtitles)==FALSE){
	if (is.logical(show.subtitles) == FALSE & is.numeric(show.subtitles)==FALSE) {
		stop("\"show.subtitles\" is not a logical flag.")
	}
  }
  
  # Xdistr
  if (!Xdistr %in% c("hist","histogram","density","none")){
    stop("\"Xdistr\" must be \"histogram\", \"density\", or \"none\".")
  }
  
  # main
  if (is.null(main)==FALSE) {
    main <- as.character(main)[1]
  }
  # Ylabel
  if (is.null(Ylabel)==TRUE) {
    Ylabel <- Y
  } else {
    if (is.character(Ylabel) == FALSE) {
      stop("Ylabel is not a string.")
    } else {
      Ylabel <- Ylabel[1]
    }   
  }
  # Dlabel  
  if (is.null(Dlabel)==TRUE) {
    Dlabel <- D   
  } else {
    if (is.character(Dlabel) == FALSE) {
      stop("Dlabel is not a string.")
    } else {
      Dlabel <- Dlabel[1]
    }   
  }
  
  #Xlabel
  if (is.null(Xlabel)==TRUE) {
    Xlabel <- X   
  } else {
    if (is.character(Xlabel) == FALSE) {
      stop("Xlabel is not a string.")
    } else {
      Xlabel <- Xlabel[1]
    }   
  }
  
  ## axis labels
  if(is.null(xlab)==FALSE){
    if (is.character(xlab) == FALSE) {
      stop("\"xlab\" is not a string.")
    }        
  }
  if(is.null(ylab)==FALSE){
    if (is.character(ylab) == FALSE) {
      stop("\"ylab\" is not a string.")
    }        
  }

  ## xlim ylim
  if (is.null(xlim)==FALSE) {
    if (is.numeric(xlim)==FALSE) {
      stop("Some element in xlim is not numeric.")
    } else {
      if (length(xlim)!=2) {
        stop("xlim must be of length 2.")
      }
    }
  }
  if (is.null(ylim)==FALSE) {
    if (is.numeric(ylim)==FALSE) {
      stop("Some element in ylim is not numeric.")
    } else {
      if (length(ylim)!=2) {
        stop("ylim must be of length 2.")
      }
    }
  }
  
  ## theme.bw
  if (is.logical(theme.bw) == FALSE & is.numeric(theme.bw)==FALSE) {
    stop("\"theme.bw\" is not a logical flag.")
  }
  
  # interval
  if (is.null(interval)==FALSE) {
	if (is.numeric(interval)==FALSE) {
      stop("Some element in interval is not numeric.")
    } 
  }
  
  ## show.grid
  if (is.logical(show.grid) == FALSE & is.numeric(show.grid)==FALSE) {
    stop("\"show.grid\" is not a logical flag.")
  }
  
  ## font size
  if (is.null(cex.main)==FALSE) {
    if (is.numeric(cex.main)==FALSE) {
      stop("\"cex.main\" is not numeric.")
    }
  }
  if (is.null(cex.sub)==FALSE) {
    if (is.numeric(cex.sub)==FALSE) {
      stop("\"cex.sub\" is not numeric.")
    }
  }
  if (is.null(cex.lab)==FALSE) {
    if (is.numeric(cex.lab)==FALSE) {
      stop("\"cex.lab\" is not numeric.")
    }
  }
  if (is.null(cex.axis)==FALSE) {
    if (is.numeric(cex.axis)==FALSE) {
      stop("\"cex.axis\" is not numeric.")
    }    
  }
  
  # file
  if (is.null(file)==FALSE) {
	if (is.character(file)==FALSE) {
      stop("Wrong file name.")
    } 
  }
  
  # ncols
  if (is.null(ncols) == FALSE) {
    if (is.numeric(ncols)==FALSE) {
      stop("\"ncols\" is not a positive integer.")
    } else {
      ncols <- ncols[1]
      if (ncols%%1 != 0 | ncols < 1) {
        stop("\"ncols\" is not a positive number.")
      }
    } 
  }
  
  ## pool
  if (is.logical(pool) == FALSE & is.numeric(pool)==FALSE) {
    stop("\"pool\" is not a logical flag.")
  }
  
  ## color
  if(is.null(color)==FALSE){
	color <- as.character(color)
	color.in <- c()
	for(char in color){
		res <- try(col2rgb(char),silent=TRUE)
		if(!"try-error"%in%class(res)){
			color.in <- c(color.in,char)
		}else{stop(paste0(char," is not a color name.\n"))}
	}
	color <- color.in
  }
  
  ## legend.title
  if (is.null(legend.title)==FALSE) {
    legend.title <- as.character(legend.title)[1]
  }
  
  ## diff.values
  if(is.null(diff.values)==FALSE){
	if(is.numeric(diff.values)==FALSE){
		stop("\"diff.values\" is not numeric.")
	}
	if(length(diff.values)!=3){
		stop("\"diff.values\" must be of length 3.")
	}
	min.XX <- min(data[,X])
	max.XX <- max(data[,X])
	for(a in diff.values){
		if(a<min.XX|a>max.XX){
			stop("Elements in \"diff.values\" should be greater than the minimum and less than the maximum of the moderator.")
		}
	}
  }
}
  
if(TRUE){ #TREAT CHECK

  show.subtitle <- show.subtitles
  subtitle <- subtitles
  
  if(is.null(treat.type)==T){
    if(is.numeric(data[,D])==T){
      if(length(unique(data[,D]))>5){
        treat.type <- 'continuous'
      }
      else{
        treat.type <- 'discrete'
      }
    }
    else{
      treat.type <- 'discrete'
    }
  }
  
  if (!treat.type %in% c("discrete","continuous") ){
    stop("\"treat.type\" must be one of the following: \"discrete\",\"continuous\".")
  }

  if (treat.type=='discrete') {
    if(length(unique(data[,D]))>5) {
      warning("More than 5 kinds of treatments")
    }
    
    if(length(unique(data[,D]))>10) {
      stop("Too many kinds of treatments, treatments may be continuous")
    }
    
    data[,D] <- as.character(data[,D])
    
    
    if(is.null(base)==TRUE) {
      base=sort(unique(data[,D]))[1]
      f=sprintf("Base group has not been specified, choose treat = %s as base group.\n",base)
      cat(f)
    }
    else {
      base <- as.character(base)
      if (!base %in% unique(data[,D])){
        stop("\"base\" must be one kind of treatment")
      }
      
      f=sprintf("Base group: treat = %s \n",base)
      cat(f)
    }
    

	
    all_treat=sort(unique(data[,D]))
	ntreat <- length(all_treat)
	names(all_treat) <- paste("Group",c(1:length(all_treat)),sep  = '.')
    other_treat <- all_treat[which(all_treat!=base)]
    other_treat <- sort(other_treat)
    nothertreat <- length(other_treat)
	
    if(is.null(order)==F){
	  
	  order <- as.character(order)
	  
	  if(length(order)!=length(unique(order))){
        stop("\"order\" should not contain repeated values.")
      }
	  
      if(length(order)!=length(other_treat)){
        stop("\"order\" should include all kinds of treatment arms except for the baseline category.")
	  }

      if(sum(!is.element(order,other_treat))!=0 | sum(!is.element(other_treat,order))!=0){
        stop("\"order\" should include all kinds of treatment arms except for the baseline category.")
	  }
	  
      other_treat <- order
	  
	  colnames.p <- c()
	  for(char in other_treat){
		colnames.p <- c(colnames.p,names(all_treat[which(all_treat==char)]))
	  }
	  names(other_treat) <- colnames.p
	}
	
	all_treat.origin <- all_treat
	other_treat.origin <- other_treat
	base.origin <- base
	data.origin <- data
	for(char in names(all_treat)){
		data[which(data[,D]==all_treat[char]),D] <- char
	}

	all_treat <- names(all_treat.origin)
	other_treat <- names(other_treat.origin)
	base <- names(all_treat.origin[which(all_treat.origin==base.origin)])
	names(all_treat) <- all_treat.origin
	names(other_treat) <- other_treat.origin
	
	
	if(is.null(subtitle)==F){
	    if(length(subtitle)!=length(other_treat)){
	      stop("The number of elements in \"subtitles\" should be m-1(m is the number of different treatment arms).")
      }
	  }
	  
	  if (is.logical(show.subtitle) == F & is.numeric(show.subtitle)==F & is.null(show.subtitle)==F) {
	    stop("\"show.subtitles\" is not a logical flag.")
	  }
    
  }
  
  if (treat.type=='continuous') {
    if (is.numeric(data[,D])==F) {
      stop("\"D\" is not a numeric variable")
    }
	#if(predict==TRUE){
	if(is.null(D.ref)==TRUE){
    D.sample <- quantile(data[,D],probs = c(0.25,0.5,0.75),na.rm=T)
	all_treat <- names(D.sample)
    ntreat <- length(D.sample)
    labelname <- c()
    for (targetD in D.sample){
      labelname <- c(labelname,paste0("D=",round(targetD,2)))
    }
    labelname <- paste0(labelname,' (',all_treat,')')
	} else{
	if (is.numeric(D.ref)==F) {
      stop("\"D.ref\" is not a numeric variable")
    }
	D.sample <- D.ref
	labelname <- c()
    for (targetD in D.sample){
      labelname <- c(labelname,paste0("D=",round(targetD,2)))
    }
	names(D.sample) <- labelname
	all_treat <- labelname
    ntreat <- length(D.sample)
	}
	data.origin <- data
	all_treat.origin <- all_treat
	#}
  }

}
  
####################################################
  
if(TRUE){ # Preprocess

  if (is.null(ncols) == FALSE) {
    if (ncols%%1 != 0) {
      stop("\"ncols\" is not a positive integer.")
    } else {
      ncols <- ncols[1]
    }
    if (ncols < 1) {
      stop("\"ncols\" is not a positive integer.")
    }
  } else{
  if(treat.type=="discrete"){
	ncols <- length(other_treat)
  }
  if(treat.type=="continuous"){
	ncols <- 1
  }
  }
 
  n<-dim(data)[1]
  # factor covariates
  to_dummy_var <- c()
  for(a in Z){
	if(is.factor(data[,a])==TRUE){
		to_dummy_var <- c(to_dummy_var,a)
	}	
  }
  if(length(to_dummy_var)>0){
	fnames <- paste("factor(", to_dummy_var, ")", sep = "")
	contr.list <- list(contr.sum, contr.sum)
	names(contr.list) <- fnames
	to_dummy_form <- as.formula(paste("~", paste(fnames, collapse = " + ")))
	suppressWarnings(
	to_dummy_mat <- model.matrix(to_dummy_form, data = data,
                          contrasts.arg = contr.list)[, -1]
	)
	to_dummy_mat <- as.matrix(to_dummy_mat)
	dummy_colnames <- c()
	for(i in 1:dim(to_dummy_mat)[2]){
		dummy_colnames <- c(dummy_colnames,paste0("Dummy.Covariate.",i))
	}
	colnames(to_dummy_mat) <- dummy_colnames
	data <- cbind(data,to_dummy_mat)
	Z <- Z[!Z %in% to_dummy_var]
	Z <- c(Z,dummy_colnames)
  }

# fully moderated model
  if(full==T){
    cat("Use a fully moderated model.\n")
    full_names <- names(data)
    names_Z <- Z
    for(char in Z){
      fm_name <- paste0(X,"_",char)
      names_Z <- c(names_Z,fm_name)
      full_names <- c(full_names,fm_name)
      data <- cbind(data,data[,X]*data[,char])
    }
    names(data) <- full_names
    Z <- names_Z
  }
  
  if(bw.adaptive==T){
    cat("Use adaptive bandwidth.\n")
  }
  
  
  ## transform moderator into uniform
  if (Xunif == TRUE) {
    x <- data[, X]
    data[, X] <- rank(x, ties.method = "average")/length(x)*100
    Xlabel <- paste(Xlabel,"(Percentile)")
  }
 
  ## change fixed effect variable to factor
  if (is.null(FE)==FALSE) {
    if (length(FE) == 1) {
      data[, FE] <- as.numeric(as.factor(data[, FE]))
    } else {
      data[, FE] <- sapply(data[,FE],function(vec){as.numeric(as.factor(vec))})
    }
  }
  
  ## labels
  if (is.null(Xlabel)==TRUE) {Xlabel = X}
  if (is.null(Ylabel)==TRUE) {Ylabel = Y}
  if (is.null(Dlabel)==TRUE) {Dlabel = D}
  
  ## X more than 5 values
  if (length(unique(data[ ,X]))< 5) {
    warning("Moderator has less than 5 values; consider a fully saturated model.")
  }
  
  if(is.null(diff.values)==TRUE){
	diff.values.plot <- NULL
	diff.values <- quantile(data[,X],probs = c(0.25,0.5,0.75))
	difference.name <- c("50% vs 25%","75% vs 50%","75% vs 25%")
  }else{
	diff.values.plot<-diff.values
	#difference.name <- c("2nd-1st","3rd-2nd","3rd-1st")
	difference.name <- c(paste0(diff.values[2]," vs ",diff.values[1]),
				   paste0(diff.values[3]," vs ",diff.values[2]),
				   paste0(diff.values[3]," vs ",diff.values[1]))
  }
  
  if(TRUE){ # kernel: preparation
## preparation for parallel computing
  if (parallel == TRUE & (CI == TRUE|is.null(bw))) {
    requireNamespace("doParallel")
    ## require(iterators)
    maxcores <- detectCores()
    cores <- min(maxcores, cores)
    pcl<-makeCluster(cores)  
    doParallel::registerDoParallel(pcl)
    cat("Parallel computing with", cores,"cores...\n") 
  }
  
  ## evaluation points
  if(quantile.eval==FALSE){
	X.eval <- seq(min(data[,X]), max(data[,X]), length.out = neval)
  }
  if(quantile.eval==TRUE){
	X.eval <- quantile(data[,X], probs=seq(0,1,length.out = neval))
  }
}
 
}
   
if(TRUE){# Bandwidth Selection: Cross-Validation
  if(is.null(bw) == TRUE){
    CV <- 1
    if (length(grid) == 1) {
      rangeX <- max(data[, X]) - min(data[, X])
      grid <- exp(seq(log(rangeX/50), log(rangeX), length.out = grid))
    }
    cv.out <- crossvalidate.new(data = data, X.eval = X.eval,kfold = kfold,treat.type = treat.type,CV.method=CV.method,
                            Y = Y, D = D, X = X, Z = Z,
                            FE = FE, cl = cl,
                            weights = weights,
                            seed = seed,
                            bw.adaptive = bw.adaptive,
                            grid = grid, metric = metric,
                            parallel = parallel)
    bw <- cv.out$opt.bw
  }  
  else {
    CV <- 0
  }
}
  
if(predict==TRUE){ #A function that returns predicted value of Y under kernel model specifications.
	
	## when there exists fixed effects, using demeaned Y to make predictions
  	data.demean <- data
	if(is.null(FE)==FALSE){
		data.touse <- as.matrix(data[,Y])
		data.fe <- as.matrix(data[,FE])
		if(is.null(weights)==FALSE){
			w <- as.matrix(data[,weights])
		}else{w <- rep(1,dim(data.touse)[1])}
		fastplm_demean <- fastplm(data = data.touse, FE = data.fe,
							      weight = w, FEcoefs = 0)
		data.demean[,Y] <- fastplm_demean$residuals
	}
  
	if (treat.type=='continuous'){
	# kernel predict
    gen_Ey_kernel <- function(data,Y,D,X,FE,weights,Z=NULL,bw,X.pred){
      n<-dim(data)[1]
	  X.eval2 <- X.pred
	  npred <- length(X.pred)
      suppressWarnings(
      coef <- coefs.new(data=data,D=D,Y=Y,Z=Z,X=X,bw=bw,FE=NULL,base = base,weights = weights,predict=TRUE,
                        treat.type = treat.type, X.eval = X.eval2,bw.adaptive = bw.adaptive)
      )
	
     key  <- c("Intercept","D")
     key <- c(key,Z)
     Knn <- as.matrix(coef[,key])

     Ey_all <- matrix(NA,nrow=npred,ncol=0)
     data_predict_start <- matrix(rep(1,npred))
      
     for(target.D in D.sample){
      data_predict <- data_predict_start
      data_predict <- cbind(data_predict,rep(target.D,npred))
      for (char in Z){
          data_predict <- cbind(data_predict,mean(data[,char]))
      }
      output<-matrix(lapply(1:length(data_predict), function(i){data_predict[i]*Knn[[i]]}),
                       nrow=nrow(data_predict), ncol=ncol(data_predict))
      output <- rowSums(matrix(unlist(output),npred))
      Ey_all <- cbind(Ey_all,output)
     }
      colnames(Ey_all) <- names(D.sample)
      return(Ey_all)
    }
  }
  
    if (treat.type=='discrete'){
	gen_Ey_kernel <- function(data,Y,D,X,FE,weights,Z=NULL,bw,X.pred){
    n<-dim(data)[1]
    all_treat.new <- sort(unique(data[,D]))
    if(length(all_treat.new)<ntreat){ # in case some kinds of treatments are not in bootstrap samples
      return(matrix(NA,nrow=npred,ncol=0))
    }
	
	#first demean
	if(is.null(FE)==FALSE){
		data.touse <- as.matrix(data[,Y])
		data.fe <- as.matrix(data[,FE])
		if(is.null(weights)==FALSE){
			w <- as.matrix(data[,weights])
		}else{w <- rep(1,dim(data.touse)[1])}
		fastplm_demean <- fastplm(data = data.touse, FE = data.fe,
                           weight = w, FEcoefs = 0)
		data[,Y] <- fastplm_demean$residuals
	}
	
	X.eval2 <- X.pred
	npred <- length(X.pred)
    
    coef <- coefs.new(data=data,D=D,Y=Y,Z=Z,X=X,bw=bw,FE=NULL,base = base,weights = weights,predict=TRUE,
                       treat.type = 'discrete',X.eval = X.eval2,bw.adaptive = bw.adaptive)
	
    #print(coef)
    coef[,paste0("D.",base)] <- rep(0,npred)

    key  <- c("Intercept")
    for(char in all_treat){
      key <- c(key,paste0("D.",char))
    }
    key <- c(key,Z)  
	Knn <- as.matrix(coef[,key]) 
    Ey_all <- matrix(NA,nrow=npred,ncol=0)
    data_predict_start <- matrix(rep(1,npred))
    
    for(target_treat in all_treat){
      data_predict <- data_predict_start
      for(char in all_treat){
        if(char==target_treat){
          data_predict <- cbind(data_predict,rep(1,npred))
        }
        else{
          data_predict <- cbind(data_predict,rep(0,npred))
        }
      }
      for (char in Z){
        data_predict <- cbind(data_predict,mean(data[,char]))
      }
      output<-matrix(lapply(1:length(data_predict), function(i){data_predict[i]*Knn[[i]]}),
                     nrow=nrow(data_predict), ncol=ncol(data_predict))
      output <- rowSums(matrix(unlist(output),npred))
	  
      Ey_all <- cbind(Ey_all,output)
    }
    colnames(Ey_all) <- all_treat
    return(Ey_all)
  }
}
  
	## predict
	output_kernel <- gen_Ey_kernel(data = data.demean,Y = Y,D = D,Z = Z,X = X,FE=NULL,
								   weights = weights,bw = bw,X.pred=X.eval)
	est_predict_kernel <- list()
  }

if(is.null(diff.values)==FALSE){ #A function that gives the difference between treatment effects at two different values of the moderator.
	if(treat.type=='discrete'){
		gen_difference <- function(data,diff.values,coef){
			n<-dim(data)[1]
			all_treat.new <- sort(unique(data[,D]))
			if(length(all_treat.new)<ntreat){ # in case some kinds of treatments are not in bootstrap samples
				return(matrix(NA,nrow=neval,ncol=0))
			}
			
			est.TE<-function(x){ ##estimate the treatment effects for an observation
				Xnew<-abs(X.eval-x)
				d1<-min(Xnew)     
				label1<-which.min(Xnew)
				Xnew[label1]<-Inf
				d2<-min(Xnew)     
				label2<-which.min(Xnew)
				if(d1==0){
					func <-  (coef[label1,paste0("D.",other_treat)])
				}  
				else if(d2==0){
					func <-  (coef[label2,paste0("D.",other_treat)])
				} 
				else{ ## weighted average
					func1 <- coef[label1,paste0("D.",other_treat)]+1*(x-X.eval[label1])* coef[label1,paste0("DX.",other_treat)]
					func2 <- coef[label2,paste0("D.",other_treat)]+1*(x-X.eval[label2])* coef[label2,paste0("DX.",other_treat)]
					func <- ((func1 * d2 + func2 * d1)/(d1 + d2))
				}
				names(func) <- paste0("TE.",other_treat)
				return(func)
			}
			
			difference1 <- est.TE(diff.values[2])-est.TE(diff.values[1])
			difference2 <- est.TE(diff.values[3])-est.TE(diff.values[2])
			difference3 <- est.TE(diff.values[3])-est.TE(diff.values[1])
			difference <- rbind(difference1,difference2,difference3)
			rownames(difference) <- difference.name
			colnames(difference) <- paste0("Delta.TE.",other_treat)
			return(difference)
	}
		
	}
		
	if(treat.type=='continuous'){
		gen_difference <- function(data,diff.values,coef){
			n<-dim(data)[1]

			est.TE<-function(x){ ##estimate the treatment effects for an observation
				Xnew<-abs(X.eval-x)
				d1<-min(Xnew)     
				label1<-which.min(Xnew)
				Xnew[label1]<-Inf
				d2<-min(Xnew)     
				label2<-which.min(Xnew)
				if(d1==0){
					func <-  coef[label1,"D"]
				}  
				else if(d2==0){
					func <-  coef[label2,"D"]
				} 
				else{ ## weighted average
					func1 <- coef[label1,"D"]+1*(x-X.eval[label1])* coef[label1,"Dx"]
					func2 <- coef[label2,"D"]+1*(x-X.eval[label2])* coef[label2,"Dx"]
					func <- ((func1 * d2 + func2 * d1)/(d1 + d2))
				}
				names(func) <-"TE"
				
				return(func)
			}
			
			difference1 <- est.TE(diff.values[2])-est.TE(diff.values[1])
			difference2 <- est.TE(diff.values[3])-est.TE(diff.values[2])
			difference3 <- est.TE(diff.values[3])-est.TE(diff.values[1])
			difference <- rbind(difference1,difference2,difference3)
			colnames(difference) <- "Delta.TE"
			rownames(difference) <- difference.name
			return(difference)	
			}
	
  }
}

if(TRUE){## Estimates
    if (CI == FALSE) {
	est.matrix <- NULL
    if(treat.type=='discrete'){
    coef.all <- coefs.new(data = data, bw = bw, Y = Y, X = X, D = D, Z = Z, 
                      treat.type = 'discrete', base=base,FE = FE, X.eval = X.eval,
                      bw.adaptive=bw.adaptive,weights = weights)
	coef <- coef.all[,c(3:(2+nothertreat))]
    
    est_discrete <- list()
    for(char in other_treat) {
      if(nothertreat>1){
        est <- data.frame(cbind("X" = X.eval, "ME" = coef[,paste0("D.",char)]),
                                "Treatment"=rep(other_treat.origin[char],length(X.eval)))
      }
      else{
        est <- data.frame(cbind("X" = X.eval, "ME" = coef),"Treatment"=rep(other_treat.origin[char],length(X.eval)))
      }
      est_discrete[[other_treat.origin[char]]] <- est
    }
    est <-est_discrete 
	
	if(predict==TRUE){
		for(char in all_treat){
			est_predict_kernel[[all_treat.origin[char]]] <- cbind.data.frame(X = X.eval, EY = output_kernel[,char],
																			 Treat=rep(all_treat.origin[char],neval))
	  }	
	}
	
	if(is.null(diff.values)==FALSE){
		difference.all <- gen_difference(data=data,diff.values=diff.values,coef=coef.all)
		
		
		diff.table.list <- list()
			for(char in other_treat){
				diff.table <- round(as.matrix(difference.all[,paste0("Delta.TE.",char)]),3)
				colnames(diff.table) <- c("Difference")
				rownames(diff.table) <- difference.name
				diff.table.list[[other_treat.origin[char]]] <- diff.table
			}
			diff.table <- diff.table.list	
	}
  }
    
    if(treat.type=='continuous'){
      coef.all <- coefs.new(data = data, bw = bw, Y = Y, X = X, D = D, Z = Z, treat.type = 'continuous',
                       FE = FE, bw.adaptive=bw.adaptive, X.eval = X.eval, weights = weights)
	  coef <- coef.all[,3]
      est<-data.frame(cbind("X" = X.eval, "ME" = coef))
	if(predict==TRUE){
		for(char in all_treat){
			est_predict_kernel[[char]] <- cbind.data.frame(X = X.eval, EY = output_kernel[,char],Treat=rep(char,neval))
	  }
	}
	if(is.null(diff.values)==FALSE){
		difference.all <- gen_difference(data=data,diff.values=diff.values,coef=coef.all)
		diff.table <- round(as.matrix(difference.all[,"Delta.TE"]),3)
		colnames(diff.table) <- c("Difference")
		rownames(diff.table) <- difference.name
	}
    }
  } 
  else { ## bootstrap
    
    if(treat.type=='discrete'){
      coef.all <- coefs.new(data = data, bw = bw, Y = Y, X = X, D = D, Z = Z, treat.type = 'discrete',base=base,
                       FE = FE, bw.adaptive=bw.adaptive, X.eval = X.eval, weights = weights)
	  coef <- coef.all[,c(3:(2+nothertreat))]
    }
    
    if(treat.type=='continuous'){
      coef.all <- coefs.new(data = data, bw = bw, Y = Y, X = X, D = D, Z = Z, treat.type = 'continuous',
                       FE = FE, bw.adaptive=bw.adaptive, X.eval = X.eval, weights = weights)
	  coef <- coef.all[,3]
    }
	
	if(is.null(diff.values)==FALSE){
		difference.all <- gen_difference(data=data,diff.values=diff.values,coef=coef.all)	
	}
	
	
    ## boostrap SEs
    if (is.null(cl)==FALSE) { ## find clusters
      clusters<-unique(data[,cl])
      id.list<-split(1:n,data[,cl])
    }
	
    oneboot.new <- function() {
      if (is.null(cl)==TRUE) {
        smp<-sample(1:n,n,replace=TRUE)
      } else { ## block bootstrap
        cluster.boot<-sample(clusters,length(clusters),replace=TRUE)
        smp<-unlist(id.list[match(cluster.boot,clusters)])
      }   
      s<-data[smp,]
      
      if(treat.type=='discrete'){
        if(length(unique(s[,D]))<length(unique(data[,D]))){
          return(matrix(NA,length(X.eval),0))
        }
      }
      
	  all.colnames <- c()
	  
      if(treat.type=='discrete'){
        out.all <- coefs.new(data = s, bw = bw, Y = Y, X = X, D = D, Z = Z, treat.type = 'discrete', base=base,
                         FE = FE, bw.adaptive=bw.adaptive, X.eval = X.eval, weights = weights)
		out <- out.all[,c(3:(2+nothertreat))]
		all.colnames <- c(all.colnames,paste0("D.",other_treat))
      }
	  
      if(treat.type=='continuous'){
        out.all <- coefs.new(data = s, bw = bw, Y = Y, X = X, D = D, Z = Z, treat.type = 'continuous',
                         FE = FE, bw.adaptive=bw.adaptive, X.eval = X.eval, weights = weights)
		out <- out.all[,3]
		all.colnames <- c(all.colnames,"out")
      }
	  

	  
	  if(is.null(diff.values)==FALSE){
		if(treat.type=='continuous'){
			temp.difference <- gen_difference(data=s,diff.values=diff.values,coef=out.all)
			output3 <- matrix(NA,neval,1)
			output3[1:3,1] <- temp.difference[1:3,"Delta.TE"]
			colnames(output3) <- 'coef.inter'
			rownames(output3) <- NULL
			out <- cbind(out,output3)
			all.colnames <- c(all.colnames,"coef.inter")
		}
		if(treat.type=='discrete'){
			temp.difference <- gen_difference(data=s,diff.values=diff.values,coef=out.all)
			output3 <- matrix(NA,neval,0)
			colnames.output3 <- c()
			for(char in other_treat){
				newcoef <- matrix(NA,neval,1)
				newcoef[1:3,1] <- temp.difference[1:3,paste0("Delta.TE.",char)]
				output3 <- cbind(output3,newcoef)
				colnames.output3 <- c(colnames.output3,paste0('coef.inter.',char))
			}
			colnames(output3) <- colnames.output3
			out <- cbind(out,output3)
			all.colnames <- c(all.colnames,colnames.output3)
		}
	  }
	  
	  if(predict==TRUE){
		ss<-data.demean[smp,]
		pred_out_kernel=gen_Ey_kernel(data = ss,Y=Y,D=D,X=X,FE=NULL,weights=weights,Z=Z,bw=bw,X.pred=X.eval)
		colnames(pred_out_kernel) <- paste0("pred.kernel.",all_treat)
		pred_out <- pred_out_kernel
		out <- cbind(out,pred_out)
		all.colnames <- c(all.colnames,paste0("pred.kernel.",all_treat))
	  }
	  
	  
	  out <- as.matrix(out)
	  colnames(out) <- all.colnames
	  
	  
	  
      return(out)
    }
    
    cat("Bootstrapping ...\n")
    if (parallel==TRUE) {
      suppressWarnings(
        bootout <- foreach (i=1:nboots, .combine=cbind,
                            .export=c("oneboot.new","coefs.new"),
                            .inorder=FALSE) %dopar% {oneboot.new()}
      ) 
      cat("\r")
    } 
    else {
        bootout<-matrix(NA,length(X.eval),0)
        for (i in 1:nboots) {
          tempdata <- oneboot.new()
          if(is.null(tempdata)==F){
            bootout<- cbind(bootout,tempdata)
          }
          if (i%%50==0) cat(i) else cat(".")
        }
      cat("\r")
      }

	

    if(treat.type=='continuous'){
		## summary
		CI.lvl <- c((1-conf.level)/2, (1-(1-conf.level)/2))
		marg.con <- matrix(NA,nrow=neval,ncol=0)
		bootout_predict_kernel <- list()
		coef.inter <- matrix(NA,nrow=3,ncol=0)
		for (char in all_treat){
				bootout_predict_kernel[[char]] <- matrix(NA,neval,0)
		}
			
		if(predict==TRUE){
			if(is.null(diff.values)==FALSE){
				kcol <- 2+ntreat
			}else{kcol <- 1+ntreat}
		}else{
			if(is.null(diff.values)==FALSE){
				kcol <- 2
			}else{kcol <- 1}
		}
		
		trueboots <- nboots/kcol
		for(k in 0:(trueboots-1)){
			start <- kcol*k+1
			end <- kcol*(k+1)
			bootout_seg <- as.matrix(bootout[,start:end])
			
			# ME
			if(dim(bootout_seg)[2]>1){
				marg.con <- cbind(marg.con,bootout_seg[,'out'])
			} else{marg.con <- cbind(marg.con,bootout_seg)}
			
			# Predict
			if(predict==TRUE){
				for(char in all_treat){
					bootout_predict_kernel[[char]] <- cbind(bootout_predict_kernel[[char]],bootout_seg[,paste0("pred.kernel.",char)])
				}
			}
			
			if(is.null(diff.values)==FALSE){
					coef.inter <- cbind(coef.inter,bootout_seg[1:3,'coef.inter'])
			}
		}
		
		if(predict==TRUE){
			for(char in all_treat){
				ci_kernel<- t(apply(bootout_predict_kernel[[char]], 1, quantile, CI.lvl))
				sd_kernel <- apply(bootout_predict_kernel[[char]],1,sd)
	
				est_predict_kernel[[char]] <- cbind.data.frame(X = X.eval, EY = output_kernel[,char],
                                         SE=sd_kernel ,Treatment=rep(char,neval),
                                         CI_lower=ci_kernel[,1], CI_upper=ci_kernel[,2])
								 
			}
		}
		
		if(is.null(diff.values)==FALSE){
		
			diff.table.start <- matrix(0,nrow=0,ncol=6)
			colnames(diff.table.start) <- c("Difference","se","Z-Score","P-value","CI-lower(95%)","CI-upper(95%)")
			for(j in 1:3){
				difference <- difference.all[j,"Delta.TE"]
				difference.sd <- sqrt(var(coef.inter[j,]))
				difference.z <- difference/difference.sd
				difference.pvalue2sided=2*pnorm(-abs(difference.z))
				difference_interval<- quantile(coef.inter[j,],c(0.025,0.975))
				difference.lbound <- difference_interval[1]
				difference.ubound <- difference_interval[2]
				diff.table <- round(c(difference,difference.sd,difference.z,difference.pvalue2sided,difference.lbound,difference.ubound),3)
				diff.table.start <- rbind(diff.table.start,diff.table)
			}
			rownames(diff.table.start) <- difference.name
			diff.table <- diff.table.start

			}
		

	
		ci<-t(apply(marg.con, 1, quantile, CI.lvl))
		est<-data.frame(cbind("X" = X.eval, "ME" = coef,
                          "SE"=apply(marg.con,1,sd),
                          "CI_lower"=ci[,1], "CI_upper"=ci[,2]))
		est.matrix <- cov(t(marg.con))				  
	
    }
    
    if(treat.type=='discrete'){
      CI.lvl <- c((1-conf.level)/2, (1-(1-conf.level)/2))
      marg.list <- list()
	  pred.list <- list()
	  coef.inter.list <- list()
	  est.matrix <- list()
	  for(char in other_treat){
		marg.list[[char]] <- matrix(NA,nrow=neval,ncol=0)
		coef.inter.list[[char]] <- matrix(NA,nrow=3,ncol=0)
	  }
	  for(char in all_treat){
			pred.list[[char]] <- matrix(NA,nrow=neval,ncol=0)
	  }
		
	  if(predict==TRUE){
			if(is.null(diff.values)==FALSE){
				kcol <- 2*length(other_treat)+length(all_treat)
			} else{kcol <- length(other_treat)+length(all_treat)}
	  }
	  if(predict==FALSE){
			if(is.null(diff.values)==FALSE){
				kcol <- 2*length(other_treat)
			} else{kcol <- length(other_treat)}
	  }
      trueboot <- dim(bootout)[2]/kcol
	  for(k in 0:(trueboot-1)){
			start <- kcol*k+1
			end <- kcol*(k+1)
			bootout_seg <- as.matrix(bootout[,start:end])
			for(char in other_treat){
				tempmarg <- marg.list[[char]]
				if(dim(bootout_seg)[2]>1){
				tempmarg <- cbind(tempmarg,bootout_seg[,paste0("D.",char)])
				} else{tempmarg <- cbind(tempmarg,bootout_seg)}
				marg.list[[char]] <- tempmarg 
			}
			
			if(predict==TRUE){
				for(char in all_treat){
					temppred <- pred.list[[char]]
					temppred <- cbind(temppred,bootout_seg[,paste0("pred.kernel.",char)])
					pred.list[[char]] <- temppred
				}
			}
			
			if(is.null(diff.values)==FALSE){
				for(char in other_treat){
					temp.coef.inter <- coef.inter.list[[char]]
					temp.coef.inter <- cbind(temp.coef.inter,bootout_seg[1:3,paste0('coef.inter.',char)])
					coef.inter.list[[char]] <- temp.coef.inter 
				}
			}
	    }
		
		bootout_predict_kernel <- pred.list
		
		if(predict==TRUE){
			for(char in all_treat){
				ci_kernel<- t(apply(bootout_predict_kernel[[char]], 1, quantile, CI.lvl))
				sd_kernel <- apply(bootout_predict_kernel[[char]],1,sd)
	
				est_predict_kernel[[all_treat.origin[char]]] <- cbind.data.frame(X = X.eval, EY = output_kernel[,char],
                                         SE=sd_kernel ,Treatment=rep(all_treat.origin[char],neval),
                                         CI_lower=ci_kernel[,1], CI_upper=ci_kernel[,2])
								 
			}
		}
		
		est_discrete <- list()
		bootout_coef <- marg.list
		for(char in other_treat) {
			ci <- t(apply(bootout_coef[[char]], 1, quantile, CI.lvl))
			if(nothertreat>1){
				est <- data.frame(cbind("X" = X.eval, "ME" = coef[,paste0("D.",char)],
                              "SE"=apply(bootout_coef[[char]],1,sd),
                              "CI_lower"=ci[,1], "CI_upper"=ci[,2]),"Treatment"=rep(other_treat.origin[char],neval))
			}
			else{
				est <- data.frame(cbind("X" = X.eval, "ME" = coef,
                                  "SE"=apply(bootout_coef[[char]],1,sd),
                                  "CI_lower"=ci[,1], "CI_upper"=ci[,2]),"Treatment"=rep(other_treat.origin[char],neval))			  
			}
        est_discrete[[other_treat.origin[char]]] <- est
		est.matrix[[other_treat.origin[char]]] <- cov(t(marg.list[[char]]))
      }
	  est <- est_discrete
	  
	  if(is.null(diff.values)==FALSE){
		diff.table.list <- list()
		for(char in other_treat){
			
			diff.table.start <- matrix(0,nrow=0,ncol=6)
			colnames(diff.table.start) <- c("Difference","se","Z-Score","P-value","CI-lower(95%)","CI-upper(95%)")
			
			for(j in 1:3){
				difference <- difference.all[j,paste0("Delta.TE.",char)]
				difference.sd <- sqrt(var(coef.inter.list[[char]][j,]))
				difference.z <- difference/difference.sd
				difference.pvalue2sided <- 2*pnorm(-abs(difference.z))
				difference_interval<- quantile(coef.inter.list[[char]][j,],c(0.025,0.975))
				difference.lbound <- difference_interval[1]
				difference.ubound <- difference_interval[2]
				diff.table <- round(c(difference,difference.sd,difference.z,difference.pvalue2sided,difference.lbound,difference.ubound),3)
				diff.table.start <- rbind(diff.table.start,diff.table)
			}
			rownames(diff.table.start) <- difference.name
			diff.table.list[[other_treat.origin[char]]] <- diff.table.start

		}
			diff.table <- diff.table.list
	} 
}
}
  
  if (parallel == TRUE & (CI == TRUE|CV==1)) {
    suppressWarnings(stopCluster(pcl))
    cat("\n") 
  }

  
  
  
  }
 
if(TRUE){ #density or histogram
  if (is.null(weights) == FALSE) {
     dataweights <- rep(1,n)
     } 
  else {
       dataweights <- data[,weights]
     }   
   
   if (treat.type=='discrete') { ## discrete D
     # density
     if (is.null(weights)==TRUE) {
       de <- density(data[,X])
     } else {
       suppressWarnings(de <- density(data[,X],weights=dataweights))
     }
     
     treat_den <- list()
     for (char in all_treat) {
       de.name <- paste0("den.",char)
       if (is.null(weights)==TRUE) {
         de.tr <- density(data[data[,D]==char,X])
       } 
       else {
         suppressWarnings(de.tr <- density(data[data[,D]==char,X],
                                           weights=dataweights[data[,D]==char]))
       }
       treat_den[[all_treat.origin[char]]] <- de.tr
     }
     
     # histogram
     if (is.null(weights)==TRUE) {
       hist.out<-hist(data[,X],breaks=80,plot=FALSE)
     } else {
       suppressWarnings(hist.out<-hist(data[,X],weights,
                                       breaks=80,plot=FALSE))
     } 
     n.hist<-length(hist.out$mids)
     
     # count the number of treated
     treat_hist <- list()
     for (char in all_treat) {
       count1<-rep(0,n.hist)
       treat_index<-which(data[,D]==char)
       for (i in 1:n.hist) {
         count1[i]<-sum(data[treat_index,X]>=hist.out$breaks[i] &
                          data[treat_index,X]<hist.out$breaks[(i+1)])
       }
       count1[n.hist]<-sum(data[treat_index,X]>=hist.out$breaks[n.hist] &
                             data[treat_index,X]<=hist.out$breaks[n.hist+1])
       
       treat_hist[[all_treat.origin[char]]] <- count1
     }    
   }  
   
   if (treat.type=='continuous') { ## continuous D
     if (is.null(weights)==TRUE) {
       de <- density(data[,X])
     } else {
       suppressWarnings(de <- density(data[,X],weights=dataweights))
     }
     if (is.null(weights)==TRUE) {
       hist.out<-hist(data[,X],breaks=80,plot=FALSE)
     } else {
       suppressWarnings(hist.out<-hist(data[,X],weights,
                                       breaks=80,plot=FALSE))
     }
     de.co <- de.tr <- NULL 
     count1 <- NULL
   }
  }
  
if(TRUE){ #Storage
  if(treat.type=='discrete'){
    output<-list(
      type = "kernel",
      bw = bw,
      est = est,
	  vcov.matrix = est.matrix,
      treat.type = treat.type, # binary treatment
      treatlevels = all_treat.origin,
	  order = order,
      base=base.origin,
      Xlabel = Xlabel,
      Dlabel = Dlabel,
      Ylabel = Ylabel,
      de = de,
      de.tr = treat_den, # density
      hist.out = hist.out,
      count.tr = treat_hist,
	  CI=CI,
	  ttest.diffs = diff.table
    )
  }
  
  if(treat.type=="continuous"){
    output<-list(
      type = "kernel",
      bw = bw,
      est = est,  
	  vcov.matrix = est.matrix,
      treat.type = treat.type,
      treatlevels= NULL,
	  order = NULL,
      base= NULL,
      Xlabel = Xlabel,
      Dlabel = Dlabel,
      Ylabel = Ylabel,
      de = de, # density
      de.tr = de.tr,
      hist.out = hist.out,
      count.tr = NULL,
	  CI=CI,
	  ttest.diffs = diff.table
    )
  }
  
  if (CV == 1) {
    output <- c(output, list(CV.out = cv.out$CV.out))
  }
  class(output) <- "interflex"
  }
      

#plot
  if (figure==TRUE & pool==FALSE){
  suppressMessages(
  graph <- plot.interflex(x = output, CI = CI, xlab = xlab, ylab = ylab,
                      Ylabel = Ylabel, Dlabel = Dlabel, Xlabel = Xlabel, order = order,
                      subtitles = subtitles,diff.values = diff.values.plot,
                      show.subtitles = show.subtitles,
                      main = main, xlim = xlim, ylim = ylim, Xdistr = Xdistr,interval = interval,color=color,
                      file = file, theme.bw = theme.bw, show.grid = show.grid,ncols = ncols,pool=pool,jitter=jitter,
                      cex.main = cex.main,cex.sub = cex.sub, cex.axis = cex.axis, cex.lab = cex.lab)
  )
  output <- c(output, list(graph = graph))
  class(output) <- "interflex"
  }
  
  if (figure==TRUE & pool==TRUE & treat.type=='discrete'){

  graph <- plot.interflex(x = output, CI = CI, xlab = xlab, ylab = ylab,
                      Ylabel = Ylabel, Dlabel = Dlabel, Xlabel = Xlabel, subtitles = subtitles,
                      show.subtitles = show.subtitles,diff.values = diff.values.plot,
                      main = main, xlim = xlim, ylim = ylim, Xdistr = Xdistr,interval = interval,color=color,
                      file = file, theme.bw = theme.bw, show.grid = show.grid,ncols = ncols,pool=pool,jitter=jitter,
                      cex.main = cex.main, cex.axis = cex.axis, cex.lab = cex.lab)

  output <- c(output, list(graph = graph))
  class(output) <- "interflex"
  }	
  
  if(predict==TRUE){
	
   if(treat.type=='discrete'){
	 labelname <- NULL
	 D.ref <- NULL
   }
  
   if(treat.type=='continuous'){
		all_treat.origin <- names(D.sample)
	}
  output <- c(output, list(est.predict = est_predict_kernel,labelname = labelname,all.treat = all_treat.origin,
						   X=X,Y=Y,D=D,D.ref=D.ref,predict=TRUE))
  class(output) <- "interflex"
  }
  
  if(predict==FALSE){
	output <- c(output, list(predict=FALSE))
	class(output) <- "interflex"
  }
  return (output)
  }



