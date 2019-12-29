inter.linear<-function(data,
                            Y, # outcome
                            D, # treatment indicator
                            X, # moderator
                            treat.type = NULL, #discrete or continuous
                            base = NULL, # base group when treatments are discrete
                            Z = NULL, # covariates
                            FE = NULL, # fixed effects
                            weights = NULL, # weighting variable
                            full.moderate = TRUE, # whether use fully moderated model
                            na.rm = FALSE,
                            Xunif = FALSE,
							CI = TRUE,
                            vartype = "robust", # variance type
                            ##  "homoscedastic" (default); "robust"; "cluster", "pcse", "bootstrap"
							nboots = 200,
							parallel = TRUE,
							cores = 4,
                            cl = NULL, # variable to be clustered on
                            time = NULL, # time variable for pcse
                            pairwise = TRUE, # pcse option
                            predict = FALSE,
							D.ref = NULL,
							figure = TRUE,
                            order = NULL,
                            subtitles = NULL,
                            show.subtitles = NULL,
                            Xdistr = "histogram", # c("density","histogram","none")
                            main = NULL,
                            Ylabel = NULL,
                            Dlabel = NULL,
                            Xlabel = NULL,
                            xlab = NULL,
                            ylab = NULL,
                            xlim = NULL,
                            ylim = NULL,
                            theme.bw = FALSE,
                            show.grid = TRUE,
                            cex.main = NULL,
                            cex.sub = NULL,
                            cex.lab = NULL,
                            cex.axis = NULL,                       
                            interval = NULL,
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
  
  
  if(TRUE){ #INPUT CHECK
  
  ## in case data is in tibble format
  if (is.data.frame(data) == FALSE) {
    data <- as.data.frame(data)
  }
  n<-dim(data)[1]
  ## check input
  if (is.character(Y) == FALSE) {
    stop("\"Y\" is not a string.")
  } else {
    Y <- Y[1]
  }
  if (is.character(D) == FALSE) {
    stop("\"D\" is not a string.")
  } else {
    D <- D[1]
  }
  
  if (is.character(X) == FALSE) {
    stop("\"X\" is not a string.")
  } else {
    X <- X[1]
  }
  
  if (is.null(Z) == FALSE) {
    for (i in 1:length(Z)) {
      if (is.character(Z[i]) == FALSE) {
        stop("Some element in Z is not a string.")
      }
    }
  }
  if (is.null(FE) == FALSE) {
    requireNamespace("lfe")
    for (i in 1:length(FE)) {
      if (is.character(FE[i]) == FALSE) {
        stop("Some element in FE is not a string.")
      }
    }
    if (vartype == "pcse") {
      vartype <- "cluster"
      warning("Fixed-effect models do not allow panel corrected standard errors; changed to clustered standard errors.")
    }
  }
  
  if (is.null(weights) == FALSE) {
    if (is.character(weights) == FALSE) {
      stop("\"weights\" is not a string.")
    } else {
      dataweights <- data[,weights]
    }   
  }
  
  if (is.logical(na.rm) == FALSE & is.numeric(na.rm)==FALSE) {
    stop("\"na.rm\" is not a logical flag.")
  }
  
  if (is.null(vartype)==TRUE) {
    vartype <- "homoscedastic"
  }
  if (!vartype %in% c("homoscedastic","robust","cluster","pcse","bootstrap")){
    stop("\"vartype\" must be one of the following: \"homoscedastic\",\"robust\",\"cluster\",\"pcse\",\"bootstrap\".")
  } else if (vartype == "cluster") {
    if (is.null(cl)==TRUE) {
      stop("\"cl\" not specified; set cl = \"varname\".")
    }
  } else if (vartype == "pcse") {
    if (is.null(cl)==TRUE | is.null(time)==TRUE) {
      stop("\"cl\" or \"time\" not specified; set cl = \"varname1\", time = \"varname2\":.")
    }
  }
  if (is.null(cl)==FALSE) {
    if (is.character(cl) == FALSE) {
      stop("cl is not a string.")
    } else {
      cl <- cl[1]
      if (vartype != "pcse" & vartype != "bootstrap") {
        vartype <- "cluster"
      }
    }
  }
  if (is.null(time)==FALSE) {
    if (is.character(time) == FALSE) {
      stop("time is not a string.")
    } else {
      time <- time[1]
    }
  }
  if (is.logical(pairwise) == FALSE & is.numeric(pairwise)==FALSE) {
    stop("pairwise is not a logical flag.")
  } 
  if (is.null(Ylabel)==TRUE) {
    Ylabel <- Y
  } else {
    if (is.character(Ylabel) == FALSE) {
      stop("Ylabel is not a string.")
    } else {
      Ylabel <- Ylabel[1]
    }   
  } 
  if (is.null(Dlabel)==TRUE) {
    Dlabel <- D   
  } else {
    if (is.character(Dlabel) == FALSE) {
      stop("Dlabel is not a string.")
    } else {
      Dlabel <- Dlabel[1]
    }   
  }
  if (is.null(Xlabel)==TRUE) {
    Xlabel <- X   
  } else {
    if (is.character(Xlabel) == FALSE) {
      stop("Xlabel is not a string.")
    } else {
      Xlabel <- Xlabel[1]
    }   
  }
  if (is.null(main)==FALSE) {
    main <- main[1]
  }
  if (!Xdistr %in% c("hist","histogram","density","none")){
    stop("\"Xdistr\" must be \"histogram\", \"density\", or \"none\".")
  }
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
  
  if (is.null(interval)==FALSE) {
	if (is.numeric(interval)==FALSE) {
      stop("Some element in interval is not numeric.")
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
  
  ## check missing values
  vars <- c(Y, D, X, Z, FE, cl, weights)
  if (na.rm == TRUE) {        
    data <- na.omit(data[,vars])
  } else {
    if (sum(is.na(data[,vars]))>0) {
      stop("Missing values. Try option na.rm = TRUE\n")
    }
  }
  
  if (is.logical(full.moderate) == FALSE & is.numeric(full.moderate)==FALSE) {
    stop("full.moderate is not a logical flag.")
  }else{
    full <- full.moderate
  } 
  
  ## font size
  if (is.null(cex.lab)==FALSE) {
    if (is.numeric(cex.lab)==FALSE) {
      stop("\"cex.lab\" is not numeric.")
    }
  }
  if (is.null(cex.main)==FALSE) {
    if (is.numeric(cex.main)==FALSE) {
      stop("\"cex.main\" is not numeric.")
    }
  }
  if (is.null(cex.sub)==FALSE) {
    if (is.numeric(cex.sub)==FALSE) {
      stop("\"cex.main\" is not numeric.")
    }
  }
  if (is.null(cex.axis)==FALSE) {
    if (is.numeric(cex.axis)==FALSE) {
      stop("\"cex.axis\" is not numeric.")
    }    
  }  

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
  if (is.logical(parallel) == FALSE & is.numeric(parallel)==FALSE) {
    stop("\"paralell\" is not a logical flag.")
  }
  if (is.numeric(cores)==FALSE) {
    stop("\"cores\" is not a positive integer.")
  } else {
    cores <- cores[1]
    if (cores%%1!= 0 | cores<=0) {
      stop("\"cores\" is not a positive integer.")
    }
  }
  
  if(is.null(diff.values)==FALSE){
	if(is.numeric(diff.values)==FALSE){
		stop("\"diff.values\" is not numeric.")
	}
	if(length(diff.values)!=2){
		stop("\"diff.values\" must be of length 2.")
	}
  }
}



  if(TRUE){ #treat.type
  
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
  
  #if treat is discrete
  if (treat.type=='discrete') {
    if(length(unique(data[,D]))>5) {
      warning("More than 5 kinds of treatments.")
    }
    
    if(length(unique(data[,D]))>20) {
      stop("Too many kinds of treatments.")
    }
    data[,D] <- as.character(data[,D])
    
    
    if(is.null(base)==TRUE) {
      base=sort(unique(data[,D]))[1]
      f=sprintf("Base group has not been specified, choose treat = %s as base group. \n",base)
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
    
	## in case there are special characters in D
    all_treat=sort(unique(data[,D]))
	names(all_treat) <- paste("Group",c(1:length(all_treat)),sep  = '.')
    other_treat <- all_treat[which(all_treat!=base)]
    other_treat <- sort(other_treat)
    
    if(is.null(order)==F){
	  order <- as.character(order)
      if(length(order)!=length(other_treat)){
        stop("\"order\" should contain all kinds of treatments except for the base group.")
      }
      if(length(order)!=length(unique(order))){
        stop("\"order\" should not contain repeated values.")
      }
      
      if(sum(!is.element(order,other_treat))!=0 | sum(!is.element(other_treat,order))!=0){
        stop("\"order\" should contain all kinds of treatments except for the base group.")
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
        stop("\"subtitle\" had a wrong length.")
      }
    }
	
    if (is.logical(show.subtitle) == F & is.numeric(show.subtitle)==F & is.null(show.subtitle)==F) {
      stop("\"show.subtitle\" is not a logical flag.")
    } 
  }
  
  if (treat.type=='continuous') {
    if (is.numeric(data[,D])==F) {
      stop("\"D\" is not a numeric variable")
    }
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
  }
  

  
  # number of columns
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
		ncols <- length(unique(data[,D]))-1
	}
	if(treat.type=="continuous"){
		ncols <- 1
	}
  }
}

  
  if(TRUE){ #PREPROCESS
  #factor cov
  panel_vars <- c()
  for(a in Z){
	if(is.factor(data[,a])==TRUE){
		panel_vars <- c(panel_vars,a)
	}	
  }
  if(length(panel_vars)>0){
	fnames <- paste("factor(", panel_vars, ")", sep = "")
	contr.list <- list(contr.sum, contr.sum)
	names(contr.list) <- fnames
	panel_form <- as.formula(paste("~", paste(fnames, collapse = " + ")))
	suppressWarnings(
	panel_mat <- model.matrix(panel_form, data = data,
                          contrasts.arg = contr.list)[, -1]
	)
	dummy_colnames <- c()
	for(i in 1:dim(panel_mat)[2]){
		dummy_colnames <- c(dummy_colnames,paste0("Dummy.Covariate.",i))
	}
	colnames(panel_mat) <- dummy_colnames
	data <- cbind(data,panel_mat)
	Z <- Z[!Z %in% panel_vars]
	Z <- c(Z,dummy_colnames)
  }
  
  
  ## fully moderated model
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
  ## fully moderated model END
  
  ## transform moderator into uniform
  if (Xunif == TRUE) {
    x <- data[, X]
    data[, X] <- rank(x, ties.method = "average")/length(x)*100
    Xlabel <- paste(Xlabel,"(Percentile)")
  }
  ## transform moderator into uniform END
  
  ## change fixed effect variable to factor
  if (is.null(FE)==FALSE) {
    if (length(FE) == 1) {
      data[, FE] <- as.numeric(as.factor(data[, FE]))
    } else {
      data[, FE] <- sapply(data[,FE],function(vec){as.numeric(as.factor(vec))})
    }
  }
  ## change fixed effect variable to factor END
  
  if(predict==FALSE){
	est.pred.linear <- NULL
  }
  
  
  }

  cat("Use a linear interaction model\n")
  
  
  if(TRUE){ #LINEAR Model FORMULA
	##  make a vector of the marginal effect of D on Y as X changes
	X.lvls<-as.numeric(quantile(data[,X], probs=seq(0,1,0.01)))
  
	## linear model formula
	if (treat.type=="discrete") {
		data[,D] <- as.factor(data[,D])
		data[,D] <- relevel(data[,D], ref=base)
	}
	mod.f<-paste0(Y,"~",D,"+",X,"+",D,"*",X)
	if (is.null(Z)==FALSE) {
		mod.f <- paste0(mod.f, "+", paste0(Z,collapse="+"))
	}
	if (is.null(FE)==FALSE) {
		mod.f <- paste0(mod.f, "|",paste0(FE, collapse = "+"))
    if (vartype=="cluster") {
      mod.f <- paste0(mod.f, "| 0 |",paste0(cl,collapse = "+"))
    }
  }
	mod.f <- as.formula(mod.f)
  }


  if(vartype!='bootstrap'){
	## a naive fit
	if (is.null(FE)==TRUE) { #OLS
		if (is.null(weights)==TRUE) {
			mod.naive<-lm(mod.f,data=data)
		} else {
			mod.naive<-lm(mod.f,data=data,weights=dataweights)
		}
	} else { # FE
		if (is.null(weights)==TRUE) {
			mod.naive<-felm(mod.f,data=data)
		} else {
			mod.naive<-felm(mod.f,data=data,weights=dataweights)
		}
	}
  
	## coefficients
	coefs<-summary(mod.naive)$coefficients[,1]
	if(treat.type=='continuous'){
	coef.D<-coefs[D]
	coef.X<-coefs[X]
	coef.DX<-coefs[paste(D,X,sep=":")] #interaction
	}
  
	if(treat.type=='discrete'){
		coef_list <- list()
		coef_inter_list <- list()
		for(char in other_treat){
			coef_list[[char]] <- coefs[paste0(D,char)]
			coef_inter_list[[char]] <- coefs[paste(paste0(D,char),X,sep=":")]
		}
	}
	

	## # variance    
	if (is.null(FE)==TRUE) { #OLS
	if (vartype=="homoscedastic") {
			v<-vcov(mod.naive)
    } else if (vartype=="robust") {
			requireNamespace("sandwich")
			v<-vcov(mod.naive,type="HC1") # White with small sample correction
    } else if (vartype=="cluster") {
      v<-vcovCluster(mod.naive,cluster = data[,cl])
    } else if (vartype=="pcse") {
      requireNamespace("pcse")
      v<-pcse(mod.naive,groupN=data[,cl],groupT=data[,time],pairwise=pairwise)$vcov
    }
  } else { # FE
    if (vartype=="homoscedastic") {
      v<-vcov(mod.naive, type = "iid")
    } else if (vartype=="robust") {
      v<-vcov(mod.naive, type="robust") 
    } else if (vartype=="cluster") {
      v<-vcov(mod.naive, type = "cluster") 
    }
  }
  
  ## get variance
  if(treat.type=='continuous'){
    if (vartype=="pcse") {
      var.D<-v[D,D]
      var.DX<-v[paste(D,X,sep="."),paste(D,X,sep=".")]
      cov.DX<-v[D,paste(D,X,sep=".")]
    } 
    else {
      var.D<-v[D,D]
      var.DX<-v[paste(D,X,sep=":"),paste(D,X,sep=":")]
      cov.DX<-v[D,paste(D,X,sep=":")]
    }
  }
  
  if(treat.type=='discrete'){
	var_list <- list()
	varinter_list <- list()
	cov_list <- list()
    if (vartype=="pcse") {
      for(char in other_treat){
		var_list[[char]] <- v[paste0(D,char),paste0(D,char)]
		varinter_list[[char]] <- v[paste(paste0(D,char),X,sep="."),paste(paste0(D,char),X,sep=".")]
		cov_list[[char]] <- v[paste0(D,char),paste(paste0(D,char),X,sep=".")]
      }
    } 
    else {
      for(char in other_treat){
		var_list[[char]] <- v[paste0(D,char),paste0(D,char)]
		varinter_list[[char]] <- v[paste(paste0(D,char),X,sep=":"),paste(paste0(D,char),X,sep=":")]
		cov_list[[char]] <- v[paste0(D,char),paste(paste0(D,char),X,sep=":")]
      }
    }
  }
  
  	if(is.null(diff.values)==FALSE){
		x.diff <- diff.values[2]-diff.values[1]
		if(treat.type=='continuous'){
			difference <- coef.DX*x.diff
			difference.sd <- abs(x.diff)*sqrt(var.DX)
			difference.z <- difference/difference.sd
			difference.pvalue2sided=2*pnorm(-abs(difference.z))
			difference.lbound <- difference-1.96*difference.sd
			difference.ubound <- difference+1.96*difference.sd
			diff.table <- round(c(difference,difference.sd,difference.z,difference.pvalue2sided,difference.lbound,difference.ubound),3)
			names(diff.table) <- c("Difference","se","Z-Score","P-value","CI-lower(95%)","CI-upper(95%)")
		}
		
		if(treat.type=='discrete'){
			diff.table.list <- list()
			for(char in other_treat){
				difference <- coef_inter_list[[char]]*x.diff
				difference.sd <- abs(x.diff)*sqrt(varinter_list[[char]])
				difference.z <- difference/difference.sd
				difference.pvalue2sided=2*pnorm(-abs(difference.z))
				difference.lbound <- difference-1.96*difference.sd
				difference.ubound <- difference+1.96*difference.sd
				diff.table <- round(c(difference,difference.sd,difference.z,difference.pvalue2sided,difference.lbound,difference.ubound),3)
				names(diff.table) <- c("Difference","se","Z-Score","P-value","CI-lower(95%)","CI-upper(95%)")
				diff.table.list[[other_treat.origin[char]]] <- diff.table
			}
			diff.table <- diff.table.list
		}
	} else{diff.table <- NULL}
	
	
  
  
  if (treat.type=='continuous'){
    marg<-coef.D + coef.DX*X.lvls
  
    ## the variance is var(B1_D) + X^2*var(B_3) + 2*inst*cov(D, X)
    se<-sqrt(var.D +  X.lvls^2*var.DX + 2*X.lvls*cov.DX)
    df<-mod.naive$df.residual
    crit<-abs(qt(.025, df=df)) # critical values
  
    ##make 95% confidence bands. 
    lb<-marg-crit*se
    ub<-marg+crit*se
    est.lin<-data.frame(X.lvls, marg, lb, ub)
  }
  
  if (treat.type=='discrete'){
    df <- mod.naive$df.residual
    crit<-abs(qt(.025, df=df))
    est.lin<-list()
    
    for(char in other_treat) {
	  marg <- coef_list[[char]] + coef_inter_list[[char]]*X.lvls
      se <- sqrt(var_list[[char]] + X.lvls^2*varinter_list[[char]]+2*X.lvls*cov_list[[char]])
	  #print(se)
	  #print(df)
      lb <- marg-crit*se
      ub <- marg+crit*se
      tempest <- data.frame(X.lvls,marg,se,lb,ub)
      tempest[,'Treatment']<- rep(other_treat.origin[char],dim(tempest)[1])
      est.lin[[other_treat.origin[char]]] <- tempest
    }
  }


	if(predict==TRUE){
		data.old <- data
		FE.old <- FE
		if(is.null(FE)==FALSE){
								data.touse <- as.matrix(data[,Y])
								data.fe <- as.matrix(data[,FE])
								if(is.null(weights)==FALSE){
									w <- as.matrix(data[,weights])
								}else{w <- rep(1,dim(data.touse)[1])}
										fastplm_demean <- fastplm(data = data.touse, FE = data.fe,
															  weight = w, FEcoefs = 0)
								data[,Y] <- fastplm_demean$residuals
								FE <- NULL
		}
	
		data_predict <- data[,c(Y,X)]
		predict_colname <- c(Y,X)
		formula <- paste0(Y,"~",X)
		if(treat.type=='discrete'){
			for(char in other_treat){
				data_predict <- cbind(data_predict,as.numeric(data[,D]==char))
				predict_colname <- c(predict_colname,paste0("D.",char))
				formula <- paste0(formula,"+",paste0("D.",char))
			}
			for(char in other_treat){
				data_predict <- cbind(data_predict,as.numeric(data[,D]==char)*data[,X])
				predict_colname <- c(predict_colname,paste0("DX.",char))
				formula <- paste0(formula,"+",paste0("DX.",char))
			}
		}
		if(treat.type=='continuous'){
				data_predict <- cbind(data_predict,data[,D])
				predict_colname <- c(predict_colname,D)
				formula <- paste0(formula,"+",D)
				data_predict <- cbind(data_predict,data[,D]*data[,X])
				predict_colname <- c(predict_colname,"DX")
				formula <- paste0(formula,"+","DX")
		}
		if(is.null(Z)==FALSE){
			for(char in Z){
				data_predict <- cbind(data_predict,data[,char])
				predict_colname <- c(predict_colname,char)
				formula <- paste0(formula,"+",char)
			}
		}
		predict_colname0 <- predict_colname
		#if (is.null(FE)==FALSE) {
		#	data_predict <- cbind(data_predict,data[,FE])
		#	formula <- paste0(formula, "|",paste0(FE, collapse = "+"))
		#	predict_colname <- c(predict_colname,FE)
		#	if (vartype=="cluster") {
		#		formula<- paste0(formula, "| 0 |",paste0(cl,collapse = "+"))
		#		if(!cl%in%FE){
		#			data_predict <- cbind(data_predict,data[,cl])
		#			predict_colname <- c(predict_colname,cl)
		#		}
		#	}
		#}
		
		colnames(data_predict) <- predict_colname
		formula <- as.formula(formula)
		
		if (is.null(FE)==TRUE) { #OLS
			if (is.null(weights)==TRUE) {
				mod.naive2<-lm(formula,data=data_predict)
			} else {
				mod.naive2<-lm(formula,data=data_predict,weights=dataweights)
			}
		} else { # FE
			if (is.null(weights)==TRUE) {
				mod.naive2<-felm(formula,data=data_predict)
			} else {
				mod.naive2<-felm(formula,data=data_predict,weights=dataweights)
			}
		}
		
		
		
		if (is.null(FE)==TRUE) { #OLS
		  if (vartype=="homoscedastic") {
				v<-vcov(mod.naive2)
		} else if (vartype=="robust") {
			requireNamespace("sandwich")
			v<-vcov(mod.naive2,type="HC1") # White with small sample correction
		} else if (vartype=="cluster") {
			v<-vcovCluster(mod.naive2,cluster = data[,cl])
		} else if (vartype=="pcse") {
			requireNamespace("pcse")
			v<-pcse(mod.naive2,groupN=data[,cl],groupT=data[,time],pairwise=pairwise)$vcov
		}
		} else { # FE
			if (vartype=="homoscedastic") {
				v<-vcov(mod.naive2, type = "iid")
			} else if (vartype=="robust") {
				v<-vcov(mod.naive2, type="robust") 
			} else if (vartype=="cluster") {
				v<-vcov(mod.naive2, type = "cluster") 
		}
		}
		v[which(is.na(v))] <- 0 
		coef2 <- as.matrix(coef(mod.naive2))
		coef2[which(is.na(coef2))] <- 0
		X.pred <- seq(min(data[,X]),max(data[,X]),length.out=101)
		est.predict.linear <- list()
		
		for(char in all_treat){
			newdata <- X.pred
			if(treat.type=='discrete'){
				for(char0 in other_treat){
					newdata <- cbind(newdata,as.numeric(char0==char))
				}
				for(char0 in other_treat){
					newdata <- cbind(newdata,as.numeric(char0==char)*X.pred)
				}
			}
			if(treat.type=='continuous'){
				newdata <- cbind(newdata,D.sample[char])
				newdata <- cbind(newdata,D.sample[char]*X.pred)
			}
			
			for(char0 in Z){
				newdata <- cbind(newdata,mean(data[,char0]))
			}
			newdata <- cbind(1,newdata)
			
			newdata <- as.data.frame(newdata)

			colnames(newdata) <- predict_colname0
			m.mat <- model.matrix(mod.naive2$terms,data=newdata)
			if(is.null(FE)==FALSE){
				m.mat <- m.mat[,-1]
			}
			
			fit <- as.vector(m.mat %*% coef2)
			se.fit <- sqrt(diag(m.mat%*%v%*%t(m.mat)))
			df2 <- mod.naive2$df
			
			CI.lvl <- c((1-0.95)/2, (1-(1-0.95)/2))
			crit<-abs(qt(CI.lvl[1], df=df2))
			lb<-fit-crit*se.fit
			ub<-fit+crit*se.fit
			if(treat.type=='discrete'){
				est.predict.linear[[all_treat.origin[char]]] <- cbind.data.frame(X = X.pred, EY = fit, 
                                           SE = se.fit ,Treatment=all_treat.origin[char],
                                           CI_lower=lb, CI_upper=ub
                                           )
			}
			if(treat.type=='continuous'){
				est.predict.linear[[char]] <- cbind.data.frame(X = X.pred, EY = fit, 
                                           SE = se.fit ,Treatment=char,
                                           CI_lower=lb, CI_upper=ub
                                           )
			}
		}
		data <- data.old
		FE <- FE.old
}

}


  if(vartype=='bootstrap'){
  
	## a function that will return marginal effects and binning coeficients
	gen_marg <- function(data,Df=FALSE){
		if(is.null(weights)==FALSE){
			dataweights <- data[,weights]
		}
		
		if (is.null(FE)==TRUE) { #OLS
			if (is.null(weights)==TRUE) {
				mod.naive<-lm(mod.f,data=data)
			} 	else {
				mod.naive<-lm(mod.f,data=data,weights=dataweights)
			}
		} else { # FE
			if (is.null(weights)==TRUE) {
				mod.naive<-felm(mod.f,data=data)
			} else {
				mod.naive<-felm(mod.f,data=data,weights=dataweights)
			}
		}
		
		## coefficients
		coefs<-summary(mod.naive)$coefficients[,1]
	
		if(treat.type=='continuous'){
			coef.D<-coefs[D]
			coef.X<-coefs[X]
			coef.DX<-coefs[paste(D,X,sep=":")] #interaction
			marg<-coef.D + coef.DX*X.lvls
			marg[which(is.na(marg))] <- 0
		}	
  
		if(treat.type=='discrete'){
			coef_list <- list()
			coef_inter_list <- list()
			marg_list <- list()
			for(char in other_treat){
				coef_list[[char]] <- coefs[paste0(D,char)]
				coef_inter_list[[char]] <- coefs[paste(paste0(D,char),X,sep=":")]
				temp_marg <- coef_list[[char]] + coef_inter_list[[char]]*X.lvls
				temp_marg[which(is.na(temp_marg))] <- 0
				marg_list[[char]] <- temp_marg
			}
			coef.DX <- coef_inter_list
		}
  
		if(treat.type=='continuous'){
			output <- matrix(NA,nrow=length(X.lvls),ncol=1)
			output[,1] <- marg
			colnames(output) <- c('ME')
		}
		if(treat.type=='discrete'){
			output <- matrix(NA,nrow=length(X.lvls),ncol=length(other_treat))
			k <- 1
			output_colname <- c()
			for(char in other_treat){
				output[,k] <- marg_list[[char]]
				output_colname <- c(output_colname,paste0("ME.",char))
				k <- k + 1
		}
		colnames(output) <- output_colname
	}
 
		if(Df==FALSE){
			return(output)
		}
		if(Df==TRUE) {
			return(list(output=output,df=df,coef.inter=coef.DX))
		}
	}
	
	## a function that will return linear prediction
	if(predict==TRUE){
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
		
		if(treat.type=='continuous'){
			gen_Ey_linear <- function(data,Y,D,X,FE,weights,Z=NULL,X.pred){
								n<-dim(data)[1]
								data_touse <- data
								if(is.null(weights)==F){
									weight_touse <- as.matrix(data[,weights])
								}
								else{
									weight_touse <- as.matrix(rep(1,n))
								}
      
								data1 <- as.data.frame(data_touse[,c(Y,X,D)])
								data1 <- cbind(data1,data1[,X]*data1[,D])
								colnames <- c("Y","X","D","D.X")
								formula <- "Y ~ X + D + D.X"
								for(char in Z){
									data1 <- cbind(data1,data_touse[,char])
									colnames <- c(colnames,char)
									formula <- paste(formula,char,sep = "+")
								}
								
								if (is.null(FE)==FALSE) {
										formula <- paste0(formula, "|",paste0(FE, collapse = "+"))
										data1 <- cbind(data1,data_touse[,FE])
										colnames <- c(colnames,FE)
								}
								colnames(data1) <- colnames
								formula <- as.formula(formula)
								# coef
								if(is.null(FE)==TRUE){
										naive_linearfit <- lm(formula=formula,data=data1,weights = weight_touse)
									}
								if(is.null(FE)==FALSE){
										suppressWarnings(naive_linearfit <- felm(formula=formula,data=data1,weights = weight_touse))
									}
								
								naive.coef <- naive_linearfit$coefficients
								naive.coef[which(is.na(naive.coef))] <- 0
								naive.coef <- as.matrix(naive.coef)

								# predict
								x_predict <- X.pred
								npred <- length(X.pred)
								data_predict_start <- matrix(rep(1,npred))
								data_predict_start <- cbind(data_predict_start,x_predict)
								Ey_all <- matrix(NA,nrow=npred,ncol=0)
								for(target.D in D.sample){
									data_predict <- data_predict_start
									data_predict <- cbind(data_predict,rep(target.D,npred))
									data_predict <- cbind(data_predict,target.D*x_predict)
								for(char in Z){
									data_predict <- cbind(data_predict,mean(data1[,char]))
								}
								output <- as.double(data_predict%*%naive.coef)
								Ey_all <- cbind(Ey_all,output)
								}
								colnames(Ey_all) <- names(D.sample)
								return(Ey_all)
			}
		}
		
		if(treat.type=='discrete'){
			ntreat <- length(all_treat)
			gen_Ey_linear <- function(data,Y,D,X,FE,weights,Z=NULL,X.pred){
									  n<-dim(data)[1]
									  all_treat.new <- sort(unique(data[,D]))
									  if(length(all_treat.new)<ntreat){ # in case some kinds of treatments are not in bootstrap samples
											return(matrix(NA,nrow=npred,ncol=0))
									  }
									  
									if(is.null(FE)==FALSE){
									data.touse <- as.matrix(data[,Y])
									data.fe <- as.matrix(data[,FE])
									if(is.null(weights)==FALSE){
										w <- as.matrix(data[,weights])
									}else{w <- rep(1,dim(data.touse)[1])}
									fastplm_demean <- fastplm(data = data.touse, FE = data.fe,
															  weight = w, FEcoefs = 0)
									data[,Y] <- fastplm_demean$residuals
									FE <- NULL
									}
									  
									  
									  
									  
									  
									  data_touse <- data
									  if(is.null(weights)==F){
											weight_touse <- as.matrix(data[,weights])
									  }
									  else{
											weight_touse <- as.matrix(rep(1,n))
									  }
           
									  data1 <- as.data.frame(data_touse[,c(Y,X)])
									  colnames <- c("Y","X")
								      formula <- "Y ~ X"
									  for(char in Z){
										data1 <- cbind(data1,data_touse[,char])
										colnames <- c(colnames,char)
										formula <- paste(formula,char,sep = "+")
									  }
									  for(char in other_treat){
										data1 <- cbind(data1,as.numeric(data_touse[,D]==char))
										colnames <- c(colnames,paste0("D.",char))
										formula <- paste(formula,paste0("D.",char),sep = "+")
									  }
									  for(char in other_treat){
										data1 <- cbind(data1,data_touse[,X]*as.numeric(data_touse[,D]==char))
										colnames <- c(colnames,paste0("DX.",char))
										formula <- paste(formula,paste0("DX.",char),sep = "+")
									  }
									  
									  if (is.null(FE)==FALSE) {
											formula <- paste0(formula, "|",paste0(FE, collapse = "+"))
											data1 <- cbind(data1,data_touse[,FE])
											colnames <- c(colnames,FE)
										}
									  colnames(data1) <- colnames
									  formula <- as.formula(formula)
									  #coef
									  if(is.null(FE)==TRUE){
										naive_linearfit <- lm(formula=formula,data=data1,weights = weight_touse)
									  }
									  if(is.null(FE)==FALSE){
										suppressWarnings(naive_linearfit <- felm(formula=formula,data=data1,weights = weight_touse))
									  }
									  naive.coef <- naive_linearfit$coefficients
								      naive.coef[which(is.na(naive.coef))] <- 0
									  naive.coef <- as.matrix(naive.coef)

									  x_predict <- X.pred
									  npred <- length(X.pred)
									  data_predict_start <- matrix(rep(1,npred))
									  data_predict_start <- cbind(data_predict_start,x_predict)								  
									  for (char in Z){
										data_predict_start <- cbind(data_predict_start,mean(data1[,char]))
									  }
      
									  Ey_all <- matrix(NA,nrow=npred,ncol=0)
      
									  for(target_treat in all_treat){
										data_predict <- data_predict_start
										for (char in other_treat){
											if(char==target_treat){
												data_predict <- cbind(data_predict,rep(1,npred))
											}
											else{
												data_predict <- cbind(data_predict,rep(0,npred))
											}
										}
									  for (char in other_treat){
										if(char==target_treat){
											data_predict <- cbind(data_predict,rep(1,npred)*x_predict)
										}
										else{
											data_predict <- cbind(data_predict,rep(0,npred))
										}
									  }
									  output <- as.double(data_predict%*%naive.coef)
								      Ey_all <- cbind(Ey_all,output)
									  }
									  colnames(Ey_all) <- all_treat
									  return(Ey_all)
									  }
		}
	
		X.pred <- seq(min(data[,X]),max(data[,X]),length.out=101)
		EY_output <- gen_Ey_linear(data=data.demean,Y=Y,D=D,X=X,FE=NULL,weights=weights,Z=Z,X.pred=X.pred)
	}
	
	output.list <- gen_marg(data=data,Df=TRUE)
	output <- output.list$output
	df.X <- output.list$df
	coef.inter.T <- output.list$coef.inter
	
	
	if (is.null(cl)==FALSE) { ## find clusters
      clusters<-unique(data[,cl])
      id.list<-split(1:n,data[,cl])
    }
 
	one.boot <- function(){
      if (is.null(cl)==TRUE) {
        smp<-sample(1:n,n,replace=TRUE)
      } else { ## block bootstrap
        cluster.boot<-sample(clusters,length(clusters),replace=TRUE)
        smp<-unlist(id.list[match(cluster.boot,clusters)])
      }   
      s<-data[smp,]
	  
	  if(treat.type=='discrete'){
        if(length(unique(s[,D]))<length(unique(data[,D]))){
          return(matrix(NA,length(X.lvls),0))
        }
      }
	  
	  runflag <- try(output.all <- gen_marg(data=s,Df=TRUE),silent=T)
	  if(class(runflag)=='try-error'){
	  #print('wrong')
	  return(matrix(NA,length(X.lvls),0))
	  }
	  output <- output.all$output
	  
	  if(predict==TRUE){
			  ss <- data.demean[smp,]
			  runflag <- try(output2 <- gen_Ey_linear(data=ss,Y=Y,D=D,X=X,FE=NULL,weights=weights,Z=Z,X.pred=X.pred),silent=T)
			  if(class(runflag)=='try-error'){
				return(matrix(NA,length(X.lvls),0))
			  }
			  colnames(output2) <- paste0("pred.",colnames(output2))
			  output <- cbind(output,output2)
	  }
	  
	  if(is.null(diff.values)==FALSE){
		output.coef <- output.all$coef.inter
		if(treat.type=='continuous'){
			output3 <- matrix(output.coef,length(X.lvls),1)
			colnames(output3) <- 'coef.inter'
			output <- cbind(output,output3)
		}
		if(treat.type=='discrete'){
		
			output3 <- matrix(NA,length(X.lvls),0)
			colnames.output3 <- c()
			for(char in other_treat){
				newcoef <- matrix(output.coef[[char]],length(X.lvls),1)
				output3 <- cbind(output3,newcoef)
				colnames.output3 <- c(colnames.output3,paste0('coef.inter.',char))
			}
			colnames(output3) <- colnames.output3
			output3 <- as.data.frame(output3)
			output <- cbind(output,output3)
		}
	  }
	  return(output)
	}
	
	
	if (parallel==TRUE) {
	   requireNamespace("doParallel")
		## require(iterators)
		maxcores <- detectCores()
		cores <- min(maxcores, cores)
		pcl<-makeCluster(cores)  
		doParallel::registerDoParallel(pcl)
		cat("Parallel computing with", cores,"cores...\n") 

      suppressWarnings(
        bootout <- foreach (i=1:nboots, .combine=cbind,
                            .export=c("one.boot"),.packages=c("lfe"),
                            .inorder=FALSE) %dopar% {one.boot()}
      ) 
	  suppressWarnings(stopCluster(pcl))
      cat("\r")
    } 
    else {
        bootout<-matrix(NA,length(X.lvls),0)
        for (i in 1:nboots) {
          tempdata <- one.boot()
          if(is.null(tempdata)==F){
            bootout<- cbind(bootout,tempdata)
          }
          if (i%%50==0) cat(i) else cat(".")
        }
      cat("\r")
	}
	
	
	if(dim(bootout)[2]==0){
		stop("Bootstraping is not stable, try another vartype.")
	}
	
	if(treat.type=='discrete'){
		marg.list <- list()
		pred.list <- list()
		coef.inter.list <- list()
		for(char in other_treat){
			marg.list[[char]] <- matrix(NA,nrow=length(X.lvls),ncol=0)
			coef.inter.list[[char]] <- c()
		}
		for(char in all_treat){
			pred.list[[char]] <- matrix(NA,nrow=length(X.lvls),ncol=0)
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
				tempmarg <- cbind(tempmarg,bootout_seg[,paste0("ME.",char)])
				} else{tempmarg <- cbind(tempmarg,bootout_seg)}
				marg.list[[char]] <- tempmarg 
			}
			if(predict==TRUE){
				for(char in all_treat){
					temppred <- pred.list[[char]]
					temppred <- cbind(temppred,bootout_seg[,paste0("pred.",char)])
					pred.list[[char]] <- temppred
				}
			}
			if(is.null(diff.values)==FALSE){
				for(char in other_treat){
					temp.coef.inter <- coef.inter.list[[char]]
					temp.coef.inter <- c(temp.coef.inter,mean(bootout_seg[,paste0('coef.inter.',char)]))
					coef.inter.list[[char]] <- temp.coef.inter 
				}
			}
			
		}
	}
	
	if(treat.type=='continuous'){
		marg.con <- matrix(NA,nrow=length(X.lvls),ncol=0)
		pred.list <- list()
		coef.inter <- c()
		for(char in all_treat){
			pred.list[[char]] <- matrix(NA,nrow=length(X.lvls),ncol=0)
		}
		
		if(predict==TRUE){
			if(is.null(diff.values)==FALSE){
				kcol <- 2+length(all_treat)
			} else{kcol <- 1+length(all_treat)}	
		}
		if(predict==FALSE){
			if(is.null(diff.values)==FALSE){
				kcol <- 2
			} else{kcol <- 1}	
		}

		trueboot <- dim(bootout)[2]/kcol

		for(k in 0:(trueboot-1)){
			start <- kcol*k+1
			end <- kcol*(k+1)
			bootout_seg <- as.matrix(bootout[,start:end])
			
			if(dim(bootout_seg)[2]>1){
				marg.con <- cbind(marg.con,bootout_seg[,'ME'])
			} else{marg.con <- cbind(marg.con,bootout_seg)}
			
			if(predict==TRUE){
				for(char in all_treat){
					temppred <- pred.list[[char]]
					temppred <- cbind(temppred,bootout_seg[,paste0("pred.",char)])
					pred.list[[char]] <- temppred
				}
			}
			
			if(is.null(diff.values)==FALSE){
					coef.inter <- c(coef.inter,mean(bootout_seg[,'coef.inter']))
			}
		}
	}
	
	if(is.null(diff.values)==FALSE){
		x.diff <- diff.values[2]-diff.values[1]
		if(treat.type=='continuous'){
			difference <- coef.inter.T*x.diff
			difference.sd <- abs(x.diff)*sqrt(var(coef.inter))
			difference.z <- difference/difference.sd
			difference.pvalue2sided=2*pnorm(-abs(difference.z))
			difference_interval<- quantile(coef.inter,c(0.025,0.975))*x.diff
			difference.lbound <- difference_interval[1]
			difference.ubound <- difference_interval[2]
			diff.table <- round(c(difference,difference.sd,difference.z,difference.pvalue2sided,difference.lbound,difference.ubound),3)
			names(diff.table) <- c("Difference","se","Z-Score","P-value","CI-lower(95%)","CI-upper(95%)")
		}
		
		if(treat.type=='discrete'){
			diff.table.list <- list()
			for(char in other_treat){
				difference <- coef.inter.T[[char]]*x.diff
				difference.sd <- abs(x.diff)*sqrt(var(coef.inter.list[[char]]))
				difference.z <- difference/difference.sd
				difference.pvalue2sided <- 2*pnorm(-abs(difference.z))
				
				difference_interval<- quantile(coef.inter.list[[char]],c(0.025,0.975))*x.diff
				difference.lbound <- difference_interval[1]
				difference.ubound <- difference_interval[2]
				diff.table <- round(c(difference,difference.sd,difference.z,difference.pvalue2sided,difference.lbound,difference.ubound),3)
				names(diff.table) <- c("Difference","se","Z-Score","P-value","CI-lower(95%)","CI-upper(95%)")
				diff.table.list[[other_treat.origin[char]]] <- diff.table
			}
			diff.table <- diff.table.list
		}
	} else{diff.table <- NULL}
	


	CI.lvl <- c((1-0.95)/2, (1-(1-0.95)/2))
	if(treat.type=='discrete'){
		est.lin <- list()
		for(char in other_treat){
			marg <- output[,paste0("ME.",char)]
			marg.ci <- t(apply(marg.list[[char]], 1, quantile, CI.lvl,na.rm=TRUE))
			se <- apply(marg.list[[char]],1,sd,na.rm=TRUE)
			lb <- marg.ci[,1]
			ub <- marg.ci[,2]
			tempest <- data.frame(X.lvls,marg,se,lb,ub)
			tempest[,'Treatment']<- rep(other_treat.origin[char],dim(tempest)[1])
			est.lin[[other_treat.origin[char]]] <- tempest
		}
	}
 
	if(treat.type=='continuous'){
		marg <- output[,"ME"]
		marg.ci <- t(apply(marg.con, 1, quantile, CI.lvl))
		se <- apply(marg.con,1,sd)
		lb <- marg.ci[,1]
		ub <- marg.ci[,2]
		tempest <- data.frame(X.lvls,marg,se,lb,ub)
		est.lin <- tempest
	}
	
	
	if(predict==TRUE){
	est.predict.linear <- list()
		for(char in all_treat){
			fit <- EY_output[,char]
			pred.ci <- t(apply(pred.list[[char]], 1, quantile, CI.lvl,na.rm=TRUE))
			se.fit <- apply(pred.list[[char]],1,sd,na.rm=TRUE)
			lb <- pred.ci[,1]
			ub <- pred.ci[,2]
			
			if(treat.type=='discrete'){
				
				est.predict.linear[[all_treat.origin[char]]] <- cbind.data.frame(X = X.pred, EY = fit, 
                                           SE = se.fit ,Treatment=rep(all_treat.origin[char],101),
                                           CI_lower=lb, CI_upper=ub
                                           )
										   
			}						   
			
			if(treat.type=='continuous'){
				est.predict.linear[[char]] <- cbind.data.frame(X = X.pred, EY = fit, 
                                           SE = se.fit ,Treatment=rep(char,101),
                                           CI_lower=lb, CI_upper=ub
                                           )
			}
			}
		}
}


  ## L kurtosis
  requireNamespace("Lmoments")
  Lkurtosis <- Lmoments(data[,X],returnobject=TRUE)$ratios[4]
  
  if(TRUE){ # density or histogram
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
      #print("hhh")
      suppressWarnings(hist.out<-hist(data[,X],weights,
                                      breaks=80,plot=FALSE))
    }
    de.co <- de.tr <- NULL 
    count1 <- NULL
  } 
}


  if(TRUE){ # Tests
	tests <- list(
			treat.type = treat.type,
			X.Lkurtosis = sprintf(Lkurtosis,fmt = '%#.3f')
			)
  }
  
  if(TRUE){ # Storage
	if(treat.type=='discrete'){
	output<-list(
    type = "linear",
    est.lin = est.lin,
    treat.type = treat.type,
    treatlevels = all_treat.origin,
	order = order,
    base = base.origin,
    Xlabel = Xlabel,
    Dlabel = Dlabel,
    Ylabel = Ylabel,
    de = de,
    de.tr = treat_den, # density
    hist.out = hist.out,
    count.tr = treat_hist,
    tests = tests,
	difference.est = diff.table,
	predict = predict
  )
  }
  
  if(treat.type=="continuous"){
    output<-list(
      type = "linear",
      est.lin = est.lin,
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
      tests = tests,
	  difference.est = diff.table,
	  predict = predict
    )
  }
  
  if(predict==TRUE){
  
	if(treat.type=='discrete'){
		labelname <- NULL
		D.ref <- NULL
	}
	if(treat.type=='continuous'){
		all_treat.origin <- names(D.sample)
	}
  
   output <- c(output,list(est.predict = est.predict.linear,labelname = labelname,all.treat = all_treat.origin,
						   X=X,Y=Y,D=D,D.ref=D.ref))
  }
  
  }
  class(output) <- "interflex"
  

  if (figure==TRUE & pool==F) {
    
    class(output) <- "interflex"
    suppressMessages(
    graph <- plot.interflex(out = output, CI = CI, xlab = xlab, ylab = ylab, color = color, order = order,
                        subtitles = subtitles,
                        show.subtitles = show.subtitles,
                        Ylabel = Ylabel, Dlabel = Dlabel, Xlabel = Xlabel, 
                        main = main, xlim = xlim, ylim = ylim, Xdistr = Xdistr,interval = interval,pool=pool,
                        file = file, theme.bw = theme.bw, show.grid = show.grid, ncols=ncols,
                        cex.main = cex.main,cex.sub = cex.sub, cex.axis = cex.axis, cex.lab = cex.lab)  
    )
    output<-c(output,list(graph=graph))
    class(output) <- "interflex"
  }
  
  if(figure==TRUE & pool==TRUE & treat.type=='discrete'){
    suppressMessages(
    graph <- plot.interflex(out = output, CI = CI, xlab = xlab, ylab = ylab,order=order,subtitles = subtitles,show.subtitles = show.subtitles,
                              Ylabel = Ylabel, Dlabel = Dlabel, Xlabel = Xlabel, 
                              main = main, xlim = xlim, ylim = ylim, Xdistr = Xdistr,interval = interval,pool=pool,color=color,
                              file = file, theme.bw = theme.bw, show.grid = show.grid,
                              cex.main = cex.main, cex.axis = cex.axis, cex.lab = cex.lab,legend.title =legend.title)
    )
    output <- c(output, list(graph = graph))
  }
  class(output) <- "interflex"
  return(output)
  }
  







