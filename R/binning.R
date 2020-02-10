inter.binning<-function(data,
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
                            nbins = 3,  # No. of X bins
                            cutoffs = NULL,
							CI = TRUE,
                            vartype = "robust", # variance type
                            ##  "homoscedastic" (default); "robust"; "cluster", "pcse", "bootstrap"
							nboots = 200,
							parallel = TRUE,
							cores = 4,
                            cl = NULL, # variable to be clustered on
                            time = NULL, # time variable for pcse
                            pairwise = TRUE, # pcse option
                            wald = TRUE,
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
                            bin.labs = TRUE,                        
                            interval = NULL,
                            file = NULL,
                            ncols = NULL,
                            pool = FALSE,
							color = NULL,
                            jitter = FALSE,
							legend.title = NULL
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
  ## FE
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
  ## weights
  if (is.null(weights) == FALSE) {
    if (is.character(weights) == FALSE) {
      stop("\"weights\" is not a string.")
    } else {
      dataweights <- data[,weights]
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
  
  ## nbins
  if (is.null(nbins) == FALSE) {
    if (nbins%%1 != 0) {
      stop("\"nbins\" is not a positive integer.")
    } else {
      nbins <- nbins[1]
    }
    if (nbins < 2) {
      stop("\"nbins\" should be a positive integer larger than 2.")
    }
  }
  
  ## cutoffs
  if (is.null(cutoffs) == FALSE) {
    if (is.numeric(cutoffs) == FALSE) {
      stop("Some element of cutoffs is not numeric.")
    } 
  }

  ## CI
  if (is.logical(CI) == FALSE & is.numeric(CI)==FALSE) {
    stop("\"CI\" is not a logical flag.")
  }
  ## vartype
  if (is.null(vartype)==TRUE) {
    vartype <- "homoscedastic"
  }
  if (!vartype %in% c("homoscedastic","robust","cluster","pcse","bootstrap")){
    stop("\"vartype\" must be one of the following: \"homoscedastic\",\"robust\",\"cluster\",\"pcse\",\"bootstrap\".")
  } 
  if (vartype == "cluster") {
    if (is.null(cl)==TRUE) {
      stop("\"cl\" not specified; set cl = \"varname\".")
    }
  } 
  if (vartype == "pcse") {
    if (is.null(cl)==TRUE | is.null(time)==TRUE) {
      stop("\"cl\" or \"time\" not specified; set cl = \"varname1\", time = \"varname2\".")
    }
  }
  
  #nboots
  if (is.null(nboots) == FALSE) {
    if (is.numeric(nboots)==FALSE) {
      stop("\"nboots\" is not a positive integer.")
    } else {
      nboots <- nboots[1]
      if (nboots%%1 != 0 | nboots < 1) {
        stop("\"nboots\" is not a positive number.")
      }
    } 
  }else{nboots <- 200}
  
  # Parallel
  if (is.logical(parallel) == FALSE & is.numeric(parallel)==FALSE) {
    stop("\"paralell\" is not a logical flag.")
  }
  
  # Cores
  if (is.numeric(cores)==FALSE) {
    stop("\"cores\" is not a positive integer.")
  } else {
    cores <- cores[1]
    if (cores%%1!= 0 | cores<=0) {
      stop("\"cores\" is not a positive integer.")
    }
  }
  
  # cl
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
  
  # time
  if (is.null(time)==FALSE) {
    if (is.character(time) == FALSE) {
      stop("time is not a string.")
    } else {
      time <- time[1]
    }
  }
  
  ## check missing values
  vars <- c(Y, D, X, Z, FE, cl,time, weights)
  if (na.rm == TRUE) {        
    data <- na.omit(data[,vars])
  } else {
    if (sum(is.na(data[,vars]))>0) {
      stop("Missing values. Try option na.rm = TRUE\n")
    }
  }
  
  #pairwise
  if (is.logical(pairwise) == FALSE & is.numeric(pairwise)==FALSE) {
    stop("\"pairwise\" is not a logical flag.")
  }
  
  #wald
  if (is.logical(wald) == FALSE & is.numeric(wald)==FALSE) {
    stop("\"wald\" is not a logical flag.")
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
  
  #main
  if (is.null(main)==FALSE) {
    main <- as.character(main)[1]
  }
  #Ylabel
  if (is.null(Ylabel)==TRUE) {
    Ylabel <- Y
  } else {
    if (is.character(Ylabel) == FALSE) {
      stop("Ylabel is not a string.")
    } else {
      Ylabel <- Ylabel[1]
    }   
  }
  #Dlabel  
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
  
  ## bin.labs
  if (is.logical(bin.labs) == FALSE & is.numeric(bin.labs)==FALSE) {
    stop("\"bin.labs\" is not a logical flag.")
  }
  
  # interval
  if (is.null(interval)==FALSE) {
	if (is.numeric(interval)==FALSE) {
      stop("Some element in interval is not numeric.")
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
  
  ## jitter
  if (is.logical(jitter) == FALSE & is.numeric(jitter)==FALSE) {
    stop("\"jitter\" is not a logical flag.")
  }
  
  ## legend.title
  if (is.null(legend.title)==FALSE) {
    legend.title <- as.character(legend.title)[1]
  }
  
}

if(TRUE){ # treat.type Check
	## treat.type check
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
## treat.type checks END
  
## number of columns in plots
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
  
if(TRUE){ #Preprocess
  n<-dim(data)[1]
  #factor
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
  }
  
##  make a vector of the marginal effect of D on Y as X changes
##X.lvls<-as.numeric(quantile(data[,X], probs=seq(0,1,0.01)))
X.lvls <- seq(min(data[,X]), max(data[,X]), length.out = 50)
 
if(TRUE){ ## linear model formula
	if (treat.type=="discrete"){
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
  
if(TRUE){## binning model cutting points
	## grouping by X
	if (is.null(cutoffs)==TRUE) {
		cuts.X<-quantile(data[,X],probs=seq(0,1,1/nbins))
		while (length(unique(cuts.X))!=nbins+1) {
			nbins<-nbins-1
			cuts.X<-quantile(data[,X],probs=seq(0,1,1/nbins))
		} 
	} else {
		cutoffs <-cutoffs[which(cutoffs>min(data[,X]) & cutoffs < max(data[,X]))]
		cuts.X<- sort(unique(c(min(data[,X]),cutoffs,max(data[,X]))))
	} 
	groupX<-cut(data[,X],breaks=cuts.X, labels = FALSE)
	groupX[which(data[,X]==min(data[,X]))]<-1
	nbins <- length(unique(groupX))
  
	## X labels
	groupX2 <- cut(data[,X],breaks=cuts.X)
	gp.lab = paste(Xlabel, ": ", levels(groupX2), sep="")
	gp.lab[1] <- paste(Xlabel, ": [", substring(levels(groupX2)[1],2), sep = "")

	## mid points
	x0<-rep(NA,nbins)
	for (i in 1:nbins) x0[i]<-median(data[which(groupX==i),X], na.rm=TRUE)
}
  
if(vartype=='bootstrap'){
	## a function that will return marginal effects and binning coeficients
	## (predict==T) generate a function that returns predicted value of Y under binning model specifications.
	## bootstrap
	## some tests
	gen_marg <- function(data,nbins,cuts.X,x0,Df=FALSE){
		if(is.null(weights)==FALSE){
			dataweights <- data[,weights]
		}
		groupX<-cut(data[,X],breaks=cuts.X, labels = FALSE)
		groupX[which(data[,X]==min(data[,X]))]<-1
		
		if(length(unique(groupX))<nbins){
			return(matrix(NA,nrow=length(X.lvls),ncol=1))
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
		}
		df <- NULL
		if (treat.type=='continuous') {
			G<-DG<-GX<-DGX<-matrix(0,n,nbins)
			for (i in 1:nbins) {
				G[which(groupX==i),i]<-1
				DG[,i]<-data[,D]*G[,i]
				GX[,i]<-G[,i]*(data[,X]-x0[i])
				DGX[,i]<-DG[,i]*(data[,X]-x0[i])
			}
  
			Gs<-GXs<-DGs<-DGXs<-c()
			for (i in 1:nbins) {
				Gs<-c(Gs,paste0("G[,",i,"]"))
				GXs<-c(GXs,paste0("GX[,",i,"]"))
				DGs<-c(DGs,paste0("DG[,",i,"]"))
				DGXs<-c(DGXs,paste0("DGX[,",i,"]"))
			}
			
			Xf<-paste0(Y,"~ -1+",paste0(DGs,collapse="+"),"+",paste0(DGXs,collapse="+"),
             "+",paste0(Gs,collapse="+"),"+",paste0(GXs,collapse="+"))
			 
			if (is.null(Z)==FALSE) {
				Xf<- paste0(Xf,"+",paste0(Z,collapse="+"))
			}
			if (is.null(FE)==FALSE) {
				Xf <- paste0(Xf, "|",paste0(FE, collapse = "+"))
			}
			mod.Xf<-as.formula(Xf)    
  
			if (is.null(FE)==TRUE) { #OLS
				if (is.null(weights)==TRUE) {
					mod.X<-lm(mod.Xf,data=data)
				} else {
					mod.X<-lm(mod.Xf,data=data,weights=dataweights)
				}
			} else { # FE
				if (is.null(weights)==TRUE) {
					mod.X<-suppressWarnings(felm(mod.Xf,data=data))
				} else {
					mod.X<-suppressWarnings(felm(mod.Xf,data=data,weights=dataweights))
				}
			}	
			df <- mod.X$df
			Xcoefs<-mod.X$coefficients[1:nbins]
			#Xcoefs[which(is.na(Xcoefs))] <- 0
			output <- matrix(NA,nrow=length(X.lvls),ncol=2)
			output[,1] <- marg
			output[1:length(Xcoefs),2] <- Xcoefs 
			colnames(output) <- c('ME','BinCoef')
			}
	
		if(treat.type=='discrete'){
			G<-GX<-matrix(0,n,nbins)
			for (i in 1:nbins) {
				G[which(groupX==i),i] <- 1
				GX[,i] <- G[,i]*(data[,X]-x0[i])
			}
			for (char in other_treat) {
				DG_name <- paste0("DG.",char)
				DGX_name <- paste0("DGX.",char)
				DG_matrix <- DGX_matrix <- matrix(0,n,nbins)
				for (i in 1:nbins) {
					DG_matrix[,i] <- (data[,D]==char)*G[,i]
					DGX_matrix[,i] <- DG_matrix[,i]*(data[,X]-x0[i])
			}
			assign(DG_name,DG_matrix)
			assign(DGX_name,DGX_matrix)
		}
		Gs<-GXs<-DGs<-DGXs<-c()
		for (i in 1:nbins)  {
			Gs<-c(Gs,paste0("G[,",i,"]"))
			GXs<-c(GXs,paste0("GX[,",i,"]"))
		for (char in other_treat) {
			DGs<-c(DGs,paste0("DG.",char,"[,",i,"]"))
			DGXs<-c(DGXs,paste0("DGX.",char,"[,",i,"]"))
		}
		}
		Xf<-paste0(Y,"~ -1+",paste0(DGs,collapse="+"),"+",paste0(DGXs,collapse="+"),
               "+",paste0(Gs,collapse="+"),"+",paste0(GXs,collapse="+"))
		if (is.null(Z)==FALSE) {
			Xf<- paste0(Xf,"+",paste0(Z,collapse="+"))
		}
		if (is.null(FE)==FALSE) {
			Xf <- paste0(Xf, "|",paste0(FE, collapse = "+"))
		}
		mod.Xf<-as.formula(Xf)    
    
		## fit
		if (is.null(FE)==TRUE) { #OLS
			if (is.null(weights)==TRUE) {
				mod.X<-lm(mod.Xf,data=data)
			} else {
				mod.X<-lm(mod.Xf,data=data,weights=dataweights)
			}
		} else { # FE
			if (is.null(weights)==TRUE) {
				mod.X<-suppressWarnings(felm(mod.Xf,data=data))
			} else {
				mod.X<-suppressWarnings(felm(mod.Xf,data=data,weights=dataweights))
		}
		}
		df <- mod.X$df
		Xcoefs_list <- list()
		for (char in other_treat) {
			DGs <- c()
			for (i in 1:nbins) {
				DGs<-c(DGs,paste0("DG.",char,"[, ",i,"]"))
			}
			if(is.null(FE)==T){
				tempcoef <- mod.X$coefficients[DGs]
			}
			else{
				tempcoef <- mod.X$coefficients[DGs,]
			}
			#tempcoef[which(is.na(tempcoef))] <- 0
			Xcoefs_list[[char]] <- tempcoef
		}
		
		output <- matrix(NA,nrow=length(X.lvls),ncol=length(other_treat)*2)
		k <- 1
		output_colname <- c()
		for(char in other_treat){
			output[,k] <- marg_list[[char]]
			output_colname <- c(output_colname,paste0("ME.",char))
			k <- k + 1
			output[1:nbins,k] <- Xcoefs_list[[char]]
			output_colname <- c(output_colname,paste0("BinCoef.",char))
			k <- k + 1
		}
		colnames(output) <- output_colname
	}

 if(Df==FALSE){
	return(output)
 }
 if(Df==TRUE) {
	return(list(output=output,df=df))
 }
 }
	output.list <- gen_marg(data=data,nbins=nbins,cuts.X=cuts.X,x0=x0,Df=TRUE)
	output <- output.list$output
	df.X <- output.list$df
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
			    gen_Ey <- function(data,Y,D,X,FE,weights,Z=NULL,nbins,cuts.X,x0){
									n<-dim(data)[1]
									data_touse <- data
									if(is.null(weights)==F){
										weight_touse <- as.matrix(data[,weights])
									}
									else{
										weight_touse <- as.matrix(rep(1,n))
									}
									groupX<-cut(data[,X],breaks=cuts.X, labels = FALSE)
									groupX[which(data[,X]==min(data[,X]))]<-1
									if(length(unique(groupX))<nbins){
										return(matrix(NA,nrow=length(X.lvls),ncol=0))
									}
									data1 <- as.data.frame(data_touse[,Y])
									colnames <- c("Y")
									formula <- "Y ~ -1"
									for(char in Z){
										data1 <- cbind(data1,data_touse[,char])
										colnames <- c(colnames,char)
										formula <- paste(formula,char,sep = "+")
									}
									for (i in 1:nbins){
										data1 <- cbind(data1,as.numeric(groupX==i))
										colnames <- c(colnames,paste0("G.",i))
										formula <- paste(formula,paste0("G.",i),sep = "+")
									}
									for (i in 1:nbins){
										data1 <- cbind(data1,as.numeric(groupX==i)*data_touse[,D])
										colnames <- c(colnames,paste0("G.",i,".D"))
										formula <- paste(formula,paste0("G.",i,".D"),sep="+")
									}
									for (i in 1:nbins){
										data1 <- cbind(data1,as.numeric(groupX==i)*(data_touse[,X]-x0[i]))
										colnames <- c(colnames,paste0("G.",i,".X"))
										formula <- paste(formula,paste0("G.",i,".X"),sep = "+")
									}
									for (i in 1:nbins){
										data1 <- cbind(data1,as.numeric(groupX==i)*data_touse[,D]*(data_touse[,X]-x0[i]))
										colnames <- c(colnames,paste0("G.",i,".DX"))
										formula <- paste(formula,paste0("G.",i,".DX"),sep="+")
									}
		
									if (is.null(FE)==FALSE) {
										formula <- paste0(formula, "|",paste0(FE, collapse = "+"))
										data1 <- cbind(data1,data_touse[,FE])
										colnames <- c(colnames,FE)
									}
		
									colnames(data1) <- colnames
									formula <- as.formula(formula)
									if(is.null(FE)==TRUE){
										binning_fit <- lm(formula=formula,data=data1,weights = weight_touse)
									}
									if(is.null(FE)==FALSE){
										suppressWarnings(binning_fit <- felm(formula=formula,data=data1,weights = weight_touse))
									}
	  
									binning_coef <- binning_fit$coefficients
									binning_coef[which(is.na(binning_coef))] <- 0
									binning_coef <- as.matrix(binning_coef)
	  
									npred <- length(X.pred)
									x_predict <- X.pred
									data_predict_start <- matrix(NA,nrow = npred,ncol = 0)
									for (char in Z){
										data_predict_start <- cbind(data_predict_start,mean(data_touse[,char]))
									}
									Ey_all <- matrix(NA,nrow=npred,ncol=0)
									groupX_predict<-cut(x_predict,breaks=cuts.X, labels = FALSE)
									groupX_predict[which(x_predict==min(x_predict))]<-1
									for(target.D in D.sample){
										data_predict <- data_predict_start
										for(i in 1:nbins){
											data_predict <- cbind(data_predict,as.numeric(groupX_predict==i))
										}
										for(i in 1:nbins){
											data_predict <- cbind(data_predict,as.numeric(groupX_predict==i)*target.D)
										}
										for(i in 1:nbins){
											data_predict <- cbind(data_predict,as.numeric(groupX_predict==i)*(x_predict-x0[i]))
										}
										for(i in 1:nbins){
											data_predict <- cbind(data_predict,as.numeric(groupX_predict==i)*target.D*(x_predict-x0[i]))
										}
										output <- as.double(data_predict%*%binning_coef)
										Ey_all <- cbind(Ey_all,output)
									}	
									colnames(Ey_all) <- names(D.sample)
									return(Ey_all)
									}
		}
		if(treat.type=='discrete'){
			gen_Ey <- function(data,Y,D,X,FE,weights,Z=NULL,nbins,cuts.X,x0){
								n <- dim(data)[1]
								

								data_touse <- data
								if(is.null(weights)==F){
									weight_touse <- as.matrix(data[,weights])
								}
								else{
									weight_touse <- as.matrix(rep(1,n))
								}
	  
								groupX<-cut(data[,X],breaks=cuts.X, labels = FALSE)
								groupX[which(data[,X]==min(data[,X]))]<-1
								if(length(unique(groupX))<nbins){
									return(matrix(NA,nrow=length(X.pred),ncol=0))
								}
  
								all_treat.new <- sort(unique(data[,D]))
								if(length(all_treat.new)<length(all_treat)){ # in case some kinds of treatments are not in bootstrap samples
									return(matrix(NA,nrow=length(X.pred),ncol=0))
								}
    
								data1 <- as.data.frame(data_touse[,Y])
								colnames <- c("Y")
								formula <- "Y ~ -1"
    
								for(char in Z){
									data1 <- cbind(data1,data_touse[,char])
									colnames <- c(colnames,char)
									formula <- paste(formula,char,sep = "+")
								}
    
								for (i in 1:nbins){
									data1 <- cbind(data1,as.numeric(groupX==i))
									colnames <- c(colnames,paste0("G.",i))
									formula <- paste(formula,paste0("G.",i),sep = "+")
								}
    
								for (i in 1:nbins){
									for(char in other_treat){
										data1 <- cbind(data1,as.numeric(groupX==i)*as.numeric(data_touse[,D]==char))
										colnames <- c(colnames,paste0("G.",i,".D.",char))
										formula <- paste(formula,paste0("G.",i,".D.",char),sep="+")
									}
								}
    
								for (i in 1:nbins){
									data1 <- cbind(data1,as.numeric(groupX==i)*(data_touse[,X]-x0[i]))
									colnames <- c(colnames,paste0("G.",i,".X"))
									formula <- paste(formula,paste0("G.",i,".X"),sep = "+")
								}
    
								for (i in 1:nbins){
									for(char in other_treat){
										data1 <- cbind(data1,as.numeric(groupX==i)*as.numeric(data_touse[,D]==char)*(data_touse[,X]-x0[i]))
										colnames <- c(colnames,paste0("G.",i,".D.",char,".X"))
										formula <- paste(formula,paste0("G.",i,".D.",char,".X"),sep="+")
									}
								}
	
								if (is.null(FE)==FALSE) {
									formula <- paste0(formula, "|",paste0(FE, collapse = "+"))
									data1 <- cbind(data1,data_touse[,FE])
									colnames <- c(colnames,FE)
								}
								colnames(data1) <- colnames
								formula <- as.formula(formula)
								if(is.null(FE)==TRUE){
									binning_fit <- lm(formula=formula,data=data1,weights = weight_touse)
								}
								if(is.null(FE)==FALSE){
									suppressWarnings(binning_fit <- felm(formula=formula,data=data1,weights = weight_touse))
								}
	
								binning_coef <- binning_fit$coefficients
								binning_coef[which(is.na(binning_coef))] <- 0
								binning_coef <- as.matrix(binning_coef)
    
								x_predict <- X.pred
								npred <- length(X.pred)
								data_predict_start <- matrix(NA,nrow = npred,ncol = 0)
								for (char in Z){
									data_predict_start <- cbind(data_predict_start,mean(data_touse[,char]))
								}
    
								Ey_all <- matrix(NA,nrow=npred,ncol=0)
								groupX_predict<-cut(x_predict,breaks=cuts.X, labels = FALSE)
								groupX_predict[which(x_predict==min(x_predict))] <- 1
    
								for(target_treat in all_treat){
									data_predict <- data_predict_start
									for(i in 1:nbins){
										data_predict <- cbind(data_predict,as.numeric(groupX_predict==i))
									}
									for(i in 1:nbins){
										for(char in other_treat){
											data_predict <- cbind(data_predict,as.numeric(groupX_predict==i)*as.numeric(char==target_treat))
										}
									}
									for(i in 1:nbins){
										data_predict <- cbind(data_predict,as.numeric(groupX_predict==i)*(x_predict-x0[i]))
									}
									for(i in 1:nbins){
										for(char in other_treat){
											data_predict <- cbind(data_predict,as.numeric(groupX_predict==i)*as.numeric(char==target_treat)*(x_predict-x0[i]))
										}
									}
									output <- as.double(data_predict%*%binning_coef)
									Ey_all <- cbind(Ey_all,output)
								}
								colnames(Ey_all) <- all_treat
								return(Ey_all)
								}		
		}
		X.pred <- seq(min(data[,X]),max(data[,X]),length.out=length(X.lvls))
		predict_Ey <- gen_Ey(data=data.demean,Y=Y,D=D,X=X,FE=NULL,weights=weights,Z=Z,nbins=nbins,cuts.X=cuts.X,x0=x0)
	}
	
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
	  
	  runflag <- try(output <- gen_marg(data=s,nbins=nbins,cuts.X=cuts.X,x0=x0),silent=T)
	  
	  if(class(runflag)=='try-error'){
		return(matrix(NA,length(X.lvls),0))
	  }
	  
	  if(predict==TRUE){
		ss<-data.demean[smp,]
		runflag2 <- try(output2 <- gen_Ey(data=ss,Y=Y,D=D,X=X,FE=NULL,Z=Z,weights=weights,
										  nbins=nbins,cuts.X=cuts.X,x0=x0),silent=T)
		if(class(runflag2)=='try-error'){
			return(matrix(NA,length(X.lvls),0))
		}
		if(dim(output2)[2]==0){
			return(matrix(NA,length(X.lvls),0))
		}
		colnames(output2) <- paste0("pred.",colnames(output2))
		output <- cbind(output,output2)
	  }
	  return(output)
	}
	
	if (parallel==TRUE){
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
	stop("Bootstrap standard error is not stable, try another vartype.")
}


 if(treat.type=='discrete'){
	marg.list <- list()
	bincoef.list <- list()
	bin.var.list <- list()
	pred.list <- list()
	for(char in other_treat){
		marg.list[[char]] <- matrix(NA,nrow=length(X.lvls),ncol=0)
		bincoef.list[[char]] <- matrix(NA,nrow=nbins,ncol=0)
	}
	for(char in all_treat){
			pred.list[[char]] <- matrix(NA,nrow=length(X.lvls),ncol=0)
	}
	if(predict==TRUE){
			kcol <- 2*length(other_treat)+length(all_treat)
	}
	if(predict==FALSE){
			kcol <- 2*length(other_treat)
	}
	trueboot <- dim(bootout)[2]/kcol
	for(k in 0:(trueboot-1)){
		start <- kcol*k+1
		end <- kcol*(k+1)
		bootout_seg <- bootout[,start:end]
		
		for(char in other_treat){
			tempmarg <- marg.list[[char]]
			tempmarg <- cbind(tempmarg,bootout_seg[,paste0("ME.",char)])
			marg.list[[char]] <- tempmarg 
			tempbincoef <- bincoef.list[[char]]
			tempbincoef <- cbind(tempbincoef,bootout_seg[,paste0("BinCoef.",char)][1:nbins])
			bincoef.list[[char]] <- tempbincoef
		}
		
		if(predict==TRUE){
			for(char in all_treat){
					temppred <- pred.list[[char]]
					temppred <- cbind(temppred,bootout_seg[,paste0("pred.",char)])
					pred.list[[char]] <- temppred
				}
		}
	}	
	for(char in other_treat){
		bin.var.list[[char]] <- var(t(bincoef.list[[char]]),na.rm=TRUE)
	}	
 }
 
 
 if(treat.type=='continuous'){
	marg.con <- matrix(NA,nrow=length(X.lvls),ncol=0)
	coef.con <- matrix(NA,nrow=nbins,ncol=0)
	pred.list <- list()
	for(char in all_treat){
			pred.list[[char]] <- matrix(NA,nrow=length(X.lvls),ncol=0)
	}
	if(predict==TRUE){
			kcol <- 2+length(all_treat)
	}
	if(predict==FALSE){
			kcol <- 2
	}
	
	trueboot <- dim(bootout)[2]/kcol
		
	for(k in 0:(trueboot-1)){
		start <- kcol*k+1
		end <- kcol*(k+1)
		bootout_seg <- bootout[,start:end]
		marg.con <- cbind(marg.con,bootout_seg[,"ME"])
		coef.con <- cbind(coef.con,bootout_seg[,"BinCoef"][1:nbins])
	
		if(predict==TRUE){
			for(char in all_treat){
				temppred <- pred.list[[char]]
				temppred <- cbind(temppred,bootout_seg[,paste0("pred.",char)])
				pred.list[[char]] <- temppred
			}
		}
	}
		bin.var <- var(t(coef.con),na.rm=TRUE)
 }
 
 CI.lvl <- c((1-0.95)/2, (1-(1-0.95)/2))
 if(treat.type=='discrete'){
	est.lin <- list()
	est.bin <- list()
	est.matrix <- list()
	for(char in other_treat){
		marg <- output[,paste0("ME.",char)]
		marg.ci <- t(apply(marg.list[[char]], 1, quantile, CI.lvl,na.rm=TRUE))
		
		se <- apply(marg.list[[char]],1,sd,na.rm=TRUE)
		lb <- marg.ci[,1]
		ub <- marg.ci[,2]
		tempest <- data.frame(X.lvls,marg,se,lb,ub)
		tempest[,'Treatment']<- rep(other_treat.origin[char],dim(tempest)[1])
		est.lin[[other_treat.origin[char]]] <- tempest
		est.matrix[[other_treat.origin[char]]] <- cov(t(marg.list[[char]]))
		
		
		bin.coef <- output[,paste0("BinCoef.",char)][1:nbins]
		bin.coef.ci <- t(apply(bincoef.list[[char]], 1, quantile, CI.lvl,na.rm=TRUE))
		bin.coef.se <- apply(bincoef.list[[char]],1,sd,na.rm=TRUE)
		bin.coef.lb <- bin.coef.ci[,1]
		bin.coef.ub <- bin.coef.ci[,2]
		tempbin <- data.frame(x0,coef=bin.coef,se=bin.coef.se,CI_lower=bin.coef.lb,CI_upper=bin.coef.ub)
		rownames(tempbin) <- gp.lab
		est.bin[[other_treat.origin[char]]] <- tempbin
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
	est.matrix <- cov(t(marg.con))
	
	bin.coef <- output[,"BinCoef"][1:nbins]
	bin.coef.ci <- t(apply(coef.con, 1, quantile, CI.lvl, na.rm=TRUE))
	bin.coef.se <- apply(coef.con,1,sd, na.rm=TRUE)
	bin.coef.lb <- bin.coef.ci[,1]
	bin.coef.ub <- bin.coef.ci[,2]
	tempbin <- data.frame(x0,coef=bin.coef,se=bin.coef.se,CI_lower=bin.coef.lb,CI_upper=bin.coef.ub)
	rownames(tempbin) <- gp.lab
	est.bin <- tempbin
		
	
 }
 
 	if(predict==TRUE){
	est.predict.binning <- list()
		for(char in all_treat){
			fit <- predict_Ey[,char]
			pred.ci <- t(apply(pred.list[[char]], 1, quantile, CI.lvl,na.rm=TRUE))
			se.fit <- apply(pred.list[[char]],1,sd,na.rm=TRUE)
			lb <- pred.ci[,1]
			ub <- pred.ci[,2]
			
			if(treat.type=='discrete'){
				
				est.predict.binning[[all_treat.origin[char]]] <- cbind.data.frame(X = X.pred, EY = fit, 
                                           SE = se.fit ,Treatment=rep(all_treat.origin[char],length(X.lvls)),
                                           CI_lower=lb, CI_upper=ub
                                           )
										   
			}						   
			
			if(treat.type=='continuous'){
				est.predict.binning[[char]] <- cbind.data.frame(X = X.pred, EY = fit, 
                                           SE = se.fit ,Treatment=rep(char,length(X.lvls)),
                                           CI_lower=lb, CI_upper=ub
                                           )
			}
			}
	
	}
 
 ################### testing  ###############################
   if(treat.type=="continuous" & nbins>1){
	varD<-c()
	for (i in 1:nbins) {
		varD<-c(varD,var(data[groupX==i,D]/mean(data[groupX==i,D])))
	}
  
	## if the signs are correct
	## nbins!=3
	correctorder<-NULL
	if(nbins==3){
		correctOrder<-ifelse(as.numeric((bin.coef[1]-bin.coef[2])*(bin.coef[2]-bin.coef[3]))>0,TRUE,FALSE) 
	}

	#print(X.v)
	## p values
	pvalue <- function(i,j){
		stat <- (bin.coef[i]-bin.coef[j])/sqrt(bin.var[i,i]+bin.var[j,j]-2*bin.var[i,j])
		p <- (1-pt(abs(stat),df.X))*2
		return(p)
	}
  
	p.twosided<-NULL
  
	if (nbins==3) {
		p.twosided<-round(c(pvalue(1,2),pvalue(2,3),pvalue(1,3)),digits=4)
		names(p.twosided)<-c("p.1v2","p.2v3","p.1v3")
		names(bin.coef)<-c("X_low","X_med","X_high")
	} else if (nbins==2) {
		p.twosided<-round(pvalue(1,2),digits=4)
		names(p.twosided)<-c("p.LvH")
		names(bin.coef)<-c("X_low","X_high")
  } else if (nbins==4) {
		names(bin.coef)<-c("X_low","X_med1","X_med2","X_high")
  }
}

if(treat.type=="discrete" & nbins > 1){
  	
    group.Xcoefs <- list()
    group.p.twosided <- list()
    correctOrder.group <- list()
    
    for(char in other_treat) {
      
	  Xcoefs <- est.bin[[other_treat.origin[char]]][,"coef"]
      X.v <- bin.var.list[[char]]
      	  
      #correct order
      correctorder<-NULL
      if(nbins==3){
        correctOrder<-ifelse(as.numeric((Xcoefs[1]-Xcoefs[2])*(Xcoefs[2]-Xcoefs[3]))>0,TRUE,FALSE) 
        correctOrder.group[[other_treat.origin[char]]] <- correctOrder
      }

      pvalue<-function(i,j){
        stat<-(Xcoefs[i]-Xcoefs[j])/sqrt(X.v[i,i]+X.v[j,j]-2*X.v[i,j])
        p<-(1-pt(abs(stat),df.X))*2
        return(p)
      }
      p.twosided<-NULL
      if (nbins==3) {
        p.twosided<-round(c(pvalue(1,2),pvalue(2,3),pvalue(1,3)),digits=4)
        names(p.twosided)<-c("p.1v2","p.2v3","p.1v3")
        names(Xcoefs)<-c("X_low","X_med","X_high")
      } else if (nbins==2) {
        p.twosided<-round(pvalue(1,2),digits=4)
        names(p.twosided)<-c("p.LvH")
        names(Xcoefs)<-c("X_low","X_high")
      } else if (nbins==4) {
        names(Xcoefs)<-c("X_low","X_med1","X_med2","X_high")
      }
      
      group.Xcoefs[[other_treat.origin[char]]] <- Xcoefs
      group.p.twosided[[other_treat.origin[char]]] <- p.twosided
    }
	}
	p.wald <- NULL
 }
 
if(vartype!='bootstrap'){
  # if vartype!='bootstrap'
  ## estimate coefficients -> estimate covariance matrix 
  ## -> (predict==T) estimate expected value of Y 
  ## some tests
  if(TRUE){ #linear estimator
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
  v[which(is.na(v))] <- 0
  
  if(treat.type=='continuous'){ ## get variance
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
  
  if(treat.type=='discrete'){ ## get variance
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
  
 
  if (treat.type=='continuous'){
    marg<-coef.D + coef.DX*X.lvls
	marg[which(is.na(marg))] <- 0
    ## the variance is var(B1_D) + X^2*var(B_3) + 2*inst*cov(D, X)
    se<-sqrt(var.D +  X.lvls^2*var.DX + 2*X.lvls*cov.DX)
    df<-mod.naive$df.residual
    crit<-abs(qt(.025, df=df)) # critical values
  
    ##make 95% confidence bands. 
    lb<-marg-crit*se
    ub<-marg+crit*se
    est.lin<-data.frame(X.lvls, marg,se, lb, ub)
	
	# var-cov matrix
	gen_matrix <- function(colvec,x0){
		output <- var.D + (x0+colvec)*cov.DX + x0*colvec*var.DX
		return(output)
	}
	cov_matrix <- as.matrix(sapply(X.lvls,function(x0){gen_matrix(X.lvls,x0)}))
	est.matrix <- cov_matrix
  }
  
  if (treat.type=='discrete'){
    df <- mod.naive$df.residual
    crit<-abs(qt(.025, df=df))
    est.lin<-list()
    est.matrix <- list()
	
    for(char in other_treat) {
	  marg <- coef_list[[char]] + coef_inter_list[[char]]*X.lvls
	  marg[which(is.na(marg))] <- 0
      se <- sqrt(var_list[[char]] + X.lvls^2*varinter_list[[char]]+2*X.lvls*cov_list[[char]])
      lb <- marg-crit*se
      ub <- marg+crit*se
      tempest <- data.frame(X.lvls,marg,se,lb,ub)
      tempest[,'Treatment']<- rep(other_treat.origin[char],dim(tempest)[1])
      est.lin[[other_treat.origin[char]]] <- tempest
	  
	  # var-cov matrix
	  gen_matrix <- function(colvec,x0){
		output <- var_list[[char]] + (x0+colvec)*cov_list[[char]] + x0*colvec*varinter_list[[char]]
		return(output)
	  }
	  cov_matrix <- as.matrix(sapply(X.lvls,function(x0){gen_matrix(X.lvls,x0)}))
	  est.matrix[[other_treat.origin[char]]] <- cov_matrix
    }
  }

  }
  

 ##################################################

  if(TRUE){ # binning estimator
  if (treat.type=='continuous') {
  G<-DG<-GX<-DGX<-matrix(0,n,nbins)
  for (i in 1:nbins) {
    G[which(groupX==i),i]<-1
    DG[,i]<-data[,D]*G[,i]
    GX[,i]<-G[,i]*(data[,X]-x0[i])
    DGX[,i]<-DG[,i]*(data[,X]-x0[i])
  }
  
  ## formula and esitmation
  Gs<-GXs<-DGs<-DGXs<-c()
  for (i in 1:nbins) {
    Gs<-c(Gs,paste0("G[,",i,"]"))
    GXs<-c(GXs,paste0("GX[,",i,"]"))
    DGs<-c(DGs,paste0("DG[,",i,"]"))
    DGXs<-c(DGXs,paste0("DGX[,",i,"]"))
  }
  Xf<-paste0(Y,"~ -1+",paste0(DGs,collapse="+"),"+",paste0(DGXs,collapse="+"),
             "+",paste0(Gs,collapse="+"),"+",paste0(GXs,collapse="+"))
  if (is.null(Z)==FALSE) {
    Xf<- paste0(Xf,"+",paste0(Z,collapse="+"))
  }
  if (is.null(FE)==FALSE) {
    Xf <- paste0(Xf, "|",paste0(FE, collapse = "+"))
    if (vartype=="cluster") {
      Xf <- paste0(Xf, "| 0 |",paste0(cl,collapse = "+"))
    }
  }
  mod.Xf<-as.formula(Xf)    
  
  ## fit
  if (is.null(FE)==TRUE) { #OLS
    if (is.null(weights)==TRUE) {
      mod.X<-lm(mod.Xf,data=data)
    } else {
      mod.X<-lm(mod.Xf,data=data,weights=dataweights)
    }
  } else { # FE
    if (is.null(weights)==TRUE) {
      mod.X<-suppressWarnings(felm(mod.Xf,data=data))
    } else {
      mod.X<-suppressWarnings(felm(mod.Xf,data=data,weights=dataweights))
    }
  }
  
  ## coefficients and CIs
  if (is.null(FE)==TRUE) { #OLS
    if (vartype=="homoscedastic") {
      X.v<-vcov(mod.X)
    } else if (vartype=="robust") {
      X.v<-vcov(mod.X,type="HC1") ## White with small sample correction
    } else if (vartype=="cluster") {
      X.v<-vcovCluster(mod.X,cluster=data[,cl])
    } else if (vartype=="pcse") {
      if (is.null(Z)==FALSE) {
        exclude<-names(which(is.na(mod.X$coefficients)==TRUE))  ## drop colinear variables
        Z.ex<-setdiff(Z,exclude)
        Xf<-paste(Y,"~ -1+",paste(DGs,collapse="+"),"+",paste(DGXs,collapse="+"),
                  "+",paste(Gs,collapse="+"),"+",paste(GXs,collapse="+"),"+",paste(Z.ex,collapse="+"),sep="")
		Z <- Z.ex
        mod.X<-lm(as.formula(Xf),data=data)
      }
      X.v<-pcse(mod.X,groupN=data[,cl],groupT=data[,time],pairwise=pairwise)$vcov
    }
  } else { # FE
    if (vartype=="homoscedastic") {
      X.v<-vcov(mod.X, type = "iid")
    } else if (vartype=="robust") {
      X.v<-vcov(mod.X, type="robust") 
    } else if (vartype=="cluster") {
      X.v<-vcov(mod.X, type = "cluster") 
    }
  }
  
  # X.v should be checked by name  
  DGs <- c()
  for (i in 1:nbins){
      DGs<-c(DGs,paste0("DG[, ",i,"]"))
  }
  if(is.null(FE)==T){
        Xcoefs <- mod.X$coefficients[DGs]
  }
  else{
        Xcoefs <- mod.X$coefficients[DGs,]
  }
  
  if(vartype=='pcse'){
    pcse_DGs <- c()
    for (i in 1:nbins) {
         pcse_DGs<-c(pcse_DGs,paste0("DG","...",i,"."))
    }
         X.v <- X.v[pcse_DGs,pcse_DGs]
  }else{
        X.v <- X.v[DGs,DGs]
  }

  X.se<-sqrt(diag(as.matrix(X.v,drop=FALSE)))
  X.se[which(is.na(Xcoefs))]<-NA
  df.X<-mod.X$df.residual
  crit.X<-abs(qt(.025, df=df.X))
  lb.X<-Xcoefs-crit.X*X.se
  ub.X<-Xcoefs+crit.X*X.se
  
  est.bin <- data.frame(x0, Xcoefs, X.se, lb.X, ub.X)
  colnames(est.bin) <- c("x0", "coef", "se", "CI_lower", "CI_upper")
  rownames(est.bin) <- gp.lab
  }
  
  if(treat.type=="discrete"){ # discrete
    G<-GX<-matrix(0,n,nbins)
    for (i in 1:nbins) {
      G[which(groupX==i),i] <- 1
      GX[,i] <- G[,i]*(data[,X]-x0[i])
    }
    
    for (char in other_treat) {
      DG_name <- paste0("DG.",char)
      DGX_name <- paste0("DGX.",char)
      DG_matrix <- DGX_matrix <- matrix(0,n,nbins)
      for (i in 1:nbins) {
        DG_matrix[,i] <- (data[,D]==char)*G[,i]
        DGX_matrix[,i] <- DG_matrix[,i]*(data[,X]-x0[i])
      }
      assign(DG_name,DG_matrix)
      assign(DGX_name,DGX_matrix)
    }
    
    Gs<-GXs<-DGs<-DGXs<-c()
    for (i in 1:nbins)  {
      Gs<-c(Gs,paste0("G[,",i,"]"))
      GXs<-c(GXs,paste0("GX[,",i,"]"))
      for (char in other_treat) {
        DGs<-c(DGs,paste0("DG.",char,"[,",i,"]"))
        DGXs<-c(DGXs,paste0("DGX.",char,"[,",i,"]"))
      }
    }
    
    Xf<-paste0(Y,"~ -1+",paste0(DGs,collapse="+"),"+",paste0(DGXs,collapse="+"),
               "+",paste0(Gs,collapse="+"),"+",paste0(GXs,collapse="+"))
    
    if (is.null(Z)==FALSE) {
      Xf<- paste0(Xf,"+",paste0(Z,collapse="+"))
    }
    if (is.null(FE)==FALSE) {
      Xf <- paste0(Xf, "|",paste0(FE, collapse = "+"))
      if (vartype=="cluster") {
        Xf <- paste0(Xf, "| 0 |",paste0(cl,collapse = "+"))
      }
    }
    mod.Xf<-as.formula(Xf)    
    
    ## fit
    if (is.null(FE)==TRUE) { #OLS
      if (is.null(weights)==TRUE) {
        mod.X<-lm(mod.Xf,data=data)
      } else {
        mod.X<-lm(mod.Xf,data=data,weights=dataweights)
      }
    } else { # FE
      if (is.null(weights)==TRUE) {
        mod.X<-suppressWarnings(felm(mod.Xf,data=data))
      } else {
        mod.X<-suppressWarnings(felm(mod.Xf,data=data,weights=dataweights))
      }
    }
    
    ## variance
    if (is.null(FE)==TRUE) { #OLS
      if (vartype=="homoscedastic") {
        X.v_all<-vcov(mod.X)
      } else if (vartype=="robust") {
        X.v_all<-vcov(mod.X,type="HC1") ## White with small sample correction
      } else if (vartype=="cluster") {
        X.v_all<-vcovCluster(mod.X,cluster=data[,cl])
      } else if (vartype=="pcse") {
        if (is.null(Z)==FALSE) {
          exclude<-names(which(is.na(mod.X$coefficients)==TRUE))  ## drop colinear variables
          Z.ex<-setdiff(Z,exclude)
          Xf<-paste(Y,"~ -1+",paste(DGs,collapse="+"),"+",paste(DGXs,collapse="+"),
                    "+",paste(Gs,collapse="+"),"+",paste(GXs,collapse="+"),"+",paste(Z.ex,collapse="+"),sep="")
		  Z <- Z.ex
          mod.X<-lm(as.formula(Xf),data=data)
        }
        X.v_all<-pcse(mod.X,groupN=data[,cl],groupT=data[,time],pairwise=pairwise)$vcov
      }
    } else { # FE
      if (vartype=="homoscedastic") {
        X.v_all<-vcov(mod.X, type = "iid")
      } else if (vartype=="robust") {
        X.v_all<-vcov(mod.X, type="robust") 
      } else if (vartype=="cluster") {
        X.v_all<-vcov(mod.X, type = "cluster") 
      }
    }
    
    est.bin=list()
    df.X <- mod.X$df.residual
    crit.X <- abs(qt(.025,df=df.X))
    
    for (char in other_treat) {
      DGs <- c()
      for (i in 1:nbins) {
        DGs<-c(DGs,paste0("DG.",char,"[, ",i,"]"))
      }
      
      if(is.null(FE)==T){
        Xcoefs <- mod.X$coefficients[DGs]
        }
      else{
        Xcoefs <- mod.X$coefficients[DGs,]
      }
      
      if(vartype=='pcse') {
        pcse_DGs <- c()
        for (i in 1:nbins) {
          pcse_DGs<-c(pcse_DGs,paste0("DG.",char,"...",i,"."))
        }
        X.v <- X.v_all[pcse_DGs,pcse_DGs]
      }
      else {
        X.v <- X.v_all[DGs,DGs]
      }
      
      X.se <- sqrt(diag(as.matrix(X.v,drop=F)))
      X.se[which(is.na(Xcoefs))] <- NA
      lb.X <- Xcoefs-crit.X*X.se
      ub.X <- Xcoefs+crit.X*X.se
      est_bin <- data.frame(x0,Xcoefs,X.se,lb.X,ub.X)
      colnames(est_bin) <- c("x0", "coef", "se", "CI_lower", "CI_upper")
      rownames(est_bin) <- gp.lab
      est.bin[[other_treat.origin[char]]] <- est_bin
    }
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
  
	X.pred <- seq(min(data[,X]),max(data[,X]),length.out=length(X.lvls))
	est.predict.binning <- list()
	if(treat.type=='continuous'){
		n<-dim(data)[1]
		data_touse <- data
		if(is.null(weights)==F){
			weight_touse <- as.matrix(data_touse[,weights])
		}
		else{
			weight_touse <- as.matrix(rep(1,n))
		}

		data1 <- as.data.frame(data_touse[,Y])
		colnames <- c("Y")
		formula <- "Y ~ -1"
		for(char in Z){
			data1 <- cbind(data1,data_touse[,char])
			colnames <- c(colnames,char)
			formula <- paste(formula,char,sep = "+")
		}
		for (i in 1:nbins){
			data1 <- cbind(data1,as.numeric(groupX==i))
			colnames <- c(colnames,paste0("G.",i))
			formula <- paste(formula,paste0("G.",i),sep = "+")
		}
		for (i in 1:nbins){
			data1 <- cbind(data1,as.numeric(groupX==i)*data_touse[,D])
			colnames <- c(colnames,paste0("G.",i,".D"))
			formula <- paste(formula,paste0("G.",i,".D"),sep="+")
		}
		for (i in 1:nbins){
			data1 <- cbind(data1,as.numeric(groupX==i)*(data_touse[,X]-x0[i]))
			colnames <- c(colnames,paste0("G.",i,".X"))
			formula <- paste(formula,paste0("G.",i,".X"),sep = "+")
		}
		for (i in 1:nbins){
			data1 <- cbind(data1,as.numeric(groupX==i)*data_touse[,D]*(data_touse[,X]-x0[i]))
			colnames <- c(colnames,paste0("G.",i,".DX"))
			formula <- paste(formula,paste0("G.",i,".DX"),sep="+")
		}
			
		if (is.null(FE)==FALSE) {
			formula <- paste0(formula, "|",paste0(FE, collapse = "+"))
			data1 <- cbind(data1,data_touse[,FE])
			colnames <- c(colnames,FE)
			if(vartype=='cluster'){
				formula <- paste0(formula, "| 0 |",paste0(cl,collapse = "+"))
				if(!cl%in%FE){
					data1 <- cbind(data1,data_touse[,cl])
					colnames <- c(colnames,cl)
				}
			}
		}
		
		colnames(data1) <- colnames
		formula <- as.formula(formula)
		if(is.null(FE)==TRUE){
			binning_fit <- lm(formula=formula,data=data1,weights = weight_touse)
		}
		if(is.null(FE)==FALSE){
			suppressWarnings(binning_fit <- felm(formula=formula,data=data1,weights = weight_touse))
		}
	  
		#coef
		binning_coef <- binning_fit$coefficients
		binning_coef[which(is.na(binning_coef))] <- 0
		binning_coef <- as.matrix(binning_coef)
		
		#var
		if (is.null(FE)==TRUE) { #OLS
		  if (vartype=="homoscedastic") {
				v<-vcov(binning_fit)
		} else if (vartype=="robust") {
			requireNamespace("sandwich")
			v<-vcov(binning_fit,type="HC1") # White with small sample correction
		} else if (vartype=="cluster") {
			v<-vcovCluster(binning_fit,cluster = data[,cl])
		} else if (vartype=="pcse") {
			requireNamespace("pcse")
			v<-pcse(binning_fit,groupN=data[,cl],groupT=data[,time],pairwise=pairwise)$vcov
		}
		} else { # FE
			if (vartype=="homoscedastic") {
				v<-vcov(binning_fit, type = "iid")
			} else if (vartype=="robust") {
				v<-vcov(binning_fit, type="robust") 
			} else if (vartype=="cluster") {
				v<-vcov(binning_fit, type = "cluster") 
		}
		}
		v[which(is.na(v))] <- 0
		#newdata
		npred <- length(X.pred)
		x_predict <- X.pred
		data_predict_start <- matrix(NA,nrow = npred,ncol = 0)
		for (char in Z){
			data_predict_start <- cbind(data_predict_start,mean(data_touse[,char]))
		}

		groupX_predict<-cut(x_predict,breaks=cuts.X, labels = FALSE)
		groupX_predict[which(x_predict==min(x_predict))]<-1
		
		
		for(char in all_treat){
			target.D <- D.sample[char]
			data_predict <- data_predict_start
			for(i in 1:nbins){
				data_predict <- cbind(data_predict,as.numeric(groupX_predict==i))
			}
			for(i in 1:nbins){
				data_predict <- cbind(data_predict,as.numeric(groupX_predict==i)*target.D)
			}
			for(i in 1:nbins){
				data_predict <- cbind(data_predict,as.numeric(groupX_predict==i)*(x_predict-x0[i]))
			}
			for(i in 1:nbins){
				data_predict <- cbind(data_predict,as.numeric(groupX_predict==i)*target.D*(x_predict-x0[i]))
			}
			
			
			m.mat <- as.matrix(data_predict)
			fit <- as.vector(m.mat %*% binning_coef)
			se.fit <- sqrt(diag(m.mat%*%v%*%t(m.mat)))
			df2 <- binning_fit$df
			CI.lvl <- c((1-0.95)/2, (1-(1-0.95)/2))
			crit<-abs(qt(CI.lvl[1], df=df2))
			lb<-fit-crit*se.fit
			ub<-fit+crit*se.fit
			est.predict.binning[[char]] <- cbind.data.frame(X = X.pred, EY = fit, 
                                           SE = se.fit ,Treatment=char,
                                           CI_lower=lb, CI_upper=ub
                                           )
		}
	}
	
	if(treat.type=='discrete'){
		n <- dim(data)[1]
		data_touse <- data
		if(is.null(weights)==F){
			weight_touse <- as.matrix(data[,weights])
		}
		else{
			weight_touse <- as.matrix(rep(1,n))
		}
	    
		data1 <- as.data.frame(data_touse[,Y])
		colnames <- c("Y")
		formula <- "Y ~ -1"
		for(char in Z){
			data1 <- cbind(data1,data_touse[,char])
									colnames <- c(colnames,char)
									formula <- paste(formula,char,sep = "+")
		}
    
		for (i in 1:nbins){
			data1 <- cbind(data1,as.numeric(groupX==i))
			colnames <- c(colnames,paste0("G.",i))
			formula <- paste(formula,paste0("G.",i),sep = "+")
		}
    
		for (i in 1:nbins){
				for(char in other_treat){
						data1 <- cbind(data1,as.numeric(groupX==i)*as.numeric(data_touse[,D]==char))
						colnames <- c(colnames,paste0("G.",i,".D.",char))
						formula <- paste(formula,paste0("G.",i,".D.",char),sep="+")
				}
			}
    
		for (i in 1:nbins){
				data1 <- cbind(data1,as.numeric(groupX==i)*(data_touse[,X]-x0[i]))
				colnames <- c(colnames,paste0("G.",i,".X"))
				formula <- paste(formula,paste0("G.",i,".X"),sep = "+")
			}
    
		for (i in 1:nbins){
				for(char in other_treat){
						data1 <- cbind(data1,as.numeric(groupX==i)*as.numeric(data_touse[,D]==char)*(data_touse[,X]-x0[i]))
						colnames <- c(colnames,paste0("G.",i,".D.",char,".X"))
						formula <- paste(formula,paste0("G.",i,".D.",char,".X"),sep="+")
				}
			}
	
		if (is.null(FE)==FALSE) {
			formula <- paste0(formula, "|",paste0(FE, collapse = "+"))
			data1 <- cbind(data1,data_touse[,FE])
			colnames <- c(colnames,FE)
			if(vartype=='cluster'){
				formula <- paste0(formula, "| 0 |",paste0(cl,collapse = "+"))
				if(!cl%in%FE){
					data1 <- cbind(data1,data_touse[,cl])
					colnames <- c(colnames,cl)
				}
			}
		}
		colnames(data1) <- colnames
		formula <- as.formula(formula)
		if(is.null(FE)==TRUE){
			binning_fit <- lm(formula=formula,data=data1,weights = weight_touse)
		}
		if(is.null(FE)==FALSE){
			suppressWarnings(binning_fit <- felm(formula=formula,data=data1,weights = weight_touse))
		}
	
		binning_coef <- binning_fit$coefficients
		binning_coef[which(is.na(binning_coef))] <- 0
		binning_coef <- as.matrix(binning_coef)
    
		#var
		if (is.null(FE)==TRUE) { #OLS
		  if (vartype=="homoscedastic") {
				v<-vcov(binning_fit)
		} else if (vartype=="robust") {
			requireNamespace("sandwich")
			v<-vcov(binning_fit,type="HC1") # White with small sample correction
		} else if (vartype=="cluster") {
			v<-vcovCluster(binning_fit,cluster = data[,cl])
		} else if (vartype=="pcse") {
			requireNamespace("pcse")
			v<-pcse(binning_fit,groupN=data[,cl],groupT=data[,time],pairwise=pairwise)$vcov
		}
		} else { # FE
			if (vartype=="homoscedastic") {
				v<-vcov(binning_fit, type = "iid")
			} else if (vartype=="robust") {
				v<-vcov(binning_fit, type="robust") 
			} else if (vartype=="cluster") {
				v<-vcov(binning_fit, type = "cluster") 
		}
		}
		v[which(is.na(v))] <- 0
	

		x_predict <- X.pred
		npred <- length(X.pred)
		data_predict_start <- matrix(NA,nrow = npred,ncol = 0)
		for (char in Z){
			data_predict_start <- cbind(data_predict_start,mean(data_touse[,char]))
		}

		groupX_predict<-cut(x_predict,breaks=cuts.X, labels = FALSE)
		groupX_predict[which(x_predict==min(x_predict))] <- 1
    
		for(target_treat in all_treat){
			data_predict <- data_predict_start
			for(i in 1:nbins){
				data_predict <- cbind(data_predict,as.numeric(groupX_predict==i))
			}
			for(i in 1:nbins){
				for(char in other_treat){
						data_predict <- cbind(data_predict,as.numeric(groupX_predict==i)*as.numeric(char==target_treat))
					}
			}
			for(i in 1:nbins){
				data_predict <- cbind(data_predict,as.numeric(groupX_predict==i)*(x_predict-x0[i]))
			}
			for(i in 1:nbins){
				for(char in other_treat){
					data_predict <- cbind(data_predict,as.numeric(groupX_predict==i)*as.numeric(char==target_treat)*(x_predict-x0[i]))
					}
			}
			
			m.mat <- as.matrix(data_predict)
			fit <- as.vector(m.mat %*% binning_coef)
			se.fit <- sqrt(diag(m.mat%*%v%*%t(m.mat)))
			df2 <- binning_fit$df
			CI.lvl <- c((1-0.95)/2, (1-(1-0.95)/2))
			crit<-abs(qt(CI.lvl[1], df=df2))
			lb<-fit-crit*se.fit
			ub<-fit+crit*se.fit
			est.predict.binning[[all_treat.origin[target_treat]]] <- cbind.data.frame(X = X.pred, EY = fit, 
                                           SE = se.fit ,Treatment=rep(all_treat.origin[target_treat],length(x_predict)),
                                           CI_lower=lb, CI_upper=ub
                                           )
			}
	}
	data <- data.old
	FE <- FE.old
  }
  
  
  ################### testing  ###############################
  


  
  ## variance of treatment in each group 
  if(treat.type=="continuous" & nbins>1){
	varD<-c()
	for (i in 1:nbins) {
		varD<-c(varD,var(data[groupX==i,D]/mean(data[groupX==i,D])))
	}
  
	## if the signs are correct
	## nbins!=3
	correctorder<-NULL
	if(nbins==3){
		correctOrder<-ifelse(as.numeric((Xcoefs[1]-Xcoefs[2])*(Xcoefs[2]-Xcoefs[3]))>0,TRUE,FALSE) 
	}

	#print(X.v)
	## p values
	pvalue<-function(i,j){
		stat<-(Xcoefs[i]-Xcoefs[j])/sqrt(X.v[i,i]+X.v[j,j]-2*X.v[i,j])
		p<-(1-pt(abs(stat),df.X))*2
		return(p)
	}
  
	p.twosided<-NULL
  
	if (nbins==3) {
		p.twosided<-round(c(pvalue(1,2),pvalue(2,3),pvalue(1,3)),digits=4)
		names(p.twosided)<-c("p.1v2","p.2v3","p.1v3")
		names(Xcoefs)<-c("X_low","X_med","X_high")
	} else if (nbins==2) {
		p.twosided<-round(pvalue(1,2),digits=4)
		names(p.twosided)<-c("p.LvH")
		names(Xcoefs)<-c("X_low","X_high")
	} else if (nbins==4) {
		names(Xcoefs)<-c("X_low","X_med1","X_med2","X_high")
  }
}
  
  if(treat.type=="discrete" & nbins > 1){
  	
    group.Xcoefs <- list()
    group.p.twosided <- list()
    correctOrder.group <- list()
    
    for(char in other_treat) {
      
      DGs <- c()
      for (i in 1:nbins) {
        DGs<-c(DGs,paste0("DG.",char,"[, ",i,"]"))
      }
      
      if(is.null(FE)==T){
        Xcoefs <- mod.X$coefficients[DGs]
      }
      else{
        Xcoefs <- mod.X$coefficients[DGs,]
      }
      
      if(vartype=='pcse') {
        pcse_DGs <- c()
        for (i in 1:nbins) {
          pcse_DGs<-c(pcse_DGs,paste0("DG.",char,"...",i,"."))
        }
        X.v <- X.v_all[pcse_DGs,pcse_DGs]
      }
      else {
        X.v <- X.v_all[DGs,DGs]
      }
      	  
      #correct order
      correctorder<-NULL
      if(nbins==3){
        correctOrder<-ifelse(as.numeric((Xcoefs[1]-Xcoefs[2])*(Xcoefs[2]-Xcoefs[3]))>0,TRUE,FALSE) 
        correctOrder.group[[other_treat.origin[char]]] <- correctOrder
      }

      pvalue<-function(i,j){
        stat<-(Xcoefs[i]-Xcoefs[j])/sqrt(X.v[i,i]+X.v[j,j]-2*X.v[i,j])
        p<-(1-pt(abs(stat),df.X))*2
        return(p)
      }
      p.twosided<-NULL
      if (nbins==3) {
        p.twosided<-round(c(pvalue(1,2),pvalue(2,3),pvalue(1,3)),digits=4)
        names(p.twosided)<-c("p.1v2","p.2v3","p.1v3")
        names(Xcoefs)<-c("X_low","X_med","X_high")
      } else if (nbins==2) {
        p.twosided<-round(pvalue(1,2),digits=4)
        names(p.twosided)<-c("p.LvH")
        names(Xcoefs)<-c("X_low","X_high")
      } else if (nbins==4) {
        names(Xcoefs)<-c("X_low","X_med1","X_med2","X_high")
      }
      
      group.Xcoefs[[other_treat.origin[char]]] <- Xcoefs
      group.p.twosided[[other_treat.origin[char]]] <- p.twosided
      
    }
	}
    

  
  
  ##############  Wald Test #####################
  
  if (wald == TRUE & treat.type=='continuous') { 
    ## formula
    formula0 <- paste(Y,"~",D,"+",X,"+",D,"*",X)
    ## create dummies for bins and interactions
    ## G -- a matrix of group dummies
    ## DG -- a matrix of interactions
    G<-DG<-GX<-DGX<-matrix(0,n,(nbins-1))
    for (i in 1:(nbins-1)) {
      G[which(groupX==(i+1)),i]<-1
      DG[,i]<-data[,D]*G[,i]
      GX[,i]<-data[,X]*G[,i]
      DGX[,i]<-data[,D]*data[,X]*G[,i]
    } 
    ## formula and esitmation
    Gs<-GXs<-DGs<-DGXs<-c()
    for (i in 2:nbins)  {
      Gs<-c(Gs,paste("G",i,sep=""))
      GXs<-c(GXs,paste("GX",i,sep=""))
      DGs<-c(DGs,paste("DG",i,sep=""))
      DGXs<-c(DGXs,paste("DGX",i,sep=""))
    }
    colnames(G) <- Gs
    colnames(DG) <- DGs
    colnames(GX) <- GXs
    colnames(DGX) <- DGXs
    data.aug <- cbind.data.frame(data, G, DG, GX, DGX)
    formula1<-paste(formula0,
                    "+",paste(Gs,collapse=" + "),
                    "+",paste(GXs,collapse=" + "),
                    "+",paste(DGs,collapse=" + "),
                    "+",paste(DGXs,collapse=" + "))
    if (is.null(Z)==FALSE) {
      formula0 <- paste0(formula0, "+",paste(Z,collapse=" + "))
      formula1 <- paste0(formula1, "+",paste(Z,collapse=" + "))
    } 
    if (is.null(FE)==FALSE) {
      formula0 <- paste0(formula0, "|",paste0(FE, collapse=" + "))
      formula1 <- paste0(formula1, "|",paste0(FE, collapse = "+"))
      if (vartype=="cluster") {
        formula0 <- paste0(formula0, "| 0 |",paste0(cl,collapse = "+"))
        formula1 <- paste0(formula1, "| 0 |",paste0(cl,collapse = "+"))
      }
    }
    mod.formula0<-as.formula(formula0)
    mod.formula1<-as.formula(formula1)    
    
    ## fit
    if (is.null(FE)==TRUE) { #OLS
      ## fit
      if (is.null(weights)==TRUE) {
        mod.re<-lm(mod.formula0,data=data.aug)
        mod.un<-lm(mod.formula1,data=data.aug)
      } else {
        mod.re<-lm(mod.formula0,data=data.aug,weights=dataweights)
        mod.un<-lm(mod.formula1,data=data.aug,weights=dataweights)
      }
      
      ## vcov
      if (is.null(vartype)==TRUE) {vartype <- "homoscedastic"}
      if (vartype=="homoscedastic") {
        v<-vcov(mod.un)
      } else if (vartype=="robust") {
        v<-vcov(mod.un,type="HC1") # White with small sample correction
      } else if (vartype=="cluster") {
        v<-vcovCluster(mod.un,cluster = data.aug[,cl])
      } else if (vartype=="pcse") {
        v<-pcse(mod.un,
                groupN=data.aug[,cl],
                groupT=data.aug[,time],
                pairwise=pairwise)$vcov
      }
	  
	  ## wald test
      requireNamespace("lmtest")

      wald.out <- tryCatch(
        p.wald <- round(lmtest::waldtest(mod.re, mod.un,test="Chisq", vcov=v)[[4]][2],4),
        error = function(e){return(NULL)}
      )   
      ## warning
      if (is.null(wald.out)==TRUE) {
        p.wald <- NULL
        warning("Var-cov matrix nearly singular in the Wald test.")
      }  
    } else { # FE
      requireNamespace("lfe")
      ## fit
      if (is.null(weights)==TRUE) {
        mod.un<-suppressWarnings(felm(mod.formula1,data=data.aug))
      } else {
        mod.un<-suppressWarnings(felm(mod.formula1,data=data.aug,weights=dataweights))
      }
      ## wald test
      constraints <- as.formula(paste0("~",paste0(c(Gs,GXs,DGs,DGXs), collapse = "|")))            
      if (vartype=="homoscedastic") {
        p.wald <- lfe::waldtest(mod.un, constraints, type = "default")[1]
      } else if (vartype=="robust") {
        p.wald <- lfe::waldtest(mod.un, constraints, type = "robust")[1]
      } else {
        p.wald <- lfe::waldtest(mod.un, constraints)[1] # clustered
      }
      names(p.wald) <- NULL            
      p.wald <- round(p.wald,4)        
    }
 }
  
  
 if(wald==TRUE & treat.type=='discrete') { #discrete wald test
    formula0 <- paste(Y,"~",D,"+",X,"+",D,"*",X)
    G<-GX<-matrix(0,n,nbins-1)
    
    for (i in 1:(nbins-1)) {
      G[which(groupX==(i+1)),i] <- 1
      GX[,i] <- G[,i]*(data[,X])
    }
    Gs<-GXs<-c()
    for (i in 1:(nbins-1))  {
      Gs<-c(Gs,paste0("G.",i))
      GXs<-c(GXs,paste0("GX.",i))
    }
    colnames(G) <- Gs
    colnames(GX) <- GXs
    
    for (char in other_treat) {
      
      DG_name <- paste0("DG.",char)
      DGX_name <- paste0("DGX.",char)
      DG_matrix <- DGX_matrix <- matrix(0,n,nbins-1)
      DGs <- c()
      DGXs <- c()
      for (i in 1:(nbins-1)) {
        DG_matrix[,i] <- (data[,D]==char)*G[,i]
        DGX_matrix[,i] <- DG_matrix[,i]*(data[,X])
        DGs <- c(DGs,paste0("DG.",char,".",i))
        DGXs <- c(DGXs,paste0("DGX.",char,".",i))
      }
      colnames(DG_matrix) <- DGs
      colnames(DGX_matrix) <- DGXs
      assign(DG_name,DG_matrix)
      assign(DGX_name,DGX_matrix)
    }
    
    
    data.aug <- cbind.data.frame(data, G, GX)
    
    for (char in other_treat) {
      DG_name <- paste0("DG.",char)
      data.aug <- cbind.data.frame(data.aug, get(DG_name)) 
    }
    
    for (char in other_treat) {
      DGX_name <- paste0("DGX.",char)
      data.aug <- cbind.data.frame(data.aug, get(DGX_name)) 
    }
    
    Gs<-GXs<-DGs<-DGXs<-c()
    for (i in 1:(nbins-1))  {
      Gs<-c(Gs,paste0("G.",i))
      GXs<-c(GXs,paste0("GX.",i))
      for (char in other_treat) {
        DGs<-c(DGs,paste0("DG.",char,".",i))
        DGXs<-c(DGXs,paste0("DGX.",char,".",i))
      }
    }
    
    
    
    formula1<-paste(formula0,
                    "+",paste(Gs,collapse=" + "),
                    "+",paste(GXs,collapse=" + "),
                    "+",paste(DGs,collapse=" + "),
                    "+",paste(DGXs,collapse=" + "))
    
    if (is.null(Z)==FALSE) {
      formula0 <- paste0(formula0, "+",paste(Z,collapse=" + "))
      formula1 <- paste0(formula1, "+",paste(Z,collapse=" + "))
    } 
    if (is.null(FE)==FALSE) {
      formula0 <- paste0(formula0, "|",paste0(FE, collapse=" + "))
      formula1 <- paste0(formula1, "|",paste0(FE, collapse = "+"))
      if (vartype=="cluster") {
        formula0 <- paste0(formula0, "| 0 |",paste0(cl,collapse = "+"))
        formula1 <- paste0(formula1, "| 0 |",paste0(cl,collapse = "+"))
      }
    }
    mod.formula0<-as.formula(formula0)
    mod.formula1<-as.formula(formula1)
    
    if (is.null(FE)==TRUE) { #OLS
      ## fit
      if (is.null(weights)==TRUE) {
        mod.re<-lm(mod.formula0,data=data.aug)
        mod.un<-lm(mod.formula1,data=data.aug)
      } else {
        mod.re<-lm(mod.formula0,data=data.aug,weights=dataweights)
        mod.un<-lm(mod.formula1,data=data.aug,weights=dataweights)
      }
      ## vcov
      if (is.null(vartype)==TRUE) {vartype <- "homoscedastic"}
      if (vartype=="homoscedastic") {
        v<-vcov(mod.un)
      } else if (vartype=="robust") {
        v<-vcov(mod.un,type="HC1") # White with small sample correction
      } else if (vartype=="cluster") {
        v<-vcovCluster(mod.un,cluster = data.aug[,cl])
      } else if (vartype=="pcse") {
        v<-pcse(mod.un,
                groupN=data.aug[,cl],
                groupT=data.aug[,time],
                pairwise=pairwise)$vcov
      }
      ## wald test
      requireNamespace("lmtest")
      wald.out <- tryCatch(
        p.wald <- round(lmtest::waldtest(mod.re, mod.un,test="Chisq", vcov=v)[[4]][2],4),
        error = function(e){return(NULL)}
      )   
      ## warning
      if (is.null(wald.out)==TRUE) {
        p.wald <- NULL
        warning("Var-cov matrix nearly singular in the Wald test.")
      }  
    } else { # FE
      requireNamespace("lfe")
      ## fit
      if (is.null(weights)==TRUE) {
        mod.un<-suppressWarnings(felm(mod.formula1,data=data.aug))
      } else {
        mod.un<-suppressWarnings(felm(mod.formula1,data=data.aug,weights=dataweights))
      }
      
      ## wald test
      constraints <- as.formula(paste0("~",paste0(c(Gs,GXs,DGs,DGXs), collapse = "|")))  
      
      if (vartype=="homoscedastic") {
        p.wald <- lfe::waldtest(mod.un, constraints, type = "default")[1]
      } else if (vartype=="robust") {
        p.wald <- lfe::waldtest(mod.un, constraints, type = "robust")[1]
      } else {
        p.wald <- lfe::waldtest(mod.un, constraints)[1] # clustered
      }
      names(p.wald) <- NULL            
      p.wald <- round(p.wald,4)        
    }
    
  }

  # end of Wald test
  
}  
    
if(TRUE){## get densities and histogram
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

if(TRUE){## Storage

  ## L kurtosis
  requireNamespace("Lmoments")
  Lkurtosis <- Lmoments(data[,X],returnobject=TRUE)$ratios[4]
  if(nbins>1){
  if(treat.type=='discrete'){
  tests <- list(
    treat.type = treat.type,
    #est.binning = est.bin,
    bin.size = sprintf(c(table(groupX)/length(groupX)),fmt = '%#.3f'),
    X.Lkurtosis = sprintf(Lkurtosis,fmt = '%#.3f')
  )
  if (nbins==3) {
    tests<-c(tests,list(correctOrder=correctOrder.group))
  }
  if (nbins%in%c(2,3)) {
    tests<-c(tests,list(p.twosided=group.p.twosided, fmt = '%#.3f')) 
  }
  if (wald == TRUE & vartype!='bootstrap') {
    tests <- c(tests, list(p.wald = sprintf(p.wald, fmt = '%#.3f')))
  }    
  }
  
  if(treat.type=='continuous'){
    tests <- list(
      treat.type = treat.type,
      #est.binning = est.bin,
      bin.size = sprintf(c(table(groupX)/length(groupX)),fmt = '%#.3f'),
      treat.variation.byGroup=sprintf(varD,fmt = '%#.3f'),
      X.Lkurtosis = sprintf(Lkurtosis,fmt = '%#.3f')
    )
    if (nbins==3) {
      tests<-c(tests,list(correctOrder=correctOrder))
    }
    if (nbins%in%c(2,3)) {
      tests<-c(tests,list(p.twosided=sprintf(p.twosided, fmt = '%#.3f'))) 
    }
    if (wald == TRUE & vartype!='bootstrap') {
      tests <- c(tests, list(p.wald = sprintf(p.wald, fmt = '%#.3f')))
    }    
  }
  }

  ## saving
  if(nbins==1){
    est.bin <- NULL
    tests <- NULL
  }


  if(treat.type=='discrete'){
  output<-list(
    type = "binning",
    nbins = nbins,
    est.lin = est.lin,
	vcov.matrix = est.matrix,
    est.bin = est.bin,
    treat.type = treat.type,
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
    tests = tests,
	predict = predict
  )
  }
  
  if(treat.type=="continuous"){
    output<-list(
      type = "binning",
      nbins = nbins,
      est.lin = est.lin,
	  vcov.matrix = est.matrix,
      est.bin = est.bin,
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
	output <- c(output,list(est.predict = est.predict.binning,labelname = labelname,all.treat = all_treat.origin,
						    X=X,Y=Y,D=D,D.ref=D.ref))
  }
  
}

## Plot
class(output) <- "interflex"
if (figure==TRUE & pool==FALSE) {
    
    class(output) <- "interflex"
    suppressMessages(
    graph <- plot.interflex(x = output, CI = CI, xlab = xlab, ylab = ylab, color = color, order = order,
                        subtitles = subtitles,
                        show.subtitles = show.subtitles,
                        Ylabel = Ylabel, Dlabel = Dlabel, Xlabel = Xlabel, 
                        main = main, xlim = xlim, ylim = ylim, Xdistr = Xdistr,interval = interval,pool=pool,
                        file = file, theme.bw = theme.bw, show.grid = show.grid, bin.labs = bin.labs, ncols=ncols,
                        cex.main = cex.main,cex.sub = cex.sub, cex.axis = cex.axis, cex.lab = cex.lab)  
    )
    output<-c(output,list(graph=graph))
    class(output) <- "interflex"
  }
  
if(figure==TRUE & pool==TRUE & treat.type=='discrete'){
    suppressMessages(
    graph <- plot.interflex(x = output, CI = CI, xlab = xlab, ylab = ylab,order=order,subtitles = subtitles,show.subtitles = show.subtitles,
                              Ylabel = Ylabel, Dlabel = Dlabel, Xlabel = Xlabel, 
                              main = main, xlim = xlim, ylim = ylim, Xdistr = Xdistr,interval = interval,pool=pool,color=color,
                              file = file, theme.bw = theme.bw, show.grid = show.grid,
                              cex.main = cex.main, cex.axis = cex.axis, cex.lab = cex.lab,jitter=jitter,legend.title =legend.title)
    )
    output <- c(output, list(graph = graph))
  }
class(output) <- "interflex"
return(output)
}


