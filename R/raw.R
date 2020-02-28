inter.raw<-function(
  data,Y,D,X,
  treat.type=NULL,
  subtitles=NULL,
  order=NULL,
  Z=NULL,
  FE=NULL,
  weights=NULL,
  full.moderate = TRUE,
  nbins=3,
  cutoffs=NULL, 
  span=NULL,
  pos=NULL,
  main = NULL,
  Ylabel=NULL,
  Dlabel= NULL,
  Xlabel=NULL,
  theme.bw = FALSE, 
  show.grid = TRUE,
  cex.main = NULL,
  cex.lab = NULL,
  cex.axis = NULL,
  ncols = NULL
){ 
  
  if (is.data.frame(data) == FALSE) {
    data <- as.data.frame(data)
  }
  
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
  # treat.type default value
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
  
  if(treat.type=='discrete'){
	all.treat <- unique(data[,D])
	if(is.null(order)==FALSE){
	  order <- as.character(order)
	  if(length(order)!=length(unique(order))){
        stop("\"order\" should not contain repeated values.")
      }
      if(length(order)!=length(all.treat)){
        stop("\"order\" should contain all kinds of treatment arms.")
      }

      if(sum(!is.element(order,all.treat))!=0 | sum(!is.element(all.treat,order))!=0){
        stop("\"order\" should contain all kinds of treatment arms.")
      }
  } else {order <- sort(all.treat)}
  
  subtitle <- subtitles
  if(is.null(subtitle)==FALSE){
	subtitle <- as.character(subtitle)
	if(length(subtitle)!=length(all.treat)){
		stop("The number of elements in \"subtitles\" should equal to the number of different treatment arms.")
	}
  }
  }
  
  if (is.character(X) == FALSE) {
    stop("\"X\" is not a string.")
  } else {
    X <- X[1]
  }
  
  if (is.null(Z) == FALSE) {
    for (i in 1:length(Z)) {
      if (is.character(Z[i]) == FALSE) {
        stop("Some element in \"Z\" is not a string.")
      }
    }
  }
  
  if (is.logical(full.moderate) == FALSE & is.numeric(full.moderate)==FALSE) {
    stop("\"full.moderate\" is not a logical flag.")
  }else{
    full <- full.moderate
  } 
  
  if (is.null(FE) == FALSE) {
    requireNamespace("lfe")
    for (i in 1:length(FE)) {
      if (is.character(FE[i]) == FALSE) {
        stop("Some element in \"FE\" is not a string.")
      }
    }
  }
  
  if (is.null(weights) == FALSE) {
    if (is.character(weights) == FALSE) {
      stop("\"weights\" is not a string.")
    } else {
      weights <- weights[1]
    }   
  }
  if (is.null(nbins) == FALSE) {
    if (nbins%%1 != 0) {
      stop("\"nbins\" is not a positive integer.")
    } else {
      nbins <- nbins[1]
    }
    if (nbins < 1) {
      stop("\"nbins\" is not a positive integer.")
    }
  }
  if (is.null(cutoffs) == FALSE) {
    if (is.numeric(cutoffs) == FALSE) {
      stop("Some element in \"cutoffs\" is not numeric.")
    } 
  }
  if (is.null(span) == FALSE) {
    if (is.numeric(span) == FALSE) {
      stop("\"span\" is not numeric.")
    } else {
      span <- span[1]
      if (span <= 0) {
        stop("\"span\" is not a positive number.")
      }
    }
  } 
  if (is.null(Ylabel)==TRUE) {
    Ylabel <- Y
  } else {
    if (is.character(Ylabel) == FALSE) {
      stop("\"Ylabel\" is not a string.")
    } else {
      Ylabel <- Ylabel[1]
    }   
  } 
  if (is.null(Dlabel)==TRUE) {
    Dlabel <- D   
  } else {
    if (is.character(Dlabel) == FALSE) {
      stop("\"Dlabel\" is not a string.")
    } else {
      Dlabel <- Dlabel[1]
    }   
  }
  if (is.null(Xlabel)==TRUE) {
    Xlabel <- X   
  } else {
    if (is.character(Xlabel) == FALSE) {
      stop("\"Xlabel\" is not a string.")
    } else {
      Xlabel <- Xlabel[1]
    }   
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
  if (is.null(cex.axis)==FALSE) {
    if (is.numeric(cex.axis)==FALSE) {
      stop("\"cex.axis\" is not numeric.")
    }    
  }  
  
  ## title
  if (is.null(main)==FALSE) {
    main <- main[1]
  } 
  ## discrete or continuous treatment
  
  if (!treat.type %in% c("discrete","continuous") ){
    stop("\"treat.type\" must be one of the following: \"discrete\",\"continuous\".")
  }
  
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
	ncols <- length(unique(data[,D]))
  }
  if(treat.type=="continuous"){
	ncols <- nbins
  }
  }
  
  ## load packages
  requireNamespace("ggplot2")
  
  ## drop missing values
  data <- na.omit(data[,c(Y, D, X, Z, FE, weights)])
  n <- dim(data)[1]
  
  #factor cov
  to_dummy_var <- c()
  for(a in Z){
	if(is.factor(data[,a])==TRUE){
		to_dummy_var <- c(to_dummy_var,a)
	}	
	if(is.character(data[,a])==TRUE){
		stop("\"Z\" should be numeric or factorial.")
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
  
  
  if(full==T){
    cat("Use fully moderated model.\n")
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
  
  ## Covariates & FE
  ## formula
  mod.f<-paste0(Y,"~","-1")
  if (is.null(Z)==FALSE) {
    mod.f <- paste0(mod.f, "+", paste0(Z,collapse="+"))
  }
  if (is.null(FE)==FALSE) {
    mod.f <- paste0(mod.f, "|",paste0(FE, collapse = "+"))
  }
  mod.f <- as.formula(mod.f)
  
  
  if (is.null(FE)==TRUE) { #OLS
    if (is.null(weights)==TRUE) {
      mod.demean<-lm(mod.f,data=data)
    } else {
      mod.demean<-lm(mod.f,data=data,weights=weights)
    }
  } else { # FE
    if (is.null(weights)==TRUE) {
      mod.demean<-felm(mod.f,data=data)
    } else {
      mod.demean<-felm(mod.f,data=data,weights=data[,weights])
    }
  }
  data[,Y] <- mod.demean$residuals
  
  ## plotting
  if (treat.type=="discrete") { ## discrte case
    if(length(unique(data[,D]))>9) {
      stop("Too many kinds of treatment arms")
    }
    all_treatment=as.character(unique(data[,D]))
    
    data.aug<-data
    data.aug[,'treat'] <- paste("Treatment =",data[,D])
	order.lab <- paste("Treatment =",order)
    
	if(is.null(subtitle)==TRUE){
		subtitle <- order.lab
	}
	
    if (is.null(pos)==TRUE) {
      box.pos <- min(data[,Y])
    } else {
      box.pos <- pos[1]
    }
    
    data.aug$qt90 <- NA
    data.aug$qt50 <- NA
    data.aug$med <- NA
    for (char in all_treatment) {
      treat <- paste0("Treat",char)
      treat_90 <- paste0(treat,"90")
      treat_50 <- paste0(treat,"50")
      treat_med <- paste0(treat,"med")
      
      assign(treat,which(data.aug[,D]==char))
      assign(treat_90,quantile(data[get(treat),X],c(0.05,0.95)))
      assign(treat_50,quantile(data[get(treat),X],c(0.25,0.75)))
      assign(treat_med,median(data[get(treat),X]))
      
      data.aug$qt90[which(data.aug[,D]==char & 
                            data.aug[,X]>=get(treat_90)[1] &
                            data.aug[,X]<=get(treat_90)[2])] <- box.pos
      
      data.aug$qt50[which(data.aug[,D]==char & 
                            data.aug[,X]>=get(treat_50)[1] &
                            data.aug[,X]<=get(treat_50)[2])] <- box.pos
      
      data.aug$med[get(treat)[which.min(abs(data.aug[get(treat),X]-get(treat_med)))]] <- box.pos
    }
    
    
    ## plotting
    if (is.null(weights)==TRUE) {
      #p1 <- ggplot(data.aug, aes_string(X, Y))
	  p1 <- ggplot(transform(data.aug,treat=factor(treat,levels=order.lab,labels=subtitle)),aes_string(X, Y))
    } else {
      #p1 <- ggplot(data.aug, aes_string(X, Y, weight=weights))
      p1 <- ggplot(transform(data.aug,treat=factor(treat,levels=order.lab,labels=subtitle)),aes_string(X, Y, weight=weights))
    }
    if (theme.bw == TRUE) {
      p1 <- p1 + theme_bw() 
    }
    if (show.grid == FALSE) {
      p1 <- p1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    }
    p1 <- p1 + geom_point() +
      geom_smooth(method = "lm", se = F, fullrange = T,
                  colour = "steelblue", size = 1)
    
    if (is.null(span)==TRUE) {
      p1 <- p1+ geom_smooth(method = "loess", formula = y ~ x,
                            se = F, colour="red") 
    } else {
      p1 <- p1 + geom_smooth(method = "loess", formula = y ~ x,
                             se = F, colour="red",span=span)
    }
    p1 <- p1 + xlab(Xlabel) + ylab(Ylabel)
    
    
    p1 <- p1 + geom_line(aes_string(X,"qt90"),
                         size=1,colour="grey50", na.rm = TRUE)
    p1 <- p1 + geom_line(aes_string(X,"qt50"),
                         size=3,colour="grey50", na.rm = TRUE)
    p1 <- p1 + geom_point(aes_string(X,"med"),size=3,
                          shape=21,fill="white",colour="red", na.rm = TRUE)
    p1 <- p1 + facet_wrap(~treat, ncol=ncols) 
  } else { # continuous case
    
    if (is.numeric(data[,D])==F) {
      stop("\"D\" is not a numeric variable")
    }
    
    ## grouping by X
    if (is.null(cutoffs)==TRUE) {
      cutoff<-quantile(data[,X],probs=seq(0,1,1/nbins))
      while (length(unique(cutoff))!=nbins+1) {
        nbins<-nbins-1
        cutoff<-quantile(data[,X],probs=seq(0,1,1/nbins))
      } 
    } else {
      cutoffs <-cutoffs[which(cutoffs>min(data[,X]) & cutoffs < max(data[,X]))]
      cutoff<- sort(unique(c(min(data[,X]),cutoffs,max(data[,X]))))
    } 
    groupID<-cut(data[,X],breaks=cutoff, labels = FALSE)
    groupID[which(data[,X]==min(data[,X]))]<-1
    
    ## X labels
    groupID2 <- cut(data[,X],breaks=cutoff)
    gp.lab = paste(Xlabel, ": ", levels(groupID2), sep="")
    gp.lab[1] <- paste(Xlabel, ": [", substring(levels(groupID2)[1],2), sep = "")
    nbins <- length(unique(groupID))
    groupID <- factor(groupID, labels=gp.lab)
    data.aug <- data
    data.aug$groupID<-groupID
    
    ## plotting
    if (is.null(weights)==TRUE) {
      p1 <- ggplot(data.aug, aes_string(D, Y))
    } else {
      p1 <- ggplot(data.aug, aes_string(D, Y,weight=weights))
    }
    if (theme.bw == TRUE) {
      p1 <- p1 + theme_bw() 
    }
    if (show.grid == FALSE) {
      p1 <- p1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    }
    p1 <- p1 + geom_point() + 
      geom_smooth(method = "lm", se = F, fullrange = T,
                  colour = "steelblue", size = 1)
    if (is.null(span)==TRUE) {
      p1 <- p1 +
        geom_smooth(method = "loess", formula = y ~ x,
                    se = F, colour="red") 
    } else {
      p1 <- p1 +
        geom_smooth(method = "loess", formula = y ~ x,
                    se = F, colour="red",span=span) 
    }        
    p1 <- p1 + xlab(Dlabel) + ylab(Ylabel) + facet_wrap(.~groupID,ncol=ncols)
    
  }
  
  ## axis labels
  if (is.null(cex.lab)==TRUE) {
    cex.lab <- 15
  } else {
    cex.lab <- 15 * cex.lab
  }
  if (is.null(cex.axis)==TRUE) {
    cex.axis <- 15
  } else {
    cex.axis <- 15 * cex.axis
  }
  p1 <- p1 + theme(axis.text = element_text(size=cex.axis), axis.title = element_text(size=cex.lab))
  
  ## title
  if (is.null(cex.main)==TRUE) {
    cex.main <- 18
  } else {
    cex.main <- 18 * cex.main
  }
  if (is.null(main)==FALSE) {
    p1<-p1 + ggtitle(main) +
      theme(plot.title = element_text(hjust = 0.5, size=cex.main,lineheight=.8, face="bold"))
  } 
  
  return(p1) 
}
