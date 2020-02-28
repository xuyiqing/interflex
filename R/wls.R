######################################
## Weighted Least-Squares
######################################

coefs.new <- function(data,bw,Y,X,D,
                  treat.type,
                  base=NULL,
                  Z=NULL,
                  FE=NULL,
                  weights=NULL, #string
                  X.eval = NULL,
                  bw.adaptive = TRUE,
                  neval = 50,
				  predict = FALSE){
  
  if (!treat.type %in% c("discrete","continuous") ){
    stop("\"treat.type\" must be one of the following: \"discrete\",\"continuous\".")
  }
  
  if(bw.adaptive==T){
    Xdensity <- density(data[,X])
  }
  
  #if treat is discrete
  if (treat.type=='discrete') {
    if(length(unique(data[,D]))>5) {
      warning("More than 5 kinds of treatments")
    }
    
    if(length(unique(data[,D]))>20) {
      stop("Too many kinds of treatments, treatments may be continuous")
    }
    
    data[,D] <- as.character(data[,D])
    
    
    if(is.null(base)==TRUE) {
      base=unique(data[,D])[1]
    }
    else {
      base <- as.character(base)
      if (!base %in% unique(data[,D])){
        stop("\"base\" must be one kind of treatment")
      }
    }
    
    all.treat <- sort(unique(data[,D]))
    other.treat <- sort(all.treat[which(all.treat!=base)])
  }
  
  
  ## evaluation points
  if (is.null(X.eval) == TRUE) {
    X.eval <- seq(min(data[,X]), max(data[,X]), length.out = neval)
  }
  
  ## conditions
  noZ <- is.null(Z)
  noFE <- is.null(FE)
  
  ## survey weights
  n <- dim(data)[1]
  if (is.null(weights)==TRUE) {
    weights <- rep(1, n)
  } else {
    weights <- data[,weights]
  }
  
  ## storage
  if(treat.type=="continuous"){
  result <- matrix(NA, length(X.eval), (5 + length(Z)))
  colnames(result)<-c("x","Intercept","D","deltax","Dx",Z)
  }
  
  if(treat.type=='discrete'){
    treat_length <- length(other.treat)
    result <- matrix(NA,length(X.eval),(1+1+treat_length+1+treat_length+length(Z)))
    
    result_colnames <- c("x","Intercept")
    for (char in other.treat) {
      result_colnames <- c(result_colnames,paste0("D.",char))
    }
    result_colnames <- c(result_colnames,'deltax')
    for (char in other.treat) {
      result_colnames <- c(result_colnames,paste0("DX.",char))
    }
    for (cov in Z) {
      result_colnames <- c(result_colnames,cov)
    }
    colnames(result) <- result_colnames #x intercept D... deltax DX...Z
  }
  
  ## main algorithm
  if (noFE) { # no fixed effects, wls
    
    if (treat.type=="continuous"){
      formula <- paste(Y,"~D+deltax+DX",sep=" ")
      if (is.null(Z)==FALSE) {
        formula <- paste(formula,"+",paste(Z,collapse = "+"))
      }
      formula <- as.formula(formula)
    }
    
    if(treat.type=='discrete') {
      formula <- paste(Y,"~",sep=" ")
      for (char in other.treat) {
        formula <- paste(formula,paste0("D.",char,"+"),sep=" ")
      }
      formula <- paste0(formula,"deltax")
      for (char in other.treat) {
        formula <- paste(formula,paste0("+ DX.",char),sep=" ")
      }
      if (is.null(Z)==FALSE) {
        formula <- paste(formula,"+",paste(Z,collapse = "+"))
      }
      formula <- as.formula(formula)
    }
    
    if(treat.type=='continuous') {
      wls<-function(x, dat, weights){
        data1<-dat
        data1[,'D'] <- data1[,D]
        data1[,'deltax'] <- data1[,X]-x
        data1[,'DX'] <- data1[,'D'] * data1[,'deltax']
        
        if(bw.adaptive==T){
          temp_density <- Xdensity$y[which.min(abs(Xdensity$x-x))]
          bw.adapt <- bw*(1+log(max(Xdensity$y)/temp_density))
          w <- dnorm(data1[,'deltax']/bw.adapt) * weights
          data1$w <- dnorm(data1[,'deltax']/bw.adapt) * weights
        }else{
		  w <- dnorm(data1[,'deltax']/bw) * weights
          data1$w <- dnorm(data1[,'deltax']/bw) * weights
        }
        reg <- lm(formula, data=data1, weights = w)
        result <- c(x, reg$coef)
        result[which(is.na(result))] <- 0 
        return(result)   
      }
    }
    
    if(treat.type=='discrete') {
      wls<-function(x, dat, weights){
        data1 <- dat
        for (char in other.treat){
          data1[,paste0("D.",char)] <- as.numeric(data[,D]==char)
        }
        data1[,'deltax'] <- data1[,X]-x
        for (char in other.treat){
          data1[,paste0("DX.",char)] <- data1[,paste0("D.",char)]*data1[,'deltax']
        }
        
        if(bw.adaptive==T){
          temp_density <- Xdensity$y[which.min(abs(Xdensity$x-x))]
          bw.adapt <- bw*(1+log(max(Xdensity$y)/temp_density))
          w <- dnorm(data1[,'deltax']/bw.adapt) * weights
          data1$w <- dnorm(data1[,'deltax']/bw.adapt) * weights
        }else{
          w <- dnorm(data1[,'deltax']/bw) * weights
          data1$w <- dnorm(data1[,'deltax']/bw) * weights
        }
        
        reg <- lm(formula, data=data1, weights = w)
        result <- c(x, reg$coef) 
        result[which(is.na(result))] <- 0
        return(result)   
      }
    }
    
    dat<-data[, c(Y, D, X, Z)]
    
    for (i in 1:length(X.eval)) {
      result[i,] <- wls(X.eval[i], dat, weights = weights)
    }
    result <- data.frame(result) 
  } 
  else{
    if(treat.type=='continuous'){
      dat <- data[, c(Y, D, X, Z)] 
      dat.FE <- as.matrix(data[,FE],drop = FALSE)
      for (i in 1:length(X.eval)) {
        xx <- dat[,X]-X.eval[i]
        dat1 <- as.matrix(cbind(dat[,Y],dat[,D],
                                xx, dat[,D] * xx, dat[,Z]))

        if(bw.adaptive==T){
          temp_density <- Xdensity$y[which.min(abs(Xdensity$x-X.eval[i]))]
          bw.adapt <- bw*(1+log(max(Xdensity$y)/temp_density))
          w <- dnorm(xx/bw.adapt) * weights
        }else{
          w<-dnorm(xx/bw)* weights
        }
		
		if(predict==TRUE){

        invisible(
			capture.output(
				fastplm_res <- fastplm(data = dat1, FE = dat.FE,
									   weight = w, FEcoefs = 1),type='message'
				)
			)
			estcoef <- fastplm_res$coefficients
			estcoef[which(is.nan(estcoef))] <- 0
			result[i,] <- c(X.eval[i],fastplm_res$mu,estcoef)
		}
		if(predict==FALSE){
			fastplm_res <- fastplm(data = dat1, FE = dat.FE,
                           weight = w, FEcoefs = 0)
			estcoef <- fastplm_res$coefficients
			estcoef[which(is.nan(estcoef))] <- 0
			result[i,] <- c(X.eval[i],0,estcoef)
		}
      }
      result <- data.frame(result)
    }
    
    if(treat.type=='discrete'){
      dat.FE <- as.matrix(data[,FE],drop = FALSE)
      for (i in 1:length(X.eval)) {
        xx <- data[,X]-X.eval[i]
        
        data1 <- as.data.frame(data[,Y])
        
        for (char in other.treat){
          data1[,paste0("D.",char)] <- as.numeric(data[,D]==char)
        }
        
        data1[,'deltax'] <- xx
        
        for (char in other.treat){
          data1[,paste0("DX.",char)] <- data1[,paste0("D.",char)]*data1[,'deltax']
        }
        data1 <- cbind(data1,data[,Z])
        data1 <- as.matrix(data1)
        
        if(bw.adaptive==T){
          temp_density <- Xdensity$y[which.min(abs(Xdensity$x-X.eval[i]))]
          bw.adapt <- bw*(1+log(max(Xdensity$y)/temp_density))
          w <- dnorm(xx/bw.adapt) * weights
        }else{
          w<-dnorm(xx/bw) * weights
        }
        if(predict==T){
            #fastplm_res <- fastplm(data = dat1, FE = dat.FE,weight = w, FEcoefs = 1)
			invisible(
				capture.output(
					fastplm_res <- fastplm(data = data1, FE = dat.FE,
									   weight = w, FEcoefs = 1),type='message'
				)
			)
			estcoef <- fastplm_res$coefficients
			estcoef[which(is.nan(estcoef))] <- 0
			result[i,] <- c(X.eval[i],fastplm_res$mu,estcoef)
		}
		if(predict==F){
			fastplm_res <- fastplm(data = data1, FE = dat.FE,
								   weight = w, FEcoefs = 0)
			estcoef <- fastplm_res$coefficients
			estcoef[which(is.nan(estcoef))] <- 0
			result[i,] <- c(X.eval[i],0,estcoef)
		}
      }
      result <- data.frame(result)
    }
  } 
  return(result)
}

