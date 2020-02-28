crossvalidate.new <- function(data, Y, D, X, CV.method='simple', FE=NULL, treat.type, Z = NULL, weights = NULL, cl = NULL, base=NULL, 
                              X.eval, grid, kfold = 10, parallel = TRUE, seed=NULL, bw.adaptive = FALSE,
                              metric = "MSPE"){
  bw <- NULL
					  
  if(is.null(seed)==FALSE){
    set.seed(seed)
  }
  requireNamespace('lfe')
  cat("Cross-validating bandwidth ... \n")
  
  if(TRUE){ ## Sample
  ## generate K random folds
  ## if CV.method=="cluster" then fold by cluster
  ## if CV.method=="stratify" then conduct stratified cross-validation
  ## if CV.method=="simple" then randomly sample
  n<-dim(data)[1]
  fold<-rep(0,n)

  if(CV.method=='simple'){
	kfold <- min(n,kfold)
    cat("#folds =",kfold)
    cat("\n")
    fold <- c(0:(n-1))%%kfold + 1
    fold<-sample(fold, n, replace = FALSE)
  }
  if(CV.method=='cluster' & is.null(cl)==TRUE){
	stop("\"cl\" is not specified.")
  }
  if(CV.method=='cluster' & is.null(cl)==FALSE){
	clusters<-unique(data[,cl])
    m <- length(clusters)
    kfold <- min(m,kfold)
	cat("Use clustered cross-validation.\n")
    cat("#folds =",kfold)
	cat("\n")
    id.list<-split(1:n,data[,cl])
    cl.fold <- c(0:(m-1))%%kfold + 1
    cl.fold <- sample(cl.fold, m, replace = FALSE)
    for (i in 1:m) {
      id.list[[i]] <- rep(cl.fold[i],length(id.list[[i]]))
    }
    fold <- unlist(id.list)
  }
  
  if(CV.method=='stratify'){
	requireNamespace("caret")
	fold <- createFolds(factor(data[,D]), k = kfold, list = FALSE)
	cat("Use stratified cross-validation.\n")
	cat("#folds =",kfold)
	cat("\n")
  }
  
  }
  
  if(TRUE){ ## TREAT.TYPE
	if(treat.type=='discrete') {
		all.treat <- unique(data[,D])
		if(is.null(base)==TRUE){
			base <- all.treat[1]
		}
    other.treat <- all.treat[which(all.treat!=base)]
  }
  }


  ##FUNCTION: calculate errors in testing set
  getError.new <- function(train, test, bw, Y, X, D, Z, FE, weights, 
                           X.eval, treat.type){
	  if(is.null(weights)==T){
        w <- rep(1,dim(train)[1])
        w2 <- rep(1,dim(test)[1])
      }
      else{
        w <- train[,weights]
        w2 <- test[,weights]
      }
      
	  # FIXED EFFECTS
	  ## demean outcome(y) and estimate fixed effects
	  ## save fixed effects in add_FE
	  ## add fixed effects to predicted values
      add_FE <- rep(0,dim(test)[1])
      if(is.null(FE)==FALSE){
		train_y <- as.matrix(train[,Y])
		train_FE <- as.matrix(train[,FE])
		invisible(
			capture.output(
			fastplm_res <- fastplm(data = train_y,FE = train_FE,weight = w, FEcoefs = 1L),type='message'
			)
		)
		FEvalues <- fastplm_res$FEvalues
		FEnumbers <- dim(fastplm_res$FEvalues)[1]
		FE_coef <- matrix(FEvalues[,3],nrow = FEnumbers, ncol = 1)
		rowname <- c()
		for(i in 1:FEnumbers){
			rowname <- c(rowname,paste0(FE[FEvalues[i,1]+1],'.',FEvalues[i,2]))
		}
		rownames(FE_coef) <- rowname
		train[,Y] <- fastplm_res$residuals
      
		#addictive FE
		add_FE <- matrix(0,nrow = dim(test)[1],ncol = length(FE))
		colnames(add_FE) <- FE
		for(fe in FE){
			add_FE[,fe] <- 0
			fe_name <- paste0(fe,".",test[,fe])
			add_FE[which(fe_name %in% rownames(FE_coef)),fe] <- FE_coef[fe_name[(which(fe_name %in% rownames(FE_coef)))],]
		}
		add_FE <- rowSums(add_FE) 
		add_FE <- add_FE + fastplm_res$mu #intercept
      }
	  
	  if(treat.type=='discrete'){
		#generate dummy variable
		test_d <- test[,c(Y,X)]
		for (char in other.treat) {
			test_d[,paste0("D.",char)] <- as.numeric(test[,D]==char)
		}
      
		#get coef
		coef<-coefs.new(data=train,bw=bw,Y=Y,X=X,D=D,Z=Z,base=base,treat.type = 'discrete',
						weights = weights, X.eval= X.eval,bw.adaptive = bw.adaptive)
		coef[is.na(coef)] <- 0
		num.Z<-length(Z)
		num.treat <- length(other.treat)
		
		esCoef<-function(x){ ##obtain the coefficients for x[i]
			Xnew<-abs(X.eval-x)
			d1<-min(Xnew)     ## distance between x[i] and the first nearest x in training set
			label1<-which.min(Xnew)
			Xnew[label1]<-Inf
			d2<-min(Xnew)     ## distance between x[i] and the second nearest x in training set
			label2<-which.min(Xnew)
			if(d1==0){
				if(is.null(Z)==T){
					func <- coef[label1,c(2:(2+num.treat))] # X.eval (1), intercept (2), d (3), xx (4), d:xx (5), z
				}
				else{
					func <- coef[label1,c(c(2:(2+num.treat)),c((4+2*num.treat):(3+2*num.treat+num.Z)))] # X.eval (1), intercept (2), d (3), xx (4), d:xx (5), z  
				}
			}	  
			else if(d2==0){
				if(is.null(Z)==T){
					func <- coef[label2,c(2:(2+num.treat))] 
				}
				else{
					func <- coef[label2,c(c(2:(2+num.treat)),c((4+2*num.treat):(3+2*num.treat+num.Z)))] 
				}
			} 
			else{ ## weighted average 
				if(is.null(Z)==T){
					func1 <- coef[label1,c(2:(2+num.treat))] 
					func2 <- coef[label2,c(2:(2+num.treat))]
				}
				else{
					func1 <- coef[label1,c(c(2:(2+num.treat)),c((4+2*num.treat):(3+2*num.treat+num.Z)))]
					func2 <- coef[label2,c(c(2:(2+num.treat)),c((4+2*num.treat):(3+2*num.treat+num.Z)))]
				}
					func <- (func1 * d2 + func2 * d1)/(d1 + d2) 
			}
			return(func)
		}
      
		Knn<-t(sapply(test[,X],esCoef)) ## coefficients for test  class==matrix
      
		## predicting
		test.Y <- test[,Y]
		test.X <- as.data.frame(rep(1,dim(test)[1]))
      
		for (char in other.treat) {
			test.X[,paste0("D.",char)] <- test_d[,paste0("D.",char)]
		}
		test.X <- cbind(test.X,as.data.frame(test[,Z]))
		test.X <- as.matrix(test.X)
		sumOfEst<-matrix(lapply(1:length(test.X), function(i){test.X[i]*Knn[[i]]}),
						 nrow=nrow(test.X), ncol=ncol(test.X))
		error <- test.Y - rowSums(matrix(unlist(sumOfEst),length(test.Y))) - add_FE  
		## weights
		error <- error*w2/mean(w2)
		return(c(mean(abs(error)),mean(error^2)))
	  }
	  
	  if(treat.type=='continuous'){
		coef<-coefs.new(data=train,bw=bw,Y=Y,X=X,D=D,Z=Z,treat.type = 'continuous',bw.adaptive = bw.adaptive,
						weights = weights, X.eval= X.eval)
		coef[is.na(coef)] <- 0
		n2<-length(Z)
		esCoef<-function(x){ ##obtain the coefficients for x[i]
			Xnew<-abs(X.eval-x)
			d1<-min(Xnew)     ## distance between x[i] and the first nearest x in training set
			label1<-which.min(Xnew)
			Xnew[label1]<-Inf
			d2<-min(Xnew)     ## distance between x[i] and the second nearest x in training set
			label2<-which.min(Xnew)
			if(d1==0){
				func <- coef[label1,-c(1,4,5)] # X.eval (1), intercept (2), d (3), xx (4), d:xx (5), z
			} else if(d2==0){
				func <- coef[label2,-c(1,4,5)]  
			} else{ ## weighted average 
				func <- (coef[label1,-c(1,4,5)] * d2 + coef[label2,-c(1,4,5)] * d1)/(d1 + d2) 
			}
			return(func)
		} 
		Knn<-t(sapply(test[,X],esCoef)) ## coefficients for test  class==matrix
      
		## predicting
		test.Y <- test[,Y]
		test.X <- as.matrix(cbind(rep(1,dim(test)[1]),test[,c(D,Z)])) 
		sumOfEst<-matrix(lapply(1:length(test.X), function(i){test.X[i]*Knn[[i]]}),
						 nrow=nrow(test.X), ncol=ncol(test.X))
		error<-test.Y - rowSums(matrix(unlist(sumOfEst),length(test.Y)))-add_FE
      
		## weights
		error <- error*w2/mean(w2)
		return(c(mean(abs(error)),mean(error^2)))
    } 
  }
  
  
  ##FUNCTION: calculate MSE
  ## grid search and k fold cross-validation
  cv.new<-function(bw){
    mse<-matrix(NA,kfold,2)
    for(j in 1:kfold){ # K-fold CV
      testid <- which(fold==j)
      train <- data[-testid,]
      test <- data[testid,]
      mse[j,] <- getError.new(train= train, test = test, treat.type = treat.type,FE=FE,
                              bw = bw, Y=Y, X=X, D=D, Z=Z, weights = weights, X.eval= X.eval)
    }
    return(c(bw, apply(mse,2,mean)))
  }
  
  
  ## Parallel Computing or Not
  if (parallel == TRUE) {
    Error<-suppressWarnings(foreach(bw = grid, .combine = rbind,
                                    .export = c("coefs.new","getError.new"),
                                    .inorder = FALSE) %dopar% {cv.new(bw)})
 
  } 
  else {
    Error <- matrix(NA, length(grid), 3)
    for (i in 1:length(grid)) {
      Error[i, ] <- cv.new(grid[i])
      cat(".")
    } 
  } 
  colnames(Error) <- c("bandwidth","MAPE","MSPE")
  rownames(Error) <- NULL
  if (metric=="MAPE") {
    opt.bw <- grid[which.min(Error[,2])]
  } 
  else {
    opt.bw <- grid[which.min(Error[,3])]
  } 
  output <- list(CV.out = round(Error,3),
                   opt.bw = opt.bw)
  cat(paste("Bandwidth =", round(output$opt.bw,3),"\n"))
  return(output)
}

