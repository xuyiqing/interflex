inter.gam<-function(data,Y,D,X,
                    Z=NULL,
					weights=NULL,
                    full.moderate=TRUE,
                    FE=NULL,
                    SE = FALSE,
                    k=10,
                    angle=c(30, 100,-30,-120),
                    Ylabel = NULL, Dlabel = NULL, Xlabel = NULL){
  
  ## in case data is in tibble format
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
  
  if (is.logical(full.moderate) == FALSE & is.numeric(full.moderate)==FALSE) {
    stop("full.moderate is not a logical flag.")
  }else{
    full <- full.moderate
  } 
  
  if (is.null(FE) == FALSE) {
    for (i in 1:length(FE)) {
      if (is.character(FE[i]) == FALSE) {
        stop("Some element in FE is not a string.")
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
  
  if (is.logical(SE) == FALSE & is.numeric(SE)==FALSE) {
    stop("\"SE\" is not a logical flag.")
  }
  if (is.null(k) == FALSE) {
    if (is.numeric(k) == FALSE) {
      stop("\"k\" is not numeric.")
    } else {
      k <- k[1]
      if (k <= 0) {
        stop("\"k\" is not a positive number.")
      }
    }
  }
  if (length(angle)>=5 | length(angle)<1) {
    stop("\"angle\" must have length 1 to 4.")
  } else {
    if (is.numeric(angle)==FALSE) {
      stop("Some element in angle is not numeric.")
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
  
  if (length(unique(data[,D]))==2) {
    warning("The treatment variable is dichotomous.")
  }
  
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
  

  
  ## drop missing values
  data <- na.omit(data[,c(Y, D, X, Z, FE)])
  
  requireNamespace("mgcv") 
  if (is.null(FE)==FALSE) {
    if (is.null(Z)==TRUE) {Z<-c()}
    for (i in 1:length(FE)) {
      Z<-c(Z,paste("as.factor(",FE[i],")",sep=""))
    } 
  }
  
  if (is.null(Z)==TRUE) { # no controls
    formula<-as.formula(paste(Y,"~","s(",D,",",X,",k=",k,")"))
  } else {
    formula<-as.formula(paste(Y,"~","s(",D,",",X,",k=",k,")+",
                              paste(Z,collapse="+")))
  }
  
  
  if (is.null(weights)==TRUE) {
    model<-gam(formula, data=data)
  } else {
    model<-gam(formula, data=data,weights=data[,weights])
  }
  par(mfrow=c(2,2),mar=c(2,2,0,0))
  if (SE==0) {
    for (i in angle) {
      vis.gam(model, view=c(D,X), type="response",cex.lab=1.5,
              xlab = Dlabel, ylab = Xlabel, zlab=Ylabel,
              ticktype="detailed",plot.type="persp",n.grid=40,too.far=0.5,theta=i,phi=20)
    }
  } else {
    for (i in angle) {
      vis.gam(model, view=c(D, X), type="response",cex.lab=1.5,
              xlab = Dlabel, ylab = Xlabel, zlab=Ylabel,
              ticktype="detailed",plot.type="persp",se=2,n.grid=40,too.far=0.5,theta=30,phi=20)
    }
  } 
}