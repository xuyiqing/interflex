#gam
interflex.gam<-function(data,
                    Y,
                    D,
                    X,
                    Z = NULL,
                    Z.ref = NULL,
					weights = NULL,
                    full.moderate = TRUE,
                    FE = NULL,
                    CI = FALSE,
                    method = 'linear',
                    k=10,
                    angle=c(30, 100,-30,-120),
                    Ylabel = NULL,
					Dlabel = NULL,
					Xlabel = NULL,
                    n.grid = 50,
                    file = NULL,
                    scale = 1.1,
                    height = 7,
                    width = 10
                    ){


	n <- dim(data)[1]
    if(is.null(FE)==TRUE){
		use_fe <- 0
	}else{
		use_fe <- 1
	}
    if(is.null(weights)==TRUE){
		w <- rep(1,n)
	}
	else{
		w <- data[,weights]
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
            stop("Some element in \"angle\" is not numeric.")
        }
    }

    
    if(full.moderate==TRUE){
        if(is.null(Z)==TRUE){
            formula <- paste0(Y,"~","s(",D,",",X,",k=",k,")")
        }
        if(is.null(Z)==FALSE){
            formula <- paste0(Y,"~","s(",D,",",X,",k=",k,")")
            for(sub.z in Z){
                formula <- paste0(Y,"~","s(",D,",",X,",k=",k,")")
                for(sub.z in Z){
                    formula <- paste0(formula,"+",sub.z)
                    formula <- paste0(formula,"+","s(I(",X,"*",sub.z,"),k=",k,")")
                }
            }
        }
    }

    if(full.moderate==FALSE){
        if(is.null(Z)==TRUE){
            formula <- paste0(Y,"~","s(",D,",",X,",k=",k,")")
        }
        if(is.null(Z)==FALSE){
            formula <- paste0(Y,"~","s(",D,",",X,",k=",k,")")
            for(sub.z in Z){
                formula <- paste0(formula,"+",sub.z)
            }
        }
    }
    
    if(use_fe==1){
        for (i in 1:length(FE)) {
            formula <- paste0(formula,paste0("as.factor(",FE[i],")"))
        }
    }

    formula <- as.formula(formula)
    requireNamespace("mgcv") 
    #print(formula)
    if(method=='linear'){
        suppressWarnings(
		    model <- try(gam(formula, data=data,weights=w))
	    )
    }
    if(method=='logit'){
        suppressWarnings(
		    model <- try(gam(formula, data=data,weights=w,family=binomial(link = 'logit')))
	    )
    }
    if(method=='probit'){
        suppressWarnings(
		    model <- try(gam(formula, data=data,weights=w,family=binomial(link = 'probit')))
	    )
    }
    if(method=='poisson'){
        suppressWarnings(
		    model <- try(gam(formula, data=data,weights=w,family=poisson(link = "log")))
	    )
    }
    if(method=='nbinom'){
        suppressWarnings(
		    model <- try(gam(formula, data=data,weights=w,family=nb()))
	    )
    }

    if ('try-error' %in% class(model)) {
        stop("A term has fewer unique covariate combinations than specified maximum degrees of freedom, try to specify a smaller gam.k.\n")
    }

    condition <- list()
    if(is.null(Z)==FALSE){
        for(i in 1:length(Z.ref)){
            condition[[names(Z.ref)[i]]] <- Z.ref[i]
        }
    }


    par(mfrow=c(2,2),mar=c(2,2,0,0))
    if (CI==FALSE) {
        for (i in angle) {
            vis.gam(model, view=c(D,X), type="response",cex.lab=1.5,
                    xlab = Dlabel, ylab = Xlabel, zlab=Ylabel,
                    cond = condition,
                    ticktype="detailed",plot.type='persp',n.grid=n.grid,too.far=0.5,theta=i,phi=20)
        }
    } else {
        for (i in angle) {
            vis.gam(model, view=c(D, X), type="response",cex.lab=1.5,
                    xlab = Dlabel, ylab = Xlabel, zlab=Ylabel,
                    cond = condition,
                    ticktype="detailed",plot.type='persp',se=2,n.grid=n.grid,too.far=0.5,theta=i,phi=20)
        }
    }
    
    ## save to file
    ##if (is.null(file)==FALSE) {
	##    ggsave(file, plot = last_plot(),scale = scale,width=width,height = height)
    ##}
    return(model)
}