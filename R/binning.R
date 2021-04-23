interflex.binning <- function(data,
								Y, # outcome
								D, # treatment indicator
								X, # moderator
								treat.info,
								Z = NULL, # covariates
								FE = NULL, # fixed effects
								IV = NULL, # instrument variables
								weights = NULL, # weighting variable
								full.moderate = FALSE, # fully moderated model
								Z.X = NULL, # fully moderated terms
								neval = 50,
								method = "linear", ## "probit"; "logit"; "poisson"; "nbinom"
								nbins = 3,
								cutoffs = NULL,
								vartype = "simu", ## variance type "simu"; "bootstrap"; "delta"
								vcov.type = "robust", ##"homoscedastic"; "robust"; "cluster"; "pcse"
								time = NULL,
								pairwise = TRUE,
								nboots = 200,
								nsimu = 1000,
								parallel = TRUE,
								cores = 4,
								cl = NULL, # variable to be clustered on
								#predict = FALSE,
								Z.ref = NULL, # same length as Z, set the value of Z when estimating marginal effects/predicted value
								wald = TRUE,
								figure = TRUE,
								CI=TRUE,
								bin.labs = TRUE,
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
								show.all = FALSE,
								scale = 1.1,
  								height = 7,
  								width = 10
){	
	WEIGHTS <- NULL
	n <- dim(data)[1]
	treat.type <- treat.info[["treat.type"]]
	if(treat.type=='discrete'){
		other.treat <- treat.info[["other.treat"]]
		other.treat.origin <- names(other.treat)
		names(other.treat.origin) <- other.treat
		all.treat <- treat.info[["all.treat"]]
		all.treat.origin <- names(all.treat)
		names(all.treat.origin) <- all.treat
		base <- treat.info[["base"]]
	}
	if(treat.type=='continuous'){
		D.sample <- treat.info[["D.sample"]]
		label.name <- names(D.sample)
		#names(label.name) <- D.sample
	}
	ncols <- treat.info[["ncols"]]  

	if(is.null(cl)==FALSE){ ## find clusters
		clusters<-unique(data[,cl])
		id.list<-split(1:n,data[,cl])
	}

	if(is.null(FE)==TRUE){
		use_fe <- 0
	}else{
		use_fe <- 1
	}

	#pcse
	if(vcov.type=="pcse"){
		data.cl <- data[,cl]
		data.time <- data[,time]
	}
  
	X.eval <- seq(min(data[,X]), max(data[,X]), length.out = neval)
	
	##----------------------------------------------------------------------
	##---------The Linear Estimator----------------------------------------
	##----------------------------------------------------------------------
	
	formula <- paste0(Y,"~",X)
	if(treat.type=="discrete"){
		for(char in other.treat){
			data[,paste0("D",".",char)] <- as.numeric(data[,D]==char)
			data[,paste0("DX",".",char)] <- as.numeric(data[,D]==char)*data[,X]
			formula <- paste0(formula,"+",paste0("D",".",char),"+",paste0("DX",".",char))
		}
	}
	else{
		data[,"DX"] <- data[,D]*data[,X]
		formula <- paste0(formula,"+",D,"+DX")
	}
  
	if(is.null(Z)==FALSE){
		formula <- paste0(formula,"+",paste0(Z,collapse="+"))
		if(full.moderate==TRUE){
			formula <- paste0(formula,"+",paste0(Z.X,collapse="+"))
		}
		formula.origin <- formula
	}
	else{
		formula.origin <- formula
	}

	if (use_fe==1) {
		formula <- paste0(formula, "|",paste0(FE, collapse = "+"))
		if (vcov.type=="cluster") {
			formula <- paste0(formula, "| 0 |",paste0(cl,collapse = "+"))
		}
	}

	#iv formula
	if(!is.null(IV)){
		if(use_fe==0){
			#mod2 <- ivreg(Y~D+X+D:X+Z1|W+X+W:X+Z1,data=s1)
			formula.iv <- X
			for(sub.iv in IV){
				data[,paste0(X,".",sub.iv)] <- data[,sub.iv]*data[,X]
				formula.iv <- paste0(formula.iv,"+",sub.iv)
				formula.iv <- paste0(formula.iv,"+",paste0(X,".",sub.iv))
			}
			if(is.null(Z)==FALSE){
				formula.iv <- paste0(formula.iv,"+",paste0(Z,collapse="+"))
				if(full.moderate==TRUE){
					formula.iv <- paste0(formula.iv,"+",paste0(Z.X,collapse="+"))
				}
			}
			formula <- paste0(formula,"|",formula.iv)
		}
		if(use_fe==1){
			#mod2 <- felm(Y~X+Z1+Z2|unit+year|(D|DX ~ W+WX)|0,data = s1)
			formula <- paste0(Y,"~",X)
			formula.iv <- ""
			name.update <- c(X)
			if(is.null(Z)==FALSE){
				formula <- paste0(formula,"+",paste0(Z,collapse="+"))
				name.update <- c(name.update,Z)
				if(full.moderate==TRUE){
					formula <- paste0(formula,"+",paste0(Z.X,collapse="+"))
					name.update <- c(name.update,Z.X)
				}
			}
			if(treat.type=="discrete"){
				for(char in other.treat){
					formula.iv <- paste0(formula.iv,"|",paste0("D",".",char),"|",paste0("DX",".",char))
					name.update <- c(name.update,paste0("D",".",char),paste0("DX",".",char))
				}
			}
			else{
				formula.iv <- paste0(formula.iv,"|",D,"|DX")
				name.update <- c(name.update,D,"DX")
			}
			formula.iv <- substring(formula.iv, 2)
			formula.iv <- paste0(formula.iv," ~ ")

			for(sub.iv in IV){
				data[,paste0(X,".",sub.iv)] <- data[,sub.iv]*data[,X]
				formula.iv <- paste0(formula.iv,"+",sub.iv)
				formula.iv <- paste0(formula.iv,"+",paste0(X,".",sub.iv))
			}

			formula <- paste0(formula, "|",paste0(FE, collapse = "+"),"|")
			formula.iv <- paste0("(",formula.iv,")")
			formula <- paste0(formula,formula.iv)
			if (vcov.type=="cluster") {
				formula <- paste0(formula, "|",paste0(cl,collapse = "+"))
			}
		}
	}
	formula <- as.formula(formula)
  
	if(is.null(weights)==TRUE){
		w <- rep(1,n)
	}
	else{
		w <- data[,weights]
	}
	data[['WEIGHTS']] <- w

	## linear part
	if(use_fe==0 & is.null(IV)){
		if(method=='linear'){
			suppressWarnings(
				model <- glm(formula,data=data,weights=WEIGHTS)
			)
		}
		if(method=='logit'){
			suppressWarnings(
				model <- glm(formula,data=data,family=binomial(link = 'logit'),weights=WEIGHTS)
			)
		}
		if(method=='probit'){
			suppressWarnings(
				model <- glm(formula,data=data,family=binomial(link = 'probit'),weights=WEIGHTS)
			)
		}
		if(method=='poisson'){
			suppressWarnings(
				model <- glm(formula,data=data,family=poisson,weights=WEIGHTS)
			)
		}
		if(method=='nbinom'){
			suppressWarnings(
				model <- glm.nb(formula,data=data,weights=WEIGHTS)
			)
		}
		if(model$converged==FALSE){
			stop("Linear estimator can't converge.")
		}
	}

	if(use_fe==TRUE & is.null(IV)){
		suppressWarnings(
			model <- felm(formula,data=data,weights=w)
		)
	}

	if(!is.null(IV)){
		if(use_fe==TRUE){
			suppressWarnings(
					model <- felm(formula,data=data,weights=w)
			)
		}
		if(use_fe==0){
			suppressWarnings(
				model <- ivreg(formula,data=data,weights=WEIGHTS)
			)
		}
	}
	
	model.coef <- coef(model)
	if(!is.null(IV)){
	  if(use_fe==1){
	    names(model.coef) <- name.update
	  }
    }

	if(use_fe==FALSE){
		if(vcov.type=="homoscedastic"){
			model.vcov <- vcov(model)
		}
		if(vcov.type=="robust"){
			model.vcov <- vcovHC(model,type="HC1")
		}
		if(vcov.type=="cluster"){
			model.vcov <- vcovCluster(model,cluster = data[,cl])
		}
		if(vcov.type=="pcse"){
			model.vcov <- pcse(model,pairwise=pairwise,groupN=data.cl,groupT=data.time)$vcov
			rownames(model.vcov)[1] <- "(Intercept)"
			colnames(model.vcov)[1] <- "(Intercept)"
		}
	}

	if(use_fe==TRUE){
		if (vcov.type=="homoscedastic") {
    		model.vcov <- vcov(model, type = "iid")
    	} 
		if (vcov.type=="robust") {
    		model.vcov <- vcov(model, type = "robust") 
    	} 
		if (vcov.type=="cluster") {
    		model.vcov <- vcov(model, type = "cluster") 
    	}
	}
	
	if(!is.null(IV)){
		if(use_fe==1){
	    	rownames(model.vcov) <- colnames(model.vcov) <- name.update
	  	}
    }
	
	model.df <- model$df.residual
	model.coef[which(is.na(model.coef))] <- 0
	model.vcov[which(is.na(model.vcov))] <- 0
	
	model.vcov.all <- matrix(0,nrow=length(model.coef),ncol=length(model.coef))
	colnames(model.vcov.all) <- names(model.coef)
	rownames(model.vcov.all) <- names(model.coef)
	for(a1 in rownames(model.vcov)){
		for(a2 in colnames(model.vcov)){
			model.vcov.all[a1,a2] <- model.vcov[a1,a2]
		}
	}
	if(isSymmetric.matrix(model.vcov.all,tol = 1e-6)==FALSE){
		warning(paste0("Option vcov.type==",vcov.type,"leads to unstable standard error in linear estimator, by default use homoscedastic standard error as an alternative.\n"))
		model.vcov.all <- vcov(model)
		model.vcov.all[which(is.na(model.vcov.all))] <- 0
	}
	model.vcov <- model.vcov.all


	##----------------------------------------------------------------------
	##---------The Binning Estimator----------------------------------------
	##----------------------------------------------------------------------
	
	if(is.null(cutoffs)==TRUE){
		cuts.X<-quantile(data[,X],probs=seq(0,1,1/nbins))
		if (length(unique(cuts.X))!=nbins+1) {
			cuts.X <- unique(cuts.X)
			nbins <- length(cuts.X)-1
		} 
	}else {
		cutoffs <- cutoffs[which(cutoffs>min(data[,X]) & cutoffs < max(data[,X]))]
		cuts.X<- sort(unique(c(min(data[,X]),cutoffs,max(data[,X]))))
	} 
	groupX <- cut(data[,X],breaks=cuts.X, labels = FALSE)
	groupX[which(data[,X]==min(data[,X]))] <- 1
	nbins <- length(unique(groupX))
  
	## X labels
	groupX2 <- cut(data[,X],breaks=cuts.X)
	gp.lab <- paste(Xlabel, ":", levels(groupX2), sep="")
	gp.lab[1] <- paste(Xlabel, ":[", substring(levels(groupX2)[1],2), sep = "")

	## mid points
	x0<-rep(NA,nbins)
	for (i in 1:nbins) x0[i] <- median(data[which(groupX==i),X], na.rm=TRUE)

	## binning formula
	formula.binning <- paste0(Y,"~-1")
	data.touse <- data
	for(i in 1:nbins){
		data.touse[,paste0('G.',i)] <- as.numeric(groupX==i)
		data.touse[,paste0('GX.',i)] <- as.numeric(groupX==i)*(data.touse[,X]-x0[i])
		formula.binning <- paste0(formula.binning,'+',paste0('G.',i),'+',paste0('GX.',i))
	}
	
	if(treat.type=='discrete'){
		for(char in other.treat){
			for(i in 1:nbins){
				data.touse[,paste0('D.',char,'.','G.',i)] <- as.numeric(groupX==i)*as.numeric(data.touse[,D]==char)
				data.touse[,paste0('D.',char,'.','GX.',i)] <- as.numeric(groupX==i)*as.numeric(data.touse[,D]==char)*(data.touse[,X]-x0[i])
				formula.binning <- paste0(formula.binning,'+',paste0('D.',char,'.','G.',i),'+',paste0('D.',char,'.','GX.',i))
			}
		}
	}
	
	if(treat.type=='continuous'){
		for(i in 1:nbins){
			data.touse[,paste0('D.G.',i)] <- as.numeric(groupX==i)*data.touse[,D]
			data.touse[,paste0('D.GX.',i)] <- as.numeric(groupX==i)*data.touse[,D]*(data.touse[,X]-x0[i])
			formula.binning <- paste0(formula.binning,'+',paste0('D.G.',i),'+',paste0('D.GX.',i))
		}
	}
	
	if(is.null(Z)==FALSE){
		if(full.moderate==FALSE){
			formula.binning <- paste0(formula.binning,"+",paste0(Z,collapse="+"))
		}
		if(full.moderate==TRUE){
			for(a in Z){
				for(i in 1:nbins){
					data.touse[,paste0(a,'.G.',i)] <- as.numeric(groupX==i)*data.touse[,a]
					data.touse[,paste0(a,'.GX.',i)] <- as.numeric(groupX==i)*data.touse[,a]*(data.touse[,X]-x0[i])
					formula.binning <- paste0(formula.binning,'+',paste0(a,'.G.',i),'+',paste0(a,'.GX.',i))
				}
			}
		}
	}

	if (use_fe==1) {
		formula.binning <- paste0(formula.binning, "|",paste0(FE, collapse = "+"))
		if (vcov.type=="cluster") {
			formula.binning <- paste0(formula.binning, "| 0 |",paste0(cl,collapse = "+"))
		}
	}

	## iv regression
	if(!is.null(IV) & use_fe==FALSE){
		# use ivreg
		#Y ~ -1 + G.1 + GX.1 + G.2 + GX.2 + G.3 + GX.3 + D*G.1 + 
        #D*GX.1 + D*G.2 + D*GX.2 + D*G.3 + D*GX.3|
		#-1 + G.1 + GX.1 + G.2 + GX.2 + G.3 + GX.3 + W*G.1 + 
        #W*GX.1 + W*G.2 + W*GX.2 + W*G.3 + W*GX.3
		formula.binning.iv <- "-1"
		for(i in 1:nbins){
			formula.binning.iv <- paste0(formula.binning.iv,'+',paste0('G.',i),'+',paste0('GX.',i))
		}
		for(sub.iv in IV){
			for(i in 1:nbins){
				data.touse[,paste0(sub.iv,'.','G.',i)] <- as.numeric(groupX==i)*as.numeric(data.touse[,sub.iv])
				data.touse[,paste0(sub.iv,'.','GX.',i)] <- as.numeric(groupX==i)*as.numeric(data.touse[,sub.iv])*(data.touse[,X]-x0[i])
				formula.binning.iv <- paste0(formula.binning.iv,'+',paste0(sub.iv,'.','G.',i),'+',paste0(sub.iv,'.','GX.',i))
			}
		}
		if(is.null(Z)==FALSE){
			if(full.moderate==FALSE){
				formula.binning.iv <- paste0(formula.binning.iv,"+",paste0(Z,collapse="+"))
			}
			if(full.moderate==TRUE){
				for(a in Z){
					for(i in 1:nbins){
						formula.binning.iv <- paste0(formula.binning.iv,'+',paste0(a,'.G.',i),'+',paste0(a,'.GX.',i))
					}
				}
			}
		}
		formula.binning <- paste0(formula.binning,"|",formula.binning.iv)
	}

	if(!is.null(IV) & use_fe==TRUE){
		# use felm
		# Y ~ -1 + G.1 + GX.1 + G.2 + GX.2 + G.3 + GX.3 + Z1 + Z2
		# |unit+year
		# |(D*G.1|D*GX.1|D*G.2|D*GX.2|D*G.3|D*GX.3 ~ 
		# W*G.1+W*GX.1+W*G.2+W*GX.2+W*G.3+W*GX.3)|cluster

		formula.binning.part1 <- "Y ~ -1"
		formula.binning.part2 <- ""
		iv.reg.name <- c()
		for(i in 1:nbins){
			formula.binning.part1 <- paste0(formula.binning.part1,'+',paste0('G.',i),'+',paste0('GX.',i))
			iv.reg.name <- c(iv.reg.name,paste0('G.',i),paste0('GX.',i))
		}
		if(is.null(Z)==FALSE){
			if(full.moderate==FALSE){
				formula.binning.part1 <- paste0(formula.binning.part1,"+",paste0(Z,collapse="+"))
				iv.reg.name <- c(iv.reg.name,Z)
			}
			if(full.moderate==TRUE){
				for(a in Z){
					for(i in 1:nbins){
						formula.binning.part1 <- paste0(formula.binning.part1,'+',paste0(a,'.G.',i),'+',paste0(a,'.GX.',i))
						iv.reg.name <- c(iv.reg.name,paste0(a,'.G.',i),paste0(a,'.GX.',i))
					}
				}
			}
		}

		if(treat.type=='discrete'){
			for(char in other.treat){
				for(i in 1:nbins){
					formula.binning.part2 <- paste0(formula.binning.part2,'|',paste0('D.',char,'.','G.',i),'|',paste0('D.',char,'.','GX.',i))
					iv.reg.name <- c(iv.reg.name,paste0('D.',char,'.','G.',i),paste0('D.',char,'.','GX.',i))
				}
			}
		}

		if(treat.type=='continuous'){
			for(i in 1:nbins){
				formula.binning.part2 <- paste0(formula.binning.part2,'|',paste0('D.G.',i),'|',paste0('D.GX.',i))
				iv.reg.name <- c(iv.reg.name,paste0('D.G.',i),paste0('D.GX.',i))
			}
		}

		formula.binning.part2 <- paste0(formula.binning.part2,"~")

		for(sub.iv in IV){
			for(i in 1:nbins){
				data.touse[,paste0(sub.iv,'.','G.',i)] <- as.numeric(groupX==i)*as.numeric(data.touse[,sub.iv])
				data.touse[,paste0(sub.iv,'.','GX.',i)] <- as.numeric(groupX==i)*as.numeric(data.touse[,sub.iv])*(data.touse[,X]-x0[i])
				formula.binning.part2 <- paste0(formula.binning.part2,'+',paste0(sub.iv,'.','G.',i),'+',paste0(sub.iv,'.','GX.',i))
			}
		}

		formula.binning.part2 <- substring(formula.binning.part2, 2)
		formula.binning.part1 <- paste0(formula.binning.part1, "|",paste0(FE, collapse = "+"),"|")
		formula.binning.part2 <- paste0("(",formula.binning.part2,")")
		formula.binning <- paste0(formula.binning.part1,formula.binning.part2)
		if (vcov.type=="cluster") {
			formula.binning <- paste0(formula.binning, "|",paste0(cl,collapse = "+"))
		}
	}

	## binning fit
	formula.binning <- as.formula(formula.binning)

	if(use_fe==FALSE & is.null(IV)){
		if(method=='linear'){
			suppressWarnings(
				model.binning <- glm(formula.binning,data=data.touse,weights=WEIGHTS)
			)
		}
		if(method=='logit'){
			suppressWarnings(
				model.binning <- glm(formula.binning,data=data.touse,family=binomial(link = 'logit'),weights=WEIGHTS)
			)
		}
		if(method=='probit'){
			suppressWarnings(
				model.binning <- glm(formula.binning,data=data.touse,family=binomial(link = 'probit'),weights=WEIGHTS)
			)
		}
		if(method=='poisson'){
			suppressWarnings(
				model.binning <- glm(formula.binning,data=data.touse,family=poisson,weights=WEIGHTS)
			)
		}
		if(method=='nbinom'){
			suppressWarnings(
				model.binning <- glm.nb(formula.binning,data=data.touse,weights=WEIGHTS)
			)
		}
		
		if(model.binning$converged==FALSE){
			stop("Binning estimator can't converge.")
		}
	}


	if(use_fe==FALSE & !is.null(IV)){
		suppressWarnings(
			model.binning <- ivreg(formula.binning,data=data.touse,weights=WEIGHTS)
		)
	}

	if(use_fe==TRUE){
		w <- data.touse[,'WEIGHTS']
		suppressWarnings(
			model.binning <- felm(formula.binning,data=data.touse,weights=w)
		)
	}
	
	model.binning.df <- model.binning$df.residual
	model.binning.coef <- coef(model.binning) #keep the NaN(s) in coefficient 
	if(!is.null(IV) & use_fe==TRUE){
		names(model.binning.coef) <- iv.reg.name
	}

	if(use_fe==FALSE){
		if(vcov.type=="homoscedastic"){
			model.binning.vcov <- vcov(model.binning)
		}
		if(vcov.type=="robust"){
			model.binning.vcov <- vcovHC(model.binning,type="HC1")
		} 
		if(vcov.type=="cluster"){
			model.binning.vcov <- vcovCluster(model.binning,cluster = data[,cl])
		}
		if(vcov.type=="pcse"){
			model.binning.vcov <- pcse(model.binning,pairwise=pairwise,groupN=data.cl,groupT=data.time)$vcov
		}
	}

	if(use_fe==TRUE){
		if (vcov.type=="homoscedastic") {
    		model.binning.vcov <- vcov(model.binning, type = "iid")
    	} 
		if (vcov.type=="robust") {
    		model.binning.vcov <- vcov(model.binning, type = "robust") 
    	} 
		if (vcov.type=="cluster") {
    		model.binning.vcov <- vcov(model.binning, type = "cluster") 
    	}
	}
	if(!is.null(IV) & use_fe==TRUE){
		rownames(model.binning.vcov) <- colnames(model.binning.vcov) <- iv.reg.name
	}
	
	model.binning.vcov[which(is.na(model.binning.vcov))] <- 0
	model.binning.vcov.all <- matrix(0,nrow=length(model.binning.coef),ncol=length(model.binning.coef))
	colnames(model.binning.vcov.all) <- names(model.binning.coef)
	rownames(model.binning.vcov.all) <- names(model.binning.coef)
	for(a1 in rownames(model.binning.vcov)){
		for(a2 in colnames(model.binning.vcov)){
			model.binning.vcov.all[a1,a2] <- model.binning.vcov[a1,a2]
		}
	}
	if(isSymmetric.matrix(model.binning.vcov.all,tol = 1e-6)==FALSE){
		warning(paste0("Option vcov.type==",vcov.type,"leads to unstable standard error in binning estimator, by default use homoscedastic standard error as an alternative.\n"))
		#return(model.binning.vcov.all)
		model.binning.vcov.all <- vcov(model.binning)
		if(!is.null(IV) & use_fe==TRUE){
			rownames(model.binning.vcov.all) <- colnames(model.binning.vcov.all) <- iv.reg.name
		}
		model.binning.vcov.all[which(is.na(model.binning.vcov.all))] <- 0
	}
	model.binning.vcov <- model.binning.vcov.all


	##Function A.1 (linear)
	#1, estimate treatment effects/marginal effects given model.coef
	#2, input: model.coef; char(discrete)/D.ref(continuous);
	#3, output: marginal effects/treatment effects
	gen.general.TE <- function(model.coef,char=NULL,D.ref=NULL){
	if(is.null(char)==TRUE){
		treat.type='continuous'
	}
	if(is.null(D.ref)==TRUE){
		treat.type='discrete'
	}
	
	gen.TE <- function(model.coef,X.eval){
		neval <- length(X.eval)
		if(treat.type=='discrete'){
			link.1 <- model.coef['(Intercept)'] + X.eval*model.coef[X] + 1*model.coef[paste0('D.',char)] + X.eval*model.coef[paste0('DX.',char)]
			link.0 <- model.coef['(Intercept)'] + X.eval*model.coef[X]
			if(is.null(Z)==FALSE){
				for(a in Z){
					target.Z <- Z.ref[a]
					link.1 <- link.1 + target.Z*model.coef[a]
					link.0 <- link.0 +target.Z*model.coef[a]
					if(full.moderate==TRUE){
						link.1 <- link.1 + target.Z*model.coef[paste0(a,'.X')]*X.eval
						link.0 <- link.0 + target.Z*model.coef[paste0(a,'.X')]*X.eval
					}
				}
			}
			if(method=='linear'){
				TE <- link.1-link.0
				E.pred <- link.1
				E.base <- link.0
			}
		
			if(method=='logit'){
				E.pred <- E.prob.1 <- exp(link.1)/(1+exp(link.1))
				E.base <- E.prob.0 <- exp(link.0)/(1+exp(link.0))
				TE <- E.prob.1-E.prob.0
			}
	
			if(method=='probit'){
				E.pred <- E.prob.1 <- pnorm(link.1,0,1)
				E.base <- E.prob.0 <- pnorm(link.0,0,1)
				TE <- E.prob.1-E.prob.0
			}
		
			if(method=='poisson' | method=="nbinom"){
				E.pred <- E.y.1 <- exp(link.1)
				E.base <- E.y.0 <- exp(link.0)
				TE <- E.y.1-E.y.0
			}
				names(TE) <- rep(paste0("TE.",char),neval)
				names(E.pred) <- rep(paste0("Predict.",char),neval)
				names(E.base) <- rep(paste0("Predict.",base),neval)
				gen.TE.output <- list(TE=TE,E.pred=E.pred,E.base=E.base)
		}

		if(treat.type=='continuous'){
			link <- model.coef["(Intercept)"] + X.eval*model.coef[X] + model.coef[D]*D.ref + model.coef["DX"]*X.eval*D.ref
			if(is.null(Z)==FALSE){
				for(a in Z){
					target.Z <- Z.ref[a]
					link <- link + target.Z*model.coef[a]
					if(full.moderate==TRUE){
						link <- link + target.Z*model.coef[paste0(a,'.X')]*X.eval
					}
				}
			}
			if(method=='logit'){
				ME <- exp(link)/(1+exp(link))^2*(model.coef[D]+model.coef["DX"]*X.eval)
				E.pred <- exp(link)/(1+exp(link))
			}
			if(method=='probit'){
				ME <- (model.coef[D]+model.coef["DX"]*X.eval)*dnorm(link)
				E.pred <- pnorm(link,0,1)
			}
			if(method=='linear'){
				ME <- model.coef[D]+model.coef["DX"]*X.eval
				E.pred <- link
			}
			if(method=='poisson'|method=='nbinom'){
				ME <- exp(link)*(model.coef[D]+model.coef["DX"]*X.eval)
				E.pred <- exp(link)
			}
				names(ME) <- rep(paste0("ME.",names(D.sample)[D.sample == D.ref]),neval)
				names(E.pred) <- rep(paste0("Predict.",names(D.sample)[D.sample == D.ref]),neval)
				gen.TE.output <- list(ME=ME,E.pred=E.pred)
		}
		return(gen.TE.output)
	}

	gen.TE.fe <- function(model.coef,X.eval){
		neval <- length(X.eval)
		if(treat.type=='discrete'){
			TE <- model.coef[paste0('D.',char)] + X.eval*model.coef[paste0('DX.',char)]
			names(TE) <- rep(paste0("TE.",char),neval)
			E.pred <- rep(0,neval) #doesn't return prediction value for fixed effects
			E.base <- rep(0,neval)
			gen.TE.output <- list(TE=TE,E.pred=E.pred,E.base=E.base)
		}
		if(treat.type=='continuous'){
			ME <- model.coef[D] + model.coef["DX"]*X.eval
			names(ME) <- rep(paste0("ME.",names(D.sample)[D.sample == D.ref]),neval)
			E.pred <- rep(0,neval)
			gen.TE.output <- list(ME=ME,E.pred=E.pred)
		}
		return(gen.TE.output)
	}

	if(use_fe==0){
		gen.TE.output <- gen.TE(model.coef=model.coef,X.eval=X.eval)
	}
	if(use_fe==1){
		gen.TE.output <- gen.TE.fe(model.coef=model.coef,X.eval=X.eval)
	}
		
	if(treat.type=='discrete'){
		return(list(TE=gen.TE.output$TE,
					E.pred=gen.TE.output$E.pred,
					E.base=gen.TE.output$E.base))
	}
	
	if(treat.type=='continuous'){
		return(list(ME=gen.TE.output$ME,
					E.pred=gen.TE.output$E.pred))
	}
  }
  

	#  #Function A.2 (linear delta)
	#1, estimate sd of treatment effects/marginal effects using delta method
	#2, input: model.coef; model.vcov; char(discrete)/D.ref(continuous);
	#3, output: sd of TE/ME;
	gen.delta.TE <- function(model.coef,model.vcov,char=NULL,D.ref=NULL){
	if(is.null(char)==TRUE){
		treat.type <- 'continuous'
		flag <- 1
	}
	if(is.null(D.ref)==TRUE){
		treat.type <- 'discrete'
		if(char==base){
			flag=0
		}else{
			flag=1
		}
	}
	
	#sd for TE/ME
	#sd for TE/ME
	gen.sd.fe <- function(x,to.diff=FALSE){
		if(treat.type=='discrete'){
			target.slice <- c(paste0('D.',char),paste0('DX.',char))
			vec.1 <- c(1,x)
			vec.0 <- c(0,0)
			vec <- vec.1-vec.0
			temp.vcov.matrix <- model.vcov[target.slice,target.slice]	
			if(to.diff==TRUE){
				return(list(vec=vec,temp.vcov.matrix=temp.vcov.matrix))
			}		
			delta.sd <- sqrt((t(vec)%*%temp.vcov.matrix%*%vec)[1,1])
			return(delta.sd)
		}
		
		if(treat.type=='continuous'){
			target.slice <- c(D,'DX')
			vec <- c(1,x)
			temp.vcov.matrix <- model.vcov[target.slice,target.slice]
			if(to.diff==TRUE){
				return(list(vec=vec,temp.vcov.matrix=temp.vcov.matrix))
			}
					
			delta.sd <- sqrt((t(vec)%*%temp.vcov.matrix%*%vec)[1,1])
			return(delta.sd)
		}
	}

	gen.sd <- function(x,to.diff=FALSE){
		if(use_fe==TRUE){
			return(gen.sd.fe(x=x,to.diff=to.diff))
		}

		if(treat.type=='discrete'){
			link.1 <- model.coef['(Intercept)'] + x*model.coef[X] + 1*model.coef[paste0('D.',char)] + x*model.coef[paste0('DX.',char)]
			link.0 <- model.coef['(Intercept)'] + x*model.coef[X]
			if(is.null(Z)==FALSE){
				for(a in Z){
						target.Z <- Z.ref[a]
						link.1 <- link.1 + target.Z*model.coef[a]
						link.0 <- link.0 +target.Z*model.coef[a]
						if(full.moderate==TRUE){
							link.1 <- link.1 + target.Z*model.coef[paste0(a,'.X')]*x
							link.0 <- link.0 + target.Z*model.coef[paste0(a,'.X')]*x
						}
					}
			}
			if(is.null(Z)==FALSE){
				if(full.moderate==FALSE){
					vec.1 <- c(1,x,1,x,Z.ref)
					vec.0 <- c(1,x,0,0,Z.ref)
					target.slice <- c('(Intercept)',X,paste0('D.',char),paste0('DX.',char),Z)
				}
				if(full.moderate==TRUE){
					vec.1 <- c(1,x,1,x,Z.ref,x*Z.ref)
					vec.0 <- c(1,x,0,0,Z.ref,x*Z.ref)
					target.slice <- c('(Intercept)',X,paste0('D.',char),paste0('DX.',char),Z,Z.X)
				}
			}
			else{
				vec.1 <- c(1,x,1,x)
				vec.0 <- c(1,x,0,0)
				target.slice <- c('(Intercept)',X,paste0('D.',char),paste0('DX.',char))
			}		
			temp.vcov.matrix <- model.vcov[target.slice,target.slice]
			if(method=='logit'){
				vec <- vec.1*exp(link.1)/(1+exp(link.1))^2 - vec.0*exp(link.0)/(1+exp(link.0))^2
			}
			if(method=='probit'){
				vec <- vec.1*dnorm(link.1) - vec.0*dnorm(link.0)
			}
			if(method=='poisson' | method=='nbinom'){
				vec <- vec.1*exp(link.1) - vec.0*exp(link.0)
			}
			if(method=='linear'){
				vec <- vec.1-vec.0
			}		
			if(to.diff==TRUE){
				return(list(vec=vec,temp.vcov.matrix=temp.vcov.matrix))
			}
					
			delta.sd <- sqrt((t(vec)%*%temp.vcov.matrix%*%vec)[1,1])
			return(delta.sd)
		}
		
		if(treat.type=='continuous'){
			link <- model.coef["(Intercept)"] + x*model.coef[X] + model.coef[D]*D.ref + model.coef["DX"]*x*D.ref
			if(is.null(Z)==FALSE){
				for(a in Z){
					target.Z <- Z.ref[a]
					link <- link + target.Z*model.coef[a]
					if(full.moderate==TRUE){
						link <- link + x*target.Z*model.coef[paste0(a,'.X')]
					}
				}
			}
				
			if(is.null(Z)==FALSE){
				if(full.moderate==FALSE){
					vec1 <- c(1,x,D.ref,D.ref*x,Z.ref)
					vec0 <- c(0,0,1,x,rep(0,length(Z)))
					target.slice <- c('(Intercept)',X,D,'DX',Z)
				}
				if(full.moderate==TRUE){
					vec1 <- c(1,x,D.ref,D.ref*x,Z.ref,Z.ref*x)
					vec0 <- c(0,0,1,x,rep(0,2*length(Z)))
					target.slice <- c('(Intercept)',X,D,'DX',Z,Z.X)
				}
			}
			else{
				vec1 <- c(1,x,D.ref,D.ref*x)
				vec0 <- c(0,0,1,x)
				target.slice <- c('(Intercept)',X,D,'DX')
			}
			temp.vcov.matrix <- model.vcov[target.slice,target.slice]

			if(method=='logit'){
				vec <- -(model.coef[D]+x*model.coef['DX'])*(exp(link)-exp(-link))/(2+exp(link)+exp(-link))^2*vec1 + exp(link)/(1+exp(link))^2*vec0
			}
			if(method=='probit'){
				vec <- dnorm(link)*vec0-(model.coef[D]+x*model.coef['DX'])*link*dnorm(link)*vec1
			}
			if(method=='poisson'|method=='nbinom'){
				vec <- (model.coef[D]+x*model.coef['DX'])*exp(link)*vec1+exp(link)*vec0
			}
			if(method=='linear'){
				vec <- vec0
			}
					
			if(to.diff==TRUE){
				return(list(vec=vec,temp.vcov.matrix=temp.vcov.matrix))
			}
					
			delta.sd <- sqrt((t(vec)%*%temp.vcov.matrix%*%vec)[1,1])
			return(delta.sd)
		}
		}
	
	if(flag==1){
		TE.sd <- c(sapply(X.eval,function(x) gen.sd(x)))
	}
	else{
		TE.sd <- NULL
	}
	
	if(treat.type=='discrete' & is.null(TE.sd)==FALSE){
		names(TE.sd) <- rep(paste0("sd.",char),neval)
	}		
	
	if(treat.type=='continuous'){
		names(TE.sd) <- rep(paste0("sd.",names(D.sample)[D.sample == D.ref]),neval)
	}
	
	#output
	if(treat.type=='discrete'){
		return(list(TE.sd=TE.sd))
	}
	if(treat.type=='continuous'){
		return(list(ME.sd=TE.sd))
	}
		
}


	##Function B.1 (bining estimator)
	#1. estimate treatment effects/marginal effects at different binning points
	#2. input: model.binning.coef(can have NaN); char(discrete)/D.ref(continuous)
	#3. output: marginal effects/treatments effects

	gen.binning.TE.fe <- function(model.binning.coef,char=NULL,D.ref=NULL){
		if(treat.type=='discrete'){
			binning.output <- c()
			for(i in 1:nbins){
				TE.bin <- model.binning.coef[paste0('D.',char,'.','G.',i)]
				TE.bin <- c(TE.bin)
				names(TE.bin) <- paste0("G.",i)
				binning.output <- c(binning.output,TE.bin) #can have NaN
			}
		}

		if(treat.type=='continuous'){
			binning.output <- c()
			for(i in 1:nbins){
				ME.bin <- model.binning.coef[paste0('D.G.',i)]
				ME.bin <- c(ME.bin)
				names(ME.bin) <- paste0("G.",i)
				binning.output <- c(binning.output,ME.bin) #can have NaN
			}
		}
		return(binning.output)
	}


	gen.binning.TE <- function(model.binning.coef,char=NULL,D.ref=NULL){
		if(is.null(char)==TRUE){
			treat.type <- 'continuous'
		}
		if(is.null(D.ref)==TRUE){
			treat.type <- 'discrete'
		}

		if(use_fe==TRUE){
			return(gen.binning.TE.fe(model.binning.coef=model.binning.coef,char=char,D.ref=D.ref))
		}
				
		if(treat.type=='discrete'){
			binning.output <- c()
			for(i in 1:nbins){
				link.bin.1 <- model.binning.coef[paste0("G.",i)] + model.binning.coef[paste0('D.',char,'.','G.',i)]
				link.bin.0 <- model.binning.coef[paste0("G.",i)]
				if(is.null(Z)==FALSE){
					for(a in Z){
						target.Z <- Z.ref[a]
						if(full.moderate==FALSE){
							if(is.na(model.binning.coef[a])==TRUE){
								model.binning.coef[a] <- 0
							}
							link.bin.1 <- link.bin.1 + target.Z*model.binning.coef[a]
							link.bin.0 <- link.bin.0 + target.Z*model.binning.coef[a]
						}
						if(full.moderate==TRUE){
							if(is.na(model.binning.coef[paste0(a,".G.",i)])==TRUE){
								model.binning.coef[paste0(a,".G.",i)] <- 0
							}
							link.bin.1 <- link.bin.1 + target.Z*model.binning.coef[paste0(a,".G.",i)]
							link.bin.0 <- link.bin.0 + target.Z*model.binning.coef[paste0(a,".G.",i)]
						}
					}
				}
			
				if(method=='linear'){
					TE.bin <- link.bin.1-link.bin.0
				}
				
				if(method=='logit'){
					E.pred <- E.prob.1 <- exp(link.bin.1)/(1+exp(link.bin.1))
					E.base <- E.prob.0 <- exp(link.bin.0)/(1+exp(link.bin.0))
					TE.bin <- E.prob.1-E.prob.0
				}
	
				if(method=='probit'){
					E.pred <- E.prob.1 <- pnorm(link.bin.1,0,1)
					E.base <- E.prob.0 <- pnorm(link.bin.0,0,1)
					TE.bin <- E.prob.1-E.prob.0
				}
		
				if(method=='poisson' | method=="nbinom"){
					E.pred <- E.y.1 <- exp(link.bin.1)
					E.base <- E.y.0 <- exp(link.bin.0)
					TE.bin <- E.y.1-E.y.0
				}
				TE.bin <- c(TE.bin)
				names(TE.bin) <- paste0("G.",i)
				binning.output <- c(binning.output,TE.bin) #can have NaN
			}
		}
		
		if(treat.type=='continuous'){
			binning.output <- c()
			for(i in 1:nbins){
				link.bin <- model.binning.coef[paste0("G.",i)] + model.binning.coef[paste0('D.G.',i)]*D.ref
				if(is.null(Z)==FALSE){
					for(a in Z){
						target.Z <- Z.ref[a]
						if(full.moderate==FALSE){
							if(is.na(model.binning.coef[a])==TRUE){
								model.binning.coef[a] <- 0
							}
							link.bin <- link.bin + target.Z*model.binning.coef[a]
						}
						if(full.moderate==TRUE){
							if(is.na(model.binning.coef[paste0(a,".G.",i)])==TRUE){
								model.binning.coef[paste0(a,".G.",i)] <- 0
							}
							link.bin <- link.bin + target.Z*model.binning.coef[paste0(a,".G.",i)]
						}
					}
				}
			
				if(method=='logit'){
					ME.bin <- exp(link.bin)/(1+exp(link.bin))^2*model.binning.coef[paste0('D.G.',i)]
				}
				if(method=='probit'){
					ME.bin <- model.binning.coef[paste0('D.G.',i)]*dnorm(link.bin)
				}
				if(method=='linear'){
					ME.bin <- model.binning.coef[paste0('D.G.',i)]
				}
				if(method=='poisson'|method=='nbinom'){
					ME.bin <- exp(link.bin)*model.binning.coef[paste0('D.G.',i)]
				}
				ME.bin <- c(ME.bin)
				names(ME.bin) <- paste0("G.",i)
				binning.output <- c(binning.output,ME.bin) #can have NaN
			}
		}
		return(binning.output)
	}
	
	
	
	##  Function B.2 (binning estimator delta)
	#1. estimate sd of binning treatment effects/marginal effects using delta method
	#2. input: model.binning.coef; model.binning.vcov; char(discrete)/D.ref(continuous)
	#3. output: sd of binning treatment effects/marginal effects using delta method
	
	gen.binning.delta.TE.fe <- function(model.binning.coef, model.binning.vcov, char=NULL, D.ref=NULL){
		model.binning.coef[which(is.na(model.binning.coef))] <- 0
		binning.sd.output <- c()
		if(treat.type=='discrete'){
			for(i in 1:nbins){
				vec.1 <- c(1,1)
				vec.0 <- c(1,0)
				vec <- vec.1-vec.0
				target.slice <- c(paste0("G.",i),paste0('D.',char,'.','G.',i))
				vcov.binning.temp <- model.binning.vcov[target.slice,target.slice]
				delta.sd.bin <- sqrt((t(vec)%*%vcov.binning.temp%*%vec)[1,1])
				delta.sd.bin <- c(delta.sd.bin)
				names(delta.sd.bin) <- paste0("sd.G.",i)
				binning.sd.output <- c(binning.sd.output,delta.sd.bin)
			}
		}
		if(treat.type=='continuous'){
			for(i in 1:nbins){
				vec <- vec0 <- c(0,1)
				target.slice <- c(paste0("G.",i),paste0('D.G.',i))
				vcov.binning.temp <- model.binning.vcov[target.slice,target.slice]
				delta.sd.bin <- sqrt((t(vec)%*%vcov.binning.temp%*%vec)[1,1])
				delta.sd.bin <- c(delta.sd.bin)
				names(delta.sd.bin) <- paste0("sd.G.",i)
				binning.sd.output <- c(binning.sd.output,delta.sd.bin)
			}
		}
		return(binning.sd.output)
	}

	gen.binning.delta.TE <- function(model.binning.coef, model.binning.vcov, char=NULL, D.ref=NULL){
		if(is.null(char)==TRUE){
			treat.type='continuous'
		}
		if(is.null(D.ref)==TRUE){
			treat.type='discrete'
		}

		if(use_fe==TRUE){
			return(gen.binning.delta.TE.fe(model.binning.coef=model.binning.coef, model.binning.vcov=model.binning.vcov, char=char, D.ref=D.ref))
		}

		model.binning.coef[which(is.na(model.binning.coef))] <- 0
		binning.sd.output <- c()
		
		if(treat.type=='discrete'){
			for(i in 1:nbins){
				link.bin.1 <- model.binning.coef[paste0("G.",i)] + model.binning.coef[paste0('D.',char,'.','G.',i)]
				link.bin.0 <- model.binning.coef[paste0("G.",i)]
				
				if(is.null(Z)==FALSE){
					vec.1 <- c(1,1,Z.ref)
					vec.0 <- c(1,0,Z.ref)
					if(full.moderate==FALSE){
						target.slice <- c(paste0("G.",i),paste0('D.',char,'.','G.',i),Z)
					}

					if(full.moderate==TRUE){
						target.slice <- c(paste0("G.",i),paste0('D.',char,'.','G.',i))
					}

					for(a in Z){
						target.Z <- Z.ref[a]
						if(full.moderate==FALSE){
							link.bin.1 <- link.bin.1 + target.Z*model.binning.coef[a]
							link.bin.0 <- link.bin.0 + target.Z*model.binning.coef[a]
						}

						if(full.moderate==TRUE){
							link.bin.1 <- link.bin.1 + target.Z*model.binning.coef[paste0(a,".G.",i)]
							link.bin.0 <- link.bin.0 + target.Z*model.binning.coef[paste0(a,".G.",i)]
							target.slice <- c(target.slice, paste0(a,".G.",i))
						}
					}
				}
				else{
					vec.1 <- c(1,1)
					vec.0 <- c(1,0)
					target.slice <- c(paste0("G.",i),paste0('D.',char,'.','G.',i))
				}

				vcov.binning.temp <- model.binning.vcov[target.slice,target.slice]
			
				if(method=='logit'){
					vec <- vec.1*exp(link.bin.1)/(1+exp(link.bin.1))^2 - vec.0*exp(link.bin.0)/(1+exp(link.bin.0))^2
				}
				if(method=='probit'){
					vec <- vec.1*dnorm(link.bin.1) - vec.0*dnorm(link.bin.0)
				}
				if(method=='poisson' | method=='nbinom'){
					vec <- vec.1*exp(link.bin.1) - vec.0*exp(link.bin.0)
				}
				if(method=='linear'){
					vec <- vec.1-vec.0
				}
			
				delta.sd.bin <- sqrt((t(vec)%*%vcov.binning.temp%*%vec)[1,1])
				delta.sd.bin <- c(delta.sd.bin)
				names(delta.sd.bin) <- paste0("sd.G.",i)
				binning.sd.output <- c(binning.sd.output,delta.sd.bin)
			}
		}
		
		if(treat.type=='continuous'){
			for(i in 1:nbins){
				link.bin <- model.binning.coef[paste0("G.",i)] + model.binning.coef[paste0('D.G.',i)]*D.ref
				if(is.null(Z)==FALSE){
					vec1 <- c(1,D.ref,Z.ref)
					vec0 <- c(0,1,rep(0,length(Z)))

					if(full.moderate==FALSE){
						target.slice <- c(paste0("G.",i),paste0('D.G.',i),Z)
					}

					if(full.moderate==TRUE){
						target.slice <- c(paste0("G.",i),paste0('D.G.',i))
					}

					for(a in Z){
						target.Z <- Z.ref[a]
						if(full.moderate==FALSE){
							link.bin <- link.bin + target.Z*model.binning.coef[a]
						}

						if(full.moderate==TRUE){
							link.bin <- link.bin + target.Z*model.binning.coef[paste0(a,".G.",i)]
							target.slice <- c(target.slice, paste0(a,".G.",i))
						}
					}
				}
				else{
					vec1 <- c(1,D.ref)
					vec0 <- c(0,1)
					target.slice <- c(paste0("G.",i),paste0('D.G.',i))
				}
				vcov.binning.temp <- model.binning.vcov[target.slice,target.slice]
			
				if(method=='logit'){
					vec <- -model.binning.coef[paste0('D.G.',i)]*(exp(link.bin)-exp(-link.bin))/(2+exp(link.bin)+exp(-link.bin))^2*vec1 + exp(link.bin)/(1+exp(link.bin))^2*vec0
				}
				if(method=='probit'){
					vec <- dnorm(link.bin)*vec0-model.binning.coef[paste0('D.G.',i)]*link.bin*dnorm(link.bin)*vec1
				}
				if(method=='poisson'|method=='nbinom'){
					vec <- model.binning.coef[paste0('D.G.',i)]*exp(link.bin)*vec1+exp(link.bin)*vec0
				}
				if(method=='linear'){
					vec <- vec0
				}
			
				delta.sd.bin <- sqrt((t(vec)%*%vcov.binning.temp%*%vec)[1,1])
				delta.sd.bin <- c(delta.sd.bin)
				names(delta.sd.bin) <- paste0("sd.G.",i)
				binning.sd.output <- c(binning.sd.output,delta.sd.bin)
			}
		}
		return(binning.sd.output)
	}
	
	## Function B.3a (link predict)
	gen.link.binning <- function(model.binning.coef, X.eval, cuts.X, x0, char=NULL, D.ref=NULL){
		if(is.null(char)==TRUE){
			treat.type='continuous'
		}
		if(is.null(D.ref)==TRUE){
			treat.type=='discrete'
		}

		if(use_fe==TRUE){
			if(treat.type=='discrete'){
				link.1 <- link.0 <- rep(0,length(X.eval))
				return(list(link.1=link.1,link.0=link.0))
			}
			if(treat.type=='continuous'){
				link <- rep(0,length(X.eval))
				return(list(link=link))
			}
		}

		model.binning.coef[which(is.na(model.binning.coef))] <- 0
		group.X.eval <- cut(X.eval,breaks=cuts.X, labels = FALSE)
		group.X.eval[which(X.eval==min(X.eval))] <- 1

		if(treat.type=='discrete'){
			link.bin.1 <- 0
			link.bin.0 <- 0

			for(i in 1:nbins){
				link.bin.1 <- link.bin.1 +
							  model.binning.coef[paste0("G.",i)]*as.numeric(group.X.eval==i) + 
							  model.binning.coef[paste0("GX.",i)]*as.numeric(group.X.eval==i)*(X.eval-x0[i]) +
							  model.binning.coef[paste0('D.',char,'.','G.',i)]*as.numeric(group.X.eval==i) +
							  model.binning.coef[paste0('D.',char,'.','GX.',i)]*as.numeric(group.X.eval==i)*(X.eval-x0[i])
					
				link.bin.0 <- link.bin.0 +
							  model.binning.coef[paste0("G.",i)]*as.numeric(group.X.eval==i) + 
							  model.binning.coef[paste0("GX.",i)]*as.numeric(group.X.eval==i)*(X.eval-x0[i])
			
			}
			
			if(is.null(Z)==FALSE){
				for(a in Z){
					target.Z <- Z.ref[a]
					if(full.moderate==FALSE){
						link.bin.1 <- link.bin.1 + target.Z*model.binning.coef[a]
						link.bin.0 <- link.bin.0 + target.Z*model.binning.coef[a]
					}
					if(full.moderate==TRUE){
						for(i in 1:nbins){
							link.bin.1 <- link.bin.1 + target.Z*model.binning.coef[paste0(a,".G.",i)]*as.numeric(group.X.eval==i) + 
										  target.Z*model.binning.coef[paste0(a,".GX.",i)]*as.numeric(group.X.eval==i)*(X.eval-x0[i])

							link.bin.0 <- link.bin.0 + target.Z*model.binning.coef[paste0(a,".G.",i)]*as.numeric(group.X.eval==i) + 
										  target.Z*model.binning.coef[paste0(a,".GX.",i)]*as.numeric(group.X.eval==i)*(X.eval-x0[i])
						}
					}
				}
			}
			return(list(link.1=link.bin.1,link.0=link.bin.0))
		}
		
		if(treat.type=='continuous'){
			link.bin <- 0
			for(i in 1:nbins){
				link.bin <- link.bin + 
							model.binning.coef[paste0("G.",i)]*as.numeric(group.X.eval==i) + 
							model.binning.coef[paste0("GX.",i)]*as.numeric(group.X.eval==i)*(X.eval-x0[i]) +
							model.binning.coef[paste0('D.G.',i)]*as.numeric(group.X.eval==i)*D.ref + 
							model.binning.coef[paste0('D.GX.',i)]*as.numeric(group.X.eval==i)*D.ref*(X.eval-x0[i])
			}
			if(is.null(Z)==FALSE){
				for(a in Z){
					target.Z <- Z.ref[a]
					if(full.moderate==FALSE){
						link.bin <- link.bin + target.Z*model.binning.coef[a]
					}
					if(full.moderate==TRUE){
						for(i in 1:nbins){
							link.bin <- link.bin + target.Z*model.binning.coef[paste0(a,".G.",i)]*as.numeric(group.X.eval==i) + 
										target.Z*model.binning.coef[paste0(a,".GX.",i)]*as.numeric(group.X.eval==i)*(X.eval-x0[i])
						}
					}
				}
			}
			return(list(link=link.bin))
		}
	}


									 	
	## Function B.3 (binning predict)
	#1. estimate predicted value of binning estimator
	#2. input: model.binning.coef, X.eval, cuts.X, x0
	#3. output: predicted value using binning estimator
	gen.pred.binning <- function(model.binning.coef, X.eval, cuts.X, x0, char=NULL, D.ref=NULL){
		if(is.null(char)==TRUE){
			treat.type='continuous'
		}
		if(is.null(D.ref)==TRUE){
			treat.type=='discrete'
		}

		if(use_fe==TRUE){
			if(treat.type=='discrete'){
				link.1 <- link.0 <- E.base <- E.pred <- rep(0,length(X.eval))
				return(list(E.pred=E.pred,E.base=E.base,link.1=link.1,link.0=link.0))
			}
		
			if(treat.type=='continuous'){
				link <- E.pred <- rep(0,length(X.eval))
				return(list(E.pred=E.pred,link=link))
			}
		}

		model.binning.coef[which(is.na(model.binning.coef))] <- 0
		group.X.eval <- cut(X.eval,breaks=cuts.X, labels = FALSE)
		group.X.eval[which(X.eval==min(X.eval))] <- 1

		if(treat.type=='discrete'){
			link.bin.1 <- 0
			link.bin.0 <- 0

			for(i in 1:nbins){
				link.bin.1 <- link.bin.1 +
							  model.binning.coef[paste0("G.",i)]*as.numeric(group.X.eval==i) + 
							  model.binning.coef[paste0("GX.",i)]*as.numeric(group.X.eval==i)*(X.eval-x0[i]) +
							  model.binning.coef[paste0('D.',char,'.','G.',i)]*as.numeric(group.X.eval==i) +
							  model.binning.coef[paste0('D.',char,'.','GX.',i)]*as.numeric(group.X.eval==i)*(X.eval-x0[i])
					
				link.bin.0 <- link.bin.0 +
							  model.binning.coef[paste0("G.",i)]*as.numeric(group.X.eval==i) + 
							  model.binning.coef[paste0("GX.",i)]*as.numeric(group.X.eval==i)*(X.eval-x0[i])
			
			}
			
			if(is.null(Z)==FALSE){
				for(a in Z){
					target.Z <- Z.ref[a]
					if(full.moderate==FALSE){
						link.bin.1 <- link.bin.1 + target.Z*model.binning.coef[a]
						link.bin.0 <- link.bin.0 + target.Z*model.binning.coef[a]
					}
					if(full.moderate==TRUE){
						for(i in 1:nbins){
							link.bin.1 <- link.bin.1 + target.Z*model.binning.coef[paste0(a,".G.",i)]*as.numeric(group.X.eval==i) + 
										  target.Z*model.binning.coef[paste0(a,".GX.",i)]*as.numeric(group.X.eval==i)*(X.eval-x0[i])

							link.bin.0 <- link.bin.0 + target.Z*model.binning.coef[paste0(a,".G.",i)]*as.numeric(group.X.eval==i) + 
										  target.Z*model.binning.coef[paste0(a,".GX.",i)]*as.numeric(group.X.eval==i)*(X.eval-x0[i])
						}
					}
				}
			}
			
			if(method=='linear'){
				E.pred <- link.bin.1
				E.base <- link.bin.0
			}
		
			if(method=='logit'){
				E.pred <- E.prob.1 <- exp(link.bin.1)/(1+exp(link.bin.1))
				E.base <- E.prob.0 <- exp(link.bin.0)/(1+exp(link.bin.0))
			}
	
			if(method=='probit'){
				E.pred <- E.prob.1 <- pnorm(link.bin.1,0,1)
				E.base <- E.prob.0 <- pnorm(link.bin.0,0,1)
			}
		
			if(method=='poisson' | method=="nbinom"){
				E.pred <- E.y.1 <- exp(link.bin.1)
				E.base <- E.y.0 <- exp(link.bin.0)
			}
			return(list(E.pred=E.pred,E.base=E.base,link.1=link.bin.1,link.0=link.bin.0))
		}
		
		if(treat.type=='continuous'){
			link.bin <- 0
			for(i in 1:nbins){
				link.bin <- link.bin + 
							model.binning.coef[paste0("G.",i)]*as.numeric(group.X.eval==i) + 
							model.binning.coef[paste0("GX.",i)]*as.numeric(group.X.eval==i)*(X.eval-x0[i]) +
							model.binning.coef[paste0('D.G.',i)]*as.numeric(group.X.eval==i)*D.ref + 
							model.binning.coef[paste0('D.GX.',i)]*as.numeric(group.X.eval==i)*D.ref*(X.eval-x0[i])
			}
			if(is.null(Z)==FALSE){
				for(a in Z){
					target.Z <- Z.ref[a]
					if(full.moderate==FALSE){
						link.bin <- link.bin + target.Z*model.binning.coef[a]
					}
					if(full.moderate==TRUE){
						for(i in 1:nbins){
							link.bin <- link.bin + target.Z*model.binning.coef[paste0(a,".G.",i)]*as.numeric(group.X.eval==i) + 
										target.Z*model.binning.coef[paste0(a,".GX.",i)]*as.numeric(group.X.eval==i)*(X.eval-x0[i])
						}
					}
				}
			}
			if(method=='logit'){
				E.pred <- exp(link.bin)/(1+exp(link.bin))
			}
			if(method=='probit'){
				E.pred <- pnorm(link.bin,0,1)
			}
			if(method=='linear'){
				E.pred <- link.bin
			}
			if(method=='poisson'|method=='nbinom'){
				E.pred <- exp(link.bin)
			}
			return(list(E.pred=E.pred,link=link.bin))
		}
	}

	## Function B.4a (binning link delta)
	gen.link.binning.delta <- function(model.binning.coef, model.binning.vcov, X.eval, cuts.X, x0, char=NULL, D.ref=NULL){
		if(is.null(char)==TRUE){
			treat.type='continuous'
		}
		if(is.null(D.ref)==TRUE){
			treat.type=='discrete'
		}
		model.binning.coef[which(is.na(model.binning.coef))] <- 0

		gen.link.binning.delta.sd <- function(x){
			if(use_fe==TRUE){
				return(0)
			}
			group.xx <- cut(x,breaks=cuts.X, labels = FALSE,include.lowest = TRUE)
			if(is.na(group.xx)==TRUE){
				return(NA)
			}
			if(treat.type=='discrete'){			
				if(char!=base){
					link.bin.1 <- model.binning.coef[paste0("G.",group.xx)]+ 
								  model.binning.coef[paste0("GX.",group.xx)]*(x-x0[group.xx]) +
								  model.binning.coef[paste0('D.',char,'.','G.',group.xx)] +
								  model.binning.coef[paste0('D.',char,'.','GX.',group.xx)]*(x-x0[group.xx])
					target.slice <- c(paste0("G.",group.xx),
									  paste0("GX.",group.xx),
									  paste0('D.',char,'.','G.',group.xx),
								      paste0('D.',char,'.','GX.',group.xx))
				}
				if(char==base){
					link.bin.1 <- model.binning.coef[paste0("G.",group.xx)]+ 
								  model.binning.coef[paste0("GX.",group.xx)]*(x-x0[group.xx]) 			  
					target.slice <- c(paste0("G.",group.xx),
									  paste0("GX.",group.xx))
				}
				
				if(is.null(Z)==FALSE){
					if(full.moderate==FALSE){
						target.slice <- c(target.slice,Z)
						if(char!=base){
							vec <- c(1,x-x0[group.xx],1,x-x0[group.xx],Z.ref)
						}
						if(char==base){
							vec <- c(1,x-x0[group.xx],Z.ref)
						}
					}
					if(full.moderate==TRUE){
						if(char!=base){
							vec <- c(1,x-x0[group.xx],1,x-x0[group.xx])
						}
						if(char==base){
							vec <- c(1,x-x0[group.xx])
						}
					}
					for(a in Z){
						target.Z <- Z.ref[a]
						if(full.moderate==FALSE){
							link.bin.1 <- link.bin.1 + target.Z*model.binning.coef[a]
						}	
						if(full.moderate==TRUE){
							link.bin.1 <- link.bin.1 + target.Z*model.binning.coef[paste0(a,".G.",group.xx)] + target.Z*model.binning.coef[paste0(a,".GX.",group.xx)]*(x-x0[group.xx])
							target.slice <- c(target.slice, paste0(a,".G.",group.xx), paste0(a,".GX.",group.xx))
							vec <- c(vec, target.Z, target.Z*(x-x0[group.xx]))
						}
					}
				}
				else{
					if(char!=base){
						vec <- c(1,x-x0[group.xx],1,x-x0[group.xx])
					}
					if(char==base){
						vec <- c(1,x-x0[group.xx])
					}
				}
								
				temp.vcov.matrix <- model.binning.vcov[target.slice,target.slice]
				predict.sd <- sqrt((t(vec)%*%temp.vcov.matrix%*%vec)[1,1])
				return(predict.sd)
			}
			
			if(treat.type=='continuous'){
				target.slice <- c(paste0("G.",group.xx),
								  paste0("GX.",group.xx),
								  paste0('D.G.',group.xx),
								  paste0('D.GX.',group.xx))
				
				link.bin <- model.binning.coef[paste0("G.",group.xx)]+ 
							  model.binning.coef[paste0("GX.",group.xx)]*(x-x0[group.xx]) +
							  model.binning.coef[paste0('D.G.',group.xx)]*D.ref +
							  model.binning.coef[paste0('D.GX.',group.xx)]*(x-x0[group.xx])*D.ref

				if(is.null(Z)==FALSE){
					if(full.moderate==FALSE){
						vec <- c(1,x-x0[group.xx],D.ref,(x-x0[group.xx])*D.ref,Z.ref)
						target.slice <- c(target.slice,Z)
					}

					if(full.moderate==TRUE){
						vec <- c(1,x-x0[group.xx],D.ref,(x-x0[group.xx])*D.ref)
					}

					for(a in Z){
						target.Z <- Z.ref[a]
						if(full.moderate==FALSE){
							link.bin <- link.bin + target.Z*model.binning.coef[a]
						}	
						if(full.moderate==TRUE){
							link.bin <- link.bin + target.Z*model.binning.coef[paste0(a,".G.",group.xx)] + target.Z*model.binning.coef[paste0(a,".GX.",group.xx)]*(x-x0[group.xx])
							target.slice <- c(target.slice,paste0(a,".G.",group.xx),paste0(a,".GX.",group.xx))
							vec <- c(vec, target.Z, target.Z*(x-x0[group.xx]))
						}
					}
				}
				else{
					vec <- c(1,x-x0[group.xx],D.ref,(x-x0[group.xx])*D.ref)
				}
				temp.vcov.matrix <- model.binning.vcov[target.slice,target.slice]
				predict.sd <- sqrt((t(vec)%*%temp.vcov.matrix%*%vec)[1,1])
				return(predict.sd)
			}
		}

		binning.link.delta.sd <- c(sapply(X.eval,function(x) gen.link.binning.delta.sd(x)))
		names(binning.link.delta.sd) <- NULL
		return(binning.link.delta.sd)
	}

	## Function B.4 (binning predict delta)
	#1. estimate sd of predicted value(binning) using delta method 
	#2. input: model.binning.coef, model.binning.vcov, X.eval, cuts.X, x0, char(discrete), D.ref(continuous)
	#3. output: sd of predicted value(binning)
	gen.pred.binning.delta <- function(model.binning.coef, model.binning.vcov, X.eval, cuts.X, x0, char=NULL, D.ref=NULL){
		if(is.null(char)==TRUE){
			treat.type='continuous'
		}
		if(is.null(D.ref)==TRUE){
			treat.type=='discrete'
		}
		model.binning.coef[which(is.na(model.binning.coef))] <- 0
		
		gen.pred.binning.delta.sd <- function(x){
			if(use_fe==TRUE){
				return(0)
			}

			group.xx <- cut(x,breaks=cuts.X, labels = FALSE,include.lowest = TRUE)
			if(is.na(group.xx)==TRUE){
				return(NA)
			}
			if(treat.type=='discrete'){
						
				if(char!=base){
					link.bin.1 <- model.binning.coef[paste0("G.",group.xx)]+ 
								  model.binning.coef[paste0("GX.",group.xx)]*(x-x0[group.xx]) +
								  model.binning.coef[paste0('D.',char,'.','G.',group.xx)] +
								  model.binning.coef[paste0('D.',char,'.','GX.',group.xx)]*(x-x0[group.xx])
					
					target.slice <- c(paste0("G.",group.xx),
									  paste0("GX.",group.xx),
									  paste0('D.',char,'.','G.',group.xx),
								      paste0('D.',char,'.','GX.',group.xx))
				}
				
				if(char==base){
					link.bin.1 <- model.binning.coef[paste0("G.",group.xx)]+ 
								  model.binning.coef[paste0("GX.",group.xx)]*(x-x0[group.xx]) 
								  
					target.slice <- c(paste0("G.",group.xx),
									  paste0("GX.",group.xx))
				
				}
				
				if(is.null(Z)==FALSE){
					if(full.moderate==FALSE){
						target.slice <- c(target.slice,Z)
						if(char!=base){
							vec <- c(1,x-x0[group.xx],1,x-x0[group.xx],Z.ref)
						}
						if(char==base){
							vec <- c(1,x-x0[group.xx],Z.ref)
						}
					}

					if(full.moderate==TRUE){
						if(char!=base){
							vec <- c(1,x-x0[group.xx],1,x-x0[group.xx])
						}
						if(char==base){
							vec <- c(1,x-x0[group.xx])
						}
					}

					for(a in Z){
						target.Z <- Z.ref[a]
						if(full.moderate==FALSE){
							link.bin.1 <- link.bin.1 + target.Z*model.binning.coef[a]
						}	
						if(full.moderate==TRUE){
							link.bin.1 <- link.bin.1 + target.Z*model.binning.coef[paste0(a,".G.",group.xx)] + target.Z*model.binning.coef[paste0(a,".GX.",group.xx)]*(x-x0[group.xx])
							target.slice <- c(target.slice, paste0(a,".G.",group.xx), paste0(a,".GX.",group.xx))
							vec <- c(vec, target.Z, target.Z*(x-x0[group.xx]))
						}
					}
				}
				else{
					if(char!=base){
						vec <- c(1,x-x0[group.xx],1,x-x0[group.xx])
					}
					if(char==base){
						vec <- c(1,x-x0[group.xx])
					}
				}
				
				if(method=='logit'){
					vec <- vec*exp(link.bin.1)/(1+exp(link.bin.1))^2 
				}
				if(method=='probit'){
					vec <- vec*dnorm(link.bin.1)
				}
				if(method=='poisson' | method=='nbinom'){
					vec <- vec*exp(link.bin.1)
				}
				if(method=='linear'){
					vec <- vec
				}
				
				temp.vcov.matrix <- model.binning.vcov[target.slice,target.slice]
				predict.sd <- sqrt((t(vec)%*%temp.vcov.matrix%*%vec)[1,1])
				return(predict.sd)
				
			}
			
			if(treat.type=='continuous'){
				target.slice <- c(paste0("G.",group.xx),
								  paste0("GX.",group.xx),
								  paste0('D.G.',group.xx),
								  paste0('D.GX.',group.xx))
				
				link.bin <- model.binning.coef[paste0("G.",group.xx)]+ 
							  model.binning.coef[paste0("GX.",group.xx)]*(x-x0[group.xx]) +
							  model.binning.coef[paste0('D.G.',group.xx)]*D.ref +
							  model.binning.coef[paste0('D.GX.',group.xx)]*(x-x0[group.xx])*D.ref

				if(is.null(Z)==FALSE){
					if(full.moderate==FALSE){
						vec <- c(1,x-x0[group.xx],D.ref,(x-x0[group.xx])*D.ref,Z.ref)
						target.slice <- c(target.slice,Z)
					}

					if(full.moderate==TRUE){
						vec <- c(1,x-x0[group.xx],D.ref,(x-x0[group.xx])*D.ref)
					}

					for(a in Z){
						target.Z <- Z.ref[a]
						if(full.moderate==FALSE){
							link.bin <- link.bin + target.Z*model.binning.coef[a]
						}	
						if(full.moderate==TRUE){
							link.bin <- link.bin + target.Z*model.binning.coef[paste0(a,".G.",group.xx)] + target.Z*model.binning.coef[paste0(a,".GX.",group.xx)]*(x-x0[group.xx])
							target.slice <- c(target.slice,paste0(a,".G.",group.xx),paste0(a,".GX.",group.xx))
							vec <- c(vec, target.Z, target.Z*(x-x0[group.xx]))
						}
					}
				}
				else{
					vec <- c(1,x-x0[group.xx],D.ref,(x-x0[group.xx])*D.ref)
				}
				
				temp.vcov.matrix <- model.binning.vcov[target.slice,target.slice]
				if(method=='logit'){
					vec <- vec*exp(link.bin)/(1+exp(link.bin))^2 
				}
				if(method=='probit'){
					vec <- vec*dnorm(link.bin)
				}
				if(method=='poisson' | method=='nbinom'){
					vec <- vec*exp(link.bin)
				}
				if(method=='linear'){
					vec <- vec
				}
				predict.sd <- sqrt((t(vec)%*%temp.vcov.matrix%*%vec)[1,1])
				return(predict.sd)
			}
		}
		
		binning.pred.delta.sd <- c(sapply(X.eval,function(x) gen.pred.binning.delta.sd(x)))
		names(binning.pred.delta.sd) <- NULL
		return(binning.pred.delta.sd)
	
	}
	
	
	#vartype: simulate
	if(vartype=='simu'){
		#simu
		M <- nsimu
		simu.coef <- rmvnorm(M, model.coef, model.vcov)
		simu.binning.coef <- rmvnorm(M, model.binning.coef, model.binning.vcov)
		
		if(treat.type=='discrete'){
			TE.output.all.list <- list()
			TE.binning.output.all.list <- list()
			pred.output.all.list <- list()
			link.output.all.list <- list()
			for(char in other.treat){
				gen.general.TE.output <- gen.general.TE(model.coef=model.coef,char=char)
				TE.output <- gen.general.TE.output$TE
				TE.binning.output <- gen.binning.TE(model.binning.coef=model.binning.coef,
													char=char)
				pred.binning.output <- gen.pred.binning(model.binning.coef=model.binning.coef, 
														X.eval=X.eval, 
														cuts.X=cuts.X, 
														x0=x0, 
														char=char)
				
				#simu
				one.simu1.output <- function(model.coef){
					one.simu.TE.output <- gen.general.TE(model.coef=model.coef,char=char)
					output <- one.simu.TE.output$TE
					names(output) <- rep("TE",length(output))
					return(output)
				}
				
				one.simu2.output <- function(model.binning.coef){
					simu.TE.binning.output <- gen.binning.TE(model.binning.coef=model.binning.coef,
															 char=char)
					simu.pred.binning.output <- gen.pred.binning(model.binning.coef=model.binning.coef, 
																 X.eval=X.eval, 
																 cuts.X=cuts.X, 
																 x0=x0, 
																 char=char)
					
					output <- c(simu.TE.binning.output,simu.pred.binning.output$E.pred,simu.pred.binning.output$E.base,
								simu.pred.binning.output$link.1,simu.pred.binning.output$link.0)

					names(output) <- c(rep("bins.est",nbins),rep("pred.bin",length(X.eval)),rep("base.bin",length(X.eval)),
									   rep("link.1",length(X.eval)),rep("link.0",length(X.eval)))

					return(output)
				}
				
				TE.simu.matrix <- apply(simu.coef, 1, function(x) one.simu1.output(model.coef=x))
				binning.simu.matrix <- apply(simu.binning.coef,1,function(x) one.simu2.output(model.binning.coef=x))
				pred.simu.matrix <- binning.simu.matrix[rownames(binning.simu.matrix)=="pred.bin",]
				base.simu.matrix <- binning.simu.matrix[rownames(binning.simu.matrix)=="base.bin",]
				link.1.simu.matrix <- binning.simu.matrix[rownames(binning.simu.matrix)=="link.1",]
				link.0.simu.matrix <- binning.simu.matrix[rownames(binning.simu.matrix)=="link.0",]

				est.simu.matrix <- binning.simu.matrix[rownames(binning.simu.matrix)=="bins.est",]
				
				TE.simu.sd <- apply(TE.simu.matrix, 1, sd, na.rm=TRUE)
				pred.simu.sd <- apply(pred.simu.matrix, 1, sd, na.rm=TRUE)
				base.simu.sd <- apply(base.simu.matrix, 1, sd, na.rm=TRUE)
				link.1.simu.sd <- apply(link.1.simu.matrix, 1, sd, na.rm=TRUE)
				link.0.simu.sd <- apply(link.0.simu.matrix, 1, sd, na.rm=TRUE)
				est.simu.sd <- apply(est.simu.matrix, 1, sd, na.rm=TRUE)
				
				TE.simu.CI <- t(apply(TE.simu.matrix, 1, quantile, c(0.025,0.975),na.rm=TRUE))
				pred.simu.CI <- t(apply(pred.simu.matrix, 1, quantile, c(0.025,0.975),na.rm=TRUE))
				base.simu.CI <- t(apply(base.simu.matrix, 1, quantile, c(0.025,0.975),na.rm=TRUE))
				link.1.simu.CI <- t(apply(link.1.simu.matrix, 1, quantile, c(0.025,0.975),na.rm=TRUE))
				link.0.simu.CI <- t(apply(link.0.simu.matrix, 1, quantile, c(0.025,0.975),na.rm=TRUE))
				est.simu.CI <- t(apply(est.simu.matrix, 1, quantile, c(0.025,0.975),na.rm=TRUE))

				TE.output.all <- cbind(X.eval,TE.output,TE.simu.sd,TE.simu.CI[,1],TE.simu.CI[,2])
				colnames(TE.output.all) <- c("X","TE","sd","lower CI(95%)","upper CI(95%)")
				rownames(TE.output.all) <- NULL
				TE.output.all.list[[other.treat.origin[char]]] <- TE.output.all
				
				TE.binning.all <- cbind(x0,TE.binning.output,est.simu.sd,est.simu.CI[,1],est.simu.CI[,2])
				colnames(TE.binning.all) <- c('x0',"coef","sd","CI.lower","CI.upper")
				rownames(TE.binning.all) <- gp.lab
				TE.binning.output.all.list[[other.treat.origin[char]]] <- TE.binning.all
				
				pred.output.all <- cbind(X.eval,pred.binning.output$E.pred,pred.simu.sd,pred.simu.CI[,1],pred.simu.CI[,2])
				colnames(pred.output.all) <- c("X","E(Y)","sd","lower CI(95%)","upper CI(95%)")
				rownames(pred.output.all) <- NULL
				pred.output.all.list[[other.treat.origin[char]]] <- pred.output.all

				link.output.all <- cbind(X.eval,pred.binning.output$link.1,link.1.simu.sd,link.1.simu.CI[,1],link.1.simu.CI[,2])
				colnames(link.output.all) <- c("X","E(Y)","sd","lower CI(95%)","upper CI(95%)")
				rownames(link.output.all) <- NULL
				link.output.all.list[[other.treat.origin[char]]] <- link.output.all
				
			}
			
			#base
			base.output.all <- cbind(X.eval,pred.binning.output$E.base,base.simu.sd,base.simu.CI[,1],base.simu.CI[,2])
			colnames(base.output.all) <- c("X","E(Y)","sd","lower CI(95%)","upper CI(95%)")
			rownames(base.output.all) <- NULL
			pred.output.all.list[[all.treat.origin[base]]] <- base.output.all	

			link.0.output.all <- cbind(X.eval,pred.binning.output$link.0,link.0.simu.sd,link.0.simu.CI[,1],link.0.simu.CI[,2])
			colnames(link.0.output.all) <- c("X","E(Y)","sd","lower CI(95%)","upper CI(95%)")
			rownames(link.0.output.all) <- NULL
			link.output.all.list[[all.treat.origin[base]]] <- link.0.output.all	
		}
		
		if(treat.type=='continuous'){
			ME.output.all.list <- list()
			ME.binning.output.all.list <- list()
			pred.output.all.list <- list()
			link.output.all.list <- list()
			k <- 1
			for(D.ref in D.sample){
				gen.general.ME.output <- gen.general.TE(model.coef=model.coef,D.ref=D.ref)
				ME.output <- gen.general.ME.output$ME
				ME.binning.output <- gen.binning.TE(model.binning.coef=model.binning.coef,
													D.ref=D.ref)
				pred.binning.output <- gen.pred.binning(model.binning.coef=model.binning.coef, 
														X.eval=X.eval, 
														cuts.X=cuts.X, 
														x0=x0, 
														D.ref=D.ref)
				#simu
				one.simu1.output <- function(model.coef){
					one.simu.ME.output <- gen.general.TE(model.coef=model.coef,D.ref=D.ref)
					output <- one.simu.ME.output$ME
					names(output) <- rep("ME",length(output))
					return(output)
				}
				
				one.simu2.output <- function(model.binning.coef){
					simu.ME.binning.output <- gen.binning.TE(model.binning.coef=model.binning.coef,
															 D.ref=D.ref)
					simu.pred.binning.output <- gen.pred.binning(model.binning.coef=model.binning.coef, 
																 X.eval=X.eval, 
																 cuts.X=cuts.X, 
																 x0=x0, 
																 D.ref=D.ref)
					output <- c(simu.ME.binning.output,simu.pred.binning.output$E.pred,simu.pred.binning.output$link)
					names(output) <- c(rep("bins.est",nbins),rep("pred.bin",length(X.eval)),rep("link",length(X.eval)))
					return(output)
				}
				
				ME.simu.matrix <- apply(simu.coef, 1, function(x) one.simu1.output(model.coef=x))
				binning.simu.matrix <- apply(simu.binning.coef,1,function(x) one.simu2.output(model.binning.coef=x))
				
				pred.simu.matrix <- binning.simu.matrix[rownames(binning.simu.matrix)=="pred.bin",]
				est.simu.matrix <- binning.simu.matrix[rownames(binning.simu.matrix)=="bins.est",]
				link.simu.matrix <- binning.simu.matrix[rownames(binning.simu.matrix)=="link",]
				
				ME.simu.sd <- apply(ME.simu.matrix, 1, sd)
				pred.simu.sd <- apply(pred.simu.matrix, 1, sd)
				link.simu.sd <- apply(link.simu.matrix, 1, sd)
				est.simu.sd <- apply(est.simu.matrix, 1, sd)
				
				ME.simu.CI <- t(apply(ME.simu.matrix, 1, quantile, c(0.025,0.975)))
				pred.simu.CI <- t(apply(pred.simu.matrix, 1, quantile, c(0.025,0.975)))
				link.simu.CI <- t(apply(link.simu.matrix, 1, quantile, c(0.025,0.975)))
				est.simu.CI <- t(apply(est.simu.matrix, 1, quantile, c(0.025,0.975)))

				ME.output.all <- cbind(X.eval,ME.output,ME.simu.sd,ME.simu.CI[,1],ME.simu.CI[,2])
				colnames(ME.output.all) <- c("X","ME","sd","lower CI(95%)","upper CI(95%)")
				rownames(ME.output.all) <- NULL
				ME.output.all.list[[label.name[k]]] <- ME.output.all
				
				ME.binning.all <- cbind(x0,ME.binning.output,est.simu.sd,est.simu.CI[,1],est.simu.CI[,2])
				colnames(ME.binning.all) <- c('x0',"coef","sd","CI.lower","CI.upper")
				rownames(ME.binning.all) <- gp.lab
				ME.binning.output.all.list[[label.name[k]]] <- ME.binning.all
				
				pred.output.all <- cbind(X.eval,pred.binning.output$E.pred,pred.simu.sd,pred.simu.CI[,1],pred.simu.CI[,2])
				colnames(pred.output.all) <- c("X","E(Y)","sd","lower CI(95%)","upper CI(95%)")
				rownames(pred.output.all) <- NULL
				pred.output.all.list[[label.name[k]]] <- pred.output.all

				link.output.all <- cbind(X.eval,pred.binning.output$link,link.simu.sd,link.simu.CI[,1],link.simu.CI[,2])
				colnames(link.output.all) <- c("X","E(Y)","sd","lower CI(95%)","upper CI(95%)")
				rownames(link.output.all) <- NULL
				link.output.all.list[[label.name[k]]] <- link.output.all

				k <- k + 1
			}			
		}	
	}
	
	#vartype: delta
	if(vartype=='delta'){
		crit.lin <- abs(qt(.025, df=model.df))
		crit.bin <- abs(qt(.025, df=model.binning.df))

		if(treat.type=='discrete'){
			TE.output.all.list <- list()
			TE.binning.output.all.list <- list()
			pred.output.all.list <- list()
			link.output.all.list <- list()
			for(char in other.treat){
				gen.general.TE.output <- gen.general.TE(model.coef=model.coef,char=char)
				TE.output <- gen.general.TE.output$TE
				TE.binning.output <- gen.binning.TE(model.binning.coef=model.binning.coef,
													char=char)
				pred.binning.output <- gen.pred.binning(model.binning.coef=model.binning.coef, 
														X.eval=X.eval, 
														cuts.X=cuts.X, 
														x0=x0, 
														char=char)
				
				TE.delta.sd <- gen.delta.TE(model.coef=model.coef,model.vcov=model.vcov,char=char)$TE.sd
				TE.binning.delta.sd <- gen.binning.delta.TE(model.binning.coef=model.binning.coef, 
															model.binning.vcov=model.binning.vcov, 
															char=char)

				pred.binning.delta.sd <- gen.pred.binning.delta(model.binning.coef=model.binning.coef, 
															    model.binning.vcov=model.binning.vcov, 
																X.eval=X.eval, 
																cuts.X=cuts.X, 
																x0=x0, 
																char=char)

				link.binning.delta.sd <- gen.link.binning.delta(model.binning.coef=model.binning.coef, 
															    model.binning.vcov=model.binning.vcov, 
																X.eval=X.eval, 
																cuts.X=cuts.X, 
																x0=x0, 
																char=char)
				
				TE.output.all <- cbind(X.eval,TE.output,TE.delta.sd,TE.output-crit.lin*TE.delta.sd,TE.output+crit.lin*TE.delta.sd)
				colnames(TE.output.all) <- c("X","TE","sd","lower CI(95%)","upper CI(95%)")
				rownames(TE.output.all) <- NULL
				TE.output.all.list[[other.treat.origin[char]]] <- TE.output.all
				
				E.pred.output <- pred.binning.output$E.pred
				pred.output.all <- cbind(X.eval, E.pred.output, pred.binning.delta.sd, E.pred.output-crit.bin*pred.binning.delta.sd, E.pred.output+crit.bin*pred.binning.delta.sd)
				colnames(pred.output.all) <- c("X","E(Y)","sd","lower CI(95%)","upper CI(95%)")
				rownames(pred.output.all) <- NULL
				pred.output.all.list[[other.treat.origin[char]]] <- pred.output.all

				link.output <- pred.binning.output$link.1
				link.output.all <- cbind(X.eval, link.output, link.binning.delta.sd, link.output-crit.bin*link.binning.delta.sd, link.output+crit.bin*link.binning.delta.sd)
				colnames(link.output.all) <- c("X","E(Y)","sd","lower CI(95%)","upper CI(95%)")
				rownames(link.output.all) <- NULL
				link.output.all.list[[other.treat.origin[char]]] <- link.output.all
				
				TE.binning.all <- cbind(x0,TE.binning.output,TE.binning.delta.sd,TE.binning.output-crit.bin*TE.binning.delta.sd,
										TE.binning.output+crit.bin*TE.binning.delta.sd)
				colnames(TE.binning.all) <- c('x0',"coef","sd","CI.lower","CI.upper")
				rownames(TE.binning.all) <- gp.lab
				TE.binning.output.all.list[[other.treat.origin[char]]] <- TE.binning.all 
			}
			
			#base
			E.base.output <- pred.binning.output$E.base
			E.base.delta.sd <- gen.pred.binning.delta(model.binning.coef=model.binning.coef, 
													  model.binning.vcov=model.binning.vcov, 
													  X.eval=X.eval, 
													  cuts.X=cuts.X, 
													  x0=x0, 
													  char=base)
			base.output.all <- cbind(X.eval, E.base.output, E.base.delta.sd, E.base.output-crit.bin*E.base.delta.sd, E.base.output+crit.bin*E.base.delta.sd)
			colnames(base.output.all) <- c("X","E(Y)","sd","lower CI(95%)","upper CI(95%)")
			rownames(base.output.all) <- NULL
			pred.output.all.list[[all.treat.origin[base]]] <- base.output.all

			link.0.output <- pred.binning.output$link.0
			link.0.delta.sd <- gen.link.binning.delta(model.binning.coef=model.binning.coef, 
													  model.binning.vcov=model.binning.vcov, 
													  X.eval=X.eval, 
													  cuts.X=cuts.X, 
													  x0=x0, 
													  char=base)
			link.0.output.all <- cbind(X.eval, link.0.output, link.0.delta.sd, link.0.output-crit.bin*link.0.delta.sd, link.0.output+crit.bin*link.0.delta.sd)
			colnames(link.0.output.all) <- c("X","E(Y)","sd","lower CI(95%)","upper CI(95%)")
			rownames(link.0.output.all) <- NULL
			link.output.all.list[[all.treat.origin[base]]] <- link.0.output.all
		}
	
		if(treat.type=='continuous'){
			ME.output.all.list <- list()
			ME.binning.output.all.list <- list()
			pred.output.all.list <- list()
			link.output.all.list <- list()
			k <- 1
			for(D.ref in D.sample){
				gen.general.ME.output <- gen.general.TE(model.coef=model.coef,D.ref=D.ref)
				ME.output <- gen.general.ME.output$ME
				ME.binning.output <- gen.binning.TE(model.binning.coef=model.binning.coef,
													D.ref=D.ref)
				pred.binning.output <- gen.pred.binning(model.binning.coef=model.binning.coef, 
														X.eval=X.eval, 
														cuts.X=cuts.X, 
														x0=x0, 
														D.ref=D.ref)
				
				ME.delta.sd <- gen.delta.TE(model.coef=model.coef,model.vcov=model.vcov,D.ref=D.ref)$ME.sd
				ME.binning.delta.sd <- gen.binning.delta.TE(model.binning.coef=model.binning.coef, 
															model.binning.vcov=model.binning.vcov, 
															D.ref=D.ref)
				pred.binning.delta.sd <- gen.pred.binning.delta(model.binning.coef=model.binning.coef, 
															    model.binning.vcov=model.binning.vcov, 
																X.eval=X.eval, 
																cuts.X=cuts.X, 
																x0=x0, 
																D.ref=D.ref)

				link.binning.delta.sd <- gen.link.binning.delta(model.binning.coef=model.binning.coef, 
															    model.binning.vcov=model.binning.vcov, 
																X.eval=X.eval, 
																cuts.X=cuts.X, 
																x0=x0, 
																D.ref=D.ref)

				
				ME.output.all <- cbind(X.eval,ME.output,ME.delta.sd, ME.output-crit.lin*ME.delta.sd, ME.output+crit.lin*ME.delta.sd)
				colnames(ME.output.all) <- c("X","ME","sd","lower CI(95%)","upper CI(95%)")
				rownames(ME.output.all) <- NULL
				ME.output.all.list[[label.name[k]]] <- ME.output.all
				
				E.pred.output <- pred.binning.output$E.pred
				pred.output.all <- cbind(X.eval, E.pred.output, pred.binning.delta.sd, E.pred.output-crit.bin*pred.binning.delta.sd, E.pred.output+crit.bin*pred.binning.delta.sd)
				colnames(pred.output.all) <- c("X","E(Y)","sd","lower CI(95%)","upper CI(95%)")
				rownames(pred.output.all) <- NULL
				pred.output.all.list[[label.name[k]]] <- pred.output.all

				link.output <- pred.binning.output$link
				link.output.all <- cbind(X.eval, link.output, link.binning.delta.sd, link.output-crit.bin*link.binning.delta.sd, link.output+crit.bin*link.binning.delta.sd)
				colnames(link.output.all) <- c("X","E(Y)","sd","lower CI(95%)","upper CI(95%)")
				rownames(link.output.all) <- NULL
				link.output.all.list[[label.name[k]]] <- link.output.all
				
				ME.binning.all <- cbind(x0,ME.binning.output,ME.binning.delta.sd,ME.binning.output-crit.bin*ME.binning.delta.sd,
										ME.binning.output+crit.bin*ME.binning.delta.sd)
				colnames(ME.binning.all) <- c('x0',"coef","sd","CI.lower","CI.upper")
				rownames(ME.binning.all) <- gp.lab
				ME.binning.output.all.list[[label.name[k]]] <- ME.binning.all
				k <- k+1
			}
		}	
	}
	
	
	#vartype: bootstrap
	#when conducting binning estimation, use cuts.X and x0 in the original estimation
	if(vartype=='bootstrap'){
		if(treat.type=='discrete'){
			all.length <- neval*length(other.treat) + #TE
						  neval*length(all.treat) + #pred
						  neval*length(all.treat) + #link
						  nbins*length(other.treat)
		}
		
		if(treat.type=='continuous'){
			all.length <- neval*length(label.name) + #ME
						  neval*length(label.name) + #pred
						  neval*length(label.name) + #link
						  nbins*length(label.name) # 
		}
		
		one.boot <- function(){
			if (is.null(cl)==TRUE){
				smp <- sample(1:n,n,replace=TRUE)
			} else{ ## block bootstrap
				cluster.boot<-sample(clusters,length(clusters),replace=TRUE)
				smp<-unlist(id.list[match(cluster.boot,clusters)])
			}   
			data.boot <- data[smp,]
			
			boot.out <- matrix(NA,nrow=all.length,ncol=0)
		
			### check input...
			if(treat.type=='discrete'){
				if(length(unique(data.boot[,D]))!=length(unique(data[,D]))){
					return(boot.out)
				}
			}
			
			if(use_fe==FALSE & is.null(IV)){
				## Linear Part
				if(method=='linear'){
					suppressWarnings(
						model.boot <- glm(formula,data=data.boot,weights=WEIGHTS)
					)
				}
				if(method=='logit'){
					suppressWarnings(
						model.boot <- glm(formula,data=data.boot,family=binomial(link = 'logit'),weights=WEIGHTS)
					)
				}
				if(method=='probit'){
					suppressWarnings(
						model.boot <- glm(formula,data=data.boot,family=binomial(link = 'probit'),weights=WEIGHTS)
					)
				}
				if(method=='poisson'){
					suppressWarnings(
						model.boot <- glm(formula,data=data.boot,family=poisson,weights=WEIGHTS)
					)
				}
				if(method=='nbinom'){
					suppressWarnings(
						model.boot <- glm.nb(formula,data=data.boot,weights=WEIGHTS)
					)
				}
			
				#### check converge...
				if(model.boot$converged==FALSE){
					return(boot.out)
				}
			}

			if(use_fe==TRUE & is.null(IV)){
				w.boot <- data.boot[,'WEIGHTS']
				suppressWarnings(
					model.boot <- felm(formula,data=data.boot,weights=w.boot)
				)
			}

			if(use_fe==FALSE & !is.null(IV)){
				suppressWarnings(
					model.boot <- ivreg(formula,data=data.boot,weights=WEIGHTS)
				)
			}

			if(use_fe==TRUE & !is.null(IV)){
				w.boot <- data.boot[,'WEIGHTS']
				suppressWarnings(
					model.boot <- felm(formula,data=data.boot,weights=w.boot)
				)
			}
			
			coef.boot <- coef(model.boot)
			if(!is.null(IV)){
	  			if(use_fe==1){
	    			names(coef.boot) <- name.update
	  			}
  			}
			
			## Binning Part
			groupX.boot <- cut(data.boot[,X],breaks=cuts.X,labels=FALSE,include.lowest = TRUE)
			
			data.boot.binning <- data.boot
			for(i in 1:nbins){
				data.boot.binning[,paste0('G.',i)] <- as.numeric(groupX.boot==i)
				data.boot.binning[,paste0('GX.',i)] <- as.numeric(groupX.boot==i)*(data.boot[,X]-x0[i])
			}
	
			if(treat.type=='discrete'){
				for(char in other.treat){
					for(i in 1:nbins){
						data.boot.binning[,paste0('D.',char,'.','G.',i)] <- as.numeric(groupX.boot==i)*as.numeric(data.boot[,D]==char)
						data.boot.binning[,paste0('D.',char,'.','GX.',i)] <- as.numeric(groupX.boot==i)*as.numeric(data.boot[,D]==char)*(data.boot[,X]-x0[i])
					}
				}
			}
	
			if(treat.type=='continuous'){
				for(i in 1:nbins){
					data.boot.binning[,paste0('D.G.',i)] <- as.numeric(groupX.boot==i)*data.boot[,D]
					data.boot.binning[,paste0('D.GX.',i)] <- as.numeric(groupX.boot==i)*data.boot[,D]*(data.boot[,X]-x0[i])
				}
			}

			if(full.moderate==TRUE){
				#formula.binning <- paste0(formula.binning,"+",paste0(Z.X,collapse="+"))
				for(a in Z){
					for(i in 1:nbins){
						data.boot.binning[,paste0(a,'.G.',i)] <- as.numeric(groupX.boot==i)*data.boot[,a]
						data.boot.binning[,paste0(a,'.GX.',i)] <- as.numeric(groupX.boot==i)*data.boot[,a]*(data.boot[,X]-x0[i])
					}
				}
			}

			if(!is.null(IV)){
				for(sub.iv in IV){
					for(i in 1:nbins){
						data.boot.binning[,paste0(sub.iv,'.','G.',i)] <- as.numeric(groupX.boot==i)*as.numeric(data.boot[,sub.iv])
						data.boot.binning[,paste0(sub.iv,'.','GX.',i)] <- as.numeric(groupX.boot==i)*as.numeric(data.boot[,sub.iv])*(data.boot[,X]-x0[i])
					}
				}
			}

			if(use_fe==FALSE & is.null(IV)){
				if(method=='linear'){
					suppressWarnings(
						model.boot.binning <- glm(formula.binning,data=data.boot.binning,weights=WEIGHTS)
					)
				}
				if(method=='logit'){
					suppressWarnings(
						model.boot.binning <- glm(formula.binning,data=data.boot.binning,family=binomial(link = 'logit'),weights=WEIGHTS)
					)
				}
				if(method=='probit'){
					suppressWarnings(
						model.boot.binning <- glm(formula.binning,data=data.boot.binning,family=binomial(link = 'probit'),weights=WEIGHTS)
					)
				}
				if(method=='poisson'){
					suppressWarnings(
						model.boot.binning <- glm(formula.binning,data=data.boot.binning,family=poisson,weights=WEIGHTS)
					)
				}
				if(method=='nbinom'){
					suppressWarnings(
						model.boot.binning <- glm.nb(formula.binning,data=data.boot.binning,weights=WEIGHTS)
					)
				}
				
				if(model.binning$converged==FALSE){
					no.converge <- 1
					model.binning.coef.boot <- NULL
				}else{
					no.converge <- 0
					model.binning.coef.boot <- coef(model.boot.binning) #keep the NaN(s) in coefficient
				}
			}
			
			if(use_fe==TRUE & is.null(IV)){
				suppressWarnings(
					model.boot.binning <- felm(formula.binning,data=data.boot.binning,weights=w.boot)
				)
				model.binning.coef.boot <- coef(model.boot.binning)
				no.converge <- 0
			}

			if(use_fe==FALSE & !is.null(IV)){
				suppressWarnings(
					model.boot.binning <- ivreg(formula.binning,data=data.boot.binning,weights=WEIGHTS)
				)
				model.binning.coef.boot <- coef(model.boot.binning)
				no.converge <- 0
			}

			if(use_fe==TRUE & !is.null(IV)){
				suppressWarnings(
					model.boot.binning <- felm(formula.binning,data=data.boot.binning,weights=w.boot)
				)
				model.binning.coef.boot <- coef(model.boot.binning)
				names(model.binning.coef.boot) <- iv.reg.name
				no.converge <- 0
			}
			
			boot.one.round <- c()
		
			if(treat.type=='discrete'){
				for(char in other.treat){
					gen.general.TE.output <- gen.general.TE(model.coef=coef.boot,char=char)
					TE.output <- gen.general.TE.output$TE
					names(TE.output) <- rep(paste0("TE.",char),neval)
					
					if(no.converge==0){
						TE.binning.output <- gen.binning.TE(model.binning.coef=model.binning.coef.boot,
															char=char)
													
						pred.binning.output.list <- gen.pred.binning(model.binning.coef=model.binning.coef.boot, 
																X.eval=X.eval, 
																cuts.X=cuts.X, 
																x0=x0, 
																char=char)
						pred.binning.output <- pred.binning.output.list$E.pred
						base.binning.output <- pred.binning.output.list$E.base
						link.1.binning.output <- pred.binning.output.list$link.1
						link.0.binning.output <- pred.binning.output.list$link.0
					}else{
						TE.binning.output <- rep(NA,nbins)
						pred.binning.output <- rep(NA,length(X.eval))
						base.binning.output <- rep(NA,length(X.eval))
						link.1.binning.output <- rep(NA,length(X.eval))
						link.0.binning.output <- rep(NA,length(X.eval))
					}
					names(TE.binning.output) <- rep(paste0("bin.est.",char),nbins)
					names(pred.binning.output) <- rep(paste0("bin.pred.",char),neval)
					names(link.1.binning.output) <- rep(paste0("link.",char),neval)
					boot.one.round <- c(boot.one.round,TE.output,TE.binning.output,pred.binning.output,link.1.binning.output)
				}
				names(base.binning.output) <- rep(paste0("bin.pred.",base),neval)
				boot.one.round <- c(boot.one.round,base.binning.output)
				names(link.0.binning.output) <- rep(paste0("link.",base),neval)
				boot.one.round <- c(boot.one.round,link.0.binning.output)
			}
			
			if(treat.type=='continuous'){
				k <- 1
				for(D.ref in D.sample){
					gen.general.ME.output <- gen.general.TE(model.coef=coef.boot,D.ref=D.ref)
					ME.output <- gen.general.ME.output$ME
					names(ME.output) <- rep(paste0("ME.",label.name[k]),neval)
					if(no.converge==0){
						ME.binning.output <- gen.binning.TE(model.binning.coef=model.binning.coef.boot,
															D.ref=D.ref)
													
						pred.binning.output.list <- gen.pred.binning(model.binning.coef=model.binning.coef.boot, 
																	 X.eval=X.eval, 
																	 cuts.X=cuts.X, 
																	 x0=x0, 
																	 D.ref=D.ref)
																	 
						pred.binning.output <- pred.binning.output.list$E.pred
						link.binning.output <- pred.binning.output.list$link
						
					}else{
						ME.binning.output <- rep(NA,nbins)
						pred.binning.output <- rep(NA,length(X.eval))
						link.binning.output <- rep(NA,length(X.eval))
					}
					names(ME.binning.output) <- rep(paste0("bin.est.",label.name[k]),nbins)
					names(pred.binning.output) <- rep(paste0("bin.pred.",label.name[k]),neval)
					names(link.binning.output) <- rep(paste0("link.",label.name[k]),neval)
					boot.one.round <- c(boot.one.round,ME.output,ME.binning.output,pred.binning.output,link.binning.output)	
					k <- k + 1
				}
			}
			boot.out <- cbind(boot.out,boot.one.round)
			rownames(boot.out) <- names(boot.one.round)
			colnames(boot.out) <- NULL
			return(boot.out)
		}
		
		
		if(parallel==TRUE){
			requireNamespace("doParallel")
			## require(iterators)
			maxcores <- detectCores()
			cores <- min(maxcores, cores)
			pcl <-future::makeClusterPSOCK(cores)  
			doParallel::registerDoParallel(pcl)
			cat("Parallel computing with", cores,"cores...\n") 

			suppressWarnings(
				bootout <- foreach (i=1:nboots, .combine=cbind,
									.export=c("one.boot"),.packages=c('lfe','AER'),
									.inorder=FALSE) %dopar% {one.boot()}
			) 
			suppressWarnings(stopCluster(pcl))
			cat("\r")
		}	 
		else{
			bootout<-matrix(NA,all.length,0)
			for(i in 1:nboots){
				tempdata <- one.boot()
				if(is.null(tempdata)==FALSE){
					bootout<- cbind(bootout,tempdata)
				}
				if (i%%50==0) cat(i) else cat(".")
			}
			cat("\r")
		}

		if(treat.type=='discrete'){
			TE.output.all.list <- list()
			TE.binning.output.all.list <- list()
			pred.output.all.list <- list()
			link.output.all.list <- list()
			for(char in other.treat){
				gen.general.TE.output <- gen.general.TE(model.coef=model.coef,char=char)
				TE.output <- gen.general.TE.output$TE
				TE.binning.output <- gen.binning.TE(model.binning.coef=model.binning.coef,
													char=char)
				pred.binning.output <- gen.pred.binning(model.binning.coef=model.binning.coef, 
														X.eval=X.eval, 
														cuts.X=cuts.X, 
														x0=x0,char=char)
				
				TE.boot.matrix <- bootout[rownames(bootout)==paste0("TE.",char),]
				TE.binning.boot.matrix <- bootout[rownames(bootout)==paste0("bin.est.",char),]
				pred.boot.matrix <- bootout[rownames(bootout)==paste0("bin.pred.",char),]
				base.boot.matrix <- bootout[rownames(bootout)==paste0("bin.pred.",base),]
				link.boot.matrix <- bootout[rownames(bootout)==paste0("link.",char),]
				link0.boot.matrix <- bootout[rownames(bootout)==paste0("link.",base),]
				
				TE.boot.sd <- apply(TE.boot.matrix, 1, sd, na.rm=TRUE)
				TE.binning.sd <- apply(TE.binning.boot.matrix, 1, sd, na.rm=TRUE)
				pred.boot.sd <- apply(pred.boot.matrix, 1, sd, na.rm=TRUE)
				base.boot.sd <- apply(base.boot.matrix, 1, sd, na.rm=TRUE)
				link.boot.sd <- apply(link.boot.matrix, 1, sd, na.rm=TRUE)
				link0.boot.sd <- apply(link0.boot.matrix, 1, sd, na.rm=TRUE)
				
				TE.boot.CI <- t(apply(TE.boot.matrix, 1, quantile, c(0.025,0.975), na.rm=TRUE))
				TE.binning.CI <- t(apply(TE.binning.boot.matrix, 1, quantile, c(0.025,0.975),na.rm=TRUE))
				pred.boot.CI <- t(apply(pred.boot.matrix, 1, quantile, c(0.025,0.975),na.rm=TRUE))
				base.boot.CI <- t(apply(base.boot.matrix, 1, quantile, c(0.025,0.975),na.rm=TRUE))
				link.boot.CI <- t(apply(link.boot.matrix, 1, quantile, c(0.025,0.975),na.rm=TRUE))
				link0.boot.CI <- t(apply(link0.boot.matrix, 1, quantile, c(0.025,0.975),na.rm=TRUE))
				
				TE.output.all <- cbind(X.eval,TE.output,TE.boot.sd,TE.boot.CI[,1],TE.boot.CI[,2])
				colnames(TE.output.all) <- c("X","TE","sd","lower CI(95%)","upper CI(95%)")
				rownames(TE.output.all) <- NULL
				TE.output.all.list[[other.treat.origin[char]]] <- TE.output.all
				
				E.pred.output <- pred.binning.output$E.pred
				pred.output.all <- cbind(X.eval, E.pred.output, pred.boot.sd, pred.boot.CI[,1], pred.boot.CI[,2])
				colnames(pred.output.all) <- c("X","E(Y)","sd","lower CI(95%)","upper CI(95%)")
				rownames(pred.output.all) <- NULL
				pred.output.all.list[[other.treat.origin[char]]] <- pred.output.all

				link.output <- pred.binning.output$link.1
				link.output.all <- cbind(X.eval, link.output, link.boot.sd, link.boot.CI[,1], link.boot.CI[,2])
				colnames(link.output.all) <- c("X","E(Y)","sd","lower CI(95%)","upper CI(95%)")
				rownames(link.output.all) <- NULL
				link.output.all.list[[other.treat.origin[char]]] <- link.output.all
				
				TE.binning.all <- cbind(x0,TE.binning.output,TE.binning.sd,TE.binning.CI[,1],TE.binning.CI[,2])
				colnames(TE.binning.all) <- c('x0',"coef","sd","CI.lower","CI.upper")
				rownames(TE.binning.all) <- gp.lab
				TE.binning.output.all.list[[other.treat.origin[char]]] <- TE.binning.all 
														
			}
			
			E.base.output <- pred.binning.output$E.base
			base.output.all <- cbind(X.eval, E.base.output, base.boot.sd, base.boot.CI[,1], base.boot.CI[,2])
			colnames(base.output.all) <- c("X","E(Y)","sd","lower CI(95%)","upper CI(95%)")
			rownames(base.output.all) <- NULL
			pred.output.all.list[[all.treat.origin[base]]] <- base.output.all

			link0.output <- pred.binning.output$link.0
			link0.output.all <- cbind(X.eval, link0.output, link0.boot.sd, link0.boot.CI[,1], link0.boot.CI[,2])
			colnames(link0.output.all) <- c("X","E(Y)","sd","lower CI(95%)","upper CI(95%)")
			rownames(link0.output.all) <- NULL
			link.output.all.list[[all.treat.origin[base]]] <- link0.output.all
		}
		
		if(treat.type=='continuous'){
			ME.output.all.list <- list()
			ME.binning.output.all.list <- list()
			pred.output.all.list <- list()
			link.output.all.list <- list()
			k <- 1
			for(D.ref in D.sample){
				label <- label.name[k]
				gen.general.ME.output <- gen.general.TE(model.coef=model.coef,D.ref=D.ref)
				ME.output <- gen.general.ME.output$ME
				ME.binning.output <- gen.binning.TE(model.binning.coef=model.binning.coef,
													D.ref=D.ref)
				pred.binning.output <- gen.pred.binning(model.binning.coef=model.binning.coef, 
														X.eval=X.eval, 
														cuts.X=cuts.X, 
														x0=x0, 
														D.ref=D.ref)
												
				ME.boot.matrix <- bootout[rownames(bootout)==paste0("ME.",label),]
				ME.binning.boot.matrix <- bootout[rownames(bootout)==paste0("bin.est.",label),]
				pred.boot.matrix <- bootout[rownames(bootout)==paste0("bin.pred.",label),]		
				link.boot.matrix <- bootout[rownames(bootout)==paste0("link.",label),]

				ME.boot.sd <- apply(ME.boot.matrix, 1, sd, na.rm=TRUE)
				ME.binning.sd <- apply(ME.binning.boot.matrix, 1, sd, na.rm=TRUE)
				pred.boot.sd <- apply(pred.boot.matrix, 1, sd, na.rm=TRUE)
				link.boot.sd <- apply(link.boot.matrix, 1, sd, na.rm=TRUE)
				
				

				ME.boot.CI <- t(apply(ME.boot.matrix, 1, quantile, c(0.025,0.975),na.rm=TRUE))
				ME.binning.CI <- t(apply(ME.binning.boot.matrix, 1, quantile, c(0.025,0.975),na.rm=TRUE))
				pred.boot.CI <- t(apply(pred.boot.matrix, 1, quantile, c(0.025,0.975),na.rm=TRUE))
				link.boot.CI <- t(apply(link.boot.matrix, 1, quantile, c(0.025,0.975),na.rm=TRUE))
					
				ME.output.all <- cbind(X.eval,ME.output,ME.boot.sd,ME.boot.CI[,1],ME.boot.CI[,2])
				colnames(ME.output.all) <- c("X","ME","sd","lower CI(95%)","upper CI(95%)")
				ME.output.all.list[[label]] <- ME.output.all
				
				E.pred.output <- pred.binning.output$E.pred
				pred.output.all <- cbind(X.eval, E.pred.output, pred.boot.sd, pred.boot.CI[,1], pred.boot.CI[,2])
				colnames(pred.output.all) <- c("X","E(Y)","sd","lower CI(95%)","upper CI(95%)")
				rownames(pred.output.all) <- NULL
				pred.output.all.list[[label]] <- pred.output.all

				link.output <- pred.binning.output$link
				link.output.all <- cbind(X.eval, link.output, link.boot.sd, link.boot.CI[,1], link.boot.CI[,2])
				colnames(link.output.all) <- c("X","E(Y)","sd","lower CI(95%)","upper CI(95%)")
				rownames(link.output.all) <- NULL
				link.output.all.list[[label]] <- link.output.all
				
				ME.binning.all <- cbind(x0,ME.binning.output,ME.binning.sd,ME.binning.CI[,1],ME.binning.CI[,2])
				colnames(ME.binning.all) <- c('x0',"coef","sd","CI.lower","CI.upper")
				rownames(ME.binning.all) <- gp.lab
				ME.binning.output.all.list[[label]] <- ME.binning.all
				k <- k+1
			}
		
		}
	}
		
	
	if(TRUE){ # density or histogram
	if (treat.type=='discrete'){ ## discrete D
    # density
    if(is.null(weights)==TRUE){
      de <- density(data[,X])
    }else {
      suppressWarnings(de <- density(data[,X],weights=data[,'WEIGHTS']))
    }
    
    treat_den <- list()
    for (char in all.treat){
      de.name <- paste0("den.",char)
      if (is.null(weights)==TRUE){
        de.tr <- density(data[data[,D]==char,X])
      } 
      else {
        suppressWarnings(de.tr <- density(data[data[,D]==char,X],
                                          weights=data[data[,D]==char,'WEIGHTS']))
      }
      treat_den[[all.treat.origin[char]]] <- de.tr
    }
    
    # histogram
    if (is.null(weights)==TRUE) {
      hist.out<-hist(data[,X],breaks=80,plot=FALSE)
    } else {
      suppressWarnings(hist.out<-hist(data[,X],data[,'WEIGHTS'],
                                      breaks=80,plot=FALSE))
    } 
    n.hist<-length(hist.out$mids)

    # count the number of treated
    treat.hist <- list()
    for (char in all.treat) {
      count1<-rep(0,n.hist)
      treat_index<-which(data[,D]==char)
      for (i in 1:n.hist) {
        count1[i]<-sum(data[treat_index,X]>=hist.out$breaks[i] &
                         data[treat_index,X]<hist.out$breaks[(i+1)])
      }
      count1[n.hist]<-sum(data[treat_index,X]>=hist.out$breaks[n.hist] &
                            data[treat_index,X]<=hist.out$breaks[n.hist+1])
      
      treat.hist[[all.treat.origin[char]]] <- count1
    }    
  }  
  
  if (treat.type=='continuous'){ ## continuous D
    if (is.null(weights)==TRUE){
      de <- density(data[,X])
    } else {
      suppressWarnings(de <- density(data[,X],weights=data[,'WEIGHTS']))
    }
    if (is.null(weights)==TRUE){
      hist.out<-hist(data[,X],breaks=80,plot=FALSE)
    } else{
      suppressWarnings(hist.out<-hist(data[,X],data[,'WEIGHTS'],
                                      breaks=80,plot=FALSE))
    }
    de.co <- de.tr <- NULL 
  } 
}

	
	tests <- NULL
	
	#wald test: doesn't allow nbinom by far
	if(wald==TRUE & use_fe==FALSE & is.null(IV)){
		data.wald <- data
		sub.test <- NULL

		if(treat.type=='continuous'){
			formula0.wald <- paste0(Y,"~",X,"+",D,"+DX")
			formula1.wald <- formula0.wald
			var.name <- c(X,D,"DX")
			for(i in 1:(nbins-1)){
				data.wald[,paste0("G.",i+1)] <- as.numeric(groupX==(i+1))
				data.wald[,paste0("G.",i+1,".X")] <- as.numeric(groupX==(i+1))*data.wald[,X]
				data.wald[,paste0("DG.",i+1)] <- as.numeric(groupX==(i+1))*data.wald[,D]
				data.wald[,paste0("DG.",i+1,".X")] <- as.numeric(groupX==(i+1))*data.wald[,D]*data.wald[,X]
				new.var.name <- c(paste0("G.",i+1),paste0("G.",i+1,".X"),paste0("DG.",i+1),paste0("DG.",i+1,".X"))
				var.name <- c(var.name,new.var.name)
				formula1.wald <- paste0(formula1.wald,"+",paste0(new.var.name,collapse="+"))
			}

			if (is.null(Z)==FALSE){
				if(full.moderate==FALSE){
					formula0.wald <- paste0(formula0.wald, "+",paste(Z,collapse=" + "))
					formula1.wald <- paste0(formula1.wald, "+",paste(Z,collapse=" + "))
				}
				if(full.moderate==TRUE){ #z
					formula0.wald <- paste0(formula0.wald, "+",paste(Z,collapse=" + "))
					formula1.wald <- paste0(formula1.wald, "+",paste(Z,collapse=" + "))
					formula0.wald <- paste0(formula0.wald, "+",paste(Z.X,collapse=" + "))
					formula1.wald <- paste0(formula1.wald, "+",paste(Z.X,collapse=" + "))
					zero.var.name <- c()
					for(a in Z){
						for(i in 1:(nbins-1)){
							data.wald[,paste0("Z.",a,".G.",i+1)] <- as.numeric(groupX==(i+1))*data.wald[,a]
							data.wald[,paste0("ZX.",a,".G.",i+1)] <- as.numeric(groupX==(i+1))*data.wald[,a]*data.wald[,X]
							zero.var.name <- c(zero.var.name,paste0("Z.",a,".G.",i+1),paste0("ZX.",a,".G.",i+1))
						}
					}
					formula0.wald <- paste0(formula0.wald,"+",paste0(zero.var.name,collapse="+"))
					formula1.wald <- paste0(formula1.wald,"+",paste0(zero.var.name,collapse="+"))
				}
			}

			formula0.wald <- as.formula(formula0.wald)
			formula1.wald <- as.formula(formula1.wald)
			
			if(method=='linear'){
				suppressWarnings(
					fit0 <- lm(formula0.wald,data=data.wald,weights=WEIGHTS)
				)
				suppressWarnings(
					fit1 <- lm(formula1.wald,data=data.wald,weights=WEIGHTS)
				)
			}
			if(method=='logit'){
				suppressWarnings(
					fit0 <- glm(formula0.wald,data=data.wald,family=binomial(link = 'logit'),weights=WEIGHTS)
				)
				suppressWarnings(
					fit1 <- glm(formula1.wald,data=data.wald,family=binomial(link = 'logit'),weights=WEIGHTS)
				)
			}
			
			if(method=='probit'){
				suppressWarnings(
					fit0 <- glm(formula0.wald,data=data.wald,family=binomial(link = 'probit'),weights=WEIGHTS)
				)
				suppressWarnings(
					fit1 <- glm(formula1.wald,data=data.wald,family=binomial(link = 'probit'),weights=WEIGHTS)
				)
			}
			if(method=='poisson'){
				suppressWarnings(
					fit0 <- glm(formula0.wald,data=data.wald,family=poisson,weights=WEIGHTS)
				)
				suppressWarnings(
					fit1 <- glm(formula1.wald,data=data.wald,family=poisson,weights=WEIGHTS)
				)
			}
			if(method=='nbinom'){
				suppressWarnings(
					fit0 <- glm.nb(formula0.wald,data=data.wald,weights=WEIGHTS)
				)
				suppressWarnings(
					fit1 <- glm.nb(formula1.wald,data=data.wald,weights=WEIGHTS)
				)
			}

			#vcov
			if(vcov.type=="homoscedastic"){
				model.wald.vcov <- vcov(fit1)
			}
			if(vcov.type=="robust"){
				model.wald.vcov <- vcovHC(fit1,type="HC1")
			} 
			if(vcov.type=="cluster"){
				model.wald.vcov <- vcovCluster(fit1,cluster = data.wald[,cl])
			}
			if(vcov.type=="pcse"){
				model.wald.vcov <- pcse(fit1,pairwise=pairwise,groupN=data.wald[,'cl'],groupT=data.wald[,'time'])$vcov
				rownames(model.wald.vcov)[1] <- "(Intercept)"
				colnames(model.wald.vcov)[1] <- "(Intercept)"
			}

			model.wald.coef <- coef(fit1)
			model.wald.vcov[which(is.na(model.wald.vcov))] <- 0
			model.wald.vcov.all <- matrix(0,nrow=length(model.wald.coef),ncol=length(model.wald.coef))

			colnames(model.wald.vcov.all) <- names(model.wald.coef)
			rownames(model.wald.vcov.all) <- names(model.wald.coef)
			for(a1 in rownames(model.wald.vcov)){
				for(a2 in colnames(model.wald.vcov)){
					model.wald.vcov.all[a1,a2] <- model.wald.vcov[a1,a2]
				}
			}
			
			#lr.test
			lr.result <- tryCatch(
				round(lmtest::lrtest(fit1,fit0)[[5]][2],3),
				error = function(e){return(NULL)}
			)
			
			#wald.test
			wald.result <- tryCatch(
				round(lmtest::waldtest(fit1,fit0,test="Chisq", vcov=model.wald.vcov.all)[[4]][2],3),
				error = function(e){return(NULL)}
			)
		}
		
		if(treat.type=='discrete'){
			formula0.wald <- formula.origin
			formula1.wald <- formula0.wald
			for(i in 1:(nbins-1)){
				data.wald[,paste0("G.",i+1)] <- as.numeric(groupX==(i+1))
				data.wald[,paste0("G.",i+1,".X")] <- as.numeric(groupX==(i+1))*data.wald[,X]
				new.var.name <- c(paste0("G.",i+1),paste0("G.",i+1,".X"))
				for(char in other.treat){
					data.wald[,paste0("D.",char,".G.",i+1)] <- as.numeric(groupX==(i+1))*data.wald[,paste0("D.",char)]
					data.wald[,paste0("DX.",char,".G.",i+1)] <- as.numeric(groupX==(i+1))*data.wald[,paste0("D.",char)]*data.wald[,X]
					new.var.name <- c(new.var.name,paste0("D.",char,".G.",i+1),paste0("DX.",char,".G.",i+1))
				}
				if(is.null(Z)==FALSE & full.moderate==TRUE){ #z
					zero.var.name <- c()
					for(a in Z){
						data.wald[,paste0("Z.",a,".G.",i+1)] <- as.numeric(groupX==(i+1))*data.wald[,a]
						data.wald[,paste0("ZX.",a,".G.",i+1)] <- as.numeric(groupX==(i+1))*data.wald[,a]*data.wald[,X]
						new.var.name <- c(new.var.name,paste0("Z.",a,".G.",i+1),paste0("ZX.",a,".G.",i+1))
						zero.var.name <- c(zero.var.name,paste0("Z.",a,".G.",i+1),paste0("ZX.",a,".G.",i+1))
					}
					formula0.wald <- paste0(formula0.wald,"+",paste0(zero.var.name,collapse="+"))
				}
				formula1.wald <- paste0(formula1.wald,"+",paste0(new.var.name,collapse="+"))
			}
			
			formula0.wald <- as.formula(formula0.wald)
			formula1.wald <- as.formula(formula1.wald)

			if(method=='linear'){
				suppressWarnings(
					fit0 <- glm(formula0.wald,data=data.wald,weights=WEIGHTS)
				)
				suppressWarnings(
					fit1 <- glm(formula1.wald,data=data.wald,weights=WEIGHTS)
				)
			}
			if(method=='logit'){
				suppressWarnings(
					fit0 <- glm(formula0.wald,data=data.wald,family=binomial(link = 'logit'),weights=WEIGHTS)
				)
				suppressWarnings(
					fit1 <- glm(formula1.wald,data=data.wald,family=binomial(link = 'logit'),weights=WEIGHTS)
				)
			}
			
			if(method=='probit'){
				suppressWarnings(
					fit0 <- glm(formula0.wald,data=data.wald,family=binomial(link = 'probit'),weights=WEIGHTS)
				)
				suppressWarnings(
					fit1 <- glm(formula1.wald,data=data.wald,family=binomial(link = 'probit'),weights=WEIGHTS)
				)
			}
			if(method=='poisson'){
				suppressWarnings(
					fit0 <- glm(formula0.wald,data=data.wald,family=poisson,weights=WEIGHTS)
				)
				suppressWarnings(
					fit1 <- glm(formula1.wald,data=data.wald,family=poisson,weights=WEIGHTS)
				)
			}
			if(method=='nbinom'){
				suppressWarnings(
					fit0 <- glm.nb(formula0.wald,data=data.wald,weights=WEIGHTS)
				)
				suppressWarnings(
					fit1 <- glm.nb(formula1.wald,data=data.wald,weights=WEIGHTS)
				)
			}
			
			#requireNamespace('mdscore')
			#vcov
			if(vcov.type=="homoscedastic"){
				model.wald.vcov <- vcov(fit1)
			}
			if(vcov.type=="robust"){
				model.wald.vcov <- vcovHC(fit1,type="HC1")
			} 
			if(vcov.type=="cluster"){
				model.wald.vcov <- vcovCluster(fit1,cluster = data.wald[,cl])
			}
			if(vcov.type=="pcse"){
				model.wald.vcov <- pcse(fit1,pairwise=pairwise,groupN=data.wald[,'cl'],groupT=data.wald[,'time'])$vcov
				rownames(model.wald.vcov)[1] <- "(Intercept)"
				colnames(model.wald.vcov)[1] <- "(Intercept)"
			}

			model.wald.coef <- coef(fit1)
			model.wald.vcov[which(is.na(model.wald.vcov))] <- 0
			model.wald.vcov.all <- matrix(0,nrow=length(model.wald.coef),ncol=length(model.wald.coef))

			colnames(model.wald.vcov.all) <- names(model.wald.coef)
			rownames(model.wald.vcov.all) <- names(model.wald.coef)
			for(a1 in rownames(model.wald.vcov)){
				for(a2 in colnames(model.wald.vcov)){
					model.wald.vcov.all[a1,a2] <- model.wald.vcov[a1,a2]
				}
			}
			
			#lr.test
			lr.result <- tryCatch(
				round(lmtest::lrtest(fit1,fit0)[[5]][2],3),
				error = function(e){return(NULL)}
			)
			
			#wald.test
			wald.result <- tryCatch(
				round(lmtest::waldtest(fit1,fit0,test="Chisq", vcov=model.wald.vcov.all)[[4]][2],3),
				error = function(e){return(NULL)}
			)

			
			
			if(length(other.treat)>1){ #more than one treated group
				sub.test <- list()
				for(target.char in other.treat){
					formula0.wald.sub <- formula.origin
					new.var.name <- c()
					for(i in 1:(nbins-1)){	
						new.var.name <- c(new.var.name,paste0("G.",i+1),paste0("G.",i+1,".X"))	
						for(char in other.treat){
							if(char!=target.char){
								new.var.name <- c(new.var.name,paste0("D.",char,".G.",i+1),paste0("DX.",char,".G.",i+1))
							}
						}
						if(is.null(Z)==FALSE & full.moderate==TRUE){ #z
							for(a in Z){
								new.var.name <- c(new.var.name,paste0("Z.",a,".G.",i+1),paste0("ZX.",a,".G.",i+1))
							}
						}
					}
					formula0.wald.sub <- paste0(formula0.wald.sub,"+",paste0(new.var.name,collapse="+"))
					
					formula0.wald.sub <- as.formula(formula0.wald.sub)
					#print(formula0.wald.sub)
					#print(formula1.wald)
					if(method=='linear'){
						suppressWarnings(
							fit0.sub <- glm(formula0.wald.sub,data=data.wald,weights=WEIGHTS)
						)
					}
					if(method=='logit'){
						suppressWarnings(
							fit0.sub <- glm(formula0.wald.sub,data=data.wald,family=binomial(link = 'logit'),weights=WEIGHTS)
						)
					}
					if(method=='probit'){
						suppressWarnings(
							fit0.sub <- glm(formula0.wald.sub,data=data.wald,family=binomial(link = 'probit'),weights=WEIGHTS)
						)
					}
					if(method=='poisson'){
						suppressWarnings(
							fit0.sub <- glm(formula0.wald.sub,data=data.wald,family=poisson,weights=WEIGHTS)
						)
					}
					if(method=='nbinom'){
						suppressWarnings(
							fit0.sub <- glm.nb(formula0.wald.sub,data=data.wald,weights=WEIGHTS)
						)
					}
					

					#lr.test
					lr.result.sub <- tryCatch(
						round(lmtest::lrtest(fit1,fit0.sub)[[5]][2],3),
						error = function(e){return(NULL)}
					)
			
					#wald.test
					wald.result.sub <- tryCatch(
						round(lmtest::waldtest(fit1,fit0.sub,test="Chisq", vcov=model.wald.vcov.all)[[4]][2],3),
						error = function(e){return(NULL)}
					)
					
					sub.test[[other.treat.origin[target.char]]] <- list(p.wald = sprintf(wald.result.sub,fmt = '%#.3f'),
																		p.lr = sprintf(lr.result.sub,fmt = '%#.3f'),
																		formula.restrict = formula0.wald.sub,
																		formula.full = formula1.wald)

				}
			
			}		
		}
	
		requireNamespace("Lmoments")
		Lkurtosis <- Lmoments(data[,X],returnobject=TRUE)$ratios[4] 
		tests <- list(
					treat.type = treat.type,
					X.Lkurtosis = sprintf(Lkurtosis,fmt = '%#.3f'),
					p.wald = sprintf(wald.result,fmt = '%#.3f'),
					p.lr = sprintf(lr.result,fmt = '%#.3f'),
					formula.restrict = formula0.wald,
					formula.full = formula1.wald,
					sub.test = sub.test
				)
	}
	
	if(wald==TRUE & use_fe==FALSE & !is.null(IV)){
		# test 2 ivreg fit
		data.wald <- data
		sub.test <- NULL

		if(treat.type=='continuous'){
			# fit0 
			formula0.wald <- paste0(Y,"~",X,"+",D,"+DX")
			formula0.iv <- X
			for(sub.iv in IV){
				data.wald[,paste0(X,".",sub.iv)] <- data.wald[,sub.iv]*data.wald[,X]
				formula0.iv <- paste0(formula0.iv,"+",sub.iv)
				formula0.iv <- paste0(formula0.iv,"+",paste0(X,".",sub.iv))
			}

			# fit1
			formula1.wald <- formula0.wald
			var.name <- c(X,D,"DX")
			for(i in 1:(nbins-1)){
				data.wald[,paste0("G.",i+1)] <- as.numeric(groupX==(i+1))
				data.wald[,paste0("G.",i+1,".X")] <- as.numeric(groupX==(i+1))*data.wald[,X]
				data.wald[,paste0("DG.",i+1)] <- as.numeric(groupX==(i+1))*data.wald[,D]
				data.wald[,paste0("DG.",i+1,".X")] <- as.numeric(groupX==(i+1))*data.wald[,D]*data.wald[,X]
				new.var.name <- c(paste0("G.",i+1),paste0("G.",i+1,".X"),paste0("DG.",i+1),paste0("DG.",i+1,".X"))
				var.name <- c(var.name,new.var.name)
				formula1.wald <- paste0(formula1.wald,"+",paste0(new.var.name,collapse="+"))
			}
			formula1.iv <- formula0.iv
				
			for(i in 1:(nbins-1)){
				new.var.name <- c(paste0("G.",i+1),paste0("G.",i+1,".X"))
				for(sub.iv in IV){
					data.wald[,paste0(sub.iv,".G.",i+1)] <- as.numeric(groupX==(i+1))*data.wald[,sub.iv]
					data.wald[,paste0(sub.iv,".G.",i+1,".X")] <- as.numeric(groupX==(i+1))*data.wald[,sub.iv]*data.wald[,X]
					new.var.name <- c(new.var.name,paste0(sub.iv,".G.",i+1),paste0(sub.iv,".G.",i+1,".X"))
				}
				formula1.iv <- paste0(formula1.iv,"+",paste0(new.var.name,collapse="+"))
			}

			# add covariates
			if (is.null(Z)==FALSE){
				if(full.moderate==FALSE){
					formula0.wald <- paste0(formula0.wald, "+",paste(Z,collapse=" + "))
					formula1.wald <- paste0(formula1.wald, "+",paste(Z,collapse=" + "))
					formula0.iv <- paste0(formula0.iv, "+",paste(Z,collapse=" + "))
					formula1.iv <- paste0(formula1.iv, "+",paste(Z,collapse=" + "))
				}
				if(full.moderate==TRUE){ #z
					formula0.wald <- paste0(formula0.wald, "+",paste(Z,collapse=" + "))
					formula1.wald <- paste0(formula1.wald, "+",paste(Z,collapse=" + "))
					formula0.wald <- paste0(formula0.wald, "+",paste(Z.X,collapse=" + "))
					formula1.wald <- paste0(formula1.wald, "+",paste(Z.X,collapse=" + "))
					formula0.iv <- paste0(formula0.iv, "+",paste(Z,collapse=" + "))
					formula1.iv <- paste0(formula1.iv, "+",paste(Z,collapse=" + "))
					formula0.iv <- paste0(formula0.iv, "+",paste(Z.X,collapse=" + "))
					formula1.iv <- paste0(formula1.iv, "+",paste(Z.X,collapse=" + "))
					zero.var.name <- c()
					for(a in Z){
						for(i in 1:(nbins-1)){
							data.wald[,paste0("Z.",a,".G.",i+1)] <- as.numeric(groupX==(i+1))*data.wald[,a]
							data.wald[,paste0("ZX.",a,".G.",i+1)] <- as.numeric(groupX==(i+1))*data.wald[,a]*data.wald[,X]
							zero.var.name <- c(zero.var.name,paste0("Z.",a,".G.",i+1),paste0("ZX.",a,".G.",i+1))
						}
					}
					formula0.wald <- paste0(formula0.wald,"+",paste0(zero.var.name,collapse="+"))
					formula1.wald <- paste0(formula1.wald,"+",paste0(zero.var.name,collapse="+"))
					formula0.iv <- paste0(formula0.iv,"+",paste0(zero.var.name,collapse="+"))
					formula1.iv <- paste0(formula1.iv,"+",paste0(zero.var.name,collapse="+"))
				}
			}

			formula0.wald <- paste0(formula0.wald,"|",formula0.iv)
			formula1.wald <- paste0(formula1.wald,"|",formula1.iv)
			formula0.wald <- as.formula(formula0.wald)
			formula1.wald <- as.formula(formula1.wald)
			
			suppressWarnings(
				fit0 <- ivreg(formula0.wald,data=data.wald,weights=WEIGHTS)
			)
			suppressWarnings(
				fit1 <- ivreg(formula1.wald,data=data.wald,weights=WEIGHTS)
			)

			#vcov
			if(vcov.type=="homoscedastic"){
				model.wald.vcov <- vcov(fit1)
			}
			if(vcov.type=="robust"){
				model.wald.vcov <- vcovHC(fit1,type="HC1")
			} 
			if(vcov.type=="cluster"){
				model.wald.vcov <- vcovCluster(fit1,cluster = data.wald[,cl])
			}
			if(vcov.type=="pcse"){
				model.wald.vcov <- pcse(fit1,pairwise=pairwise,groupN=data.wald[,'cl'],groupT=data.wald[,'time'])$vcov
				rownames(model.wald.vcov)[1] <- "(Intercept)"
				colnames(model.wald.vcov)[1] <- "(Intercept)"
			}

			model.wald.coef <- coef(fit1)
			model.wald.vcov[which(is.na(model.wald.vcov))] <- 0
			model.wald.vcov.all <- matrix(0,nrow=length(model.wald.coef),ncol=length(model.wald.coef))

			colnames(model.wald.vcov.all) <- names(model.wald.coef)
			rownames(model.wald.vcov.all) <- names(model.wald.coef)
			for(a1 in rownames(model.wald.vcov)){
				for(a2 in colnames(model.wald.vcov)){
					model.wald.vcov.all[a1,a2] <- model.wald.vcov[a1,a2]
				}
			}
			
			#lr.test
			lr.result <- tryCatch(
				round(lmtest::lrtest(fit1,fit0)[[5]][2],3),
				error = function(e){return(NULL)}
			)
			
			#wald.test
			wald.result <- tryCatch(
				round(lmtest::waldtest(fit1,fit0,test="Chisq", vcov=model.wald.vcov.all)[[4]][2],3),
				error = function(e){return(NULL)}
			)
		}

		if(treat.type=='discrete'){
			#overall test
			# fit0
			formula0.wald <- paste0(Y,"~",X)
			for(char in other.treat){
				formula0.wald <- paste0(formula0.wald,"+",paste0("D",".",char),"+",paste0("DX",".",char))
			}
			formula0.iv <- X
			for(sub.iv in IV){
				data.wald[,paste0(X,".",sub.iv)] <- data.wald[,sub.iv]*data.wald[,X]
				formula0.iv <- paste0(formula0.iv,"+",sub.iv)
				formula0.iv <- paste0(formula0.iv,"+",paste0(X,".",sub.iv))
			}
			#fit1
			formula1.wald <- formula0.wald
			for(i in 1:(nbins-1)){
				data.wald[,paste0("G.",i+1)] <- as.numeric(groupX==(i+1))
				data.wald[,paste0("G.",i+1,".X")] <- as.numeric(groupX==(i+1))*data.wald[,X]
				new.var.name <- c(paste0("G.",i+1),paste0("G.",i+1,".X"))
				for(char in other.treat){
					data.wald[,paste0("D.",char,".G.",i+1)] <- as.numeric(groupX==(i+1))*data.wald[,paste0("D.",char)]
					data.wald[,paste0("DX.",char,".G.",i+1)] <- as.numeric(groupX==(i+1))*data.wald[,paste0("D.",char)]*data.wald[,X]
					new.var.name <- c(new.var.name,paste0("D.",char,".G.",i+1),paste0("DX.",char,".G.",i+1))
				}
				formula1.wald <- paste0(formula1.wald,"+",paste0(new.var.name,collapse="+"))
			}

			formula1.iv <- formula0.iv	
			for(i in 1:(nbins-1)){
				new.var.name <- c(paste0("G.",i+1),paste0("G.",i+1,".X"))
				for(sub.iv in IV){
					data.wald[,paste0(sub.iv,".G.",i+1)] <- as.numeric(groupX==(i+1))*data.wald[,sub.iv]
					data.wald[,paste0(sub.iv,".G.",i+1,".X")] <- as.numeric(groupX==(i+1))*data.wald[,sub.iv]*data.wald[,X]
					new.var.name <- c(new.var.name,paste0(sub.iv,".G.",i+1),paste0(sub.iv,".G.",i+1,".X"))
				}
				formula1.iv <- paste0(formula1.iv,"+",paste0(new.var.name,collapse="+"))
			}

			# add covariates
			if (is.null(Z)==FALSE){
				if(full.moderate==FALSE){
					formula0.wald <- paste0(formula0.wald, "+",paste(Z,collapse=" + "))
					formula1.wald <- paste0(formula1.wald, "+",paste(Z,collapse=" + "))
					formula0.iv <- paste0(formula0.iv, "+",paste(Z,collapse=" + "))
					formula1.iv <- paste0(formula1.iv, "+",paste(Z,collapse=" + "))
				}
				if(full.moderate==TRUE){ #z
					formula0.wald <- paste0(formula0.wald, "+",paste(Z,collapse=" + "))
					formula1.wald <- paste0(formula1.wald, "+",paste(Z,collapse=" + "))
					formula0.wald <- paste0(formula0.wald, "+",paste(Z.X,collapse=" + "))
					formula1.wald <- paste0(formula1.wald, "+",paste(Z.X,collapse=" + "))
					formula0.iv <- paste0(formula0.iv, "+",paste(Z,collapse=" + "))
					formula1.iv <- paste0(formula1.iv, "+",paste(Z,collapse=" + "))
					formula0.iv <- paste0(formula0.iv, "+",paste(Z.X,collapse=" + "))
					formula1.iv <- paste0(formula1.iv, "+",paste(Z.X,collapse=" + "))
					zero.var.name <- c()
					for(a in Z){
						for(i in 1:(nbins-1)){
							data.wald[,paste0("Z.",a,".G.",i+1)] <- as.numeric(groupX==(i+1))*data.wald[,a]
							data.wald[,paste0("ZX.",a,".G.",i+1)] <- as.numeric(groupX==(i+1))*data.wald[,a]*data.wald[,X]
							zero.var.name <- c(zero.var.name,paste0("Z.",a,".G.",i+1),paste0("ZX.",a,".G.",i+1))
						}
					}
					formula0.wald <- paste0(formula0.wald,"+",paste0(zero.var.name,collapse="+"))
					formula1.wald <- paste0(formula1.wald,"+",paste0(zero.var.name,collapse="+"))
					formula0.iv <- paste0(formula0.iv,"+",paste0(zero.var.name,collapse="+"))
					formula1.iv <- paste0(formula1.iv,"+",paste0(zero.var.name,collapse="+"))
				}
			}

			formula0.wald <- paste0(formula0.wald,"|",formula0.iv)
			formula1.wald <- paste0(formula1.wald,"|",formula1.iv)
			formula0.wald <- as.formula(formula0.wald)
			formula1.wald <- as.formula(formula1.wald)

			suppressWarnings(
				fit0 <- ivreg(formula0.wald,data=data.wald,weights=WEIGHTS)
			)
			suppressWarnings(
				fit1 <- ivreg(formula1.wald,data=data.wald,weights=WEIGHTS)
			)

			#vcov
			if(vcov.type=="homoscedastic"){
				model.wald.vcov <- vcov(fit1)
			}
			if(vcov.type=="robust"){
				model.wald.vcov <- vcovHC(fit1,type="HC1")
			} 
			if(vcov.type=="cluster"){
				model.wald.vcov <- vcovCluster(fit1,cluster = data.wald[,cl])
			}
			if(vcov.type=="pcse"){
				model.wald.vcov <- pcse(fit1,pairwise=pairwise,groupN=data.wald[,'cl'],groupT=data.wald[,'time'])$vcov
				rownames(model.wald.vcov)[1] <- "(Intercept)"
				colnames(model.wald.vcov)[1] <- "(Intercept)"
			}

			model.wald.coef <- coef(fit1)
			model.wald.vcov[which(is.na(model.wald.vcov))] <- 0
			model.wald.vcov.all <- matrix(0,nrow=length(model.wald.coef),ncol=length(model.wald.coef))

			colnames(model.wald.vcov.all) <- names(model.wald.coef)
			rownames(model.wald.vcov.all) <- names(model.wald.coef)
			for(a1 in rownames(model.wald.vcov)){
				for(a2 in colnames(model.wald.vcov)){
					model.wald.vcov.all[a1,a2] <- model.wald.vcov[a1,a2]
				}
			}
			
			#lr.test
			lr.result <- tryCatch(
				round(lmtest::lrtest(fit1,fit0)[[5]][2],3),
				error = function(e){return(NULL)}
			)
			
			#wald.test
			wald.result <- tryCatch(
				round(lmtest::waldtest(fit1,fit0,test="Chisq", vcov=model.wald.vcov.all)[[4]][2],3),
				error = function(e){return(NULL)}
			)

			#sub-test
			if(length(other.treat)>1){ #more than one treated group
				sub.test <- list()
				for(target.char in other.treat){
					formula0.sub.wald <- paste0(Y,"~",X)
					for(char in other.treat){
						formula0.sub.wald  <- paste0(formula0.sub.wald ,"+",paste0("D",".",char),"+",paste0("DX",".",char))
					}
					formula0.sub.iv <- formula1.iv

					for(i in 1:(nbins-1)){
						new.var.name <- c(paste0("G.",i+1),paste0("G.",i+1,".X"))
						for(char in other.treat){
							if(char != target.char){
								new.var.name <- c(new.var.name,paste0("D.",char,".G.",i+1),paste0("DX.",char,".G.",i+1))
							}
						}
						formula0.sub.wald <- paste0(formula0.sub.wald,"+",paste0(new.var.name,collapse="+"))
					}

					# add covariates
					if (is.null(Z)==FALSE){
						if(full.moderate==FALSE){
							formula0.sub.wald <- paste0(formula0.sub.wald, "+",paste(Z,collapse=" + "))
						}
						if(full.moderate==TRUE){ #z
							formula0.sub.wald <- paste0(formula0.sub.wald, "+",paste(Z,collapse=" + "))
							formula0.sub.wald <- paste0(formula0.sub.wald, "+",paste(Z.X,collapse=" + "))
							zero.var.name <- c()
							for(a in Z){
								for(i in 1:(nbins-1)){
									zero.var.name <- c(zero.var.name,paste0("Z.",a,".G.",i+1),paste0("ZX.",a,".G.",i+1))
								}
							}
							formula0.sub.wald <- paste0(formula0.sub.wald,"+",paste0(zero.var.name,collapse="+"))
						}
					}

					formula0.sub.wald <- paste0(formula0.sub.wald,"|",formula0.sub.iv)
					formula0.sub.wald <- as.formula(formula0.sub.wald)

					suppressWarnings(
						fit0.sub <- ivreg(formula0.sub.wald,data=data.wald,weights=WEIGHTS)
					)

					#lr.test
					lr.result.sub <- tryCatch(
						round(lmtest::lrtest(fit1,fit0.sub)[[5]][2],3),
						error = function(e){return(NULL)}
					)
			
					#wald.test
					wald.result.sub <- tryCatch(
						round(lmtest::waldtest(fit1,fit0.sub,test="Chisq", vcov=model.wald.vcov.all)[[4]][2],3),
						error = function(e){return(NULL)}
					)
					sub.test[[other.treat.origin[target.char]]] <- list(p.wald = sprintf(wald.result.sub,fmt = '%#.3f'),
																		p.lr = sprintf(lr.result.sub,fmt = '%#.3f'),
																		formula.restrict = formula0.sub.wald,
																		formula.full = formula1.wald)
				}
			}
		}

		requireNamespace("Lmoments")
		Lkurtosis <- Lmoments(data[,X],returnobject=TRUE)$ratios[4] 
		tests <- list(
					treat.type = treat.type,
					X.Lkurtosis = sprintf(Lkurtosis,fmt = '%#.3f'),
					p.wald = sprintf(wald.result,fmt = '%#.3f'),
					p.lr = sprintf(lr.result,fmt = '%#.3f'),
					formula.restrict = formula0.wald,
					formula.full = formula1.wald,
					sub.test = sub.test
		)
	}

	if(wald==TRUE & use_fe==TRUE & is.null(IV)){
		data.wald <- data
		sub.test <- NULL
		if(treat.type=='continuous'){
			formula1.wald <- paste0(Y,"~",X,"+",D,"+DX")
			var.name <- c(X,D,"DX")
			constrain.terms <- c()
			for(i in 1:(nbins-1)){
				data.wald[,paste0("G.",i+1)] <- as.numeric(groupX==(i+1))
				data.wald[,paste0("G.",i+1,".X")] <- as.numeric(groupX==(i+1))*data.wald[,X]
				data.wald[,paste0("DG.",i+1)] <- as.numeric(groupX==(i+1))*data.wald[,D]
				data.wald[,paste0("DG.",i+1,".X")] <- as.numeric(groupX==(i+1))*data.wald[,D]*data.wald[,X]
				new.var.name <- c(paste0("G.",i+1),paste0("G.",i+1,".X"),paste0("DG.",i+1),paste0("DG.",i+1,".X"))
				var.name <- c(var.name,new.var.name)
				formula1.wald <- paste0(formula1.wald,"+",paste0(new.var.name,collapse="+"))
				constrain.terms <- c(constrain.terms, new.var.name)
			}

			constraints <- as.formula(paste0("~",paste0(constrain.terms, collapse = "|")))

			if (is.null(Z)==FALSE){
				if(full.moderate==FALSE){
					formula1.wald <- paste0(formula1.wald, "+",paste(Z,collapse=" + "))
				}
				if(full.moderate==TRUE){ #z
					formula1.wald <- paste0(formula1.wald, "+",paste(Z,collapse=" + "))
					formula1.wald <- paste0(formula1.wald, "+",paste(Z.X,collapse=" + "))
					zero.var.name <- c()
					for(a in Z){
						for(i in 1:(nbins-1)){
							data.wald[,paste0("Z.",a,".G.",i+1)] <- as.numeric(groupX==(i+1))*data.wald[,a]
							data.wald[,paste0("ZX.",a,".G.",i+1)] <- as.numeric(groupX==(i+1))*data.wald[,a]*data.wald[,X]
							zero.var.name <- c(zero.var.name,paste0("Z.",a,".G.",i+1),paste0("ZX.",a,".G.",i+1))
						}
					}
					formula1.wald <- paste0(formula1.wald,"+",paste0(zero.var.name,collapse="+"))
				}
			}

      		formula1.wald <- paste0(formula1.wald, "|",paste0(FE, collapse = "+"))
      		if (vartype=="cluster") {
        		formula1.wald <- paste0(formula1.wald, "| 0 |",paste0(cl,collapse = "+"))
      		}

			formula1.wald <- as.formula(formula1.wald)
			w.wald <- data.wald[,'WEIGHTS']
			suppressWarnings(
				mod.un <- felm(formula1.wald,data=data.wald,weights=w.wald)
			)

			if(vartype=="homoscedastic") {
        		p.wald <- lfe::waldtest(mod.un, constraints, type = "default")[1]
     		}else if (vartype=="robust") {
        		p.wald <- lfe::waldtest(mod.un, constraints, type = "robust")[1]
      		}else {
        		p.wald <- lfe::waldtest(mod.un, constraints)[1] # clustered
      		}
      		names(p.wald) <- NULL            
     		p.wald <- round(p.wald,4)
		}

		if(treat.type=='discrete'){
			formula1.wald <- formula.origin
			constrain.terms <- c()
			for(i in 1:(nbins-1)){
				data.wald[,paste0("G.",i+1)] <- as.numeric(groupX==(i+1))
				data.wald[,paste0("G.",i+1,".X")] <- as.numeric(groupX==(i+1))*data.wald[,X]
				new.var.name <- c(paste0("G.",i+1),paste0("G.",i+1,".X"))
				for(char in other.treat){
					data.wald[,paste0("D.",char,".G.",i+1)] <- as.numeric(groupX==(i+1))*data.wald[,paste0("D.",char)]
					data.wald[,paste0("DX.",char,".G.",i+1)] <- as.numeric(groupX==(i+1))*data.wald[,paste0("D.",char)]*data.wald[,X]
					new.var.name <- c(new.var.name,paste0("D.",char,".G.",i+1),paste0("DX.",char,".G.",i+1))
				}

				constrain.terms <- c(constrain.terms,new.var.name)

				if(is.null(Z)==FALSE & full.moderate==TRUE){ #z
					zero.var.name <- c()
					for(a in Z){
						data.wald[,paste0("Z.",a,".G.",i+1)] <- as.numeric(groupX==(i+1))*data.wald[,a]
						data.wald[,paste0("ZX.",a,".G.",i+1)] <- as.numeric(groupX==(i+1))*data.wald[,a]*data.wald[,X]
						new.var.name <- c(new.var.name,paste0("Z.",a,".G.",i+1),paste0("ZX.",a,".G.",i+1))
						zero.var.name <- c(zero.var.name,paste0("Z.",a,".G.",i+1),paste0("ZX.",a,".G.",i+1))
					}
				}
				formula1.wald <- paste0(formula1.wald,"+",paste0(new.var.name,collapse="+"))
			}
			
			formula1.wald <- paste0(formula1.wald, "|",paste0(FE, collapse = "+"))
      		if (vartype=="cluster") {
        		formula1.wald <- paste0(formula1.wald, "| 0 |",paste0(cl,collapse = "+"))
      		}

			formula1.wald <- as.formula(formula1.wald)
			constraints <- as.formula(paste0("~",paste0(constrain.terms, collapse = "|")))

			w.wald <- data.wald[,'WEIGHTS']
			suppressWarnings(
				mod.un <- felm(formula1.wald,data=data.wald,weights=w.wald)
			)

			if(vartype=="homoscedastic") {
        		p.wald <- lfe::waldtest(mod.un, constraints, type = "default")[1]
     		}else if (vartype=="robust") {
        		p.wald <- lfe::waldtest(mod.un, constraints, type = "robust")[1]
      		}else {
        		p.wald <- lfe::waldtest(mod.un, constraints)[1] # clustered
      		}
      		names(p.wald) <- NULL            
     		p.wald <- round(p.wald,4)
			
			if(length(other.treat)>1){ #more than one treated group
				sub.test <- list()
				for(target.char in other.treat){
					constrain.terms.sub <- c()
					new.var.name <- c()
					for(i in 1:(nbins-1)){		
						constrain.terms.sub <- c(constrain.terms.sub,paste0("D.",target.char,".G.",i+1),paste0("DX.",target.char,".G.",i+1))
					}
					
					constraints.sub <- as.formula(paste0("~",paste0(constrain.terms.sub, collapse = "|")))
					
					if(vartype=="homoscedastic") {
        				p.wald.sub <- lfe::waldtest(mod.un, constraints.sub, type = "default")[1]
     				}else if (vartype=="robust") {
        				p.wald.sub <- lfe::waldtest(mod.un, constraints.sub, type = "robust")[1]
      				}else {
        				p.wald.sub <- lfe::waldtest(mod.un, constraints.sub)[1] # clustered
      				}
      				names(p.wald.sub) <- NULL            
     				p.wald.sub <- round(p.wald.sub,4)
					sub.test[[other.treat.origin[target.char]]] <- list(p.wald = sprintf(p.wald.sub,fmt = '%#.3f'))
				}
			}	
		}

		requireNamespace("Lmoments")
		Lkurtosis <- Lmoments(data[,X],returnobject=TRUE)$ratios[4] 
		tests <- list(
					treat.type = treat.type,
					X.Lkurtosis = sprintf(Lkurtosis,fmt = '%#.3f'),
					p.wald = sprintf(p.wald,fmt = '%#.3f'),
					formula.restrict = formula0.wald,
					formula.full = formula1.wald,
					sub.test = sub.test
				)
	}

	if(wald==TRUE & use_fe==TRUE & !is.null(IV)){
		#use the wald test in felm
		data.wald <- data
		sub.test <- NULL
		tests <- NULL
	}


	if(treat.type=='discrete'){
		final.output <- list(treat.info=treat.info,
							 est.lin=TE.output.all.list,
							 pred.bin=pred.output.all.list,
							 link.bin=link.output.all.list,
							 est.bin=TE.binning.output.all.list,
							 nbins=nbins,
							 Xlabel = Xlabel,
							 Dlabel = Dlabel,
							 Ylabel = Ylabel,
							 de = de,
							 de.tr = treat_den, # density
							 hist.out = hist.out,
							 count.tr = treat.hist,
							 tests = tests,
							 estimator = "binning",
							 model.linear = model,
							 model.binning = model.binning,
							 use.fe = use_fe
							)
	}
  
	if(treat.type=='continuous'){
		final.output <- list(treat.info=treat.info,
							est.lin=ME.output.all.list,
							pred.bin=pred.output.all.list,
							link.bin=link.output.all.list,
							est.bin=ME.binning.output.all.list,
							nbins=nbins,
							Xlabel = Xlabel,
							Dlabel = Dlabel,
							Ylabel = Ylabel,
							de = de, # density
							de.tr = de.tr,
							hist.out = hist.out,
							count.tr = NULL,
							tests = tests,
							estimator = "binning",
							model.linear = model,
							model.binning = model.binning,
							use.fe = use_fe
							)
	}
	
	#Plot
  if(figure==TRUE){
	class(final.output) <- "interflex"
	figure.output <- plot.interflex(	x=final.output,
									order = order,
									subtitles = subtitles,
									show.subtitles = show.subtitles,
									CI = CI,
									diff.values = NULL,
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
									show.grid = show.grid,     
									cex.main = cex.main,
									cex.sub = cex.sub,
									cex.lab = cex.lab,
									cex.axis = cex.axis,
									bin.labs = bin.labs, # bin labels    
									interval = interval, # interval in replicated papers
									file = file,
									ncols = ncols,
									#pool plot
									pool = pool,
									legend.title = legend.title,
									color = color,
									show.all = show.all,
									scale = scale,
  									height = height,
  									width = width
								)

	final.output <- c(final.output,list(figure=figure.output))
  
  }
  
	
	
	
	
	
	class(final.output) <- "interflex"
	return(final.output)

	}
	
	



