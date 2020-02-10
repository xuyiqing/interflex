#' @rdname ttest.interflex
#' @export
ttest <- function(
	out,
	diff.values,
	percentile=FALSE,
	k=15
){
	if(!class(out) %in% c("interflex")){
		stop("Not an \"interflex\" object.")
	}
	
	requireNamespace("mgcv")
	
	if(is.numeric(k)==FALSE){
		stop("\"k\" is not numeric.")
	}
	
	if(is.numeric(diff.values)==FALSE){
		stop("\"diff.values\" is not numeric.")
	}
	if(length(diff.values)!=3 & length(diff.values)!=2){
		stop("\"diff.values\" must be of length 2 or 3.")
	}
	
	if(is.logical(percentile) == FALSE & is.numeric(percentile)==FALSE) {
		stop("\"percentile\" is not a logical flag.")
	}

	if (k%%1 != 0) {
		stop("\"k\" is not a positive integer.")
	} else {
		k <- k[1]
	}
	if (k < 1) {
		stop("\"k\" should be a positive integer larger than 1.")
	}
	
	if(percentile==TRUE){
		for(a in diff.values){
			if(a<0|a>1){
				stop("Elements in \"diff.values\" should be between 0 and 1 when percentile==TRUE.")
			}
		}
	}
	
	treat.type <- out$treat.type
	type <- out$type
	base <- out$base
	treat.type <- out$treat.type
	all_treat <- out$treatlevels
	
	if(type=='kernel'){
		if(out$CI==FALSE){
			stop("Covariance matrix was not estimated.")
		}
	}
	
	if(type=='linear'){
		type <- 'binning'
	}
	
	if(treat.type=='discrete'){
		other_treat <- sort(all_treat[which(all_treat!=base)])
	}
	
	if(treat.type=='discrete' & type=='binning'){
		tempxx <- out$est.lin[[other_treat[1]]][,'X.lvls']
	}
	if(treat.type=='discrete' & type=='kernel'){
		tempxx <- out$est[[other_treat[1]]][,'X']
	}
	if(treat.type=='continuous' & type=='binning'){
		tempxx <- out$est.lin[,'X.lvls']
	}
	if(treat.type=='continuous' & type=='kernel'){
		tempxx <- out$est[,'X']
	}
	min.XX <- min(tempxx)
	max.XX <- max(tempxx)
	
	if(percentile==TRUE){
		diff.pc <- diff.values
		diff.values <- quantile(tempxx,probs=diff.values)
	}
	
	for(a in diff.values){
		if(a<min.XX|a>max.XX){
			stop("Elements in \"diff.values\" should be greater than the minimum and less than the maximum of the moderator.")
		}
	}
	
	if(type=='binning'|type=='linear'){
		est <- out$est.lin
	}
	if(type=='kernel'){
		est <- out$est
	}
	vcov.matrix <- out$vcov.matrix
	
	if(length(diff.values)==3){
		if(percentile==FALSE){
				diff.name <- c(paste0(diff.values[2]," vs ",diff.values[1]),
								   paste0(diff.values[3]," vs ",diff.values[2]),
								   paste0(diff.values[3]," vs ",diff.values[1]))
		}
		if(percentile==TRUE){
				diff.name <- c(paste0(round(100*diff.pc[2],3),'%',' vs ',round(100*diff.pc[1],3),'%'),
								   paste0(round(100*diff.pc[3],3),'%',' vs ',round(100*diff.pc[2],3),'%'),
								   paste0(round(100*diff.pc[3],3),'%',' vs ',round(100*diff.pc[1],3),'%'))
		}
	}
	if(length(diff.values)==2){
		if(percentile==FALSE){
				diff.name <- c(paste0(diff.values[2]," vs ",diff.values[1]))
			}
		if(percentile==TRUE){
				diff.name <- c(paste0(round(100*diff.pc[2],3),'%',' vs ',round(100*diff.pc[1],3),'%'))
			}
	}
	
	get_stat <- function(point1,point2){
			var1 <- predict(model_use,newdata = cbind.data.frame(x1=point1,x2=point1))
			var2 <- predict(model_use,newdata = cbind.data.frame(x1=point2,x2=point2))
			cov12 <- predict(model_use,newdata = cbind.data.frame(x1=point1,x2=point2))
			marg1 <- predict(model_use2,newdata = cbind.data.frame(X=point1))
			marg2 <- predict(model_use2,newdata = cbind.data.frame(X=point2))
			diff <- marg2-marg1
			var.diff <- var2+var1-2*cov12
			se.diff <- sqrt(var.diff)
			z.value <- diff/se.diff
			pvalue2sided <- 2*pnorm(-abs(z.value))
			lbound <- diff-1.96*se.diff
			ubound <- diff+1.96*se.diff
			diff.table <- round(c(diff,se.diff,z.value,pvalue2sided,lbound,ubound),3)
			return(diff.table)
	}
	
	if(treat.type=='continuous'){
			x1 <- rep(tempxx,length(tempxx))
			x2 <- rep(tempxx,each=length(tempxx))
			cov_base <- c(vcov.matrix)
			
			data_touse <- cbind.data.frame(x1=x1,x2=x2,y=cov_base)
			data_touse2 <- cbind.data.frame(X=est[,1],ME=est[,2])
			
			model_use <- gam(y~te(x1,x2,k=k),data = data_touse)
			model_use2 <- gam(ME~s(X,k = k),data=data_touse2)
			
			diff.table <- matrix(0,nrow=0,ncol=6)
			colnames(diff.table) <- c("Difference","se","Z-Score","P-value","CI-lower(95%)","CI-upper(95%)")
			
			if(length(diff.values)==3){
				diff.table <- rbind(diff.table,get_stat(diff.values[1],diff.values[2]))
				diff.table <- rbind(diff.table,get_stat(diff.values[2],diff.values[3]))
				diff.table <- rbind(diff.table,get_stat(diff.values[1],diff.values[3]))
				rownames(diff.table) <- diff.name
			}
			
			if(length(diff.values)==2){
				diff.table <- rbind(diff.table,get_stat(diff.values[1],diff.values[2]))
				rownames(diff.table) <- diff.name
			}
	}
	
	if(treat.type=='discrete'){
		x1 <- rep(tempxx,length(tempxx))
		x2 <- rep(tempxx,each=length(tempxx))
		diff.table.list <- list()
		
		for(char in other_treat){
			cov_base <- c(vcov.matrix[[char]])
			
			data_touse <- cbind.data.frame(x1=x1,x2=x2,y=cov_base)
			data_touse2 <- cbind.data.frame(X=est[[char]][,1],ME=est[[char]][,2])
			
			model_use <- gam(y~te(x1,x2,k=k),data = data_touse)
			model_use2 <- gam(ME~s(X,k = k),data=data_touse2)
			
			diff.table <- matrix(0,nrow=0,ncol=6)
			colnames(diff.table) <- c("Difference","se","Z-Score","P-value","CI-lower(95%)","CI-upper(95%)")
			
			if(length(diff.values)==3){
				diff.table <- rbind(diff.table,get_stat(diff.values[1],diff.values[2]))
				diff.table <- rbind(diff.table,get_stat(diff.values[2],diff.values[3]))
				diff.table <- rbind(diff.table,get_stat(diff.values[1],diff.values[3]))
				rownames(diff.table) <- diff.name
			}
			
			if(length(diff.values)==2){
				diff.table <- rbind(diff.table,get_stat(diff.values[1],diff.values[2]))
				rownames(diff.table) <- diff.name
			}
			
			diff.table.list[[char]] <- diff.table		
		}
		
		diff.table <- diff.table.list
	}
	
	return(diff.table)
}