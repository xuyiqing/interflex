vcovCluster <- function(
    model,
    cluster
)
{
    requireNamespace("sandwich")
    
    if(nrow(model.matrix(model))!=length(cluster)){
        stop("check your data: cluster variable has different N than model")
    }
    M <- length(unique(cluster))
    N <- length(cluster)           
    K <- model$rank   
	dfc <- (M/(M - 1)) * ((N - 1)/(N - K))

    est.matrix <- estfun(model)
	est.colname <- colnames(est.matrix)
	est.ori.length <- length(est.colname)
	coef.colname <- names(model$coefficients)
	diff.name <- coef.colname[-which(coef.colname %in% est.colname)]
	if(length(diff.name)>0){
		for(char in diff.name){
		est.matrix <- cbind(est.matrix,rep(0,N))
		est.colname <- c(est.colname,char)
	}
		colnames(est.matrix) <- est.colname
	}

	uj  <- apply(est.matrix, 2, function(x) tapply(x, cluster, sum))
	bread.empty <- matrix(0,nrow = length(coef.colname),ncol = length(coef.colname))
	colnames(bread.empty) <- est.colname
	rownames(bread.empty) <- est.colname
	
	
	bread.empty[1:est.ori.length,1:est.ori.length] <- bread(model)
	rcse.cov <- dfc * sandwich(model,bread. = bread.empty, meat. = crossprod(uj)/N)

    return(rcse.cov)
}