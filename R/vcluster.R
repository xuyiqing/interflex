vcovCluster <- function(
    model,
    cluster
)
{
    requireNamespace("sandwich")
    ## require(lmtest)
    
    if(nrow(model.matrix(model))!=length(cluster)){
        stop("check your data: cluster variable has different N than model")
    }
    M <- length(unique(cluster))
    N <- length(cluster)           
    K <- model$rank   
    ## if(M<50){
    ##     warning("Fewer than 50 clusters, variances may be unreliable (could try block bootstrap instead).")
    ## }
    dfc <- (M/(M - 1)) * ((N - 1)/(N - K))
    uj  <- apply(estfun(model), 2, function(x) tapply(x, cluster, sum));
    rcse.cov <- dfc * sandwich(model, meat. = crossprod(uj)/N)
    ##colnames(rcse.cov)<-rownames(rcse.cov)<-names(model$coefficients)
    return(rcse.cov)
}