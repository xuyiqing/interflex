# PLR Estimation and Parallel Bootstrap for Continuous D & Discrete X
# Requires: splines, glmnet, doParallel, foreach

estimateGATE_PLR <- function(
  data,
  Y,                # outcome variable
  D,                # continuous treatment
  X,                # discrete moderator
  Z = NULL,         # covariates
  FE = NULL,        # fixed effects
  XZ_design = NULL, # optional prebuilt design

  outcome_model_type   = "lasso",
  treatment_model_type = "lasso",

  basis_type           = c("polynomial","bspline","none"),
  include_interactions = TRUE,
  poly_degree          = 2,
  spline_df            = 4,
  spline_degree        = 2,

  lambda_cv            = NULL,   # list(outcome=...,treatment=...)
  lambda_seq           = NULL,

  verbose              = TRUE
) {
  # 0. Setup
  if (verbose) cat("PLR Step 0: Checking inputs...\n")
  basis_type <- match.arg(basis_type)
  stopifnot(is.data.frame(data))
  for(nm in c(Y,D,X)) if(!(nm%in%names(data))) stop(nm," missing.")
  Yv <- data[[Y]]; Dv <- data[[D]]; Xv <- data[[X]]

  # 1. Build/validate design
  if (verbose) cat("PLR Step 1: Building design matrix...\n")
  if(!is.null(XZ_design)){
    if(!is.matrix(XZ_design)||nrow(XZ_design)!=nrow(data)) stop("Bad XZ_design.")
  } else {
    Zm <- if(!is.null(Z)) data[,Z,drop=FALSE] else NULL
    # bs helper
    bs_helper <- function(vec,name){
      m <- splines::bs(vec,df=spline_df,degree=spline_degree)
      colnames(m) <- paste0(name,"_bs",seq_len(ncol(m))); m
    }
    # X dummies
    fac <- factor(Xv)
    Xd  <- model.matrix(~fac)[,-1,drop=FALSE]
    colnames(Xd) <- paste0(X,"_",levels(fac)[-1])
    # Z expansion
    Ze <- NULL
    if(!is.null(Zm)) for(j in seq_len(ncol(Zm))){
      v <- Zm[[j]]; nm <- colnames(Zm)[j]
      m <- switch(basis_type,
        none = matrix(v,ncol=1,dimnames=list(NULL,nm)),
        polynomial = { tmp<-stats::poly(v,degree=poly_degree,raw=TRUE); colnames(tmp)<-paste0(nm,"_poly",seq_len(ncol(tmp))); tmp },
        bspline    = bs_helper(v,nm)
      )
      Ze <- if(is.null(Ze)) m else cbind(Ze,m)
    }
    # FE dummies
    FEd <- NULL
    if(!is.null(FE)) for(fe in FE){ f<-factor(data[[fe]]); m<-model.matrix(~f)[,-1,drop=FALSE]; colnames(m)<-paste0(fe,'_',levels(f)[-1]); FEd<-if(is.null(FEd))m else cbind(FEd,m) }
    # interactions
    IntXZ <- NULL
    if(include_interactions && !is.null(Ze)){
      for(xn in colnames(Xd)) for(zn in colnames(Ze)){
        tmp <- Xd[,xn]*Ze[,zn]; IntXZ<-cbind(IntXZ,tmp)
      }
      colnames(IntXZ) <- as.vector(outer(colnames(Xd),colnames(Ze),paste,sep="."))
    }
    XZ_design <- cbind(Xd,Ze,FEd,IntXZ)
    if(verbose) cat("  → design with",ncol(XZ_design),"columns.\n")
  }

  # 2. Penalty factors
  pf_out <- pf_tr <- rep(1,ncol(XZ_design))
  if(!is.null(FE)){
    pat <- paste0('^(',paste(FE,collapse='|'),')_')
    idx <- grep(pat,colnames(XZ_design)); pf_out[idx]<-0; pf_tr[idx]<-0
  }

  # helper fit
  fit_glmnet <- function(y,Xm,type,lam,pf){
    if(type=='linear') return(list(fit=lm(y~.,data=as.data.frame(Xm)),type='lm'))
    a <- if(type=='lasso')1 else 0
    cv <- if(is.null(lam)) glmnet::cv.glmnet(Xm,y,alpha=a,lambda=lambda_seq,penalty.factor=pf)
          else          glmnet::glmnet(Xm,y,alpha=a,lambda=lam,penalty.factor=pf)
    list(fit=cv,type='glmnet',lambda=if(is.null(lam))cv$lambda.min else lam)
  }

  # 3. Fit
  if(verbose) cat("PLR Step 3: Fitting outcome & treatment...\n")
  out_fit <- fit_glmnet(Yv,XZ_design,outcome_model_type,   lambda_cv$outcome,   pf_out)
  trt_fit<- fit_glmnet(Dv,XZ_design,treatment_model_type, lambda_cv$treatment, pf_tr)

  # 4. Post-selection
  selected <- NULL
  if(outcome_model_type!='linear' || treatment_model_type!='linear'){
    active <- function(obj){ co<-coef(obj$fit,s=obj$lambda); which(as.numeric(co)!=0)-1 }
    ao<-active(out_fit); at<-active(trt_fit)
    uni<-unique(c(ao,at)); selected<-list(outcome=ao,treatment=at)
    if(length(uni)>0){
      subXZ <- XZ_design[,uni,drop=FALSE]
      out_fit<-list(fit=lm(Yv~.,data=as.data.frame(subXZ)),type='lm')
      trt_fit<-list(fit=lm(Dv~.,data=as.data.frame(subXZ)),type='lm')
    }
  }

  # 5. Residuals
  pred <- function(obj,Xm){
    if(obj$type=='lm')    as.numeric(predict(obj$fit,newdata=as.data.frame(Xm)))
    else                  as.numeric(predict(obj$fit,newx=Xm,s='lambda.min'))
  }
  yhat<-pred(out_fit,XZ_design); dhat<-pred(trt_fit,XZ_design)
  ytilde<-Yv-yhat; dtilde<-Dv-dhat

  # 6. Group regression
  if(verbose) cat("PLR Step 6: Computing GATE by group...\n")
  fac <- factor(Xv); lvls<-levels(fac)
  gate <- sapply(lvls,function(lev){
    idx <- which(fac==lev)
    coef(lm(ytilde[idx]~dtilde[idx]))[2]
  })
  gate_df <- data.frame(X=lvls,GATE=gate,stringsAsFactors=FALSE)

  if(verbose) cat("PLR Step 7: Done.\n")
  list(XZ_design=XZ_design, gate_df=gate_df, selected=selected)
}

bootstrapGATE_PLR <- function(
  data, Y, D, X,
  Z = NULL, FE = NULL,
  B = 200, alpha = 0.05,
  outcome_model_type   = "lasso",
  treatment_model_type = "lasso",
  basis_type           = c("polynomial","bspline","none"),
  include_interactions = TRUE,
  poly_degree          = 2,
  spline_df            = 4,
  spline_degree        = 2,
  lambda_cv            = NULL,
  lambda_seq           = NULL,
  verbose              = TRUE
) {
  if(verbose) message("BootPLR Step 1: Fitting full-sample PLR...")
  full <- estimateGATE_PLR(
    data, Y, D, X, Z, FE, NULL,
    outcome_model_type, treatment_model_type,
    basis_type, include_interactions,
    poly_degree, spline_df, spline_degree,
    lambda_cv, lambda_seq, verbose
  )
  lvls <- full$gate_df$X
  g_full <- full$gate_df$GATE
  n      <- nrow(data); K <- length(lvls)
  XZ0   <- full$XZ_design

  if(verbose) message("BootPLR Step 2: Setting up parallel backend...")
  if(!requireNamespace("doParallel",quietly=TRUE)) stop("doParallel needed")
  cl <- parallel::makeCluster(parallel::detectCores())
  doParallel::registerDoParallel(cl)
  `%dopar%` <- foreach::`%dopar%`

  if(verbose) message("BootPLR Step 3: Running",B,"bootstraps...")
  res_mat <- foreach::foreach(
    b = seq_len(B), .combine="rbind",
    .export="estimateGATE_PLR",
    .packages=c("splines","glmnet")
  ) %dopar% {
    set.seed(1000 + b)
    idx <- sample(n, n, replace=TRUE)
    db  <- data[idx, , drop=FALSE]
    XZb <- XZ0[idx, , drop=FALSE]
    fb  <- estimateGATE_PLR(
      db, Y, D, X, NULL, FE, XZb,
      outcome_model_type, treatment_model_type,
      basis_type, include_interactions,
      poly_degree, spline_df, spline_degree,
      lambda_cv, lambda_seq, FALSE
    )
    sapply(lvls, function(lv){ i<-which(fb$gate_df$X==lv); if(length(i)) fb$gate_df$GATE[i] else NA_real_ })
  }
  parallel::stopCluster(cl)
  if(verbose) message("BootPLR Step 4: Computing stats...")
  se  <- apply(res_mat,2,sd,na.rm=TRUE)
  cil <- apply(res_mat,2,quantile,probs=alpha/2,na.rm=TRUE)
  cih <- apply(res_mat,2,quantile,probs=1-alpha/2,na.rm=TRUE)
  uni <- calculate_uniform_quantiles(t(res_mat),alpha)
  results <- data.frame(
    X=lvls, GATE=g_full,
    SE=se, CI.lower=cil, CI.upper=cih,
    CI.lower.uni=uni$Q_j[,1], CI.upper.uni=uni$Q_j[,2],
    stringsAsFactors=FALSE
  )
  if(verbose) message("BootPLR Step 5: Done.")
  list(results=results, boot = res_mat, coverage=uni$coverage, full=full)
}
