#' @title Estimate Group Treatment Effects (GTE) for a Discrete Moderator
#' @description This function estimates Group Treatment Effects (the CME for a discrete
#'   moderator) using outcome, IPW, or AIPW signals. It is designed for a binary
#'   treatment and supports basis expansion for covariates, fixed effects, and
#'   post-Lasso model refitting.
#'
#' @param data A data.frame containing the variables for the analysis.
#' @param Y Character string, the name of the outcome variable in `data`.
#' @param D Character string, the name of the binary (0/1) treatment variable in `data`.
#' @param X Character string, the name of the discrete moderator variable in `data`.
#' @param Z Character vector, optional names of additional non-fixed-effect covariates.
#' @param FE Character vector, optional names of fixed effect variables. These will be
#'   included as unpenalized dummies in outcome models but excluded from the PS model.
#' @param estimand Character, the target estimand. Must be one of 'ATE' (Average
#'   Treatment Effect) or 'ATT' (Average Treatment Effect on the Treated).
#' @param signal Character, the type of statistical signal to construct. Must be one
#'   of "outcome", "ipw", or "aipw".
#' @param basis_type Character, the type of basis expansion for continuous covariates in Z.
#'   Must be one of "polynomial", "bspline", or "none".
#' @param include_interactions Logical, whether to include interaction terms between X, Z,
#'   and within Z in the design matrix.
#' @param poly_degree Integer, the degree for polynomial basis expansion.
#' @param spline_df Integer, degrees of freedom for B-spline basis expansion.
#' @param spline_degree Integer, the degree for B-spline basis expansion.
#' @param XZ_design A pre-built design matrix for covariates. If NULL, it is created
#'   internally based on X, Z, FE, and basis expansion arguments.
#' @param outcome_model_type Character, the model type for outcome regression. Must be
#'   one of "linear", "ridge", or "lasso".
#' @param ps_model_type Character, the model type for the propensity score model. Must be
#'   one of "linear", "ridge", or "lasso".
#' @param lambda_cv A list with pre-specified lambda values (e.g., from a full-sample
#'   fit) to be used instead of cross-validation. Elements can be `outcome1`,
#'   `outcome0`, `treatment`. Used during bootstrapping.
#' @param lambda_seq A custom numeric vector of lambda values for the glmnet grid search.
#' @param verbose Logical, whether to print progress messages to the console.
#'
#' @return A list containing detailed results of the estimation, including the final
#'   GTE data.frame, the constructed design matrix, fitted model objects, nuisance
#'   predictions (mu0, mu1, pi), the individual-level signal, and the set of
#'   covariates selected by Lasso.
#'
estimateGTE <- function(
  data,
  Y,              # outcome name
  D,              # binary treatment name
  X,              # discrete moderator name
  Z = NULL,       # extra covariates
  FE = NULL,      # fixed-effect vars
  estimand   = c("ATE","ATT"),
  signal     = c("outcome","ipw","aipw"),
  basis_type = c("polynomial","bspline","none"),
  include_interactions = FALSE,
  poly_degree  = 2,
  spline_df    = 5,
  spline_degree = 3,
  XZ_design    = NULL,    # optional prebuilt matrix
  outcome_model_type = "lasso",  # "linear","ridge","lasso"
  ps_model_type      = "lasso",  # "linear","ridge","lasso"
  lambda_cv    = NULL,    # list(outcome1, outcome0, treatment)
  lambda_seq   = NULL,    
  verbose      = TRUE
) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for penalized fits. Please install it.")
  }
  # 1) args & basic checks 
  estimand   <- match.arg(estimand)
  signal     <- match.arg(signal)
  basis_type <- match.arg(basis_type)
  stopifnot(is.data.frame(data),
            Y %in% names(data),
            D %in% names(data),
            X %in% names(data))
  if (!is.null(lambda_cv)) {
    allowed <- c("outcome1","outcome0","treatment")
    extra   <- setdiff(names(lambda_cv), allowed)
    if (length(extra) > 0) {
      stop("lambda_cv contains invalid entry(ies): ", paste(extra, collapse = ", "))
    }
  }
  Yvec <- data[[Y]]
  Dvec <- data[[D]]
  Xvec <- data[[X]]
  if (!all(Dvec %in% c(0,1))) stop("D must be 0/1.")

  # 2) prepare Z_nonFE & FE_dummies 
  Z_nonFE <- if (!is.null(Z)) {
    stopifnot(all(Z %in% names(data)))
    data[, Z, drop=FALSE]
  } else NULL

  FE_dummies <- NULL
  if (!is.null(FE)) {
    stopifnot(all(FE %in% names(data)))
    FE_mat <- data[, FE, drop=FALSE]
    FE_dummies <- do.call(cbind, lapply(names(FE_mat), function(v) {
      fac <- factor(FE_mat[[v]])
      mm  <- model.matrix(~ fac)[,-1,drop=FALSE]
      colnames(mm) <- paste0(v, "_", levels(fac)[-1])
      mm
    }))
    if(verbose) cat("Fixed effects -> dummies\n")
  }

  # 3) build design matrix if NULL 
  if (is.null(XZ_design)) {
    build_design_matrix <- function(X_vec, Z_nonFE, FE_dummies) {
      # 3A) discrete X -> k-1 dummies
      fX    <- factor(X_vec)
      X_mat <- model.matrix(~ fX)[,-1,drop=FALSE]
      colnames(X_mat) <- paste0(deparse(substitute(X)), "_", levels(fX)[-1])

      # 3B) expand Z_nonFE
      Z_expanded <- NULL
      if (!is.null(Z_nonFE) && ncol(Z_nonFE)>0) {
        Z_expanded <- do.call(cbind, lapply(seq_len(ncol(Z_nonFE)), function(j) {
          vec   <- Z_nonFE[[j]]
          name  <- colnames(Z_nonFE)[j]
          if (basis_type=="none") {
            m <- matrix(vec,ncol=1); colnames(m)<-name; m
          } else if (basis_type=="polynomial") {
            p <- poly(vec,degree=poly_degree,raw=TRUE)
            colnames(p) <- paste0(name,"_poly",seq_len(ncol(p))); p
          } else {
            b <- splines::bs(vec,df=spline_df,degree=spline_degree)
            colnames(b) <- paste0(name,"_bs",seq_len(ncol(b))); b
          }
        }))
      }

      # 3C) stack Z_expanded + FE_dummies
      Z_block <- if (!is.null(FE_dummies)) {
        if (is.null(Z_expanded)) FE_dummies else cbind(Z_expanded, FE_dummies)
      } else Z_expanded

      # 3D) X:Z interactions
      int_XZ <- NULL
      if (include_interactions && !is.null(Z_expanded)) {
        int_XZ <- do.call(cbind, lapply(colnames(X_mat), function(xn) {
          sweep(Z_expanded,1,X_mat[,xn],"*")
        }))
        colnames(int_XZ) <- as.vector(outer(colnames(X_mat),
                                            colnames(Z_expanded),
                                            paste,sep="."))
      }

      int_ZZ <- NULL
      if (include_interactions && !is.null(Z_expanded) && ncol(Z_expanded) > 1) {
        combs <- combn(colnames(Z_expanded), 2, simplify = FALSE)
        zz_list <- lapply(combs, function(p) Z_expanded[,p[1]] * Z_expanded[,p[2]])
        int_ZZ <- do.call(cbind, zz_list)
        colnames(int_ZZ) <- vapply(combs, paste, collapse = ".", FUN.VALUE = "")
      }

      cbind(X_mat, Z_block, int_XZ, int_ZZ)
    }
    XZ_design <- build_design_matrix(Xvec, Z_nonFE, FE_dummies)
    if(verbose) cat("Design matrix built\n")
  }

  n      <- nrow(data)
  p      <- ncol(XZ_design)
  fe_cols <- if(!is.null(FE_dummies)) {
    grep(paste0("^(", paste(FE, collapse="|"), ")_"), colnames(XZ_design))
  } else integer(0)

  # 4) fit outcome & ps w/ post-selection 
  requireNamespace("glmnet")
  selected <- list(outcome1=NULL, outcome0=NULL, ps=NULL)
  lambda_used <- list(outcome1=NULL, outcome0=NULL, treatment=NULL)

  # helper for outcome fits
  fit_outcome <- function(y, Xm, type, lam, pf) {
    if (type == "linear") {
      df_lin <- data.frame(y = y, Xm)
      colnames(df_lin)[-1] <- colnames(Xm)
      lm0 <- lm(y ~ ., data = df_lin)
      return(list(
        refit  = lm0,
        active = seq_len(ncol(Xm)),
        lambda = NULL
      ))
    }

    alpha <- ifelse(type == "lasso", 1, 0)
    if (!is.null(lam)) {
      fit0 <- glmnet::glmnet(Xm, y, alpha = alpha,
                             lambda = lam, penalty.factor = pf)
      lam0 <- lam
    } else {
      cv0  <- glmnet::cv.glmnet(Xm, y, alpha = alpha,
                                lambda = lambda_seq,
                                penalty.factor = pf)
      fit0 <- cv0
      lam0 <- cv0$lambda.min
    }
    co  <- as.numeric(coef(fit0, s = lam0))[-1]
    act <- which(co != 0)

    if (length(act) > 0) {
      df_sel <- data.frame(y = y, Xm[, act, drop = FALSE])
      lm1    <- lm(y ~ ., data = df_sel)
    } else {
      lm1    <- lm(y ~ 1, data = data.frame(y = y))
    }

    list(
      refit  = lm1,
      active = act,
      lambda = lam0
    )
  }

  # helper for PS fits
  fit_ps <- function(d, Xm, type, lam) {
    if (type == "linear") {
      df_lin <- data.frame(d = d, Xm)
      colnames(df_lin)[-1] <- colnames(Xm)
      gl0 <- glm(d ~ ., data = df_lin, family = binomial)
      return(list(
        refit  = gl0,
        active = seq_len(ncol(Xm)),
        lambda = NULL
      ))
    }

    alpha <- ifelse(type == "lasso", 1, 0)
    if (!is.null(lam)) {
      fit0 <- glmnet::glmnet(Xm, d, alpha = alpha,
                             lambda = lam, family = "binomial")
      lam0 <- lam
    } else {
      cv0  <- glmnet::cv.glmnet(Xm, d, alpha = alpha,
                                lambda = lambda_seq,
                                family = "binomial")
      fit0 <- cv0
      lam0 <- cv0$lambda.min
    }

    co  <- as.numeric(coef(fit0, s = lam0))[-1]
    act <- which(co != 0)

    if (length(act) > 0) {
      df_sel <- data.frame(d = d, Xm[, act, drop = FALSE])
      gl1    <- glm(d ~ ., data = df_sel, family = binomial)
    } else {
      gl1    <- glm(d ~ 1, data = data.frame(d = d), family = binomial)
    }

    list(
      refit  = gl1,
      active = act,
      lambda = lam0
    )
  }

  # 4A) outcome
  mu1_hat <- mu0_hat <- NULL
  if (signal %in% c("outcome","aipw")) {
    if(verbose) cat("Fitting outcome models...\n")
    idx1 <- which(Dvec==1); idx0 <- which(Dvec==0)
    X1   <- XZ_design[idx1,,drop=FALSE]; y1 <- Yvec[idx1]
    X0   <- XZ_design[idx0,,drop=FALSE]; y0 <- Yvec[idx0]
    pf <- rep(1, p); pf[fe_cols] <- 0

    o1 <- fit_outcome(y1, X1, outcome_model_type,
                      if(!is.null(lambda_cv)) lambda_cv$outcome1 else NULL, pf)
    o0 <- fit_outcome(y0, X0, outcome_model_type,
                      if(!is.null(lambda_cv)) lambda_cv$outcome0 else NULL, pf)

    mu1_hat <- predict(o1$refit, newdata=as.data.frame(XZ_design))
    mu0_hat <- predict(o0$refit, newdata=as.data.frame(XZ_design))

    selected$outcome1 <- o1$active
    selected$outcome0 <- o0$active
    lambda_used$outcome1 <- o1$lambda
    lambda_used$outcome0 <- o0$lambda
  }

  # 4B) propensity
  p_hat <- NULL; ps_fit <- NULL
  if (signal %in% c("ipw","aipw")) {
    if(verbose) cat("Fitting propensity model...\n")
    Xps <- if(length(fe_cols)>0) XZ_design[,-fe_cols,drop=FALSE] else XZ_design
    p1 <- fit_ps(Dvec, Xps, ps_model_type,
                 if(!is.null(lambda_cv)) lambda_cv$treatment else NULL)
    p_hat <- predict(p1$refit, newdata=as.data.frame(Xps), type="response")
    p_hat <- pmin(pmax(p_hat,1e-2),1-1e-2)

    selected$ps <- p1$active
    lambda_used$treatment <- p1$lambda
    ps_fit <- p1$refit
  }

  # 5) for ATT, group-mean p(D=1|X) 
  p_hat_X <- NULL
  if (estimand=="ATT") {
    grp <- tapply(Dvec, Xvec, mean)
    p_hat_X <- grp[as.character(Xvec)]
  }

  # 6) build individual signal 
  if (estimand=="ATE") {
    est_signal <- switch(signal,
      outcome = mu1_hat - mu0_hat,
      ipw     = Dvec*Yvec/p_hat - (1-Dvec)*Yvec/(1-p_hat),
      aipw    = (mu1_hat - mu0_hat) +
                  Dvec*(Yvec-mu1_hat)/p_hat - (1-Dvec)*(Yvec-mu0_hat)/(1-p_hat)
    )
  } else {
    est_signal <- switch(signal,
      outcome = (Yvec - mu0_hat)*Dvec/p_hat_X,
      ipw     = Yvec*(Dvec - p_hat)/((1-p_hat)*p_hat_X),
      aipw    = ((Yvec - mu0_hat)*(Dvec - (1-Dvec)*p_hat/(1-p_hat))) / p_hat_X
    )
  }

  # 7) collapse by X 
  gte_vals <- tapply(est_signal, Xvec, mean)
  gte_df   <- data.frame(
    X   = names(gte_vals),
    GTE = as.numeric(gte_vals),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  # 8) return everything 
  ret <- list(
    XZ_design   = XZ_design,
    mu1_hat     = mu1_hat,
    mu0_hat     = mu0_hat,
    p_hat       = p_hat,
    p_hat_X   = p_hat_X,
    est_signal  = est_signal,
    gte_df      = gte_df,
    selected    = selected,
    lambda_used = lambda_used,
    signal      = signal,
    estimand    = estimand
  )

  if (signal %in% c("outcome","aipw")) {
    ret$outcome_fit <- list(
      fit1 = o1$refit,
      fit0 = o0$refit
    )
  }

  if (signal %in% c("ipw","aipw")) {
    ret$ps_fit <- ps_fit
  }

  return(ret)
}

#' @title Bootstrap Confidence Intervals for Group Treatment Effects (GTE)
#' @description Performs a nonparametric bootstrap to compute pointwise and uniform
#'   confidence intervals for the GTE across levels of a discrete moderator.
#'
#' @param data A data.frame containing the variables for the analysis.
#' @param Y Character string, the name of the outcome variable in `data`.
#' @param D Character string, the name of the binary (0/1) treatment variable in `data`.
#' @param X Character string, the name of the discrete moderator variable in `data`.
#' @param Z Character vector, optional names of additional non-fixed-effect covariates.
#' @param FE Character vector, optional names of fixed effect variables.
#' @param estimand Character, the target estimand: 'ATE' or 'ATT'.
#' @param signal Character, the type of statistical signal: "outcome", "ipw", or "aipw".
#' @param B Integer, the number of bootstrap replications.
#' @param alpha Numeric, the significance level for confidence intervals (e.g., 0.05 for 95% CIs).
#' @param ... Additional arguments passed down to `estimateGTE`.
#'
#' @return A list containing the final results data.frame with point estimates and CIs,
#'   the matrix of bootstrap replications, the uniform coverage level, and the full-sample
#'   estimation object.
#'
bootstrapGTE <- function(
  data, Y, D, X, Z=NULL, FE=NULL,
  estimand   = c("ATE","ATT"),
  signal     = c("outcome","ipw","aipw"),
  B          = 200, alpha = 0.05,
  outcome_model_type   = "lasso",
  ps_model_type        = "lasso",
  basis_type           = c("polynomial","bspline","none"),
  include_interactions = FALSE,
  poly_degree          = 2,
  spline_df            = 4,
  spline_degree        = 2,
  lambda_seq           = NULL,
  CI = TRUE,
  verbose              = TRUE
) {
  signal     <- match.arg(signal)
  estimand   <- match.arg(estimand)
  basis_type <- match.arg(basis_type)

  if(verbose) message("1) Fit full-sample GTE...")
  fit_full <- estimateGTE(
    data                 = data,
    Y                    = Y,
    D                    = D,
    X                    = X,
    Z                    = Z,
    FE                   = FE,
    estimand             = estimand,
    signal               = signal,
    basis_type           = basis_type,
    include_interactions = include_interactions,
    poly_degree          = poly_degree,
    spline_df            = spline_df,
    spline_degree        = spline_degree,
    lambda_seq           = lambda_seq,
    verbose              = verbose
  )

  gte_full    <- fit_full$gte_df$GTE
  X_group     <- fit_full$gte_df$X
  lambda_used <- fit_full$lambda_used

  n    <- nrow(data)
  K    <- length(X_group)
  XZ0  <- fit_full$XZ_design

  if(isTRUE(CI)){
    if(verbose) message("2) Bootstrapping...")
    cl <- parallel::makeCluster(parallel::detectCores())
    doParallel::registerDoParallel(cl)
    `%dopar%` <- foreach::`%dopar%`
    idx <- seq_len(n)

    res_mat <- foreach::foreach(
      b = seq_len(B),
      .combine  = "rbind",
      .export = c("estimateGTE"),
      .packages = c("glmnet","splines")
    ) %dopar% {
      set.seed(1000 + b)
      ids    <- sample(idx, n, replace=TRUE)
      dat_b  <- data[ids,   , drop=FALSE]
      XZ_b   <- XZ0[ ids,   , drop=FALSE]

      fit_b <- estimateGTE(
        data                 = dat_b,
        Y                    = Y,
        D                    = D,
        X                    = X,
        Z                    = NULL,        # design matrix already built
        FE                   = FE,
        estimand             = estimand,
        signal               = signal,
        basis_type           = "none",      # design fixed
        include_interactions = FALSE,
        poly_degree          = poly_degree,
        spline_df            = spline_df,
        spline_degree        = spline_degree,
        XZ_design            = XZ_b,        # pass in bootstrap design
        outcome_model_type   = outcome_model_type,
        ps_model_type        = ps_model_type,
        lambda_cv            = lambda_used, 
        lambda_seq           = lambda_seq,
        verbose              = FALSE
      )
      fit_b$gte_df$GTE
    }
    parallel::stopCluster(cl)

    if(verbose) message("3) Computing SE & CIs...")
    se   <- apply(res_mat, 2, sd, na.rm=TRUE)
    ci_l <- apply(res_mat, 2, quantile, probs=alpha/2, na.rm=TRUE)
    ci_u <- apply(res_mat, 2, quantile, probs=1-alpha/2, na.rm=TRUE)
    
    uni  <- calculate_uniform_quantiles(t(res_mat), alpha)


    results <- data.frame(
      X                = as.numeric(X_group),
      GTE              = gte_full,
      SE               = se,
      CI.lower         = ci_l,
      CI.upper         = ci_u,
      CI.lower.uniform = uni$Q_j[,1],
      CI.upper.uniform = uni$Q_j[,2],
      stringsAsFactors = FALSE
    )

    if(verbose) message("Done.")
    list(
      results      = results,
      boot_results = res_mat,
      coverage     = uni$coverage,
      fit_full     = fit_full
    )    
  }
  else{
    results <- data.frame(
      X                = as.numeric(X_group),
      GTE              = gte_full,
      stringsAsFactors = FALSE
    )
    list(
      results      = results,
      fit_full     = fit_full
    )  
  }


}

