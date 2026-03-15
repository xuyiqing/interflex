
#' @title Estimate Conditional Marginal Effects (CME) using AIPW-Lasso
#' @description This function estimates the CME using outcome, IPW, or AIPW signals.
#' It supports basis expansion, fixed effects, and post-Lasso model refitting.
#'
#' @param data A data.frame containing the variables.
#' @param Y Character string, the name of the outcome variable.
#' @param D Character string, the name of the binary treatment variable.
#' @param X Character string, the name of the moderator variable.
#' @param Z Character vector, names of additional covariates.
#' @param FE Character vector, names of fixed effect variables.
#' @param estimand Character, either 'ATE' or 'ATT'.
#' @param signal Character, one of "outcome", "ipw", or "aipw".
#' @param neval Integer, number of points for the evaluation grid of X.
#' @param XZ_design A pre-built design matrix for covariates. If NULL, it's created internally.
#' @param outcome_model_type Character, one of "linear", "ridge", or "lasso".
#' @param ps_model_type Character, one of "linear", "ridge", or "lasso".
#' @param basis_type Character, one of "polynomial", "bspline", or "none".
#' @param include_interactions Logical, whether to include interactions in the design matrix.
#' @param poly_degree Integer, degree for polynomial basis expansion.
#' @param spline_df Integer, degrees of freedom for B-spline basis.
#' @param spline_degree Integer, degree for B-spline basis.
#' @param lambda_cv A list of pre-specified lambda values for glmnet, used during bootstrapping.
#' @param lambda_seq A custom sequence of lambda values for glmnet cross-validation.
#' @param reduce.dimension Character, method for smoothing the final signal: "bspline" or "kernel".
#' @param best_span A pre-specified span for LOESS smoothing, used during bootstrapping.
#' @param x.eval A numeric vector for the evaluation grid of X.
#' @param verbose Logical, whether to print progress messages.
#' @return A list containing the estimation results.
#'
estimateCME <- function(
    data,
    Y,          # name of outcome variable in 'data'
    D,          # name of treatment variable in 'data'
    X,          # name of focal variable in 'data'
    Z = NULL,   # character vector of additional covariates
    FE = NULL,  # character vector of fixed effect variable names (if any)
    
    # Which estimator(s) to compute?
    estimand = c('ATE','ATT'),
    signal = c("outcome","ipw","aipw"),
    
    neval = 100,
    # Optionally supply a pre-built design matrix for (X,Z):
    XZ_design = NULL,  
    # If NULL, we'll build the design matrix internally. 
    # If not NULL, we skip internal basis expansions and use it directly.
    
    outcome_model_type = "lasso", # can be "linear", "ridge", or "lasso"
    ps_model_type      = "lasso", # can be "linear", "ridge", or "lasso"
    
    # polynomial or B-spline expansion in the design matrix
    basis_type         = c("polynomial", "bspline", "none"),
    include_interactions = FALSE,
    poly_degree        = 2,      # only used if basis_type="polynomial"
    
    # B-spline parameters (used if basis_type="bspline", or for final CME fit)
    spline_df          = 5,
    spline_degree      = 3,
    
    lambda_cv = NULL, # stored penalty parameters
    lambda_seq        = NULL, # optional custom lambda sequence for glmnet
    reduce.dimension  = c("bspline","kernel"), 
    best_span         = NULL,
    x.eval            = NULL, # grid of X values for final CME curve
    # --- CHANGE MADE HERE ---
    # The `selected_covars` argument has been removed from this function signature.
    # It was not being used as an input, so removing it improves clarity.
    # The function still calculates and returns the selected covariates as an output.
    verbose           = TRUE
) {
  
  estimand <- match.arg(estimand)
  out_fit1 <- NULL
  out_fit0 <- NULL
  ps_fit <- NULL
  mu1_hat <- NULL
  mu0_hat <- NULL
  p_hat <- NULL
  
  if (verbose) cat("Step 1: Matching arguments and setting defaults...\n")
  # 1. Match arguments:
  signal           <- match.arg(signal)
  basis_type       <- match.arg(basis_type)
  reduce.dimension <- match.arg(reduce.dimension)
  
  if (verbose) cat("Step 2: Performing basic checks on input data...\n")
  # 2. Basic checks
  stopifnot(is.data.frame(data))
  if (!(Y %in% names(data))) stop("Variable name for Y not found in data.")
  if (!(D %in% names(data))) stop("Variable name for D not found in data.")
  if (!(X %in% names(data))) stop("Variable name for X not found in data.")
  
  Yvec <- data[[Y]]
  Dvec <- data[[D]]
  Xvec <- data[[X]]
  
  # Separate Z into non-fixed-effect covariates.
  if (!is.null(Z)) {
    if (!all(Z %in% names(data))) stop("Some variables in Z not found in data.")
    Z_nonFE <- data[, Z, drop=FALSE]
  } else {
    Z_nonFE <- NULL
  }
  
  if (verbose) cat("Step 3: Processing fixed effects (if provided)...\n")
  # Process fixed effects if provided: convert each FE variable to factor dummies,
  # dropping the first level to avoid collinearity.
  if (!is.null(FE)) {
    if (!all(FE %in% names(data))) stop("Some fixed effects variables not found in data.")
    FE_mat <- data[, FE, drop=FALSE]
    FE_dummies <- lapply(names(FE_mat), function(varname) {
      fac <- as.factor(FE_mat[[varname]])
      mm <- model.matrix(~ fac)[, -1, drop=FALSE]  # drop the intercept (first column)
      colnames(mm) <- paste0(varname, "_", levels(fac)[-1])
      mm
    })
    FE_dummies <- do.call(cbind, FE_dummies)
    if (verbose) cat("Fixed effects are included as dummies with no basis expansion.\n")
  } else {
    FE_dummies <- NULL
  }
  
  n <- nrow(data)
  if (!all(Dvec %in% c(0,1))) stop("Treatment variable D must be binary 0/1.")
  
  if (verbose) cat("Step 4: Setting up evaluation grid for X (if not provided)...\n")
  # If x.eval is not provided, build a default grid
  if (is.null(x.eval)) {
    x.eval <- seq(min(Xvec), max(Xvec), length.out = neval)
  }
  
  if (verbose) cat("Step 5: Building design matrix...\n")
  # 3. Optionally build XZ_design if not supplied
  if (is.null(XZ_design)) {
    build_design_matrix <- function(X_vec, Z_nonFE, FE_dummies) {
      # helper for B-spline expansion
      build_bspline_cols <- function(vec, varname, df, degree) {
        m <- splines::bs(vec, df=df, degree=degree)
        colnames(m) <- paste0(varname,"_bs", seq_len(ncol(m)))
        m
      }
      
      # Expand X
      if (basis_type == "none") {
        X_mat <- matrix(X_vec, ncol=1)
        colnames(X_mat) <- "X"
      } else if (basis_type == "polynomial") {
        if (poly_degree > 1) {
          poly_mat <- poly(X_vec, degree=poly_degree, raw=TRUE)
          colnames(poly_mat) <- paste0("X_poly", seq_len(poly_degree))
          X_mat <- poly_mat
        } else {
          X_mat <- matrix(X_vec, ncol=1, dimnames=list(NULL,"X_poly1"))
        }
      } else {  # bspline
        X_mat <- build_bspline_cols(X_vec, "X", spline_df, spline_degree)
      }
      
      # Expand Z_nonFE using basis expansion if needed
      Z_expanded <- NULL
      if (!is.null(Z_nonFE) && ncol(Z_nonFE) > 0) {
        Z_expanded_list <- lapply(seq_len(ncol(Z_nonFE)), function(j) {
          zj <- Z_nonFE[, j]
          zname <- colnames(Z_nonFE)[j]
          if (basis_type == "none") {
            matz <- matrix(zj, ncol=1)
            colnames(matz) <- zname
            matz
          } else if (basis_type=="polynomial") {
            if (poly_degree > 1) {
              zpoly <- poly(zj, degree=poly_degree, raw=TRUE)
              colnames(zpoly) <- paste0(zname,"_poly", seq_len(poly_degree))
              zpoly
            } else {
              matrix(zj, ncol=1, dimnames=list(NULL, zname))
            }
          } else {  # bspline
            build_bspline_cols(zj, zname, spline_df, spline_degree)
          }
        })
        Z_expanded <- do.call(cbind, Z_expanded_list)
      }
      
      # Combine the expanded Z and FE dummies (the FE dummies are NOT expanded)
      if (!is.null(FE_dummies)) {
        if (is.null(Z_expanded)) {
          Z_final <- FE_dummies
        } else {
          Z_final <- cbind(Z_expanded, FE_dummies)
        }
      } else {
        Z_final <- Z_expanded
      }
      
      # Add interactions between X and Z_expanded
      int_XZ <- NULL
      if (include_interactions && !is.null(Z_expanded) && ncol(Z_expanded) > 0) {
        for (i in seq_len(ncol(X_mat))) {
          for (j in seq_len(ncol(Z_expanded))) {
            tmp <- X_mat[, i] * Z_expanded[, j]
            int_XZ <- cbind(int_XZ, tmp)
          }
        }
        if (!is.null(int_XZ)) {
          colnames(int_XZ) <- apply(expand.grid(colnames(X_mat), colnames(Z_expanded)), 1,
                                      paste, collapse=".")
        }
      }
      
      # Add interactions among Z_expanded variables (only between different Z's)
      int_ZZ <- NULL
      if (include_interactions && !is.null(Z_expanded) && ncol(Z_expanded) > 1) {
        # Use combinations of two distinct columns
        combn_names <- combn(colnames(Z_expanded), 2, simplify = FALSE)
        int_Z_list <- lapply(combn_names, function(pair) {
          Z_expanded[, pair[1]] * Z_expanded[, pair[2]]
        })
        int_ZZ <- do.call(cbind, int_Z_list)
        colnames(int_ZZ) <- sapply(combn_names, function(pair) paste(pair, collapse="."))
      }
      
      # Combine: X_mat, Z_final (which includes both Z_expanded and FE dummies), 
      # interactions between X and Z_expanded, and interactions among Z_expanded.
      cbind(X_mat, Z_final, int_XZ, int_ZZ)
    }
    XZ_design <- build_design_matrix(Xvec, Z_nonFE, FE_dummies)
    if (verbose) cat("Design matrix successfully built.\n")
  }
  
  # Identify FE columns in XZ_design (if any)
  if (!is.null(FE)) {
    fe_pattern <- paste0("^(", paste0(FE, collapse="|"), ")_")
    fe_cols <- grep(fe_pattern, colnames(XZ_design))
    # Create penalty factor for outcome: 0 for FE columns, 1 for all others
    pf_outcome <- rep(1, ncol(XZ_design))
    pf_outcome[fe_cols] <- 0
  } else {
    fe_cols <- integer(0)
    pf_outcome <- NULL
  }
  
  if (verbose) cat("Step 6: Setting up model types based on signal argument...\n")
  # 4. If signal="aipw", we do need both outcome and PS models. 
  do_outcome <- (signal %in% c("outcome","aipw"))
  do_ps      <- (signal %in% c("ipw","aipw"))
  
  outcome_lasso <- (outcome_model_type!="linear")
  ps_lasso      <- (ps_model_type!="linear")
  
  if (verbose) cat("Step 7: Fitting outcome and propensity score models...\n")
  # 5. Fit the relevant models:
  
  # a) outcome models
  out_fit1 <- out_fit0 <- NULL
  mu1_hat  <- mu0_hat  <- rep(NA_real_, n)

  if (do_outcome) {
    if (verbose) cat("   Fitting outcome model for treated (D=1) units...\n")
    idx1 <- which(Dvec==1)
    idx0 <- which(Dvec==0)
    Y1   <- Yvec[idx1]
    Y0   <- Yvec[idx0]
    
    X1   <- XZ_design[idx1,,drop=FALSE]
    X0   <- XZ_design[idx0,,drop=FALSE]
    
    do_single_fit <- function(y_sub, x_sub, model_type, lambda_use, pf = NULL) {
      if (model_type=="linear") {
        df_temp <- data.frame(y=y_sub, x_sub)
        fit_lm  <- lm(y ~ ., data=df_temp)
        list(fit=fit_lm, type="lm", lambda = NULL)
      } 
      else if (model_type=="ridge") {
        if (is.null(lambda_use)) {
          if (is.null(pf)) {
            fit_cv <- glmnet::cv.glmnet(x=x_sub, y=y_sub, alpha=0, lambda=lambda_seq)
          } else {
            fit_cv <- glmnet::cv.glmnet(x=x_sub, y=y_sub, alpha=0, lambda=lambda_seq,
                                      penalty.factor = pf)
          }
          list(fit=fit_cv, type="glmnet", lambda = fit_cv$lambda.min)        
        } else {
          if (is.null(pf)) {
            fit <- glmnet::glmnet(x=x_sub, y=y_sub, alpha=0, lambda=lambda_use)
          } else {
            fit <- glmnet::glmnet(x=x_sub, y=y_sub, alpha=0, lambda=lambda_use,
                                penalty.factor = pf)
          }
          list(fit=fit, type="glmnet", lambda = lambda_use)  
        }
      } 
      else if (model_type=="lasso") {
        if (is.null(lambda_use)) {
          if (is.null(pf)) {
            fit_cv <- glmnet::cv.glmnet(x=x_sub, y=y_sub, alpha=1, lambda=lambda_seq)
          } else {
            fit_cv <- glmnet::cv.glmnet(x=x_sub, y=y_sub, alpha=1, lambda=lambda_seq,
                                      penalty.factor = pf)
          }
          list(fit=fit_cv, type="glmnet", lambda = fit_cv$lambda.min)        
        } else {
          if (is.null(pf)) {
            fit <- glmnet::glmnet(x=x_sub, y=y_sub, alpha=1, lambda=lambda_use)
          } else {
            fit <- glmnet::glmnet(x=x_sub, y=y_sub, alpha=1, lambda=lambda_use,
                                penalty.factor = pf)
          }
          list(fit=fit, type="glmnet", lambda = lambda_use)  
        }
      } 
      else {
        stop("Unsupported outcome_model_type: ", model_type)
      }
    }
    
    if (outcome_model_type=='linear') {
      out_fit1 <- do_single_fit(Y1, X1, outcome_model_type, NULL, NULL)
      out_fit0 <- do_single_fit(Y0, X0, outcome_model_type, NULL, NULL)
    } else {
      if (is.null(lambda_cv)) { # cross-validation
        out_fit1 <- do_single_fit(Y1, X1, outcome_model_type, NULL, pf_outcome)
        out_fit0 <- do_single_fit(Y0, X0, outcome_model_type, NULL, pf_outcome)
      } else {
        lambda_outcome1 <- lambda_cv[['outcome1']]
        lambda_outcome0 <- lambda_cv[['outcome0']]
        out_fit1 <- do_single_fit(Y1, X1, outcome_model_type, lambda_outcome1, pf_outcome)
        out_fit0 <- do_single_fit(Y0, X0, outcome_model_type, lambda_outcome0, pf_outcome)
      }
    }
  }
  
  # b) propensity score model
  ps_fit <- rep(NA_real_, n)
  if (do_ps) {
    if (verbose) cat("   Fitting propensity score model (no FE)...\n")
    
    # Exclude FE columns from PS design
    if (!is.null(FE)) {
      ps_cols <- setdiff(seq_len(ncol(XZ_design)), fe_cols)
      XZ_design_ps <- XZ_design[, ps_cols, drop=FALSE]
    } else {
      XZ_design_ps <- XZ_design
    }
    
    do_ps_fit <- function(D, Xmat, model_type, lambda_use) {
      if (model_type=="linear") {
        df_temp <- data.frame(d=D, Xmat)
        fit_glm <- glm(d ~ ., data=df_temp, family=binomial("logit"))
        list(fit=fit_glm, type="glm", lambda = NULL)
      } 
      else if (model_type=="ridge") {
        if (is.null(lambda_use)) {
          fit_cv <- glmnet::cv.glmnet(x=Xmat, y=D, alpha=0, family="binomial", lambda=lambda_seq)
          list(fit=fit_cv, type="glmnet", lambda = fit_cv$lambda.min)        
        } else {
          fit <- glmnet::glmnet(x=Xmat, y=D, alpha=0, family="binomial", lambda=lambda_use)
          list(fit=fit, type="glmnet", lambda = lambda_use)          
        }
      } 
      else if (model_type=="lasso") {
        if (is.null(lambda_use)) {
          fit_cv <- glmnet::cv.glmnet(x=Xmat, y=D, alpha=1, family="binomial", lambda=lambda_seq)
          list(fit=fit_cv, type="glmnet", lambda = fit_cv$lambda.min)        
        } else {
          fit <- glmnet::glmnet(x=Xmat, y=D, alpha=1, family="binomial", lambda=lambda_use)
          list(fit=fit, type="glmnet", lambda = lambda_use)          
        }
      } 
      else {
        stop("Unsupported ps_model_type: ", model_type)
      }
    }
    
    if (ps_model_type=='linear') {
      ps_fit <- do_ps_fit(Dvec, XZ_design_ps, ps_model_type, NULL)      
    } else {
      if (is.null(lambda_cv)) { # cross-validation
        ps_fit <- do_ps_fit(Dvec, XZ_design_ps, ps_model_type, NULL)
      } else {
        lambda_treatment <- lambda_cv[['treatment']]
        ps_fit <- do_ps_fit(Dvec, XZ_design_ps, ps_model_type, lambda_treatment)
      }
    }
  }
  
  if (verbose) cat("Step 8: Performing post-lasso variable selection (if applicable)...\n")
  # 6. Post-lasso selection if needed:
  selected_covars <- NULL
  lambda_cv_save <- NULL
  active_vars <- function(cv_obj, x_mat) {
    coefs <- coef(cv_obj$fit, s=cv_obj$fit$lambda.min)
    nz    <- which(as.numeric(coefs) != 0)
    setdiff(nz,1) - 1  # subtract 1 to align with col indices in x_mat
  }

  loss.save <- list()
  
  if (signal=="outcome" && outcome_lasso) {
    lambda_cv_save <- list(outcome1 = out_fit1$lambda, outcome0 = out_fit0$lambda)
    idx1 <- active_vars(list(fit=out_fit1$fit), XZ_design)
    idx0 <- active_vars(list(fit=out_fit0$fit), XZ_design)
    union_active <- unique(c(idx1, idx0))
    selected_covars <- list(outcome1 = idx1, outcome0 = idx0)
    if (length(union_active) > 0) {
      if (verbose) cat("   Post-lasso: refitting outcome models on selected variables...\n")
      idx1_full <- which(Dvec==1)
      idx0_full <- which(Dvec==0)
      X1_sub <- XZ_design[idx1_full, idx1, drop=FALSE]
      X0_sub <- XZ_design[idx0_full, idx0, drop=FALSE]
      df1_sub <- data.frame(y=Yvec[idx1_full], X1_sub)
      df0_sub <- data.frame(y=Yvec[idx0_full], X0_sub)
      lm1 <- lm(y ~ ., data=df1_sub)
      lm0 <- lm(y ~ ., data=df0_sub)
      out_fit1 <- list(fit=lm1, type="lm")
      out_fit0 <- list(fit=lm0, type="lm")
      loo_rmse_lm <- function(fit) {
        e   <- residuals(fit)
        h   <- hatvalues(fit)
        ecv <- e / (1 - h)          # LOO residuals
        sqrt(mean(ecv^2))
      }
      loss.save[['outcome1']] <- loo_rmse_lm(lm1)
      loss.save[['outcome0']] <- loo_rmse_lm(lm0)
    }
  } else if (signal=="ipw" && ps_lasso) {
    lambda_cv_save <- list(treatment = ps_fit$lambda)
    idxp <- active_vars(list(fit=ps_fit$fit), XZ_design_ps)
    selected_covars <- list(ps = idxp)
    if (length(idxp) > 0) {
      if (verbose) cat("   Post-lasso: refitting propensity score model on selected variables...\n")
      XZ_sub <- XZ_design_ps[, idxp, drop=FALSE]
      dfp_sub <- data.frame(d=Dvec, XZ_sub)
      final_logit <- glm(d ~ ., data=dfp_sub, family=binomial("logit"))
      ps_fit <- list(fit=final_logit, type="glm")

      loo_logloss_logit_fast <- function(fit, eps = 1e-15) {
        y   <- fit$y
        mu  <- fitted(fit)
        eta <- fit$linear.predictors
        h   <- hatvalues(fit) 
        adj     <- (h / (1 - h)) * (y - mu) / (mu * (1 - mu))
        p_loo   <- plogis(eta - adj)
        p_loo <- pmin(pmax(p_loo, eps), 1 - eps)
        -mean(y * log(p_loo) + (1 - y) * log(1 - p_loo))
      }
      #logloss_mean <- -as.numeric(logLik(final_logit)) / nobs(final_logit)
      loss.save[['treatment']] <- loo_logloss_logit_fast(final_logit) 
    }
    
  } else if (signal=="aipw" && outcome_lasso && ps_lasso) {
    lambda_cv_save <- list(outcome1 = out_fit1$lambda, outcome0 = out_fit0$lambda, treatment = ps_fit$lambda)
    idx1 <- active_vars(list(fit=out_fit1$fit), XZ_design)
    idx0 <- active_vars(list(fit=out_fit0$fit), XZ_design)
    idxp <- active_vars(list(fit=ps_fit$fit),  XZ_design_ps)
    union_active <- unique(c(idx1, idx0, idxp))
    selected_covars <- list(outcome1 = idx1, outcome0 = idx0, ps = idxp)
    if (length(union_active) > 0) {
      if (verbose) cat("   Post-lasso: refitting outcome and propensity score models on selected variables...\n")
      idx1_full <- which(Dvec==1)
      idx0_full <- which(Dvec==0)
      X1_sub <- XZ_design[idx1_full, idx1, drop=FALSE]
      X0_sub <- XZ_design[idx0_full, idx0, drop=FALSE]
      df1_sub <- data.frame(y=Yvec[idx1_full], X1_sub)
      df0_sub <- data.frame(y=Yvec[idx0_full], X0_sub)
      lm1 <- lm(y ~ ., data=df1_sub)
      lm0 <- lm(y ~ ., data=df0_sub)
      out_fit1 <- list(fit=lm1, type="lm")
      out_fit0 <- list(fit=lm0, type="lm")
      
      XZ_sub <- XZ_design_ps[, idxp, drop=FALSE]
      dfp_sub <- data.frame(d=Dvec, XZ_sub)
      final_logit <- glm(d ~ ., data=dfp_sub, family=binomial("logit"))
      ps_fit <- list(fit=final_logit, type="glm")

      loo_rmse_lm <- function(fit) {
        e   <- residuals(fit)
        h   <- hatvalues(fit)
        ecv <- e / (1 - h)          # LOO residuals
        sqrt(mean(ecv^2))
      }
      loss.save[['outcome1']] <- loo_rmse_lm(lm1)
      loss.save[['outcome0']] <- loo_rmse_lm(lm0)
      loo_logloss_logit_fast <- function(fit, eps = 1e-15) {
        y   <- fit$y
        mu  <- fitted(fit)
        eta <- fit$linear.predictors
        h   <- hatvalues(fit) 
        adj     <- (h / (1 - h)) * (y - mu) / (mu * (1 - mu))
        p_loo   <- plogis(eta - adj)
        p_loo <- pmin(pmax(p_loo, eps), 1 - eps)
        -mean(y * log(p_loo) + (1 - y) * log(1 - p_loo))
      }
      #logloss_mean <- -as.numeric(logLik(final_logit)) / nobs(final_logit)
      loss.save[['treatment']] <- loo_logloss_logit_fast(final_logit) 
    }
  }
  
  if (verbose) cat("Step 9: Generating predictions from fitted models...\n")
  # 7. Predictions
  predict_outcome <- function(subfit, newX) {
    if (subfit$type=="lm") {
      vars_used <- names(subfit$fit$coefficients)[-1]
      df_temp   <- data.frame(newX[, vars_used, drop=FALSE])
      colnames(df_temp) <- vars_used
      as.numeric(predict(subfit$fit, newdata=df_temp))
    } else if (subfit$type=="glmnet") {
      as.numeric(predict(subfit$fit, newx=newX, s="lambda.min"))
    } else {
      stop("Unknown subfit$type for outcome model.")
    }
  }
  predict_ps <- function(subfit, newX) {
    if (subfit$type=="glm") {
      vars_used <- names(subfit$fit$coefficients)[-1]
      df_temp   <- data.frame(newX[, vars_used, drop=FALSE])
      colnames(df_temp) <- vars_used
      as.numeric(predict(subfit$fit, newdata=df_temp, type="response"))
    } else if (subfit$type=="glmnet") {
      as.numeric(predict(subfit$fit, newx=newX, s="lambda.min", type="response"))
    } else {
      stop("Unknown subfit$type for ps model.")
    }
  }
  
  if (do_outcome) {
    mu1_hat <- predict_outcome(out_fit1, XZ_design)
    mu0_hat <- predict_outcome(out_fit0, XZ_design)
  }
  if (do_ps) {
    p_hat <- predict_ps(ps_fit, XZ_design_ps)
    p_hat <- pmin(pmax(p_hat, 1e-2), 1-1e-2)
  }
  
  if (verbose) cat("Step 10: Computing signals and performing dimension reduction...\n")
  # 8. Compute signals & do dimension reduction 
  cme_out_eval <- rep(NA_real_, length(x.eval))
  cme_ipw_eval <- rep(NA_real_, length(x.eval))
  cme_dr_eval  <- rep(NA_real_, length(x.eval))
  do_reduce_dim <- function(y_signal, x_vec, x_eval, method,
                            spline_df, spline_degree,
                            best_span) {
    if (method=="bspline") {
      fit_lm <- lm(y_signal ~ splines::bs(x_vec, df=spline_df, degree=spline_degree))
      preds  <- predict(fit_lm, newdata=data.frame(x_vec=x_eval))
      return(list(preds=preds, best_span = NULL))
    } else {
      if (is.null(best_span)) {
        cand_spans <- seq(0.1, 1, by=0.05)
        cv_loess_fold <- function(span, x, y, K=10, degree=1) {
          n <- length(y)
          folds <- sample(rep(1:K, length.out=n))
          mse_vec <- numeric(K)
          for (k in seq_len(K)) {
            tr <- which(folds != k)
            te <- which(folds == k)
            model_k <- loess(y[tr] ~ x[tr], span=span, degree=degree)
            pred_k  <- predict(model_k, newdata=x[te])
            mse_vec[k] <- mean((y[te]-pred_k)^2, na.rm=TRUE)
          }
          mean(mse_vec, na.rm=TRUE)
        }

        cv_errs <- sapply(cand_spans, cv_loess_fold, x=x_vec, y=y_signal)
        best_span <- cand_spans[which.min(cv_errs)]
      }
      fit_loess <- loess(y_signal ~ x_vec, span=best_span, degree=1)
      preds     <- predict(fit_loess, newdata=data.frame(x_vec=x_eval))
      return(list(preds=preds, best_span = best_span))
    }
  }
  if(estimand == 'ATE'){
    p_hat_gam <- NULL
    if (signal == "outcome") {
      est.signal <- out_signal <- mu1_hat - mu0_hat
      cme_out_eval <- do_reduce_dim(out_signal, Xvec, x.eval, reduce.dimension,
                                    spline_df, spline_degree, best_span)
      cme_df <- data.frame(X.eval = x.eval, CME = cme_out_eval$preds)
      best_span <- cme_out_eval$best_span
    }
    if (signal == "ipw") {
      est.signal <- ipw_signal <- Dvec * Yvec / p_hat - (1 - Dvec) * Yvec / (1 - p_hat)
      cme_ipw_eval <- do_reduce_dim(ipw_signal, Xvec, x.eval, reduce.dimension,
                                    spline_df, spline_degree, best_span)
      cme_df <- data.frame(X.eval = x.eval, CME = cme_ipw_eval$preds)
      best_span <- cme_ipw_eval$best_span
    }
    if (signal=="aipw") {
      est.signal <- dr_signal <- (mu1_hat - mu0_hat) +
        Dvec * (Yvec - mu1_hat) / p_hat - (1 - Dvec) * (Yvec - mu0_hat) / (1 - p_hat)

      cme_dr_eval <- do_reduce_dim(dr_signal, Xvec, x.eval, reduce.dimension,
                                     spline_df, spline_degree, best_span)
      cme_df <- data.frame(X.eval = x.eval, CME = cme_dr_eval$preds)
      best_span <- cme_dr_eval$best_span
    }    
  }
  else if(estimand == 'ATT'){
    fit_gam <- mgcv::gam(Dvec ~ s(Xvec), family = binomial)
    p_hat_gam <- predict(fit_gam, type = "response")
    
    if (signal == "outcome") {
      est.signal <- out_signal <- (Yvec - mu0_hat)*Dvec/p_hat_gam
      cme_out_eval <- do_reduce_dim(out_signal, Xvec, x.eval, reduce.dimension,
                                    spline_df, spline_degree, best_span)
      cme_df <- data.frame(X.eval = x.eval, CME = cme_out_eval$preds)
      best_span <- cme_out_eval$best_span
    }
    if (signal == "ipw") {
      est.signal <- ipw_signal <- Yvec*(Dvec - p_hat)/((1 - p_hat)*p_hat_gam)
      cme_ipw_eval <- do_reduce_dim(ipw_signal, Xvec, x.eval, reduce.dimension,
                                    spline_df, spline_degree, best_span)
      cme_df <- data.frame(X.eval = x.eval, CME = cme_ipw_eval$preds)
      best_span <- cme_ipw_eval$best_span
    }
    if (signal=="aipw") {
      est.signal <- dr_signal <- ((Yvec - mu0_hat)*(Dvec - (1-Dvec)*p_hat/(1-p_hat)))/p_hat_gam
        
      cme_dr_eval <- do_reduce_dim(dr_signal, Xvec, x.eval, reduce.dimension,
                                     spline_df, spline_degree, best_span)
      cme_df <- data.frame(X.eval = x.eval, CME = cme_dr_eval$preds)
      best_span <- cme_dr_eval$best_span
    }     
  }

  
  if (verbose) cat("Step 11: Estimation complete. Returning results.\n")
  out <- list(
    XZ_design   = XZ_design,
    mu1_fit     = out_fit1,
    mu0_fit     = out_fit0,
    ps_fit      = ps_fit,
    p_hat_gam   = p_hat_gam,
    mu1_hat     = mu1_hat,
    mu0_hat     = mu0_hat,
    p_hat       = p_hat,
    est.signal  = est.signal,
    cme_df      = cme_df,
    signal      = signal,
    estimand    = estimand,
    reduce.dimension = reduce.dimension,
    best_span   = best_span,
    lambda_cv = lambda_cv_save,
    loss = loss.save,
    selected_covars = selected_covars
  )
  return(out)
}



#' @title Bootstrap Confidence Intervals for Conditional Marginal Effects (CME)
#' @description This function performs a nonparametric bootstrap to compute pointwise and
#' uniform confidence intervals for the CME curve.
#'
#' @param data A data.frame containing the variables.
#' @param Y Character string, the name of the outcome variable.
#' @param D Character string, the name of the binary treatment variable.
#' @param X Character string, the name of the moderator variable.
#' @param Z Character vector, names of additional covariates.
#' @param FE Character vector, names of fixed effect variables.
#' @param estimand Character, either 'ATE' or 'ATT'.
#' @param signal Character, one of "outcome", "ipw", or "aipw".
#' @param B Integer, number of bootstrap replications.
#' @param alpha Numeric, significance level for confidence intervals.
#' @param ... Additional arguments passed to `estimateCME`.
#' @return A list with bootstrap results, including a data.frame with the final CME estimates and CIs.
#'
bootstrapCME <- function(
    data,
    Y, D, X, Z = NULL, FE = NULL,
    
    # Which single estimator to use?
    estimand = 'ATE',
    signal = c("outcome","ipw","aipw"),
    
    # Number of bootstrap draws and alpha
    B = 200,
    alpha = 0.05,
    
    # Model types
    outcome_model_type = "lasso",
    ps_model_type      = "lasso",
    
    # Polynomial or B-spline expansions
    basis_type         = c("polynomial","bspline","none"),
    include_interactions = TRUE,
    poly_degree        = 2,
    
    # B-spline parameters
    spline_df          = 4,
    spline_degree      = 2,
    
    # LASSO/GLMNET parameters
    lambda_seq         = NULL,
    
    # Dimension-reduction method for final CME curve
    reduce.dimension   = c("bspline","kernel"),
    best_span          = NULL,
    neval = 100, 
    x.eval             = NULL,
    CI = TRUE,
    verbose = TRUE
) {
  if (verbose) message("BootstrapCME Step 1: Running baseline CME estimation on full data...")
  signal           <- match.arg(signal)
  basis_type       <- match.arg(basis_type)
  reduce.dimension <- match.arg(reduce.dimension)
  
  # 1. Fit once on the full dataset
  fit_full <- estimateCME(
    data               = data,
    Y                  = Y,
    D                  = D,
    X                  = X,
    Z                  = Z,
    FE                 = FE,
    estimand           = estimand,
    signal             = signal,
    outcome_model_type = outcome_model_type,
    ps_model_type      = ps_model_type,
    basis_type         = basis_type,
    include_interactions= include_interactions,
    poly_degree        = poly_degree,
    spline_df          = spline_df,
    spline_degree      = spline_degree,
    lambda_seq         = lambda_seq,
    reduce.dimension   = reduce.dimension,
    best_span          = best_span,
    neval = neval,
    x.eval             = x.eval,
    verbose            = verbose
  )
  if (verbose) message("Baseline CME estimation complete.")
  
  cme_full <- fit_full$cme_df$CME
  x.eval   <- fit_full$cme_df$X.eval
  nEval    <- length(x.eval)
  n        <- nrow(data)
  
  XZ_full        <- fit_full$XZ_design
  lambda_cv      <- fit_full$lambda_cv
  best_span_full <- fit_full$best_span

  if(isTRUE(CI)){
    
    if (verbose) message("BootstrapCME Step 2: Setting up bootstrap storage and parallel cluster...")
    cme_mat_bs <- matrix(NA, nrow = B, ncol = nEval)
    idx_seq <- seq_len(n)
    
    if (!requireNamespace("doParallel", quietly = TRUE)) {
      stop("Package 'doParallel' not installed. Please install or remove parallel usage.")
    }
    nCores <- parallel::detectCores()
    cl <- parallel::makeCluster(nCores)
    doParallel::registerDoParallel(cl)
    `%dopar%` <- foreach::`%dopar%`
    
    if (verbose) message("BootstrapCME Step 3: Starting bootstrap loop...")
    res_list <- foreach::foreach(
      b = 1:B,
      .combine = "rbind",
      .export  = "estimateCME",
      .packages = c("splines","glmnet","mgcv")
    ) %dopar% {
      set.seed(1000 + b)
      idx_b <- sample(idx_seq, size = n, replace = TRUE)
      data_b <- data[idx_b, , drop = FALSE]
      XZ_b <- XZ_full[idx_b, , drop = FALSE]
      use_out_model <- outcome_model_type
      use_ps_model  <- ps_model_type
      
      fit_b <- estimateCME(
        data               = data_b,
        Y                  = Y,
        D                  = D,
        X                  = X,
        Z                  = NULL, # design matrix already provided
        FE                 = FE,   
        estimand           = estimand,
        signal             = signal,
        XZ_design          = XZ_b,
        outcome_model_type = use_out_model,
        ps_model_type      = use_ps_model,
        basis_type         = "none",  # not needed if XZ_design is given
        include_interactions= FALSE,
        spline_df          = spline_df,
        spline_degree      = spline_degree,
        lambda_seq         = lambda_seq,
        reduce.dimension   = reduce.dimension,
        best_span          = best_span_full,
        x.eval             = x.eval,
        # --- CHANGE MADE HERE ---
        # The `selected_covars` argument has been removed from this function call.
        # The logic remains correct, as each bootstrap run now correctly performs its own
        # variable selection internally.
        lambda_cv = lambda_cv,
        verbose            = FALSE
      )
      fit_b$cme_df$CME
    }
    parallel::stopCluster(cl)
    if (verbose) message("BootstrapCME Step 4: Bootstrap loop complete.")
    
    cme_mat_bs[,] <- res_list
    
    if (verbose) message("BootstrapCME Step 5: Computing bootstrap SE and confidence intervals...")
    alpha_lower <- alpha / 2
    alpha_upper <- 1 - alpha_lower
    
    cme_se <- apply(cme_mat_bs, 2, sd, na.rm = TRUE)
    cme_ci_l <- apply(cme_mat_bs, 2, quantile, probs = alpha_lower, na.rm = TRUE)
    cme_ci_u <- apply(cme_mat_bs, 2, quantile, probs = alpha_upper, na.rm = TRUE)
    
    theta_mat <- t(cme_mat_bs)
    uni_res   <- calculate_uniform_quantiles(theta_mat, alpha)

    
    cme_ci_l_uni <- uni_res$Q_j[,1]
    cme_ci_u_uni <- uni_res$Q_j[,2]
    coverage_uni <- uni_res$coverage
    
    final_df <- data.frame(
      X                = x.eval,
      CME              = cme_full,
      SE               = cme_se,
      CI.lower         = cme_ci_l,
      CI.upper         = cme_ci_u,
      CI.lower.uniform = cme_ci_l_uni,
      CI.upper.uniform = cme_ci_u_uni
    )
    
    if (verbose) message("BootstrapCME Step 6: Finalizing and returning bootstrap results.")
    out_list <- list(
      results    = final_df,
      boot_results = cme_mat_bs,
      coverage = coverage_uni,
      fit_full = fit_full
    )
    
    return(out_list)    
  }
  else{
    final_df <- data.frame(
      X                = x.eval,
      CME              = cme_full
    )
    out_list <- list(
      results    = final_df,
      fit_full = fit_full
    )
    return(out_list)    
  }
}

