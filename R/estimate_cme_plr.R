#' @title Estimate Conditional Marginal Effects (CME) using PO-Lasso for PLRM
#' @description This function implements the "partialling out" (PO-Lasso) estimator
#' for a Partially Linear Regression Model (PLRM) with a continuous treatment.
#'
#' @param data A data.frame containing the variables.
#' @param Y Character string, the name of the outcome variable.
#' @param D Character string, the name of the continuous treatment variable.
#' @param X Character string, the name of the moderator variable.
#' @param Z Character vector, names of additional covariates.
#' @param FE Character vector, names of fixed effect variables.
#' @param XZ_design A pre-built design matrix. If NULL, it's created internally.
#' @param outcome_model_type Character, one of "linear", "ridge", or "lasso".
#' @param treatment_model_type Character, one of "linear", "ridge", or "lasso".
#' @param basis_type Character, one of "polynomial", "bspline", or "none".
#' @param include_interactions Logical, whether to include interactions.
#' @param poly_degree Integer, degree for polynomial basis expansion.
#' @param spline_df Integer, degrees of freedom for B-spline basis.
#' @param spline_degree Integer, degree for B-spline basis.
#' @param lambda_cv List of pre-specified lambda values for glmnet.
#' @param lambda_seq Custom sequence of lambda values for glmnet cross-validation.
#' @param reduce.dimension Method for final smoothing: "bspline" or "kernel".
#' @param bw Pre-specified bandwidth for kernel smoothing.
#' @param x.eval Numeric vector for the evaluation grid of X.
#' @param verbose Logical, whether to print progress messages.
#' @return A list containing the estimation results.
#'
estimateCME_PLR <- function(
    data,
    Y,          # name of outcome variable in 'data'
    D,          # name of treatment variable in 'data'
    X,          # name of focal variable in 'data'
    Z = NULL,   # character vector of additional covariates
    FE = NULL,  # character vector of fixed-effect variable names (if any)
    XZ_design = NULL,  # Optionally supply a pre-built design matrix for (X,Z,FE)
                       # If NULL, we build it internally (including FE dummies);
                       # otherwise we skip internal basis expansions and use it directly.

    outcome_model_type   = "lasso",  # can be "linear", "ridge", or "lasso"
    treatment_model_type = "lasso",  # can be "linear", "ridge", or "lasso"

    # polynomial or B-spline expansion in the design matrix
    basis_type           = c("polynomial", "bspline", "none"),
    include_interactions = FALSE,
    poly_degree          = 2,       # only used if basis_type="polynomial"

    # B-spline parameters (used if basis_type="bspline", or for final CME fit)
    spline_df            = 4,
    spline_degree        = 2,

    lambda_cv            = NULL,    # stored penalty parameters (list with
                                    #   elements 'outcome' and 'treatment' if pre-tuned)
    lambda_seq           = NULL,    # optional custom lambda sequence for glmnet
    reduce.dimension     = c("bspline","kernel"),
    bw                   = NULL,
    x.eval               = NULL,    # grid of X values for final CME curve
    neval = 100,
    # --- CHANGE MADE HERE ---
    # The `selected_covars` argument has been removed from this function signature.
    # It was not being used as an input. The function still correctly calculates
    # and returns the selected covariates as an output.
    verbose              = TRUE
) {

  basis_type       <- match.arg(basis_type)
  reduce.dimension <- match.arg(reduce.dimension)

  ##############################################################################
  # 0. Basic checks and extract variables
  ##############################################################################
  if (verbose) cat("Step 0: Checking inputs and extracting variables...\n")
  stopifnot(is.data.frame(data))
  if (!is.character(Y) || !is.character(D) || !is.character(X)) {
    stop("Y, D, X must be character strings corresponding to variable names in 'data'.")
  }
  if (!(Y %in% names(data))) stop("Variable name for Y not found in data.")
  if (!(D %in% names(data))) stop("Variable name for D not found in data.")
  if (!(X %in% names(data))) stop("Variable name for X not found in data.")
  if (!is.null(Z) && !all(Z %in% names(data))) {
    stop("Some variables in Z not found in data.")
  }
  if (!is.null(FE) && !all(FE %in% names(data))) {
    stop("Some variables in FE not found in data.")
  }

  Yvec <- data[[Y]]
  Dvec <- data[[D]]
  Xvec <- data[[X]]

  if (length(Yvec) != nrow(data) ||
      length(Dvec) != nrow(data) ||
      length(Xvec) != nrow(data)) {
    stop("Lengths of Y, D, X do not match number of rows in data.")
  }

  if (!is.null(Z)) {
    Zmat <- data[, Z, drop = FALSE]
  } else {
    Zmat <- NULL
  }

  n <- nrow(data)
  if (is.null(x.eval)) {
    x.eval <- seq(min(Xvec), max(Xvec), length.out = neval)
  }

  ##############################################################################
  # 1. Build or validate design matrix (XZ_design), including FE dummies
  ##############################################################################
  if (verbose) cat("Step 1: Building or validating design matrix (including FE)...\n")

  if (!is.null(XZ_design)) {
    # User supplied a design matrix.  We just check dimensions.
    if (!is.matrix(XZ_design)) {
      stop("'XZ_design' must be a matrix (or coercible to matrix).")
    }
    if (nrow(XZ_design) != n) {
      stop("XZ_design must have the same number of rows as 'data'.")
    }
  } else {
    # We need splines if basis_type == "bspline"
    if (basis_type == "bspline" && !requireNamespace("splines", quietly = TRUE)) {
      stop("Package 'splines' is required for B-spline expansions.")
    }

    build_design_matrix <- function(X_vec, Z_mat, FE_vars) {
      # Helper for B-spline expansions
      build_bspline_cols <- function(vec, varname, df, degree) {
        bs_mat <- splines::bs(vec, df = df, degree = degree)
        colnames(bs_mat) <- paste0(varname, "_bs", seq_len(ncol(bs_mat)))
        bs_mat
      }

      # 1a) Expand X
      if (basis_type == "none") {
        X_mat <- matrix(X_vec, ncol = 1)
        colnames(X_mat) <- "X"
      } else if (basis_type == "polynomial") {
        if (poly_degree > 1) {
          poly_mat <- stats::poly(X_vec, degree = poly_degree, raw = TRUE)
          colnames(poly_mat) <- paste0("X_poly", seq_len(poly_degree))
          X_mat <- poly_mat
        } else {
          X_mat <- matrix(X_vec, ncol=1)
          colnames(X_mat) <- "X_poly1"
        }
      } else {  # "bspline"
        X_mat <- build_bspline_cols(X_vec, "X", spline_df, spline_degree)
      }

      # 1b) Expand Z (non-FE) if provided
      Z_expanded <- NULL
      if (!is.null(Z_mat) && ncol(Z_mat) > 0) {
        Z_expanded_list <- lapply(seq_len(ncol(Z_mat)), function(j) {
          zj <- Z_mat[, j]
          zname <- colnames(Z_mat)[j]
          if (basis_type == "none") {
            m <- matrix(zj, ncol = 1)
            colnames(m) <- zname
            m
          } else if (basis_type == "polynomial") {
            if (poly_degree > 1) {
              zpoly <- stats::poly(zj, degree = poly_degree, raw = TRUE)
              colnames(zpoly) <- paste0(zname, "_poly", seq_len(poly_degree))
              zpoly
            } else {
              m <- matrix(zj, ncol=1)
              colnames(m) <- zname
              m
            }
          } else {  # "bspline"
            build_bspline_cols(zj, zname, spline_df, spline_degree)
          }
        })
        Z_expanded <- do.call(cbind, Z_expanded_list)
      }

      # 1c) Create FE dummies (no basis expansion, just 0/1 flags)
      FE_dummies <- NULL
      if (!is.null(FE_vars) && length(FE_vars) > 0) {
        FE_list <- lapply(FE_vars, function(varname) {
          fac   <- as.factor(data[[varname]])
          mm    <- model.matrix(~ fac)[, -1, drop=FALSE]  # drop reference level
          colnames(mm) <- paste0(varname, "_", levels(fac)[-1])
          mm
        })
        FE_dummies <- do.call(cbind, FE_list)
      }

      # 1d) Interactions between X and Z_expanded
      int_XZ <- NULL
      if (include_interactions && !is.null(Z_expanded) && ncol(Z_expanded) > 0) {
        for (i in seq_len(ncol(X_mat))) {
          for (j in seq_len(ncol(Z_expanded))) {
            tmp <- X_mat[, i] * Z_expanded[, j]
            int_XZ <- cbind(int_XZ, tmp)
          }
        }
        if (!is.null(int_XZ)) {
          colnames(int_XZ) <- apply(
            expand.grid(colnames(X_mat), colnames(Z_expanded)), 1,
            paste, collapse="."
          )
        }
      }

      # 1e) Interactions among Z_expanded
      int_ZZ <- NULL
      if (include_interactions && !is.null(Z_expanded) && ncol(Z_expanded) > 1) {
        combn_names <- combn(colnames(Z_expanded), 2, simplify = FALSE)
        int_list <- lapply(combn_names, function(pair) {
          Z_expanded[, pair[1]] * Z_expanded[, pair[2]]
        })
        int_ZZ <- do.call(cbind, int_list)
        colnames(int_ZZ) <- sapply(combn_names, function(p) paste(p, collapse = "."))
      }

      # Combine columns: X_mat, Z_expanded, FE_dummies, int_XZ, int_ZZ
      cbind(X_mat, Z_expanded, FE_dummies, int_XZ, int_ZZ)
    }

    XZ_design <- build_design_matrix(Xvec, Zmat, FE)
    if (verbose) cat("  -> Design matrix built with",
                      ncol(XZ_design), "columns (including FE dummies).\n")
  }

  ##############################################################################
  # 2. Identify which columns are FE dummies, build penalty-factor vectors
  ##############################################################################
  if (verbose) cat("Step 2: Identifying FE columns and building penalty factors...\n")
  if (!is.null(FE)) {
    # Look for column names starting with any FE prefix + "_"
    fe_pattern <- paste0("^(", paste0(FE, collapse="|"), ")_")
    fe_cols    <- grep(fe_pattern, colnames(XZ_design))
    if (length(fe_cols) == 0) {
      stop("No FE dummy columns found in XZ_design. Check FE argument.")
    }
    # Outcome penalty factor: 0 for FE columns, 1 for others
    pf_outcome <- rep(1, ncol(XZ_design))
    pf_outcome[fe_cols] <- 0
    # Treatment penalty factor: same (we include FE in treatment regression unpenalized)
    pf_treatment <- pf_outcome
    if (verbose) cat("  -> Found", length(fe_cols), "FE-dummy columns; penalty.factor=0 for those.\n")
  } else {
    fe_cols      <- integer(0)
    pf_outcome   <- NULL
    pf_treatment <- NULL
    if (verbose) cat("  -> No FE specified; penalty factors default to NULL (full penalty).\n")
  }

  ##############################################################################
  # 3. Fit outcome and treatment models (with or without LASSO/Ridge)
  ##############################################################################
  if (verbose) cat("Step 3: Fitting outcome and treatment models...\n")
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for penalized regression.")
  }
  library(glmnet)

  # Generic helper to fit linear, ridge, or lasso, *with optional penalty.factor*
  do_single_fit <- function(y_sub, x_sub, model_type, lambda_use, pf = NULL) {
    if (model_type == "linear") {
      df_temp <- data.frame(y = y_sub, x_sub)
      fit_lm  <- lm(y ~ ., data = df_temp)
      return(list(fit = fit_lm, type = "lm", lambda = NULL))

    } else if (model_type == "ridge") {
      if (is.null(lambda_use)) {
        if (is.null(pf)) {
          fit_cv <- glmnet::cv.glmnet(x = x_sub, y = y_sub, alpha = 0,
                                      lambda = lambda_seq)
        } else {
          fit_cv <- glmnet::cv.glmnet(x = x_sub, y = y_sub, alpha = 0,
                                      lambda = lambda_seq,
                                      penalty.factor = pf)
        }
        return(list(fit = fit_cv, type = "glmnet", lambda = fit_cv$lambda.min))
      } else {
        if (is.null(pf)) {
          fit <- glmnet::glmnet(x = x_sub, y = y_sub, alpha = 0,
                                lambda = lambda_use)
        } else {
          fit <- glmnet::glmnet(x = x_sub, y = y_sub, alpha = 0,
                                lambda = lambda_use,
                                penalty.factor = pf)
        }
        return(list(fit = fit, type = "glmnet", lambda = lambda_use))
      }

    } else if (model_type == "lasso") {
      if (is.null(lambda_use)) {
        if (is.null(pf)) {
          fit_cv <- glmnet::cv.glmnet(x = x_sub, y = y_sub, alpha = 1,
                                      lambda = lambda_seq)
        } else {
          fit_cv <- glmnet::cv.glmnet(x = x_sub, y = y_sub, alpha = 1,
                                      lambda = lambda_seq,
                                      penalty.factor = pf)
        }
        return(list(fit = fit_cv, type = "glmnet", lambda = fit_cv$lambda.min))
      } else {
        if (is.null(pf)) {
          fit <- glmnet::glmnet(x = x_sub, y = y_sub, alpha = 1,
                                lambda = lambda_use)
        } else {
          fit <- glmnet::glmnet(x = x_sub, y = y_sub, alpha = 1,
                                lambda = lambda_use,
                                penalty.factor = pf)
        }
        return(list(fit = fit, type = "glmnet", lambda = lambda_use))
      }

    } else {
      stop("Unsupported model_type: ", model_type)
    }
  }

  # (A) Fit outcome regression: Y ~ XZ_design + FE (if any)
  if (verbose) cat("  -> Fitting outcome model (type =", outcome_model_type, ")...\n")
  if (is.null(outcome_model_type)) stop("Specify an outcome_model_type.")
  if (outcome_model_type == "linear") {
    outcome_fit <- do_single_fit(Yvec, XZ_design, "linear", NULL, NULL)
  } else {
    if (is.null(lambda_cv)) {
      outcome_fit <- do_single_fit(Yvec, XZ_design, outcome_model_type,
                                   NULL, pf_outcome)
    } else {
      lam_o <- lambda_cv[["outcome"]]
      outcome_fit <- do_single_fit(Yvec, XZ_design, outcome_model_type,
                                   lam_o, pf_outcome)
    }
  }

  # (B) Fit treatment regression: D ~ XZ_design + FE (if any)
  if (verbose) cat("  -> Fitting treatment model (type =", treatment_model_type, ")...\n")
  if (treatment_model_type == "linear") {
    treatment_fit <- do_single_fit(Dvec, XZ_design, "linear", NULL, NULL)
  } else {
    if (is.null(lambda_cv)) {
      treatment_fit <- do_single_fit(Dvec, XZ_design, treatment_model_type,
                                     NULL, pf_treatment)
    } else {
      lam_t <- lambda_cv[["treatment"]]
      treatment_fit <- do_single_fit(Dvec, XZ_design, treatment_model_type,
                                     lam_t, pf_treatment)
    }
  }

  ##############################################################################
  # 4. (Optional) Post-lasso / post-ridge variable selection and refit
  ##############################################################################
  need_tune <- (outcome_model_type != "linear" || treatment_model_type != "linear")
  if (need_tune && verbose) cat("Step 4: Performing post-selection if LASSO/Ridge was used...\n")

  selected_covars_list <- NULL
  lambda_cv_save       <- NULL
  if (need_tune) {
    # (i) Record lambdas
    lambda_cv_save <- list(
      outcome   = outcome_fit$lambda,
      treatment = treatment_fit$lambda
    )
    # (ii) Extract active sets from each penalized fit
    active_vars <- function(cv_obj, x_mat) {
      coefs <- coef(cv_obj$fit, s = cv_obj$fit$lambda.min)
      nz    <- which(as.numeric(coefs) != 0)
      setdiff(nz, 1) - 1  # subtract 1 for intercept
    }
    act_o <- active_vars(list(fit = outcome_fit$fit),   XZ_design)
    act_t <- active_vars(list(fit = treatment_fit$fit), XZ_design)
    selected_union <- unique(c(act_o, act_t))
    selected_covars_list <- list(outcome = act_o, treatment = act_t)

    if (length(selected_union) > 0) {
      # Refit outcome on act_o
      X_sub_outcome <- XZ_design[, act_o, drop = FALSE]
      df_out_sub    <- data.frame(y = Yvec, X_sub_outcome)
      final_out_lm  <- lm(y ~ ., data = df_out_sub)
      outcome_fit   <- list(fit = final_out_lm, type = "lm")

      # Refit treatment on act_t
      X_sub_treat   <- XZ_design[, act_t, drop = FALSE]
      df_tr_sub     <- data.frame(d = Dvec, X_sub_treat)
      final_tr_lm   <- lm(d ~ ., data = df_tr_sub)
      treatment_fit <- list(fit = final_tr_lm, type = "lm")
    }
  }

  ##############################################################################
  # 5. Generate in-sample predictions and residuals (signals)
  ##############################################################################
  if (verbose) cat("Step 5: Generating predictions and residual signals...\n")

  predict_model <- function(subfit, newX) {
    if (subfit$type == "lm") {
      vars_used <- names(subfit$fit$coefficients)[-1]
      df_temp   <- data.frame(newX[, vars_used, drop = FALSE])
      colnames(df_temp) <- vars_used
      return(as.numeric(predict(subfit$fit, newdata = df_temp)))
    } else if (subfit$type == "glmnet") {
      return(as.numeric(predict(subfit$fit, newx = newX, s = "lambda.min")))
    } else {
      stop("Unknown subfit$type in predict_model")
    }
  }

  outcome_hat   <- predict_model(outcome_fit,   XZ_design)
  treatment_hat <- predict_model(treatment_fit, XZ_design)

  outcome_signal   <- Yvec - outcome_hat
  treatment_signal <- Dvec - treatment_hat

  ##############################################################################
  # 6. Dimension-reduction and estimate CME curve
  ##############################################################################
  if (verbose) cat("Step 6: Estimating CME curve via", reduce.dimension, "...\n")

  if (reduce.dimension == "bspline") {
    if (verbose) cat("  -> Fitting bivariate spline of (tilde Y ~ tilde D * bs(X))...\n")
    fit_spline <- function(Y_tilde, D_tilde, X_vec) {
      df_fit <- data.frame(Y = Y_tilde, D = D_tilde, X = X_vec)
      lm(Y ~ D * splines::bs(X, df = spline_df, degree = spline_degree),
         data = df_fit)
    }
    fit_spline_model <- fit_spline(outcome_signal, treatment_signal, Xvec)
    cme_1 <- predict(fit_spline_model, newdata = data.frame(D = 1, X = x.eval))
    cme_0 <- predict(fit_spline_model, newdata = data.frame(D = 0, X = x.eval))
    cme_fit <- cme_1 - cme_0

  } else if (reduce.dimension == "kernel") {
    if (verbose) cat("  -> Running kernel-based conditional ATE via interflex::interflex()\n")
    data_k <- data.frame(
      Y = outcome_signal,
      D = treatment_signal,
      X = Xvec
    )
    if (is.null(bw)) {
      sol_k <- interflex::interflex(
        estimator = "kernel",
        Y = "Y", D = "D", X = "X",
        data = data_k,
        X.eval = x.eval,
        CV     = TRUE,
        parallel = TRUE,
        cores    = parallel::detectCores()
      )
    } else {
      sol_k <- interflex::interflex(
        estimator = "kernel",
        Y = "Y", D = "D", X = "X",
        data = data_k,
        X.eval = x.eval,
        CV     = FALSE,
        bw     = bw
      )
    }
    cme_fit <- sol_k$est.kernel[[1]][, 2]  # second column is estimate for CME

    if (!is.null(sol_k$bw) && verbose) {
      cat("  -> Selected bandwidth =", sol_k$bw, "\n")
    }
  } else {
    stop("Unsupported reduce.dimension; choose 'bspline' or 'kernel'.")
  }

  cme_df <- data.frame(
    X.eval = x.eval,
    CME_fit = cme_fit
  )
  if (verbose) cat("  -> CME estimation on grid complete.\n")

  ##############################################################################
  # 7. Package results and return
  ##############################################################################
  if (verbose) cat("Step 7: Returning results.\n")
  out_list <- list(
    # Final design matrix used
    XZ_design     = XZ_design,

    # Fitted outcome and treatment models
    outcome_fit   = outcome_fit,
    treatment_fit = treatment_fit,

    # In-sample predictions
    outcome_hat   = outcome_hat,
    treatment_hat = treatment_hat,

    # Residual signals
    outcome_signal   = outcome_signal,
    treatment_signal = treatment_signal,

    # CME curve on x.eval
    cme_df           = cme_df,
    reduce.dimension = reduce.dimension,

    # If double-selection was used
    selected_covars = selected_covars_list,
    lambda_cv       = lambda_cv_save
  )

  if (reduce.dimension == "bspline") {
    out_list$spline_fit <- fit_spline_model
  }
  if (reduce.dimension == "kernel") {
    out_list$kernel_fit <- sol_k
    out_list$bw         <- if (!is.null(sol_k$bw)) sol_k$bw else NULL
  }

  return(out_list)
}


#' @title Bootstrap Confidence Intervals for the PO-Lasso Estimator
#' @description This function performs a nonparametric bootstrap to compute pointwise and
#' uniform confidence intervals for the CME curve from a PLRM.
#'
#' @param data A data.frame containing the variables.
#' @param Y Character string, the name of the outcome variable.
#' @param D Character string, the name of the continuous treatment variable.
#' @param X Character string, the name of the moderator variable.
#' @param Z Character vector, names of additional covariates.
#' @param FE Character vector, names of fixed effect variables.
#' @param B Integer, number of bootstrap replications.
#' @param alpha Numeric, significance level for confidence intervals.
#' @param ... Additional arguments passed to `estimateCME_PLR`.
#' @return A list with bootstrap results, including a data.frame with the final CME estimates and CIs.
#'
bootstrapCME_PLR <- function(
    data,
    Y, D, X, Z = NULL, FE = NULL,

    B = 200,
    alpha = 0.05,

    outcome_model_type   = "lasso",  # can be "linear", "ridge", or "lasso"
    treatment_model_type = "lasso",  # can be "linear", "ridge", or "lasso"

    # polynomial or B-spline expansion in the design matrix
    basis_type           = c("polynomial","bspline","none"),
    include_interactions = TRUE,
    poly_degree          = 2,        # only used if basis_type="polynomial"

    # B-spline parameters (used if basis_type="bspline", or for final CME fit)
    spline_df            = 4,
    spline_degree        = 2,

    lambda_seq           = NULL,     # optional custom lambda sequence for glmnet
    reduce.dimension     = c("bspline","kernel"),
    bw                   = NULL,
    x.eval               = NULL,     # grid of X values for final CME curve
    neval = 100,
    CI = TRUE,
    verbose              = TRUE
) {
  basis_type       <- match.arg(basis_type)
  reduce.dimension <- match.arg(reduce.dimension)

  ###########################################################################
  # 1. Fit once on the full dataset (including FE)
  ###########################################################################
  if (verbose) cat("BootstrapPLR Step 1: Fitting full-sample PLR model...\n")
  fit_full <- estimateCME_PLR(
    data                 = data,
    Y                    = Y,
    D                    = D,
    X                    = X,
    Z                    = Z,
    FE                   = FE,
    outcome_model_type   = outcome_model_type,
    treatment_model_type = treatment_model_type,
    basis_type           = basis_type,
    include_interactions = include_interactions,
    poly_degree          = poly_degree,
    spline_df            = spline_df,
    spline_degree        = spline_degree,
    lambda_cv            = NULL,
    lambda_seq           = lambda_seq,
    reduce.dimension     = reduce.dimension,
    bw                   = bw,
    x.eval               = x.eval,
    neval = neval, 
    verbose              = verbose
  )
  if (verbose) cat("  -> Full-sample fit complete.\n")

  # Extract full-sample CME
  cme_out <- fit_full$cme_df$CME_fit
  x.eval  <- fit_full$cme_df$X.eval
  nEval   <- length(x.eval)
  n       <- nrow(data)

  # If kernel dimension reduction, retrieve bandwidth from full fit
  if (fit_full$reduce.dimension == "kernel") {
    bw <- fit_full$bw
    if (verbose) cat("  -> Using bandwidth =", bw, "\n")
  } else {
    bw <- NULL
  }

  # Check if double-selection (outcome + treatment) was used
  selected_covars <- fit_full$selected_covars
  lambda_cv       <- fit_full$lambda_cv
  need_tune       <- !is.null(selected_covars)

  if(CI == TRUE){
    ###########################################################################
    # 2. Prepare storage for bootstrap draws
    ###########################################################################
    if (verbose) cat("BootstrapPLR Step 2: Preparing for", B, "bootstrap draws...\n")
    fit_mat_bs <- matrix(NA, nrow = B, ncol = nEval)
    idx_seq    <- seq_len(n)

    ###########################################################################
    # 3. Set up parallel backend
    ###########################################################################
    if (verbose) cat("BootstrapPLR Step 3: Setting up parallel backend...\n")
    if (!requireNamespace("doParallel", quietly = TRUE)) {
      stop("Package 'doParallel' is required for parallel bootstrap.")
    }
    nCores <- parallel::detectCores()
    cl     <- parallel::makeCluster(nCores)
    doParallel::registerDoParallel(cl)

    ###########################################################################
    # 4. Parallel bootstrap loop
    ###########################################################################
    if (verbose) cat("BootstrapPLR Step 4: Running bootstrap loop...\n")
    `%dopar%` <- foreach::`%dopar%`

    res_list <- foreach::foreach(
      b = 1:B,
      .combine  = "rbind",
      .export   = "estimateCME_PLR",
      .packages = c("splines","glmnet","interflex")
    ) %dopar% {
      set.seed(b + 1234)

      # (a) Resample indices and data
      idx_b  <- sample(idx_seq, size = n, replace = TRUE)
      data_b <- data[idx_b, , drop = FALSE]

      # (b) Grab the corresponding rows of XZ_design (including FE dummies)
      XZ_b <- fit_full$XZ_design[idx_b, , drop = FALSE]

      # (c) Re-fit PLR on bootstrap sample, passing FE and precomputed XZ_b
      fit_b <- estimateCME_PLR(
        data                 = data_b,
        Y                    = Y,
        D                    = D,
        X                    = X,
        Z                    = NULL,
        FE                   = FE,
        XZ_design            = XZ_b,
        outcome_model_type   = outcome_model_type,
        treatment_model_type = treatment_model_type,
        basis_type           = basis_type,
        include_interactions = include_interactions,
        poly_degree          = poly_degree,
        spline_df            = spline_df,
        spline_degree        = spline_degree,
        lambda_cv            = if (need_tune) lambda_cv else NULL,
        lambda_seq           = lambda_seq,
        reduce.dimension     = reduce.dimension,
        bw                   = bw,
        x.eval               = x.eval,
        verbose              = FALSE
      )

      # Return the CME curve for this bootstrap as a single row
      c(fit_b$cme_df$CME_fit)
    }

    # Stop parallel cluster
    parallel::stopCluster(cl)
    if (verbose) cat("  -> Bootstrap loop complete.\n")

    # Fill fit_mat_bs with bootstrap results
    for (b in seq_len(B)) {
      fit_mat_bs[b, ] <- res_list[b, 1:nEval]
    }

    ###########################################################################
    # 5. Pointwise normal-based intervals
    ###########################################################################
    if (verbose) cat("BootstrapPLR Step 5: Computing pointwise SEs and CIs...\n")
    fit_se    <- apply(fit_mat_bs, 2, sd, na.rm = TRUE)
    zcrit     <- qnorm(1 - alpha/2)

    alpha_lower <- alpha / 2
    alpha_upper <- 1 - alpha_lower

    fit_ci_l <- apply(fit_mat_bs, 2, quantile, probs = alpha_lower, na.rm = TRUE)
    fit_ci_u <- apply(fit_mat_bs, 2, quantile, probs = alpha_upper, na.rm = TRUE)

    ###########################################################################
    # 6. Uniform confidence intervals (if desired)
    ###########################################################################
    if (verbose) cat("BootstrapPLR Step 6: Computing uniform CIs...\n")
    theta_matrix <- t(fit_mat_bs)  # nEval x B

    uni_res      <- calculate_uniform_quantiles(theta_matrix, alpha)
    Q_j          <- uni_res$Q_j      # nEval x 2
    coverage     <- uni_res$coverage

    fit_ci_l_uni <- Q_j[, 1]
    fit_ci_u_uni <- Q_j[, 2]

    ###########################################################################
    # 7. Assemble final results
    ###########################################################################
    if (verbose) cat("BootstrapPLR Step 7: Assembling output...\n")
    fit_df <- data.frame(
      X                = x.eval,
      CME              = cme_out,
      SE               = fit_se,
      CI.lower         = fit_ci_l,
      CI.upper         = fit_ci_u,
      CI.lower.uniform = fit_ci_l_uni,
      CI.upper.uniform = fit_ci_u_uni
    )

    ###########################################################################
    # 8. Return list
    ###########################################################################
    return(list(
      results         = fit_df,
      selected_covars = selected_covars,
      fit_full        = fit_full,
      coverage        = coverage
    ))
  }
  else{
    fit_df <- data.frame(
      X                = x.eval,
      CME              = cme_out
    )
    return(list(
      results         = fit_df,
      selected_covars = selected_covars,
      fit_full        = fit_full
    ))
  }

}



