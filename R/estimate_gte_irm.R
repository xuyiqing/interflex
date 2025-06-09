# Revised functions for discrete moderator X
# Requires packages: splines, glmnet

# Function to estimate group treatment effect (GTE) when X is discrete
estimateGTE <- function(
  data,
  Y,           # outcome variable name
  D,           # binary treatment variable name (0/1)
  X,           # discrete moderator variable name
  Z = NULL,    # vector of additional covariate names
  FE = NULL,   # vector of fixed-effect variable names
  estimand = c("ATE", "ATT"),
  signal = c("outcome", "ipw", "aipw"),
  basis_type = c("polynomial", "bspline", "none"),
  include_interactions = TRUE,
  poly_degree = 2,
  spline_df = 4,
  spline_degree = 2,
  XZ_design = NULL,    # optional prebuilt design matrix
  outcome_model_type = "lasso",  # "linear", "ridge", or "lasso"
  ps_model_type      = "lasso",  # "linear", "ridge", or "lasso"
  lambda_cv = NULL,
  lambda_seq = NULL,
  verbose = TRUE
) {
  # 1. Argument matching and input checks
  estimand   <- match.arg(estimand)
  signal     <- match.arg(signal)
  basis_type <- match.arg(basis_type)
  stopifnot(is.data.frame(data))
  if (!(Y %in% names(data))) stop("Y not found in data.")
  if (!(D %in% names(data))) stop("D not found in data.")
  if (!(X %in% names(data))) stop("X not found in data.")
  Yvec <- data[[Y]]; Dvec <- data[[D]]; Xvec <- data[[X]]
  if (!all(Dvec %in% c(0,1))) stop("D must be binary 0/1.")

  # 2. Prepare Z (non-fixed effects) and FE dummies
  Z_nonFE <- if (!is.null(Z)) {
    if (!all(Z %in% names(data))) stop("Some Z not in data.")
    data[, Z, drop=FALSE]
  } else NULL
  FE_dummies <- NULL
  if (!is.null(FE)) {
    if (!all(FE %in% names(data))) stop("Some FE not in data.")
    FE_mat <- data[, FE, drop=FALSE]
    FE_dummies <- do.call(cbind, lapply(names(FE_mat), function(var) {
      fac <- factor(FE_mat[[var]])
      mm <- model.matrix(~ fac)[, -1, drop=FALSE]
      colnames(mm) <- paste0(var, "_", levels(fac)[-1])
      mm
    }))
    if (verbose) cat("Fixed effects converted to dummies.\n")
  }

  # 3. Build design matrix if not supplied
  if (is.null(XZ_design)) {
    build_design_matrix <- function(X_vec, Z_nonFE, FE_dummies) {
      # Encode discrete X as factor dummies
      fac <- factor(X_vec)
      X_mat <- model.matrix(~ fac)[, -1, drop=FALSE]
      colnames(X_mat) <- paste0(deparse(substitute(X)), "_", levels(fac)[-1])

      # Expand Z if needed
      Z_expanded <- NULL
      if (!is.null(Z_nonFE) && ncol(Z_nonFE)>0) {
        Z_expanded <- do.call(cbind, lapply(seq_len(ncol(Z_nonFE)), function(j) {
          vec <- Z_nonFE[[j]]; name <- colnames(Z_nonFE)[j]
          if (basis_type=="none") {
            mat <- matrix(vec, ncol=1)
            colnames(mat) <- name; mat
          } else if (basis_type=="polynomial") {
            pm <- poly(vec, degree=poly_degree, raw=TRUE)
            colnames(pm) <- paste0(name, "_poly", seq_len(ncol(pm))); pm
          } else {
            bspl <- splines::bs(vec, df=spline_df, degree=spline_degree)
            colnames(bspl) <- paste0(name, "_bs", seq_len(ncol(bspl))); bspl
          }
        }))
      }
      # Combine Z_expanded and FE
      Z_block <- if (!is.null(FE_dummies)) {
        if (is.null(Z_expanded)) FE_dummies else cbind(Z_expanded, FE_dummies)
      } else Z_expanded

      # Interaction terms between X and Z_expanded
      int_XZ <- NULL
      if (include_interactions && !is.null(Z_expanded)) {
        int_XZ <- do.call(cbind, lapply(colnames(X_mat), function(xn) {
          X_mat[, xn] * Z_expanded
        }))
        colnames(int_XZ) <- as.vector(outer(colnames(X_mat), colnames(Z_expanded), paste, sep="."))
      }

      cbind(X_mat, Z_block, int_XZ)
    }
    XZ_design <- build_design_matrix(Xvec, Z_nonFE, FE_dummies)
    if (verbose) cat("Design matrix built.\n")
  }

  # 4. Fit outcome and propensity score models
  do_outcome <- signal %in% c("outcome","aipw")
  do_ps      <- signal %in% c("ipw","aipw")
  out_fit1 <- out_fit0 <- ps_fit <- NULL
  mu1_hat <- mu0_hat <- p_hat <- NULL

  # Outcome
  if (do_outcome) {
    idx1 <- which(Dvec==1); idx0 <- which(Dvec==0)
    fit_model <- function(y, Xmat, type, lambda, pf) {
      if (type=="linear") {
        df <- data.frame(y=y, Xmat); lm(y~., data=df)
      } else {
        alpha <- ifelse(type=="lasso",1,0)
        cv <- glmnet::cv.glmnet(Xmat, y, alpha=alpha, lambda=lambda_seq, penalty.factor=pf)
        cv
      }
    }
    pf_out <- rep(1, ncol(XZ_design))
    if (!is.null(FE_dummies)) {
      fe_idx <- grep(paste0("^(", paste(FE, collapse="|"), ")_"), colnames(XZ_design))
      pf_out[fe_idx] <- 0
    }
    X1 <- XZ_design[idx1,,drop=FALSE]; X0 <- XZ_design[idx0,,drop=FALSE]
    fit1 <- fit_model(Yvec[idx1], X1, outcome_model_type, lambda_cv$outcome1, pf_out)
    fit0 <- fit_model(Yvec[idx0], X0, outcome_model_type, lambda_cv$outcome0, pf_out)
    pred_out <- function(fit, Xmat) {
      if (inherits(fit, "cv.glmnet")) as.numeric(predict(fit, newx=Xmat, s="lambda.min")) else as.numeric(predict(fit, newdata=as.data.frame(Xmat)))
    }
    mu1_hat <- pred_out(fit1, XZ_design)
    mu0_hat <- pred_out(fit0, XZ_design)
    out_fit1 <- fit1; out_fit0 <- fit0
  }

  # Propensity score
  if (do_ps) {
    X_ps <- XZ_design
    if (!is.null(FE_dummies)) X_ps <- X_ps[, -grep(paste0("^", paste(FE,collapse="|"),"_"), colnames(XZ_design)), drop=FALSE]
    fit_ps <- function(d, Xmat, type) {
      if (type=="linear") glm(d~., data=data.frame(d=d,Xmat), family=binomial())
      else glmnet::cv.glmnet(Xmat, d, alpha=ifelse(type=="lasso",1,0), family="binomial", lambda=lambda_seq)
    }
    fitp <- fit_ps(Dvec, X_ps, ps_model_type)
    pred_ps <- function(fit, Xmat) {
      if (inherits(fit, "cv.glmnet")) as.numeric(predict(fit, newx=Xmat, s="lambda.min", type="response"))
      else as.numeric(predict(fit, newdata=as.data.frame(Xmat), type="response"))
    }
    p_hat <- pred_ps(fitp, X_ps)
    p_hat <- pmin(pmax(p_hat, 1e-2), 1-1e-2)
    ps_fit <- fitp
  }

  # 5. Compute group-level p_hat_gam for ATT
  if (estimand=="ATT") {
    fac <- factor(Xvec)
    p_hat_gam_group <- tapply(Dvec, fac, mean)
    p_hat_gam <- as.numeric(p_hat_gam_group[as.character(fac)])
  }

  # 6. Compute individual signals based on estimand
  if (estimand=="ATE") {
    est_signal <- switch(signal,
      outcome = mu1_hat - mu0_hat,
      ipw     = Dvec*Yvec/p_hat - (1-Dvec)*Yvec/(1-p_hat),
      aipw    = (mu1_hat - mu0_hat) + Dvec*(Yvec-mu1_hat)/p_hat - (1-Dvec)*(Yvec-mu0_hat)/(1-p_hat)
    )
  } else { # ATT
    est_signal <- switch(signal,
      outcome = (Yvec - mu0_hat) * Dvec / p_hat_gam,
      ipw     = Yvec * (Dvec - p_hat) / ((1 - p_hat) * p_hat_gam),
      aipw    = ((Yvec - mu0_hat) * (Dvec - (1-Dvec)*p_hat/(1-p_hat))) / p_hat_gam
    )
  }

  # 7. Group average by X
  fac <- factor(Xvec)
  levels_X <- levels(fac)
  gte_vals <- tapply(est_signal, fac, mean)
  gte_df <- data.frame(X = levels_X, GTE = as.numeric(gte_vals), row.names = NULL)

  # 8. Return
  list(
    XZ_design = XZ_design,
    mu1_hat   = mu1_hat,
    mu0_hat   = mu0_hat,
    p_hat     = p_hat,
    est_signal= est_signal,
    gte_df    = gte_df,
    signal    = signal,
    estimand  = estimand
  )
}

# Bootstrap function for GTE (works for both ATE and ATT)
bootstrapGTE <- function(
  data, Y, D, X,
  Z = NULL, FE = NULL,
  estimand = c("ATE","ATT"),
  signal   = c("outcome","ipw","aipw"),
  B        = 200, alpha = 0.05,
  basis_type = c("polynomial","bspline","none"),
  include_interactions = TRUE,
  poly_degree          = 2,
  spline_df            = 4,
  spline_degree        = 2,
  outcome_model_type   = "lasso",
  ps_model_type        = "lasso",
  lambda_seq           = NULL,
  verbose              = TRUE
) {
  estimand   <- match.arg(estimand)
  signal     <- match.arg(signal)
  basis_type <- match.arg(basis_type)

  if(verbose) message("BootstrapGTE Step 1: Fitting full-sample GTE...")
  full <- estimateGTE(
    data, Y, D, X, Z, FE,
    estimand, signal,
    basis_type, include_interactions,
    poly_degree, spline_df, spline_degree,
    NULL, outcome_model_type, ps_model_type,
    NULL, lambda_seq, verbose
  )
  if(verbose) message("  → Done.")

  Xlev    <- full$gte_df$X
  gte_full<- full$gte_df$GTE
  n       <- nrow(data)
  K       <- length(Xlev)
  XZ_full <- full$XZ_design

  if(verbose) message("BootstrapGTE Step 2: Initializing cluster and storage...")
  gte_mat <- matrix(NA_real_, nrow=B, ncol=K)
  idx_seq <- seq_len(n)
  if(!requireNamespace("doParallel",quietly=TRUE)) stop("doParallel needed")
  cl <- parallel::makeCluster(parallel::detectCores())
  doParallel::registerDoParallel(cl)
  `%dopar%` <- foreach::`%dopar%`

  if(verbose) message("BootstrapGTE Step 3: Starting bootstraps...")
  res_mat <- foreach::foreach(
    b = seq_len(B),
    .combine  = "rbind",
    .export   = "estimateGTE",
    .packages = c("splines","glmnet")
  ) %dopar% {
    set.seed(1000 + b)
    idx_b <- sample(idx_seq, size=n, replace=TRUE)
    db    <- data[idx_b, , drop=FALSE]
    XZb   <- XZ_full[idx_b, , drop=FALSE]
    fb    <- estimateGTE(
      db, Y, D, X, NULL, FE,
      estimand, signal,
      basis_type, FALSE,
      poly_degree, spline_df, spline_degree,
      XZb, outcome_model_type, ps_model_type,
      NULL, lambda_seq, FALSE
    )
    fb$gte_df$GTE
  }
  parallel::stopCluster(cl)
  if(verbose) message("  → Bootstraps complete.")

  gte_mat[,] <- res_mat
  if(verbose) message("BootstrapGTE Step 4: Computing SE and CIs...")
  se   <- apply(gte_mat, 2, sd, na.rm=TRUE)
  ci_l <- apply(gte_mat, 2, quantile, probs=alpha/2, na.rm=TRUE)
  ci_u <- apply(gte_mat, 2, quantile, probs=1-alpha/2, na.rm=TRUE)

  # Uniform confidence bands
  theta_mat    <- t(gte_mat)
  uni_res      <- calculate_uniform_quantiles(theta_mat, alpha)
  ci_l_uni     <- uni_res$Q_j[,1]
  ci_u_uni     <- uni_res$Q_j[,2]
  coverage_uni <- uni_res$coverage

  results <- data.frame(
    X                = Xlev,
    GTE              = gte_full,
    SE               = se,
    CI.lower         = ci_l,
    CI.upper         = ci_u,
    CI.lower.uniform = ci_l_uni,
    CI.upper.uniform = ci_u_uni,
    stringsAsFactors = FALSE
  )

  if(verbose) message("BootstrapGTE Step 5: Returning results.")
  list(
    results      = results,
    boot_results = gte_mat,
    coverage     = coverage_uni,
    fit_full     = full
  )
}