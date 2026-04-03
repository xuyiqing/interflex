#' @title Estimate Group Average Treatment Effects (GATE) for a PLRM
#' @description Implements the PO-Lasso estimator for a Partially Linear Regression
#'   Model (PLRM) with a continuous treatment and a discrete moderator. The function
#'   estimates the Group Average Treatment Effect (GATE) for each level of the
#'   moderator.
#'
#' @param data A data.frame containing the variables for the analysis.
#' @param Y Character string, the name of the outcome variable.
#' @param D Character string, the name of the continuous treatment variable.
#' @param X Character string, the name of the discrete moderator variable.
#' @param Z Character vector, optional names of additional continuous covariates.
#' @param FE Character vector, optional names of fixed effect variables.
#' @param XZ_design A pre-built design matrix. If NULL, it is created internally.
#' @param outcome_model_type Character, the model for the outcome nuisance function:
#'   "linear", "ridge", or "lasso".
#' @param treatment_model_type Character, the model for the treatment nuisance function:
#'   "linear", "ridge", or "lasso".
#' @param basis_type Character, the basis expansion for continuous covariates in Z:
#'   "polynomial", "bspline", or "none".
#' @param include_interactions Logical, whether to include interaction terms.
#' @param poly_degree Integer, the degree for polynomial basis expansion.
#' @param spline_df Integer, degrees of freedom for B-spline basis expansion.
#' @param spline_degree Integer, the degree for B-spline basis expansion.
#' @param lambda_cv A named list with pre-specified lambda values (e.g., `list(outcome=..., treatment=...)`)
#'   to be used instead of cross-validation. Typically used during bootstrapping.
#' @param lambda_seq A custom numeric vector of lambda values for the glmnet grid search.
#' @param verbose Logical, whether to print progress messages.
#'
#' @return A list containing detailed results, including the final GATE data.frame,
#'   the design matrix, fitted model objects, selected covariates, and lambdas used.
#'
estimateGATE_PLR <- function(
  data,
  Y,                      # outcome variable name
  D,                      # continuous treatment name
  X,                      # discrete moderator name
  Z = NULL,               # additional covariates
  FE = NULL,              # fixed-effect vars
  XZ_design = NULL,     # optional prebuilt design matrix

  outcome_model_type   = "lasso",  # "linear", "ridge", "lasso"
  treatment_model_type = "lasso",  # "linear", "ridge", "lasso"
  basis_type           = c("polynomial","bspline","none"),
  include_interactions = FALSE,
  poly_degree          = 2,
  spline_df            = 5,
  spline_degree        = 3,

  lambda_cv  = NULL,  # list(outcome=..., treatment=...)
  lambda_seq = NULL,  # lambda grid for glmnet

  verbose    = TRUE
) {
  # 0) prerequisites & lambda_cv check 
  if (!requireNamespace("glmnet", quietly=TRUE)) {
    stop("Package 'glmnet' is required.")
  }
  if (!requireNamespace("splines", quietly=TRUE)) {
    stop("Package 'splines' is required.")
  }
  basis_type <- match.arg(basis_type)
  if (!is.null(lambda_cv)) {
    allowed <- c("outcome","treatment")
    extra   <- setdiff(names(lambda_cv), allowed)
    if (length(extra) > 0) {
      stop("lambda_cv contains invalid entry(ies): ",
           paste(extra, collapse=", "))
    }
  }
  stopifnot(is.data.frame(data),
            all(c(Y,D,X) %in% names(data)))
  Yv <- data[[Y]]
  Dv <- data[[D]]
  Xv <- data[[X]]

  # 1) build or validate XZ_design 
  if (!is.null(XZ_design)) {
    if (!is.matrix(XZ_design) || nrow(XZ_design)!=nrow(data)) {
      stop("Provided XZ_design must be a matrix with same #rows as data.")
    }
  } else {
    # 1A) prepare Z
    Zm <- if (!is.null(Z)) {
      stopifnot(all(Z %in% names(data)))
      data[, Z, drop=FALSE]
    } else NULL
    # 1B) factor-dummy for X
    fac  <- factor(Xv)
    Xd   <- model.matrix(~ fac)[,-1,drop=FALSE]
    colnames(Xd) <- paste0(X, "_", levels(fac)[-1])

    # 1C) expand Z
    Ze <- NULL
    if (!is.null(Zm)) {
      Ze_list <- lapply(seq_len(ncol(Zm)), function(j) {
        vec  <- Zm[[j]]
        name <- colnames(Zm)[j]
        if (basis_type=="none") {
          m <- matrix(vec, ncol=1); colnames(m)<-name; m
        } else if (basis_type=="polynomial") {
          p <- poly(vec, degree=poly_degree, raw=TRUE)
          colnames(p) <- paste0(name,"_poly", seq_len(ncol(p))); p
        } else {
          b <- splines::bs(vec, df=spline_df, degree=spline_degree)
          colnames(b) <- paste0(name,"_bs", seq_len(ncol(b))); b
        }
      })
      Ze <- do.call(cbind, Ze_list)
    }

    # 1D) FE dummies
    FEd <- NULL
    if (!is.null(FE)) {
      stopifnot(all(FE %in% names(data)))
      FEd_list <- lapply(FE, function(v) {
        f <- factor(data[[v]])
        m <- model.matrix(~ f)[,-1,drop=FALSE]
        colnames(m) <- paste0(v, "_", levels(f)[-1])
        m
      })
      FEd <- do.call(cbind, FEd_list)
    }

    # 1E) X:Z interactions
    IntXZ <- NULL
    if (include_interactions && !is.null(Ze)) {
      IntXZ <- do.call(cbind, lapply(colnames(Xd), function(xn) {
        sweep(Ze, 1, Xd[,xn], "*")
      }))
      colnames(IntXZ) <- as.vector(outer(
        colnames(Xd), colnames(Ze),
        FUN = function(a,b) paste(a,b,sep=".")
      ))
    }

    # 1F) Z:Z interactions
    IntZZ <- NULL
    if (include_interactions && !is.null(Ze) && ncol(Ze)>1) {
      combs   <- combn(colnames(Ze), 2, simplify=FALSE)
      zz_list <- lapply(combs, function(p) Ze[,p[1]] * Ze[,p[2]])
      IntZZ   <- do.call(cbind, zz_list)
      colnames(IntZZ) <- vapply(combs, paste, collapse=".", FUN.VALUE="")
    }

    XZ_design <- cbind(Xd, Ze, FEd, IntXZ, IntZZ)
    if (verbose) cat("Design matrix built with", ncol(XZ_design), "columns.\n")
  }

  # 2) penalty factors for FE 
  p        <- ncol(XZ_design)
  pf_out   <- pf_tr <- rep(1, p)
  if (!is.null(FE)) {
    pat     <- paste0("^(", paste(FE, collapse="|"), ")_")
    idx_fe  <- grep(pat, colnames(XZ_design))
    pf_out[idx_fe] <- 0
    pf_tr[idx_fe]  <- 0
  }

  # 3) glmnet helper 
  fit_glmnet <- function(y, Xm, type, lam, pf) {
    if (type == "linear") {
      df <- data.frame(y=y, Xm)
      colnames(df)[-1] <- colnames(Xm)
      lm0 <- lm(y ~ ., data=df)
      return(list(
        fit    = lm0,
        type   = "lm",
        lambda = NULL,
        active = seq_len(ncol(Xm))
      ))
    }
    alpha <- if (type=="lasso") 1 else 0
    if (!is.null(lam)) {
      fit0 <- glmnet::glmnet(Xm, y, alpha=alpha, lambda=lam, penalty.factor=pf)
      lam0 <- lam
    } else {
      cv0  <- glmnet::cv.glmnet(Xm, y, alpha=alpha,
                                lambda=lambda_seq,
                                penalty.factor=pf)
      fit0 <- cv0
      lam0 <- cv0$lambda.min
    }
    co    <- as.numeric(coef(fit0, s=lam0))[-1]
    act   <- which(co != 0)
    list(
      fit    = fit0,
      type   = "glmnet",
      lambda = lam0,
      active = act
    )
  }

  # 4) initial penalized fits 
  if (verbose) cat("Fitting outcome & treatment models...\n")
  out0 <- fit_glmnet(Yv, XZ_design, outcome_model_type,
                     lambda_cv$outcome, pf_out)
  tr0  <- fit_glmnet(Dv, XZ_design, treatment_model_type,
                     lambda_cv$treatment, pf_tr)

  # 5) collect lambda & active sets 
  lambda_used <- list(outcome=out0$lambda, treatment=tr0$lambda)
  selected    <- list(outcome=out0$active, treatment=tr0$active)

  # 6) post-selection refits if any penalized 
  if (outcome_model_type != "linear") {
    ao <- out0$active
    if (length(ao) > 0) {
      sub_out <- XZ_design[, ao, drop = FALSE]
      df_out  <- data.frame(y = Yv, sub_out)
      out0$fit  <- lm(y ~ ., data = df_out)
      out0$type <- "lm"
    }
  }

  if (treatment_model_type != "linear") {
    at <- tr0$active
    if (length(at) > 0) {
      sub_tr <- XZ_design[, at, drop = FALSE]
      df_tr  <- data.frame(d = Dv, sub_tr)
      tr0$fit  <- lm(d ~ ., data = df_tr)
      tr0$type  <- "lm"
    }
  }

  # 7) residualize 
  predict_helper <- function(obj, Xm) {
    if (obj$type=="lm") {
      predict(obj$fit, newdata=as.data.frame(Xm))
    } else {
      as.numeric(predict(obj$fit, newx=Xm, s="lambda.min"))
    }
  }
  yhat   <- predict_helper(out0, XZ_design)
  dhat   <- predict_helper(tr0, XZ_design)
  ytilde <- Yv - yhat
  dtilde <- Dv - dhat

  # 8) group-by-X regression for GATE 
  if (verbose) cat("Computing GATE by group...\n")
  fac   <- factor(Xv)
  lvls  <- levels(fac)
  gate  <- sapply(lvls, function(lv) {
    idx <- which(fac == lv)
    coef(lm(ytilde[idx] ~ dtilde[idx]))[2]
  })
  gate_df <- data.frame(X=lvls, GATE=gate, stringsAsFactors=FALSE)

  # 9) return everything 
  ret <- list(
    XZ_design     = XZ_design,
    gate_df       = gate_df,
    selected      = selected,      # active sets before refit
    lambda_used   = lambda_used,   # penalty params
    outcome_fit   = out0$fit,      # final refit
    treatment_fit = tr0$fit        # final refit
  )
  return(ret)
}

#' @title Bootstrap Confidence Intervals for GATE from a PLRM
#' @description Performs a nonparametric bootstrap to compute pointwise and uniform
#'   confidence intervals for the Group Average Treatment Effect (GATE) across
#'   levels of a discrete moderator in a PLRM.
#'
#' @param data A data.frame containing the variables for the analysis.
#' @param Y Character string, the name of the outcome variable.
#' @param D Character string, the name of the continuous treatment variable.
#' @param X Character string, the name of the discrete moderator variable.
#' @param Z Character vector, optional names of additional continuous covariates.
#' @param FE Character vector, optional names of fixed effect variables.
#' @param B Integer, the number of bootstrap replications.
#' @param alpha Numeric, the significance level for confidence intervals (e.g., 0.05 for 95% CIs).
#' @param ... Additional arguments passed down to `estimateGATE_PLR`.
#'
#' @return A list containing the final results data.frame with point estimates and CIs,
#'   the matrix of bootstrap replications, uniform coverage level, and the full-sample
#'   estimation object.
#'
bootstrapGATE_PLR <- function(
  data, Y, D, X,
  Z = NULL, FE = NULL,
  B = 200, alpha = 0.05,
  outcome_model_type   = "lasso",
  treatment_model_type = "lasso",
  basis_type           = c("polynomial","bspline","none"),
  include_interactions = FALSE,
  poly_degree          = 2,
  spline_df            = 4,
  spline_degree        = 2,
  lambda_seq           = NULL,
  CI = TRUE,
  cores = 8,
  verbose              = TRUE
) {
  # 1) Prep & full-sample fit ------------------------------------------------
  basis_type <- match.arg(basis_type)
  if (verbose) message("BootstrapGATE_PLR: fitting full-sample PLR...")
  full <- estimateGATE_PLR(
    data                 = data,
    Y                    = Y,
    D                    = D,
    X                    = X,
    Z                    = Z,
    FE                   = FE,
    XZ_design            = NULL,
    outcome_model_type   = outcome_model_type,
    treatment_model_type = treatment_model_type,
    basis_type           = basis_type,
    include_interactions = include_interactions,
    poly_degree          = poly_degree,
    spline_df            = spline_df,
    spline_degree        = spline_degree,
    lambda_cv            = NULL,      # let full-sample CV select lambda's
    lambda_seq           = lambda_seq,
    verbose              = verbose
  )

  lvls        <- full$gate_df$X
  g_full      <- full$gate_df$GATE
  n           <- nrow(data)
  XZ0         <- full$XZ_design
  lambda_used <- full$lambda_used    # list(outcome=..., treatment=...)

  if(isTRUE(CI)){
    # 2) Parallel bootstrap ----------------------------------------------------
    if (verbose) message("BootstrapGATE_PLR: launching cluster...")
    if (!requireNamespace("doFuture", quietly=TRUE)) {
      stop("Package 'doFuture' required for parallel bootstrap.")
    }
    doFuture::registerDoFuture()
    future::plan(future::multisession, workers = cores)
    on.exit(future::plan(future::sequential), add = TRUE)

    if (verbose) message("BootstrapGATE_PLR: running ", B, " bootstrap draws...")
    res_mat <- foreach::foreach(
      b = seq_len(B),
      .combine  = "rbind",
      .export   = "estimateGATE_PLR",
      .packages = c("glmnet","splines"),
      .options.future = list(seed = TRUE)
    ) %dopar% {
      idx <- sample(n, n, replace = TRUE)
      db  <- data[idx, , drop = FALSE]
      XZb <- XZ0[idx, , drop = FALSE]

      fb <- estimateGATE_PLR(
        data                 = db,
        Y                    = Y,
        D                    = D,
        X                    = X,
        Z                    = NULL,        # design already in XZb
        FE                   = FE,
        XZ_design            = XZb,
        outcome_model_type   = outcome_model_type,
        treatment_model_type = treatment_model_type,
        basis_type           = basis_type,
        include_interactions = include_interactions,
        poly_degree          = poly_degree,
        spline_df            = spline_df,
        spline_degree        = spline_degree,
        lambda_cv            = lambda_used,  # reuse full-sample lambda's
        lambda_seq           = lambda_seq,
        verbose              = FALSE
      )

      # extract GATEs in same X-order
      sapply(lvls, function(lv) {
        pos <- which(fb$gate_df$X == lv)
        if (length(pos)==1) fb$gate_df$GATE[pos] else NA_real_
      })
    }

    # 3) Compute SE & CIs -------------------------------------------------------
    if (verbose) message("BootstrapGATE_PLR: summarizing draws...")
    se    <- apply(res_mat, 2, sd, na.rm=TRUE)
    cil   <- apply(res_mat, 2, quantile, probs=alpha/2,   na.rm=TRUE)
    cih   <- apply(res_mat, 2, quantile, probs=1-alpha/2, na.rm=TRUE)
    
    uni   <- calculate_uniform_quantiles(t(res_mat), alpha)


    results <- data.frame(
      X                = as.numeric(lvls),
      GATE             = g_full,
      SE               = se,
      CI.lower         = cil,
      CI.upper         = cih,
      CI.lower.uniform = uni$Q_j[,1],
      CI.upper.uniform = uni$Q_j[,2],
      stringsAsFactors = FALSE
    )

    if (verbose) message("BootstrapGATE_PLR: done.")
    list(
      results      = results,
      boot_results = res_mat,
      coverage     = uni$coverage,
      fit_full     = full
    )    
  }
  else{
    results <- data.frame(
      X                = as.numeric(lvls),
      GATE             = g_full,
      stringsAsFactors = FALSE
    )
    list(
      results      = results,
      fit_full     = full
    )     
  }


}


