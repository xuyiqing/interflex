estimateCME_PLR <- function(data,
                            Y,          # name of outcome variable in 'data'
                            D,          # name of treatment variable in 'data'
                            X,          # name of focal variable in 'data'
                            Z = NULL,   # character vector of additional covariates
                            # (NEW) Optionally supply a pre-built design matrix for (X,Z):
                            XZ_design = NULL,  
                            # If NULL, we'll build the design matrix internally. 
                            # If not NULL, we skip internal basis expansions and use it directly.
                            
                            outcome_model_type = "lasso", # can be "linear", "ridge", or "lasso"
                            treatment_model_type      = "lasso", # can be "linear", "ridge", or "lasso"
                            
                            # polynomial or B-spline expansion in the design matrix
                            basis_type         = c("polynomial", "bspline", "none"),
                            include_interactions = TRUE,
                            poly_degree        = 2,    # only used if basis_type="polynomial"
                            
                            # B-spline parameters (used if basis_type="bspline", or for final CME fit)
                            spline_df          = 4,
                            spline_degree      = 2,
                            
                            lambda_seq         = NULL, # optional custom lambda sequence for glmnet
                            reduce.dimension = c("bspline","kernel"), 
                            bw = NULL,
                            x.eval             = NULL, # grid of X values for final CME curve
                            selected_covars = NULL,
                            verbose            = TRUE){
  
  basis_type <- match.arg(basis_type)
  reduce.dimension <- match.arg(reduce.dimension)
  ##############################################################################
  # 0. Extract variables & basic checks
  ##############################################################################
  stopifnot(is.data.frame(data))
  if (!is.character(Y) || !is.character(D) || !is.character(X)) {
    stop("Y, D, X must be character strings corresponding to variable names in 'data'.")
  }
  if (!(Y %in% names(data))) stop("Variable name for Y not found in data.")
  if (!(D %in% names(data))) stop("Variable name for D not found in data.")
  if (!(X %in% names(data))) stop("Variable name for X not found in data.")
  
  Yvec <- data[[Y]]
  Dvec <- data[[D]]
  Xvec <- data[[X]]

  if (!is.null(Z)) {
    # check existence
    if (!all(Z %in% names(data))) {
      stop("Some variables in Z not found in data.")
    }
    Zmat <- data[, Z, drop = FALSE]
  } else {
    Zmat <- NULL
  }
  
  n <- nrow(data)
  if (length(Yvec) != n || length(Dvec) != n || length(Xvec) != n) {
    stop("Lengths of Y, D, X do not match.")
  }
  # Default grid for X if none provided
  if (is.null(x.eval)) {
    x.eval <- seq(min(Xvec), max(Xvec), length.out = 100)
  }
  
  if (!is.null(XZ_design)) {
    # user has supplied a design matrix
    if (!is.matrix(XZ_design)) {
      stop("'XZ_design' must be a matrix (or coerceable to matrix).")
    }
    if (nrow(XZ_design) != n) {
      stop("XZ_design must have the same number of rows as 'data'.")
    }
  } else {
    # Build the design matrix internally
    # We only need 'splines' if basis_type == "bspline", but we check anyway:
    if (!requireNamespace("splines", quietly = TRUE)) {
      stop("Package 'splines' is required (for B-spline expansions).")
    }
    
    build_design_matrix <- function(X_vec, Z_mat) {
      
      # Helper for B-spline expansions
      build_bspline_cols <- function(vec, varname, df, degree) {
        bs_mat <- bs(vec, df = df, degree = degree)
        colnames(bs_mat) <- paste0(varname, "_bs", seq_len(ncol(bs_mat)))
        bs_mat
      }
      
      #-----------------------------
      # 1) Expand X
      #-----------------------------
      if (basis_type == "none") {
        # Keep X as-is, single column named "X"
        X_mat <- matrix(X_vec, ncol=1)
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
        
      } else {
        # basis_type == "bspline"
        X_mat <- build_bspline_cols(X_vec, "X", spline_df, spline_degree)
      }
      
      #-----------------------------
      # 2) Expand each Z
      #-----------------------------
      Z_expanded <- NULL
      if (!is.null(Z_mat) && ncol(Z_mat) > 0) {
        Z_expanded_list <- lapply(seq_len(ncol(Z_mat)), function(j) {
          zj <- Z_mat[, j]
          zname <- colnames(Z_mat)[j]
          
          if (basis_type == "none") {
            # Keep Zj as-is
            mat_z <- matrix(zj, ncol=1)
            colnames(mat_z) <- zname
            mat_z
            
          } else if (basis_type == "polynomial") {
            if (poly_degree > 1) {
              zj_poly <- stats::poly(zj, degree = poly_degree, raw = TRUE)
              colnames(zj_poly) <- paste0(zname, "_poly", seq_len(poly_degree))
              zj_poly
            } else {
              matrix(zj, ncol=1, dimnames=list(NULL, zname))
            }
            
          } else {
            # basis_type == "bspline"
            build_bspline_cols(zj, zname, spline_df, spline_degree)
          }
        })
        Z_expanded <- do.call(cbind, Z_expanded_list)
      }
      
      #-----------------------------
      # 3) Interactions
      #-----------------------------
      int_mat <- NULL
      if (include_interactions && !is.null(Z_expanded)) {
        for (i_col in seq_len(ncol(X_mat))) {
          for (j_col in seq_len(ncol(Z_expanded))) {
            tmp <- X_mat[, i_col] * Z_expanded[, j_col]
            int_mat <- cbind(int_mat, tmp)
          }
        }
        if (!is.null(int_mat)) {
          colnames(int_mat) <- apply(
            expand.grid(colnames(X_mat), colnames(Z_expanded)), 1,
            paste, collapse="."
          )
        }
      }
      
      # Combine everything
      design <- cbind(X_mat, Z_expanded, int_mat)
      return(as.matrix(design))
    }
    
    XZ_design <- build_design_matrix(Xvec, Zmat)
  }
  
  both_lasso <- (outcome_model_type == "lasso" && treatment_model_type == "lasso")
  exactly_one_lasso <- (outcome_model_type == "lasso") != (treatment_model_type == "lasso")
  
  if (exactly_one_lasso) {
    stop("Error: If one model uses 'lasso', the other must also be 'lasso'.")
  }
  if (both_lasso && verbose) {
    message("Both outcome_model_type and treatment_model_type are 'lasso': ",
            "post-selection / double-selection LASSO is utilized.")
  }
  
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for penalized regression.")
  }
  library(glmnet)
  
  
  do_single_fit <- function(y_sub, x_sub, model_type) {
    if (model_type == "linear") {
      df_temp <- data.frame(y = y_sub, x_sub)
      fit_lm  <- lm(y ~ ., data = df_temp)
      return(list(fit = fit_lm, type = "lm"))
    } else if (model_type == "ridge") {
      fit_cv <- glmnet::cv.glmnet(x = x_sub, y = y_sub, alpha = 0,
                                  lambda = lambda_seq)
      return(list(fit = fit_cv, type = "glmnet"))
    } else if (model_type == "lasso") {
      fit_cv <- glmnet::cv.glmnet(x = x_sub, y = y_sub, alpha = 1,
                                  lambda = lambda_seq)
      return(list(fit = fit_cv, type = "glmnet"))
    } else {
      stop("Unsupported outcome_model_type: ", model_type)
    }
  }

  
  if(!is.null(selected_covars)){
    act.outcome <- selected_covars[['outcome']]
    act.treatment <- selected_covars[['treatment']]
    XZ_sub_outcome <- XZ_design[, act.outcome, drop=FALSE]
    XZ_sub_treatment <- XZ_design[, act.treatment, drop=FALSE]
    outcome_fit <- do_single_fit(Yvec, XZ_sub_outcome, outcome_model_type)
    treatment_fit <- do_single_fit(Dvec, XZ_sub_treatment, treatment_model_type)
  }
  else{
    outcome_fit <- do_single_fit(Yvec, XZ_design, outcome_model_type)
    treatment_fit <- do_single_fit(Dvec, XZ_design, treatment_model_type)    
  }
  
  selected_covars <- NULL
  
  active_vars <- function(cv_fit, x_mat) {
    # coefs includes intercept => index=1 is intercept
    coefs <- coef(cv_fit$fit, s = cv_fit$fit$lambda.min)
    which_nz <- which(as.numeric(coefs) != 0)
    # skip intercept
    active_idx <- setdiff(which_nz, 1) - 1
    # '-1' because the remaining indices map to columns in x_mat
    return(active_idx)
  }
  
  if (both_lasso) {
    outcome_fit_lasso <- outcome_fit
    treatment_fit_lasso <- treatment_fit
    
    # Gather active sets
    act.outcome <- active_vars(list(fit = outcome_fit$fit), XZ_design)
    act.treatment <- active_vars(list(fit = treatment_fit$fit), XZ_design)
    
    union_active <- unique(c(act.outcome, act.treatment))
    selected_covars <- union_active
    
    # Refit final OLS/logit on union
    selected_covars <- list(outcome = act.outcome,treatment = act.treatment)
    # outcome
    XZ_sub_outcome <- XZ_design[, act.outcome, drop=FALSE]
    XZ_sub_treatment <- XZ_design[, act.treatment, drop=FALSE]
    outcome_sub <- data.frame(y = Yvec, XZ_sub_outcome)
    final_out_lm <- lm(y ~ ., data = outcome_sub)

    treatment_sub <- data.frame(d = Dvec, XZ_sub_treatment)
    final_treatment_lm <- lm(d ~ ., data = treatment_sub)
    
    # Overwrite the lasso fits with final unpenalized fits
    outcome_fit <- list(fit = final_out_lm, type="lm")
    treatment_fit   <- list(fit = final_treatment_lm, type="lm")
  }
  
  
  predict_outcome <- function(subfit, newX) {
    if (subfit$type == "lm") {
      vars_used <- names(subfit$fit$coefficients)[-1]
      df_temp   <- data.frame(newX[, vars_used, drop=FALSE])
      colnames(df_temp) <- vars_used
      preds <- predict(subfit$fit, newdata = df_temp)
      return(as.numeric(preds))
    } else if (subfit$type == "glmnet") {
      preds <- as.numeric(predict(subfit$fit, newx = newX, s="lambda.min"))
      return(preds)
    } else {
      stop("Unknown subfit$type in predict_outcome")
    }
  }
  
  outcome_hat <- predict_outcome(outcome_fit, XZ_design)
  treatment_hat <- predict_outcome(treatment_fit, XZ_design)
  
  outcome_signal <- Yvec - outcome_hat
  treatment_signal <- Dvec - treatment_hat
  
  if(reduce.dimension == 'bspline'){
    fit_spline <- function(Y, D, X_vec) {
      data_fit <- data.frame(
        Y = Y,
        D = D,
        X = X_vec
      )
      lm(
        formula = Y ~ D*bs(X, df = spline_df, degree = spline_degree),
        data    = data_fit
      )
    }
    
    fit_spline_model <- fit_spline(
      Y    = outcome_signal,      # Your tilde{Y}
      D    = treatment_signal,    # Your tilde{D}
      X_vec= Xvec                # The covariate
    )
    
    cme_fit_1 <- predict(
      fit_spline_model,
      newdata = data.frame(
        D = 1,
        X = x.eval
      )
    )
    
    cme_fit_0 <- predict(
      fit_spline_model,
      newdata = data.frame(
        D = 0,
        X = x.eval
      )
    )
    
    cme_fit <- cme_fit_1 - cme_fit_0
  }
  
  if (reduce.dimension == "kernel") {
    data.k <- cbind.data.frame(
      Y    = outcome_signal,      # Your tilde{Y}
      D    = treatment_signal,    # Your tilde{D}
      X = Xvec                # The covariate
    )
    if(is.null(bw)){
      sol.k <- interflex::interflex(estimator = 'kernel',Y='Y',D='D',X='X',
                                    data =data.k, X.eval = x.eval,CV = TRUE,parallel = T,cores = 31)      
    }
    else{
      sol.k <- interflex::interflex(estimator = 'kernel',Y='Y',D='D',X='X',
                                    data =data.k, X.eval = x.eval,CV = FALSE, bw = bw)   
    }

    cme_fit <- sol.k$est.kernel[[1]][,2]
  }
  
  cme_df <- data.frame(
    X.eval      = x.eval,
    CME_fit = cme_fit
  )
  
  if (verbose) {
    message("Estimation complete.")
  }
  
  #selected_covars_names <- NULL
  #if (!is.null(selected_covars)) {
  #  selected_covars_names <- colnames(XZ_design)[selected_covars]
  #}
  
  out_list <- list(
    # The final design matrix used
    XZ_design = XZ_design,
    
    # Fitted outcome models
    outcome_fit   = outcome_fit,
    treatment_fit    = treatment_fit,
    
    # Observed-sample predictions
    outcome_hat   = outcome_hat,
    treatment_hat   = treatment_hat,
    
    # Signals
    outcome_signal = outcome_signal,
    treatment_signal = treatment_signal,
    
    # CME curve on x.eval
    cme_df = cme_df,
    reduce.dimension = reduce.dimension,
    
    # If double-selection was used
    selected_covars = selected_covars
  )
  
  if(reduce.dimension == 'bspline'){
    out_list$spline_fit <- fit_spline_model
  }
  
  if(reduce.dimension == 'kernel'){
    out_list$kernel_fit <- sol.k
    out_list$bw <- sol.k$bw
  }
  
  if (both_lasso) {
    # also store the initial lasso fits (before final OLS/logit)
    out_list$outcome_fit_lasso <- outcome_fit_lasso
    out_list$treatment_fit_lasso <- treatment_fit_lasso
  }
  
  return(out_list)
}

interflex.aipw <- function(data,
                          Y,
                          D,
                          X,
                          treat.info,
                          diff.info,
                          Z = NULL,
                          weights = NULL,
                          B = 200,
                          alpha = 0.05,
                          outcome_model_type = "lasso",
                          treatment_model_type = "lasso",
                          basis_type         = c("polynomial", "bspline", "none"),
                          include_interactions = TRUE,
                          poly_degree        = 2,   
                          spline_df          = 4,
                          spline_degree      = 2,
                          lambda_seq         = NULL,
                          reduce.dimension   = c("bspline","kernel"), 
                          bw                 = NULL,
                          x.eval             = NULL,
                          verbose            = TRUE
                          figure = TRUE,
                          CI = CI,
                          order = NULL,
                          subtitles = NULL,
                          show.subtitles = NULL,
                          Xdistr = "histogram",
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
                          width = 10) {
    covariates <- c(X, Z)
    length.covariates <- length(covariates)

    diff.values.plot <- diff.info[["diff.values.plot"]]
    treat.type <- treat.info[["treat.type"]]
    if (treat.type == "discrete") {
        other.treat <- treat.info[["other.treat"]]
        other.treat.origin <- names(other.treat)
        names(other.treat.origin) <- other.treat
        all.treat <- treat.info[["all.treat"]]
        all.treat.origin <- names(all.treat)
        names(all.treat.origin) <- all.treat
    }
    if (treat.type == "continuous") {
        D.sample <- treat.info[["D.sample"]]
        label.name <- names(D.sample)
        # names(label.name) <- D.sample
    }
    if (TRUE) {
        if (treat.type == "discrete") {
            if (is.null(weights) == TRUE) {
                de <- density(data[, X])
            } else {
                suppressWarnings(de <- density(data[, X], weights = data[, "WEIGHTS"]))
            }

            treat_den <- list()
            for (char in all.treat) {
                if (is.null(weights) == TRUE) {
                    de.tr <- density(data[data[, D] == char, X])
                } else {
                    suppressWarnings(de.tr <- density(data[data[, D] == char, X], weights = data[data[, D] == char, "WEIGHTS"]))
                }
                treat_den[[all.treat.origin[char]]] <- de.tr
            }

            if (is.null(weights) == TRUE) {
                hist.out <- hist(data[, X], breaks = 80, plot = FALSE)
            } else {
                suppressWarnings(hist.out <- hist(data[, X], data[, "WEIGHTS"],
                    breaks = 80, plot = FALSE
                ))
            }
            n.hist <- length(hist.out$mids)

            treat.hist <- list()
            for (char in all.treat) {
                count1 <- rep(0, n.hist)
                treat_index <- which(data[, D] == char)
                for (i in 1:n.hist) {
                    count1[i] <- sum(data[treat_index, X] >= hist.out$breaks[i] &
                        data[treat_index, X] < hist.out$breaks[(i + 1)])
                }
                count1[n.hist] <- sum(data[treat_index, X] >= hist.out$breaks[n.hist] &
                    data[treat_index, X] <= hist.out$breaks[n.hist + 1])

                treat.hist[[all.treat.origin[char]]] <- count1
            }
        }
        if (treat.type == "continuous") { ## continuous D
            if (is.null(weights) == TRUE) {
                de <- density(data[, X])
            } else {
                suppressWarnings(de <- density(data[, X], weights = data[, "WEIGHTS"]))
            }
            if (is.null(weights) == TRUE) {
                hist.out <- hist(data[, X], breaks = 80, plot = FALSE)
            } else {
                suppressWarnings(hist.out <- hist(data[, X], data[, "WEIGHTS"],
                    breaks = 80, plot = FALSE
                ))
            }
            de.tr <- NULL
        }
    }

    treat.base <- treat.info[["base"]]

    TE.output.all.list <- list()
    if (treat.type == "discrete") {
        for (char in other.treat) {
            data_part <- data[data[[D]] %in% c(treat.base, char), ]
            data_part[data_part[[D]] == treat.base, D] <- 0L
            data_part[data_part[[D]] == char, D] <- 1L
            data_part$D <- as.numeric(data_part$D)
            causal.forest <- causal_forest(data_part[covariates], data_part[[Y]], data_part[[D]], num.trees = num.trees)
            X.test <- matrix(0, 50, length.covariates)
            X.test[, 1] <- seq(min(data_part[[X]]), max(data_part[[X]]), length.out = 50)
            causal.forest.hat <- predict(causal.forest, X.test, estimate.variance = TRUE)
            causal.forest.hat.pred <- causal.forest.hat$predictions
            causal.forest.hat.sigma <- sqrt(causal.forest.hat$variance.estimates)

            TE.output.all <- data.frame(
                "X" = X.test[, 1],
                "ME" = causal.forest.hat.pred,
                "sd" = causal.forest.hat.sigma,
                "lower CI(95%)" = causal.forest.hat.pred - 1.96 * causal.forest.hat.sigma,
                "upper CI(95%)" = causal.forest.hat.pred + 1.96 * causal.forest.hat.sigma
            )
            TE.output.all.list[[other.treat.origin[char]]] <- TE.output.all
        }
    }
    if (treat.type == "discrete") {
        final.output <- list(
            diff.info = diff.info,
            treat.info = treat.info,
            est.grf = TE.output.all.list,
            Xlabel = Xlabel,
            Dlabel = Dlabel,
            Ylabel = Ylabel,
            de = de,
            hist.out = hist.out,
            de.tr = treat_den,
            count.tr = treat.hist,
            estimator = "grf"
        )
    }

    # Plot
    if (figure == TRUE) {
        class(final.output) <- "interflex"
        figure.output <- plot.interflex(
            x = final.output,
            order = order,
            subtitles = subtitles,
            show.subtitles = show.subtitles,
            CI = CI,
            diff.values = diff.values.plot,
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
            # bin.labs = bin.labs, # bin labels
            interval = interval, # interval in replicated papers
            file = file,
            ncols = ncols,
            # pool plot
            pool = pool,
            legend.title = legend.title,
            color = color,
            show.all = show.all,
            scale = scale,
            height = height,
            width = width
        )
        final.output <- c(final.output, list(figure = figure.output))
    }

    class(final.output) <- "interflex"
    return(final.output)
}
