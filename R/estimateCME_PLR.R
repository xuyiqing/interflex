estimateCME_PLR <- function(data,
                            Y, # name of outcome variable in 'data'
                            D, # name of treatment variable in 'data'
                            X, # name of focal variable in 'data'
                            Z = NULL, # character vector of additional covariates
                            # (NEW) Optionally supply a pre-built design matrix for (X,Z):
                            XZ_design = NULL,
                            # If NULL, we'll build the design matrix internally.
                            # If not NULL, we skip internal basis expansions and use it directly.
                            outcome_model_type = "lasso", # can be "linear", "ridge", or "lasso"
                            treatment_model_type = "lasso", # can be "linear", "ridge", or "lasso"
                            # polynomial or B-spline expansion in the design matrix
                            basis_type = c("polynomial", "bspline", "none"),
                            include_interactions = TRUE,
                            poly_degree = 2, # only used if basis_type="polynomial"
                            # B-spline parameters (used if basis_type="bspline", or for final CME fit)
                            spline_df = 4,
                            spline_degree = 2,
                            lambda_seq = NULL, # optional custom lambda sequence for glmnet
                            reduce.dimension = c("bspline", "kernel"),
                            bw = NULL,
                            x.eval = NULL, # grid of X values for final CME curve
                            selected_covars = NULL,
                            verbose = TRUE) {
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
                X_mat <- matrix(X_vec, ncol = 1)
                colnames(X_mat) <- "X"
            } else if (basis_type == "polynomial") {
                if (poly_degree > 1) {
                    poly_mat <- stats::poly(X_vec, degree = poly_degree, raw = TRUE)
                    colnames(poly_mat) <- paste0("X_poly", seq_len(poly_degree))
                    X_mat <- poly_mat
                } else {
                    X_mat <- matrix(X_vec, ncol = 1)
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
                        mat_z <- matrix(zj, ncol = 1)
                        colnames(mat_z) <- zname
                        mat_z
                    } else if (basis_type == "polynomial") {
                        if (poly_degree > 1) {
                            zj_poly <- stats::poly(zj, degree = poly_degree, raw = TRUE)
                            colnames(zj_poly) <- paste0(zname, "_poly", seq_len(poly_degree))
                            zj_poly
                        } else {
                            matrix(zj, ncol = 1, dimnames = list(NULL, zname))
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
                        paste,
                        collapse = "."
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
        message(
            "Both outcome_model_type and treatment_model_type are 'lasso': ",
            "post-selection / double-selection LASSO is utilized."
        )
    }

    if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("Package 'glmnet' is required for penalized regression.")
    }
    library(glmnet)


    do_single_fit <- function(y_sub, x_sub, model_type) {
        if (model_type == "linear") {
            df_temp <- data.frame(y = y_sub, x_sub)
            fit_lm <- lm(y ~ ., data = df_temp)
            return(list(fit = fit_lm, type = "lm"))
        } else if (model_type == "ridge") {
            fit_cv <- cv.glmnet(
                x = x_sub, y = y_sub, alpha = 0,
                lambda = lambda_seq
            )
            return(list(fit = fit_cv, type = "glmnet"))
        } else if (model_type == "lasso") {
            fit_cv <- cv.glmnet(
                x = x_sub, y = y_sub, alpha = 1,
                lambda = lambda_seq
            )
            return(list(fit = fit_cv, type = "glmnet"))
        } else {
            stop("Unsupported outcome_model_type: ", model_type)
        }
    }


    if (!is.null(selected_covars)) {
        act.outcome <- selected_covars[["outcome"]]
        act.treatment <- selected_covars[["treatment"]]
        XZ_sub_outcome <- XZ_design[, act.outcome, drop = FALSE]
        XZ_sub_treatment <- XZ_design[, act.treatment, drop = FALSE]
        outcome_fit <- do_single_fit(Yvec, XZ_sub_outcome, outcome_model_type)
        treatment_fit <- do_single_fit(Dvec, XZ_sub_treatment, treatment_model_type)
    } else {
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
        selected_covars <- list(outcome = act.outcome, treatment = act.treatment)
        # outcome
        XZ_sub_outcome <- XZ_design[, act.outcome, drop = FALSE]
        XZ_sub_treatment <- XZ_design[, act.treatment, drop = FALSE]
        outcome_sub <- data.frame(y = Yvec, XZ_sub_outcome)
        final_out_lm <- lm(y ~ ., data = outcome_sub)

        treatment_sub <- data.frame(d = Dvec, XZ_sub_treatment)
        final_treatment_lm <- lm(d ~ ., data = treatment_sub)

        # Overwrite the lasso fits with final unpenalized fits
        outcome_fit <- list(fit = final_out_lm, type = "lm")
        treatment_fit <- list(fit = final_treatment_lm, type = "lm")
    }


    predict_outcome <- function(subfit, newX) {
        if (subfit$type == "lm") {
            vars_used <- names(subfit$fit$coefficients)[-1]
            df_temp <- data.frame(newX[, vars_used, drop = FALSE])
            colnames(df_temp) <- vars_used
            preds <- predict(subfit$fit, newdata = df_temp)
            return(as.numeric(preds))
        } else if (subfit$type == "glmnet") {
            preds <- as.numeric(predict(subfit$fit, newx = newX, s = "lambda.min"))
            return(preds)
        } else {
            stop("Unknown subfit$type in predict_outcome")
        }
    }

    outcome_hat <- predict_outcome(outcome_fit, XZ_design)
    treatment_hat <- predict_outcome(treatment_fit, XZ_design)

    outcome_signal <- Yvec - outcome_hat
    treatment_signal <- Dvec - treatment_hat

    if (reduce.dimension == "bspline") {
        fit_spline <- function(Y, D, X_vec) {
            data_fit <- data.frame(
                Y = Y,
                D = D,
                X = X_vec
            )
            lm(
                formula = Y ~ D * bs(X, df = spline_df, degree = spline_degree),
                data    = data_fit
            )
        }

        fit_spline_model <- fit_spline(
            Y = outcome_signal, # Your tilde{Y}
            D = treatment_signal, # Your tilde{D}
            X_vec = Xvec # The covariate
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
            Y = outcome_signal, # Your tilde{Y}
            D = treatment_signal, # Your tilde{D}
            X = Xvec # The covariate
        )
        if (is.null(bw)) {
            sol.k <- interflex::interflex(
                estimator = "kernel", Y = "Y", D = "D", X = "X",
                data = data.k, X.eval = x.eval, CV = TRUE, parallel = T, cores = 31
            )
        } else {
            sol.k <- interflex::interflex(
                estimator = "kernel", Y = "Y", D = "D", X = "X",
                data = data.k, X.eval = x.eval, CV = FALSE, bw = bw
            )
        }

        cme_fit <- sol.k$est.kernel[[1]][, 2]
    }

    cme_df <- data.frame(
        X.eval = x.eval,
        CME_fit = cme_fit
    )

    if (verbose) {
        message("Estimation complete.")
    }

    # selected_covars_names <- NULL
    # if (!is.null(selected_covars)) {
    #  selected_covars_names <- colnames(XZ_design)[selected_covars]
    # }

    out_list <- list(
        # The final design matrix used
        XZ_design = XZ_design,

        # Fitted outcome models
        outcome_fit = outcome_fit,
        treatment_fit = treatment_fit,

        # Observed-sample predictions
        outcome_hat = outcome_hat,
        treatment_hat = treatment_hat,

        # Signals
        outcome_signal = outcome_signal,
        treatment_signal = treatment_signal,

        # CME curve on x.eval
        cme_df = cme_df,
        reduce.dimension = reduce.dimension,

        # If double-selection was used
        selected_covars = selected_covars
    )

    if (reduce.dimension == "bspline") {
        out_list$spline_fit <- fit_spline_model
    }

    if (reduce.dimension == "kernel") {
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

bootstrapCME_PLR <- function(
    data,
    Y,
    D,
    X,
    Z = NULL,
    B = 200,
    alpha = 0.05,
    outcome_model_type = "lasso", # can be "linear", "ridge", or "lasso"
    treatment_model_type = "lasso", # can be "linear", "ridge", or "lasso"
    # polynomial or B-spline expansion in the design matrix
    basis_type = c("polynomial", "bspline", "none"),
    include_interactions = TRUE,
    poly_degree = 2, # only used if basis_type="polynomial"
    # B-spline parameters (used if basis_type="bspline", or for final CME fit)
    spline_df = 4,
    spline_degree = 2,
    lambda_seq = NULL, # optional custom lambda sequence for glmnet
    reduce.dimension = c("bspline", "kernel"),
    bw = NULL,
    x.eval = NULL, # grid of X values for final CME curve
    verbose = TRUE) {
    ###########################################################################
    # 1. Fit once on the full dataset
    ###########################################################################
    fit_full <- estimateCME_PLR(
        data = data,
        Y = Y, D = D, X = X, Z = Z,
        outcome_model_type = outcome_model_type,
        treatment_model_type = treatment_model_type,
        basis_type = basis_type,
        include_interactions = include_interactions,
        poly_degree = poly_degree,
        spline_df = spline_df,
        spline_degree = spline_degree,
        lambda_seq = lambda_seq,
        reduce.dimension = reduce.dimension,
        bw = bw,
        x.eval = x.eval,
        verbose = verbose
    )
    message("Baseline Finished.\n")

    # Full-sample CME curves
    cme_out <- fit_full$cme_df$CME_fit
    x.eval <- fit_full$cme_df$X.eval
    nEval <- length(x.eval)
    n <- nrow(data)

    # If kernel dimension reduction, retrieve bandwidth
    if (fit_full$reduce.dimension == "kernel") {
        bw <- fit_full$bw
    } else {
        bw <- NULL
    }

    # Check if both-lasso was used
    selected_covars <- fit_full$selected_covars
    both_lasso <- !is.null(selected_covars) # if union_covars is non-NULL => both-lasso

    ###########################################################################
    # 2. Prepare storage for bootstrap
    ###########################################################################
    fit_mat_bs <- matrix(NA, nrow = B, ncol = nEval) # B x (# of eval points)
    idx_seq <- seq_len(n)

    ###########################################################################
    # 3. Set up parallel backend
    ###########################################################################
    if (!requireNamespace("doParallel", quietly = TRUE)) {
        stop("doParallel not available. Please install doParallel or remove parallel usage.")
    }
    nCores <- parallel::detectCores()
    cl <- parallel::makeCluster(nCores)
    doParallel::registerDoParallel(cl)

    ###########################################################################
    # 4. Parallel bootstrap loop
    ###########################################################################
    `%dopar%` <- foreach::`%dopar%`

    res_list <- foreach::foreach(
        b = 1:B,
        .combine = "rbind",
        .export = "estimateCME_PLR",
        .packages = c("splines", "glmnet", "interflex")
    ) %dopar% {
        set.seed(b + 1234) # optional reproducibility

        # resample indices
        idx_b <- sample(idx_seq, size = n, replace = TRUE)
        data_b <- data[idx_b, , drop = FALSE]

        # For both-lasso, use post-selection columns
        # and force linear as the final regression
        if (both_lasso) {
            XZ_b <- fit_full$XZ_design[idx_b, , drop = FALSE]
            out_model_b <- "linear"
            treatment_model_b <- "linear"
        } else {
            # otherwise use the user-specified model types
            XZ_b <- fit_full$XZ_design[idx_b, , drop = FALSE]
            out_model_b <- outcome_model_type
            treatment_model_b <- treatment_model_type
        }

        # Re-run estimateCME on the bootstrap sample
        fit_b <- estimateCME_PLR(
            data = data_b,
            Y = Y,
            D = D,
            X = X,
            Z = NULL,
            XZ_design = XZ_b,
            outcome_model_type = out_model_b,
            treatment_model_type = treatment_model_b,
            bw = bw,
            x.eval = x.eval,
            spline_df = spline_df,
            spline_degree = spline_degree,
            lambda_seq = lambda_seq,
            reduce.dimension = reduce.dimension,
            selected_covars = selected_covars
        )

        # Return CME estimates as a single row
        c(fit_b$cme_df$CME_fit)
    }

    # Stop cluster
    parallel::stopCluster(cl)

    # Store bootstrap CME results in fit_mat_bs
    for (b in seq_len(B)) {
        fit_mat_bs[b, ] <- res_list[b, 1:nEval]
    }

    ###########################################################################
    # 5. Pointwise normal-based intervals
    ###########################################################################
    fit_se <- apply(fit_mat_bs, 2, sd, na.rm = TRUE)
    zcrit <- qnorm(1 - alpha / 2)

    alpha_lower <- alpha / 2
    alpha_upper <- 1 - alpha / 2

    fit_ci_l <- apply(fit_mat_bs, 2, quantile, probs = alpha_lower, na.rm = TRUE)
    fit_ci_u <- apply(fit_mat_bs, 2, quantile, probs = alpha_upper, na.rm = TRUE)

    ###########################################################################
    # 6. Compute uniform confidence intervals using 'calculate_uniform_quantiles'
    ###########################################################################
    # The function expects dimension k x N:
    #   k = # of parameters (here, # of x.eval points),
    #   N = # of bootstrap draws.
    # Our 'fit_mat_bs' is B x nEval, so we transpose:
    theta_matrix <- t(fit_mat_bs) # now nEval x B
    uni_res <- calculate_uniform_quantiles(theta_matrix, alpha)
    Q_j <- uni_res$Q_j # nEval x 2
    coverage <- uni_res$coverage # coverage estimate

    # Q_j[,1] is the uniform lower bound for each X point
    # Q_j[,2] is the uniform upper bound


    # For simplicity, let's store them directly from Q_j:
    fit_ci_l_uni <- Q_j[, 1]
    fit_ci_u_uni <- Q_j[, 2]

    ###########################################################################
    # 7. Combine results into a final data frame
    ###########################################################################
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
    # 8. Return
    ###########################################################################
    return(list(
        results = fit_df,
        selected_covars = selected_covars,
        fit_full = fit_full,
        coverage = coverage # coverage of the uniform band
    ))
}
