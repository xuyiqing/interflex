# ============================================================================
# interflex DML estimator — Pure R implementation using DoubleML + mlr3
# ============================================================================
# Replaces the former Python/reticulate backend with R-native code.
# Uses the R DoubleML package (R6 interface) for core DML estimation,
# mlr3/mlr3learners for ML models, and manual BLP for CATE/GATE computation.
# ============================================================================

# --------------------------------------------------------------------------
# Helper: Map model name to mlr3 learner
# --------------------------------------------------------------------------
.set_dml_learner <- function(model, param, discrete_outcome) {
    # Ensure mlr3learners is loaded so learners (ranger, nnet, xgboost, etc.)
    # are registered in mlr3's dictionary
    if (requireNamespace("mlr3learners", quietly = TRUE)) {
        # Loading the namespace registers the learners
        loadNamespace("mlr3learners")
    }

    model_lower <- tolower(model)
    if (is.null(param)) param <- list()

    # --- linear / logistic / default ---
    if (model_lower %in% c("default", "linear", "logistic", "l", "d")) {
        if (discrete_outcome) {
            learner <- mlr3::lrn("classif.log_reg", predict_type = "prob")
        } else {
            learner <- mlr3::lrn("regr.lm")
        }
        return(learner)
    }

    # --- regularization (glmnet) ---
    if (model_lower %in% c("regularization", "r")) {
        if (discrete_outcome) {
            learner <- do.call(mlr3::lrn, c(list("classif.glmnet", predict_type = "prob"), param))
        } else {
            learner <- do.call(mlr3::lrn, c(list("regr.glmnet"), param))
        }
        return(learner)
    }

    # --- random forest (ranger) ---
    # Default params match sklearn RandomForest: max_features=1.0, min_samples_leaf=1
    if (model_lower %in% c("randomforest", "random forest", "random_forest",
                            "rf", "forest")) {
        mapped <- .map_rf_params(param)
        # Set sklearn-matching defaults if user didn't specify
        # sklearn RandomForest: n_estimators=100, min_samples_leaf=1
        if (is.null(mapped[["min.node.size"]])) mapped[["min.node.size"]] <- 1L
        if (is.null(mapped[["num.trees"]])) mapped[["num.trees"]] <- 100L
        if (discrete_outcome) {
            learner <- do.call(mlr3::lrn, c(list("classif.ranger", predict_type = "prob"), mapped))
        } else {
            # For regression, sklearn uses max_features=1.0 (all features)
            # ranger defaults to mtry=p/3 which produces much noisier nuisance estimates
            # We override to use all features unless user specified
            if (is.null(mapped[["mtry"]])) mapped[["mtry.ratio"]] <- 1.0
            learner <- do.call(mlr3::lrn, c(list("regr.ranger"), mapped))
        }
        return(learner)
    }

    # --- boosting (lightgbm / xgboost fallback) ---
    # lightgbm is algorithmically closest to sklearn HistGradientBoosting.
    # Falls back to xgboost if lightgbm is not available.
    if (model_lower %in% c("boosting", "gradient_boosting", "gradient boosting",
                            "hist_gradient_boosting", "hist gradient boosting",
                            "boost", "gradient_boost", "gradient boost",
                            "hist_gradient_boost", "hist gradient boost",
                            "b", "gb", "hgb")) {
        mapped <- .map_boosting_params(param)

        # Try lightgbm first (closest to sklearn HistGradientBoosting)
        has_lgb <- requireNamespace("mlr3extralearners", quietly = TRUE) &&
                   requireNamespace("lightgbm", quietly = TRUE)
        if (has_lgb) {
            loadNamespace("mlr3extralearners")
            # sklearn HGB defaults: learning_rate=0.1, max_iter=100, max_leaf_nodes=31,
            # min_samples_leaf=20.  sklearn also uses early_stopping by default
            # (n_iter_no_change=10, validation_fraction=0.1), so effective iterations
            # are typically 30-60.  We use num_iterations=50 to approximate this.
            if (is.null(mapped[["num_iterations"]])) mapped[["num_iterations"]] <- 50L
            if (is.null(mapped[["learning_rate"]])) mapped[["learning_rate"]] <- 0.1
            if (is.null(mapped[["num_leaves"]])) mapped[["num_leaves"]] <- 31L
            if (is.null(mapped[["min_data_in_leaf"]])) mapped[["min_data_in_leaf"]] <- 20L
            if (is.null(mapped[["verbose"]])) mapped[["verbose"]] <- -1L
            if (discrete_outcome) {
                learner <- do.call(mlr3::lrn, c(list("classif.lightgbm", predict_type = "prob"), mapped))
            } else {
                learner <- do.call(mlr3::lrn, c(list("regr.lightgbm"), mapped))
            }
        } else {
            # Fallback to xgboost
            if (is.null(mapped[["nrounds"]])) mapped[["nrounds"]] <- 100L
            if (is.null(mapped[["eta"]])) mapped[["eta"]] <- 0.1
            if (is.null(mapped[["max_depth"]])) mapped[["max_depth"]] <- 5L
            if (is.null(mapped[["min_child_weight"]])) mapped[["min_child_weight"]] <- 20L
            if (is.null(mapped[["subsample"]])) mapped[["subsample"]] <- 0.8
            if (is.null(mapped[["colsample_bytree"]])) mapped[["colsample_bytree"]] <- 0.8
            if (is.null(mapped[["verbose"]])) mapped[["verbose"]] <- 0L
            if (discrete_outcome) {
                learner <- do.call(mlr3::lrn, c(list("classif.xgboost", predict_type = "prob"), mapped))
            } else {
                learner <- do.call(mlr3::lrn, c(list("regr.xgboost"), mapped))
            }
        }
        return(learner)
    }

    # --- neural network (nnet) ---
    if (model_lower %in% c("network", "neural_network", "neural network", "nn")) {
        mapped <- .map_nn_params(param)
        if (is.null(mapped[["size"]])) mapped[["size"]] <- 20L
        if (is.null(mapped[["decay"]])) mapped[["decay"]] <- 0.01
        if (is.null(mapped[["trace"]])) mapped[["trace"]] <- FALSE
        if (discrete_outcome) {
            learner <- do.call(mlr3::lrn, c(list("classif.nnet", predict_type = "prob"), mapped))
        } else {
            learner <- do.call(mlr3::lrn, c(list("regr.nnet"), mapped))
        }
        return(learner)
    }

    stop("'model_y' and 'model_t' should be one of 'linear', 'regularization', 'rf', 'hgb', and 'nn'.")
}

# --------------------------------------------------------------------------
# Helper: Map sklearn RandomForest params -> ranger params
# --------------------------------------------------------------------------
.map_rf_params <- function(param) {
    if (is.null(param) || length(param) == 0L) return(list())

    mapping <- c(
        n_estimators    = "num.trees",
        max_depth       = "max.depth",
        min_samples_leaf = "min.node.size",
        max_features    = "mtry",
        n_jobs          = "num.threads",
        random_state    = "seed"
    )

    out <- list()
    for (nm in names(param)) {
        mapped_name <- if (nm %in% names(mapping)) mapping[[nm]] else nm
        out[[mapped_name]] <- param[[nm]]
    }
    out
}

# --------------------------------------------------------------------------
# Helper: Map sklearn HistGradientBoosting params -> xgboost params
# --------------------------------------------------------------------------
.map_boosting_params <- function(param) {
    if (is.null(param) || length(param) == 0L) return(list(verbose = 0L))

    mapping <- c(
        max_iter          = "nrounds",
        n_estimators      = "nrounds",
        max_depth         = "max_depth",
        learning_rate     = "eta",
        min_samples_leaf  = "min_child_weight",
        max_leaf_nodes    = "max_leaves",
        l2_regularization = "lambda",
        max_features      = "colsample_bytree"
    )

    out <- list(verbose = 0L)
    for (nm in names(param)) {
        mapped_name <- if (nm %in% names(mapping)) mapping[[nm]] else nm
        out[[mapped_name]] <- param[[nm]]
    }
    out
}

# --------------------------------------------------------------------------
# Helper: Map sklearn MLP params -> nnet params
# --------------------------------------------------------------------------
.map_nn_params <- function(param) {
    if (is.null(param) || length(param) == 0L) return(list(maxit = 3000L))

    mapping <- c(
        hidden_layer_sizes = "size",
        max_iter           = "maxit"
    )

    out <- list()
    for (nm in names(param)) {
        mapped_name <- if (nm %in% names(mapping)) mapping[[nm]] else nm
        val <- param[[nm]]
        # nnet supports only a single hidden layer
        if (mapped_name == "size" && length(val) > 1L) {
            warning("nnet supports only a single hidden layer; using the first element of hidden_layer_sizes.")
            val <- val[1L]
        }
        out[[mapped_name]] <- val
    }
    # default maxit if not set
    if (is.null(out[["maxit"]])) out[["maxit"]] <- 3000L
    out
}

# --------------------------------------------------------------------------
# Helper: Build paradox ParamSet from user param_grid
# --------------------------------------------------------------------------
.build_dml_tuning_grid <- function(param_grid, learner_type) {
    if (is.null(param_grid) || length(param_grid) == 0L) return(NULL)

    # Determine the right parameter name mapping based on learner_type
    model_lower <- tolower(learner_type)
    if (model_lower %in% c("randomforest", "random forest", "random_forest",
                            "rf", "forest")) {
        mapping <- c(
            n_estimators    = "num.trees",
            max_depth       = "max.depth",
            min_samples_leaf = "min.node.size",
            max_features    = "mtry"
        )
    } else if (model_lower %in% c("boosting", "gradient_boosting", "gradient boosting",
                                    "hist_gradient_boosting", "hist gradient boosting",
                                    "boost", "gradient_boost", "gradient boost",
                                    "hist_gradient_boost", "hist gradient boost",
                                    "b", "gb", "hgb")) {
        mapping <- c(
            max_iter          = "nrounds",
            n_estimators      = "nrounds",
            max_depth         = "max_depth",
            learning_rate     = "eta",
            min_samples_leaf  = "min_child_weight",
            max_leaf_nodes    = "max_leaves",
            l2_regularization = "lambda",
            max_features      = "colsample_bytree"
        )
    } else if (model_lower %in% c("network", "neural_network", "neural network", "nn")) {
        mapping <- c(
            hidden_layer_sizes = "size",
            max_iter           = "maxit"
        )
    } else {
        mapping <- character(0)
    }

    param_list <- list()
    for (nm in names(param_grid)) {
        mapped_name <- if (nm %in% names(mapping)) mapping[[nm]] else nm
        vals <- param_grid[[nm]]
        if (is.integer(vals) || all(vals == floor(vals))) {
            param_list[[mapped_name]] <- paradox::p_int(
                lower = as.integer(min(vals)),
                upper = as.integer(max(vals))
            )
        } else {
            param_list[[mapped_name]] <- paradox::p_dbl(
                lower = min(vals),
                upper = max(vals)
            )
        }
    }
    do.call(paradox::ps, param_list)
}

# --------------------------------------------------------------------------
# Helper: Compute CATE via BLP of pseudo-outcomes onto B-spline basis
# --------------------------------------------------------------------------
.compute_cate_blp <- function(dml_model, data_X, X_name, n_grid = 50L) {
    n <- length(data_X)

    # 1. Extract pseudo-outcomes (influence function values)
    psi <- .extract_psi(dml_model, n)

    # 2. B-spline basis for projection (df=5, degree=2 per spec)
    #    Include intercept to match Python patsy.dmatrix() behaviour
    B_spline <- splines::bs(data_X, df = 5, degree = 2)
    B_train <- cbind(1, B_spline)

    # 3. Evaluation grid
    x_grid <- seq(min(data_X), max(data_X), length.out = n_grid)
    B_grid <- cbind(1, predict(B_spline, newx = x_grid))

    # 4. OLS projection: beta_hat = (B'B)^{-1} B' psi
    BtB <- crossprod(B_train)
    BtB_inv <- .safe_solve(BtB)
    beta_hat <- BtB_inv %*% crossprod(B_train, psi)

    # 5. Predicted CATE on grid
    cate_hat <- as.numeric(B_grid %*% beta_hat)

    # 6. HC0 sandwich variance: Omega = (B'B)^{-1} B' diag(e^2) B (B'B)^{-1}
    #    Matches Python statsmodels model.cov_HC0 used by DoubleML BLP
    #    crossprod(B*e) = t(B*e) %*% (B*e) = sum_i e_i^2 * B_i B_i' = B'diag(e^2)B
    residuals <- as.numeric(psi - B_train %*% beta_hat)
    meat <- crossprod(B_train * residuals)
    Omega <- BtB_inv %*% meat %*% BtB_inv

    # 7. Variance-covariance on the grid
    Sigma_grid <- B_grid %*% Omega %*% t(B_grid)
    se_grid <- sqrt(pmax(diag(Sigma_grid), 0))

    # 8. Pointwise CIs
    z_alpha <- stats::qnorm(0.025)
    lower_pw <- cate_hat + z_alpha * se_grid
    upper_pw <- cate_hat - z_alpha * se_grid

    # 9. Uniform CIs via calculate_delta_uniformCI from uniform.R
    c_alpha <- tryCatch(
        calculate_delta_uniformCI(Sigma_grid, alpha = 0.05, N = 2000L),
        error = function(e) stats::qnorm(1 - 0.05 / (2 * n_grid))
    )
    lower_unif <- cate_hat - c_alpha * se_grid
    upper_unif <- cate_hat + c_alpha * se_grid

    data.frame(
        `X`                       = x_grid,
        `ME`                      = cate_hat,
        `sd`                      = se_grid,
        `lower CI(95%)`           = lower_pw,
        `upper CI(95%)`           = upper_pw,
        `lower uniform CI(95%)`   = lower_unif,
        `upper uniform CI(95%)`   = upper_unif,
        check.names = FALSE
    )
}

# --------------------------------------------------------------------------
# Helper: Compute GATE via BLP of pseudo-outcomes onto group dummies
# --------------------------------------------------------------------------
.compute_gate_blp <- function(dml_model, data_X, X_name) {
    n <- length(data_X)

    # 1. Extract pseudo-outcomes
    psi <- .extract_psi(dml_model, n)

    # 2. Group dummy basis
    groups <- sort(unique(data_X))
    K <- length(groups)
    B_dummy <- matrix(0, nrow = n, ncol = K)
    for (k in seq_len(K)) {
        B_dummy[, k] <- as.numeric(data_X == groups[k])
    }

    # 3. OLS projection
    BtB <- crossprod(B_dummy)
    BtB_inv <- .safe_solve(BtB)
    beta_hat <- as.numeric(BtB_inv %*% crossprod(B_dummy, psi))

    # 4. HC0 sandwich variance (matching Python statsmodels)
    #    crossprod(B*e) = B'diag(e^2)B
    residuals <- as.numeric(psi - B_dummy %*% beta_hat)
    meat <- crossprod(B_dummy * residuals)
    Omega <- BtB_inv %*% meat %*% BtB_inv
    se <- sqrt(pmax(diag(Omega), 0))

    # 5. Pointwise CIs
    z_alpha <- stats::qnorm(0.025)
    lower_pw <- beta_hat + z_alpha * se
    upper_pw <- beta_hat - z_alpha * se

    # 6. Uniform CIs
    c_alpha <- tryCatch(
        calculate_delta_uniformCI(Omega, alpha = 0.05, N = 2000L),
        error = function(e) stats::qnorm(1 - 0.05 / (2 * K))
    )
    lower_unif <- beta_hat - c_alpha * se
    upper_unif <- beta_hat + c_alpha * se

    data.frame(
        `X`                       = groups,
        `ME`                      = beta_hat,
        `sd`                      = se,
        `lower CI(95%)`           = lower_pw,
        `upper CI(95%)`           = upper_pw,
        `lower uniform CI(95%)`   = lower_unif,
        `upper uniform CI(95%)`   = upper_unif,
        check.names = FALSE
    )
}

# --------------------------------------------------------------------------
# Helper: Extract per-observation influence function values from DML model
# --------------------------------------------------------------------------
.extract_psi <- function(dml_model, n) {
    # Primary: use psi_b from the DoubleML model
    psi <- tryCatch({
        psi_raw <- dml_model$psi_b
        if (is.array(psi_raw) && length(dim(psi_raw)) == 3L) {
            psi_raw[, 1, 1]
        } else if (is.matrix(psi_raw)) {
            psi_raw[, 1]
        } else {
            as.numeric(psi_raw)
        }
    }, error = function(e) NULL)

    if (!is.null(psi) && length(psi) == n) return(psi)

    # Fallback: reconstruct from stored predictions and coefficients
    # For PLR: psi_b = D_tilde * Y_tilde where Y_tilde = Y - l_hat, D_tilde = D - m_hat
    # For IRM: use the AIPW score
    psi2 <- tryCatch({
        theta_hat <- dml_model$coef
        # Access all_psi which is psi_a * theta + psi_b
        all_psi <- dml_model$all_psi
        if (is.array(all_psi) && length(dim(all_psi)) == 3L) {
            # all_psi[i, 1, 1] = psi_a[i] * theta + psi_b[i], which is the centered score
            # The pseudo-outcome is theta + psi / (-E[psi_a])
            all_psi[, 1, 1] + theta_hat
        } else {
            NULL
        }
    }, error = function(e) NULL)

    if (!is.null(psi2) && length(psi2) == n) return(psi2)

    # Last resort: use the coefficient as a constant (ATE for everyone)
    warning("Could not extract individual pseudo-outcomes from DML model. Using constant ATE.")
    rep(dml_model$coef, n)
}

# --------------------------------------------------------------------------
# Helper: Safe matrix inversion with ridge regularisation fallback
# --------------------------------------------------------------------------
.safe_solve <- function(mat, tol = 1e-8) {
    tryCatch(
        solve(mat),
        error = function(e) {
            solve(mat + tol * diag(nrow(mat)))
        }
    )
}

# --------------------------------------------------------------------------
# Core worker: Run a single DML estimation
# --------------------------------------------------------------------------
.run_dml_estimation <- function(data, Y, D, X, Z,
                                model.y, param.y, param.grid.y, scoring.y,
                                model.t, param.t, param.grid.t, scoring.t,
                                CV, n.folds, cf.n.folds, cf.n.rep, gate) {

    # 1. Detect outcome / treatment type
    y_vals <- sort(unique(data[[Y]]))
    discrete_outcome <- length(y_vals) == 2L && all(y_vals == c(0, 1))
    discrete_treatment <- length(unique(data[[D]])) <= 5L

    # 2. Covariates
    if (is.null(Z)) Z <- character(0)
    if (is.character(Z) && length(Z) == 1L && Z == "") Z <- character(0)
    covariates <- c(Z, X)

    # 3. Learners
    ml_y <- .set_dml_learner(model.y, param.y, discrete_outcome)
    ml_t <- .set_dml_learner(model.t, param.t, discrete_treatment)

    # 4. B-spline expansion for linear models
    is_linear_y <- tolower(model.y) %in% c("linear", "logistic", "l", "d", "regularization", "r")
    is_linear_t <- tolower(model.t) %in% c("linear", "logistic", "l", "d", "regularization", "r")

    if (is_linear_y && is_linear_t) {
        bs_mat <- splines::bs(data[[X]], degree = 3, df = 5)
        bs_colnames <- paste0("X_bs_", seq_len(ncol(bs_mat)))
        bs_df <- as.data.frame(bs_mat)
        colnames(bs_df) <- bs_colnames
        data_for_dml <- cbind(data[, c(Y, D, X, Z), drop = FALSE], bs_df)
        covariates <- c(covariates, bs_colnames)
    } else {
        data_for_dml <- data[, c(Y, D, X, Z), drop = FALSE]
    }

    # Ensure treatment is numeric
    data_for_dml[[D]] <- as.numeric(data_for_dml[[D]])

    # 4b. Standardise covariates for models sensitive to scale (nn, hgb).
    #     Python sklearn's MLPRegressor and HistGradientBoosting do this
    #     internally; R's nnet and lightgbm do NOT.
    is_scale_sensitive <- tolower(model.y) %in% c("network","neural_network","neural network","nn",
                                                   "boosting","gradient_boosting","gradient boosting",
                                                   "hist_gradient_boosting","hist gradient boosting",
                                                   "boost","gradient_boost","gradient boost",
                                                   "hist_gradient_boost","hist gradient boost",
                                                   "b","gb","hgb") ||
                          tolower(model.t) %in% c("network","neural_network","neural network","nn",
                                                   "boosting","gradient_boosting","gradient boosting",
                                                   "hist_gradient_boosting","hist gradient boosting",
                                                   "boost","gradient_boost","gradient boost",
                                                   "hist_gradient_boost","hist gradient boost",
                                                   "b","gb","hgb")
    if (is_scale_sensitive) {
        # Standardise Y and all covariates (not D)
        scale_info <- list()
        for (col in c(Y, covariates)) {
            mu <- mean(data_for_dml[[col]], na.rm = TRUE)
            s  <- sd(data_for_dml[[col]], na.rm = TRUE)
            if (s < 1e-12) s <- 1  # avoid division by zero for constant columns
            scale_info[[col]] <- list(mean = mu, sd = s)
            data_for_dml[[col]] <- (data_for_dml[[col]] - mu) / s
        }
    }

    # 5. Create DoubleMLData
    dml_data <- DoubleML::DoubleMLData$new(
        data = data.table::as.data.table(data_for_dml),
        y_col = Y,
        d_cols = D,
        x_cols = covariates
    )

    # 6. Create DML model
    if (discrete_treatment) {
        # R ranger can produce propensity scores of exactly 0 or 1 (unlike sklearn RF
        # which outputs vote proportions that are naturally bounded away from 0/1).
        # Trimming is essential to prevent IPW blow-up in psi_b.
        dml_model <- DoubleML::DoubleMLIRM$new(
            data = dml_data,
            ml_g = ml_y,
            ml_m = ml_t,
            n_folds = as.integer(cf.n.folds),
            n_rep = as.integer(cf.n.rep),
            trimming_threshold = 0.01
        )
    } else {
        dml_model <- DoubleML::DoubleMLPLR$new(
            data = dml_data,
            ml_l = ml_y,
            ml_m = ml_t,
            n_folds = as.integer(cf.n.folds),
            n_rep = as.integer(cf.n.rep)
        )
    }

    # 7. Tuning (if CV = TRUE and grids provided)
    params <- list("model.y" = NULL, "model.t" = NULL)
    if (isTRUE(CV)) {
        tune_param_set <- list()

        if (!is.null(param.grid.y) && length(param.grid.y) > 0L) {
            learner_y_name <- if (discrete_treatment) "ml_g" else "ml_l"
            ps_y <- .build_dml_tuning_grid(param.grid.y, model.y)
            if (!is.null(ps_y)) tune_param_set[[learner_y_name]] <- ps_y
        }
        if (!is.null(param.grid.t) && length(param.grid.t) > 0L) {
            ps_t <- .build_dml_tuning_grid(param.grid.t, model.t)
            if (!is.null(ps_t)) tune_param_set[["ml_m"]] <- ps_t
        }

        if (length(tune_param_set) > 0L) {
            tune_settings <- list(
                terminator = mlr3tuning::trm("evals", n_evals = 100L),
                algorithm  = mlr3tuning::tnr("grid_search"),
                rsmp_tune  = mlr3::rsmp("cv", folds = as.integer(n.folds))
            )
            tryCatch(
                dml_model$tune(param_set = tune_param_set, tune_settings = tune_settings),
                error = function(e) {
                    warning(paste0("DML tuning failed, using default hyperparameters. Error: ",
                                   conditionMessage(e)))
                }
            )
        }

        # Record tuned params
        tryCatch({
            if (discrete_treatment) {
                params[["model.y"]] <- dml_model$params[["ml_g"]]
                params[["model.t"]] <- dml_model$params[["ml_m"]]
            } else {
                params[["model.y"]] <- dml_model$params[["ml_l"]]
                params[["model.t"]] <- dml_model$params[["ml_m"]]
            }
        }, error = function(e) NULL)
    }

    # 8. Fit
    tryCatch(
        dml_model$fit(store_predictions = TRUE),
        error = function(e) {
            stop(paste0(
                "DML estimation failed. ",
                "Ensure required packages (DoubleML, mlr3, mlr3learners) are installed. ",
                "Original error: ", conditionMessage(e)
            ))
        }
    )

    # 9. Extract nuisance losses
    losses <- tryCatch(dml_model$nuisance_loss, error = function(e) NULL)

    # 10. Compute CATE (always)
    #     Use ORIGINAL (unscaled) X for the BLP grid, since the plot is in original units
    cate_df <- .compute_cate_blp(dml_model, data[[X]], X)

    #     If Y was scaled, rescale psi-derived quantities back to original units
    if (is_scale_sensitive && !is.null(scale_info[[Y]])) {
        y_sd <- scale_info[[Y]]$sd
        cate_df[["ME"]]                    <- cate_df[["ME"]]                    * y_sd
        cate_df[["sd"]]                    <- cate_df[["sd"]]                    * y_sd
        cate_df[["lower CI(95%)"]]         <- cate_df[["lower CI(95%)"]]         * y_sd
        cate_df[["upper CI(95%)"]]         <- cate_df[["upper CI(95%)"]]         * y_sd
        cate_df[["lower uniform CI(95%)"]] <- cate_df[["lower uniform CI(95%)"]] * y_sd
        cate_df[["upper uniform CI(95%)"]] <- cate_df[["upper uniform CI(95%)"]] * y_sd
    }

    # 11. Compute GATE (conditional)
    gate_df <- data.frame()
    if (isTRUE(gate)) {
        gate_df <- .compute_gate_blp(dml_model, data[[X]], X)
    }

    # 12. Return
    list(
        cate_df = cate_df,
        gate_df = gate_df,
        params  = params,
        losses  = losses
    )
}

# ============================================================================
# Main function: interflex.dml
# ============================================================================
interflex.dml <- function(data,
                          Y, # outcome
                          D, # treatment indicator
                          X, # moderator
                          treat.info,
                          diff.info,
                          Z = NULL, # covariates
                          weights = NULL, # weighting variable
                          model.y = "rf",
                          param.y = NULL,
                          param.grid.y = NULL,
                          scoring.y = "neg_mean_squared_error",
                          model.t = "rf",
                          param.t = NULL,
                          param.grid.t = NULL,
                          scoring.t = "neg_mean_squared_error",
                          CV = FALSE,
                          n.folds = 10,
                          n.jobs = -1,
                          cf.n.folds = 5,
                          cf.n.rep = 1,
                          gate = FALSE,
                          figure = TRUE,
                          show.uniform.CI = TRUE,
                          CI = CI,
                          order = NULL,
                          subtitles = NULL,
                          show.subtitles = NULL,
                          Xdistr = "histogram", # ("density","histogram","none")
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

    # ---- PREAMBLE (preserve existing logic) ----
    diff.values.plot <- diff.info[["diff.values.plot"]]
    ti <- .extract_treat_info(treat.info)
    treat.type <- ti$treat.type
    all.treat <- all.treat.origin <- NULL
    if (treat.type == "discrete") {
        other.treat <- ti$other.treat
        other.treat.origin <- ti$other.treat.origin
        all.treat <- ti$all.treat
        all.treat.origin <- ti$all.treat.origin
    }
    if (treat.type == "continuous") {
        D.sample <- ti$D.sample
        label.name <- ti$label.name
    }

    n <- dim(data)[1]
    if (is.null(weights)) {
        w <- rep(1, n)
    } else {
        w <- data[, weights]
    }
    data[["WEIGHTS"]] <- w

    dens <- .compute_density(data, X, D, weights, treat.type, all.treat, all.treat.origin)
    de <- dens$de
    treat_den <- dens$treat_den
    hists <- .compute_histograms(data, X, D, weights, treat.type, all.treat, all.treat.origin)
    hist.out <- hists$hist.out
    treat.hist <- hists$treat.hist
    if (treat.type == "continuous") {
        de.tr <- NULL
    }

    treat.base <- treat.info[["base"]]

    # ---- DML ESTIMATION ----
    TE.output.all.list <- list()
    TE.G.output.all.list <- list()
    dml_params <- NULL
    dml_losses <- NULL

    if (treat.type == "discrete") {
        for (char in other.treat) {
            # Subset and binarise treatment
            data_part <- data[data[[D]] %in% c(treat.base, char), ]
            data_part[[D]] <- ifelse(data_part[[D]] == treat.base, 0L, 1L)

            result <- .run_dml_estimation(
                data_part, Y, D, X, Z,
                model.y, param.y, param.grid.y, scoring.y,
                model.t, param.t, param.grid.t, scoring.t,
                CV, n.folds, cf.n.folds, cf.n.rep, gate
            )

            TE.output.all.list[[other.treat.origin[char]]] <- result$cate_df
            TE.G.output.all.list[[other.treat.origin[char]]] <- result$gate_df
            dml_params <- result$params
            dml_losses <- result$losses
        }
    } else if (treat.type == "continuous") {
        result <- .run_dml_estimation(
            data, Y, D, X, Z,
            model.y, param.y, param.grid.y, scoring.y,
            model.t, param.t, param.grid.t, scoring.t,
            CV, n.folds, cf.n.folds, cf.n.rep, gate
        )

        for (k in seq_along(D.sample)) {
            label <- label.name[k]
            TE.output.all.list[[label]] <- result$cate_df
            TE.G.output.all.list[[label]] <- result$gate_df
        }
        dml_params <- result$params
        dml_losses <- result$losses
    }

    # ---- ASSEMBLE OUTPUT ----
    if (treat.type == "discrete") {
        final.output <- list(
            diff.info  = diff.info,
            treat.info = treat.info,
            est.dml    = TE.output.all.list,
            g.est.dml  = TE.G.output.all.list,
            Xlabel     = Xlabel,
            Dlabel     = Dlabel,
            Ylabel     = Ylabel,
            de         = de,
            hist.out   = hist.out,
            de.tr      = treat_den,
            count.tr   = treat.hist,
            dml.models = dml_params,
            dml.losses = dml_losses,
            estimator  = "dml"
        )
    } else if (treat.type == "continuous") {
        final.output <- list(
            diff.info  = diff.info,
            treat.info = treat.info,
            est.dml    = TE.output.all.list,
            g.est.dml  = TE.G.output.all.list,
            Xlabel     = Xlabel,
            Dlabel     = Dlabel,
            Ylabel     = Ylabel,
            de         = de,
            hist.out   = hist.out,
            de.tr      = de.tr,
            count.tr   = NULL,
            dml.models = dml_params,
            dml.losses = dml_losses,
            estimator  = "dml"
        )
    }

    # ---- PLOT ----
    if (figure) {
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
            show.uniform.CI = show.uniform.CI,
            cex.main = cex.main,
            cex.sub = cex.sub,
            cex.lab = cex.lab,
            cex.axis = cex.axis,
            interval = interval,
            file = file,
            ncols = ncols,
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
