interflex.DML <- function(data,
                          Y, # outcome
                          D, # treatment indicator
                          X, # moderator
                          treat.info,
                          diff.info,
                          Z = NULL, # covariates
                          weights = NULL, # weighting variable
                          ml_method = "rf",
                          trimming_threshold = 0.01,
                          n_estimators = 500,
                          solver = "adam",
                          max_iter = 10000,
                          alpha = 1e-5,
                          hidden_layer_sizes = c(5, 3, 2),
                          random_state = 1,
                          dml_method = "default",
                          poly_degree = 3,
                          lasso_alpha = 0.1,
                          casual_forest_criterion = "mse",
                          casual_forest_n_estimators = 1000,
                          casual_forest_in_impurity_decrease = 0.001,
                          CV_y = FALSE,
                          param_grid_y = NULL,
                          n_folds_y = 10,
                          scoring_y = "neg_mean_squared_error",
                          CV_t = FALSE,
                          param_grid_t = NULL,
                          n_folds_t = 10,
                          scoring_t = "neg_mean_squared_error",
                          CV_f = FALSE,
                          param_grid_f = NULL,
                          n_folds_f = 10,
                          scoring_f = "neg_mean_squared_error",
                          n_jobs = -1,
                          figure = TRUE,
                          CI = CI,
                          order = NULL,
                          subtitles = NULL,
                          show.subtitles = NULL,
                          Xdistr = "histogram", # c("density","histogram","none")
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

    # reticulate::use_virtualenv("r-reticulate", required = TRUE)
    reticulate::use_condaenv(condaenv = "r-reticulate")

    TE.output.all.list <- list()
    if (treat.type == "discrete") {
        python_script_path <- system.file("python/dml_treatment.py", package = "interflex")
        reticulate::source_python(python_script_path)
        for (char in other.treat) {
            data_part <- data[data[[D]] %in% c(treat.base, char), ]
            data_part[data[[D]] == treat.base, D] <- 0
            data_part[data[[D]] == char, D] <- 1
            TE.output.all.python <- marginal_effect_for_treatment(data_part,
                ml_method = ml_method, Y = Y, D = D, X = X, Z = Z, d_ref = 1,
                n_estimators = n_estimators,
                solver = solver, max_iter = max_iter, alpha = alpha, hidden_layer_sizes = hidden_layer_sizes, random_state = random_state,
                dml_method = dml_method,
                poly_degree = poly_degree, lasso_alpha = lasso_alpha,
                casual_forest_criterion = casual_forest_criterion,
                casual_forest_n_estimators = casual_forest_n_estimators,
                casual_forest_in_impurity_decrease = casual_forest_in_impurity_decrease,
                CV_y=CV_y, param_grid_y=param_grid_y, n_folds_y=n_folds_y, scoring_y=scoring_y,
                CV_t=CV_t, param_grid_t=param_grid_t, n_folds_t=n_folds_t, scoring_t=scoring_t,
                CV_f=CV_f, param_grid_f=param_grid_f, n_folds_f=n_folds_f, scoring_f=scoring_f,
                n_jobs=n_jobs
            )
            TE.output.all <- data.frame(TE.output.all.python, check.names = FALSE)
            TE.output.all.list[[other.treat.origin[char]]] <- TE.output.all
        }
    } else if (treat.type == "continuous") {
        python_script_path <- system.file("python/dml_treatment.py", package = "interflex")
        reticulate::source_python(python_script_path)
        k <- 1
        for (d_ref in D.sample) {
            label <- label.name[k]
            TE.output.all.python <- marginal_effect_for_treatment(data,
                ml_method = ml_method, Y = Y, D = D, X = X, Z = Z, d_ref = d_ref,
                n_estimators = n_estimators,
                solver = solver, max_iter = max_iter, alpha = alpha, hidden_layer_sizes = hidden_layer_sizes, random_state = random_state,
                dml_method = dml_method,
                poly_degree = poly_degree, lasso_alpha = lasso_alpha,
                casual_forest_criterion = casual_forest_criterion,
                casual_forest_n_estimators = casual_forest_n_estimators,
                casual_forest_in_impurity_decrease = casual_forest_in_impurity_decrease,
                CV_y=CV_y, param_grid_y=param_grid_y, n_folds_y=n_folds_y, scoring_y=scoring_y,
                CV_t=CV_t, param_grid_t=param_grid_t, n_folds_t=n_folds_t, scoring_t=scoring_t,
                CV_f=CV_f, param_grid_f=param_grid_f, n_folds_f=n_folds_f, scoring_f=scoring_f,
                n_jobs=n_jobs
            )
            TE.output.all <- data.frame(TE.output.all.python, check.names = FALSE)
            TE.output.all.list[[label]] <- TE.output.all
            k <- k + 1
        }
    }
    if (treat.type == "discrete") {
        final.output <- list(
            diff.info = diff.info,
            treat.info = treat.info,
            est.dml = TE.output.all.list,
            Xlabel = Xlabel,
            Dlabel = Dlabel,
            Ylabel = Ylabel,
            de = de,
            hist.out = hist.out,
            de.tr = treat_den,
            count.tr = treat.hist,
            estimator = "DML"
        )
    } else if (treat.type == "continuous") {
        final.output <- list(
            diff.info = diff.info,
            treat.info = treat.info,
            est.dml = TE.output.all.list,
            Xlabel = Xlabel,
            Dlabel = Dlabel,
            Ylabel = Ylabel,
            de = de,
            hist.out = hist.out,
            de.tr = de.tr,
            count.tr = NULL,
            estimator = "DML"
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
