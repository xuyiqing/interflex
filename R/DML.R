interflex.DML <- function(data,
                          Y, # outcome
                          D, # treatment indicator
                          X, # moderator
                          treat.info,
                          diff.info,
                          Z = NULL, # covariates
                          weights = NULL, # weighting variable
                          CV.first = 2,
                          model.y = "rf",
                          param.y = NULL,
                          CV.y = FALSE,
                          param.grid.y = NULL,
                          n.folds.y = 10,
                          scoring.y = "neg_mean_squared_error",
                          model.t = "rf",
                          param.t = NULL,
                          CV.t = FALSE,
                          param.grid.t = NULL,
                          n.folds.t = 10,
                          scoring.t = "neg_mean_squared_error",
                          model.final = "linear",
                          param.final = NULL,
                          CV.final = FALSE,
                          param.grid.final = NULL,
                          n.folds.final = 10,
                          scoring.final = "neg_mean_squared_error",
                          n.jobs = -1,
                          figure = TRUE,
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
        python_script_path <- system.file("python/dml.py", package = "interflex")
        reticulate::source_python(python_script_path)
        for (char in other.treat) {
            data_part <- data[data[[D]] %in% c(treat.base, char), ]
            data_part[data_part[[D]] == treat.base, D] <- 0
            data_part[data_part[[D]] == char, D] <- 1
            result <- marginal_effect_for_treatment(data_part,
                Y = Y, D = D, X = X, Z = Z, d_ref = 1,
                CV_first = CV.first,
                model_y = model.y,
                param_y = param.y,
                CV_y = CV.y,
                param_grid_y = param.grid.y,
                n_folds_y = n.folds.y,
                scoring_y = scoring.y,
                model_t = model.t,
                param_t = param.t,
                CV_t = CV.t,
                param_grid_t = param.grid.t,
                n_folds_t = n.folds.t,
                scoring_t = scoring.t,
                model_final = model.final,
                param_final = param.final,
                CV_final = CV.final,
                param_grid_final = param.grid.final,
                n_folds_final = n.folds.final,
                scoring_final = scoring.final,
                n_jobs = n.jobs
            )
            TE.output.all <- data.frame(result[1], check.names = FALSE)
            TE.output.all.list[[other.treat.origin[char]]] <- TE.output.all
        }
    } else if (treat.type == "continuous") {
        python_script_path <- system.file("python/dml.py", package = "interflex")
        reticulate::source_python(python_script_path)
        k <- 1
        for (d_ref in D.sample) {
            label <- label.name[k]
            result <- marginal_effect_for_treatment(data,
                Y = Y, D = D, X = X, Z = Z, d_ref = d_ref,
                CV_first = CV.first,
                model_y = model.y,
                param_y = dict(param.y),
                CV_y = CV.y,
                param_grid_y = dict(param.grid.y),
                n_folds_y = n.folds.y,
                scoring_y = scoring.y,
                model_t = model.t,
                param_t = dict(param.t),
                CV_t = CV.t,
                param_grid_t = dict(param.grid.t),
                n_folds_t = n.folds.t,
                scoring_t = scoring.t,
                model_final = model.final,
                param_final = dict(param.final),
                CV_final = CV.final,
                param_grid_final = dict(param.grid.final),
                n_folds_final = n.folds.final,
                scoring_final = scoring.final,
                n_jobs = n.jobs
            )
            TE.output.all <- data.frame(result[1], check.names = FALSE)
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
            model.y = result[2][[1]],
            model.t = result[3][[1]],
            model.final = result[4][[1]],
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
            model.y = result[2][[1]],
            model.t = result[3][[1]],
            model.final = result[4][[1]],
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
