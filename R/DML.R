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

    # reticulate::use_virtualenv("r-reticulate", required = TRUE)
    # use_condaenv(condaenv = "r-reticulate", required = TRUE)

    TE.output.all.list <- list()
    TE.G.output.all.list <- list()
    python_script_path <- system.file("python/dml.py", package = "interflex")
    if (!nzchar(python_script_path)) {
        stop("DML Python script not found. The file 'python/dml.py' is missing from the interflex package installation.")
    }
    tryCatch(
        source_python(python_script_path),
        error = function(e) {
            stop(paste0(
                "Failed to load DML Python script. ",
                "Ensure Python and required packages (econml, sklearn) are installed. ",
                "Original error: ", conditionMessage(e)
            ))
        }
    )
    if (treat.type == "discrete") {
        for (char in other.treat) {
            data_part <- data[data[[D]] %in% c(treat.base, char), ]
            data_part[data_part[[D]] == treat.base, D] <- 0L
            data_part[data_part[[D]] == char, D] <- 1L
            result <- marginal_effect_for_treatment(data_part,
                Y = Y, D = D, X = X, Z = Z,
                model_y = model.y,
                param_y = dict(param.y),
                param_grid_y = dict(param.grid.y),
                scoring_y = scoring.y,
                model_t = model.t,
                param_t = dict(param.t),
                param_grid_t = dict(param.grid.t),
                scoring_t = scoring.t,
                CV = CV,
                n_folds = n.folds,
                n_jobs = n.jobs,
                cf_n_folds = cf.n.folds,
                cf_n_rep = cf.n.rep,
                gate = gate
            )
            TE.output.all <- data.frame(result[1], check.names = FALSE)
            TE.G.output.all <- data.frame(result[2], check.names = FALSE)
            TE.output.all.list[[other.treat.origin[char]]] <- TE.output.all
            TE.G.output.all.list[[other.treat.origin[char]]] <- TE.G.output.all
        }
    } else if (treat.type == "continuous") {
        result <- marginal_effect_for_treatment(data,
            Y = Y, D = D, X = X, Z = Z,
            model_y = model.y,
            param_y = dict(param.y),
            param_grid_y = dict(param.grid.y),
            scoring_y = scoring.y,
            model_t = model.t,
            param_t = dict(param.t),
            param_grid_t = dict(param.grid.t),
            scoring_t = scoring.t,
            CV = CV,
            n_folds = n.folds,
            n_jobs = n.jobs,
            cf_n_folds = cf.n.folds,
            cf_n_rep = cf.n.rep,
            gate = gate
        )
        TE.output.all <- data.frame(result[1], check.names = FALSE)
        TE.G.output.all <- data.frame(result[2], check.names = FALSE)

        k <- 1
        for (d_ref in D.sample) {
            label <- label.name[k]
            TE.output.all.list[[label]] <- TE.output.all
            TE.G.output.all.list[[label]] <- TE.G.output.all
            k <- k + 1
        }
    }
    if (treat.type == "discrete") {
        final.output <- list(
            diff.info = diff.info,
            treat.info = treat.info,
            est.dml = TE.output.all.list,
            g.est.dml = TE.G.output.all.list,
            Xlabel = Xlabel,
            Dlabel = Dlabel,
            Ylabel = Ylabel,
            de = de,
            hist.out = hist.out,
            de.tr = treat_den,
            count.tr = treat.hist,
            dml.models = result[3][[1]],
            dml.losses = result[4],
            estimator = "dml"
        )
    } else if (treat.type == "continuous") {
        final.output <- list(
            diff.info = diff.info,
            treat.info = treat.info,
            est.dml = TE.output.all.list,
            g.est.dml = TE.G.output.all.list,
            Xlabel = Xlabel,
            Dlabel = Dlabel,
            Ylabel = Ylabel,
            de = de,
            hist.out = hist.out,
            de.tr = de.tr,
            count.tr = NULL,
            dml.models = result[3][[1]],
            dml.losses = result[4],
            estimator = "dml"
        )
    }

    # Plot
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
