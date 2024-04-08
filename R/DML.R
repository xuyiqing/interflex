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
                          n_folds = 5,
                          n_estimators = 500,
                          solver = "adam",
                          max_iter = 10000,
                          alpha = 1e-5,
                          hidden_layer_sizes = c(5, 3, 2),
                          random_state = 1,
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
        warning("continuous")
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
        if (treat.type == "continuous") {
            warnings("continuous")
        }
    }

    treat.base <- treat.info[["base"]]

    reticulate::use_virtualenv("r-reticulate", required = TRUE)
    reticulate::source_python(paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/../Python/btco_demo.py"))

    TE.output.all.list <- list()
    for (char in other.treat) {
        data_part <- data[data[[D]] %in% c(treat.base, char), ]
        TE.output.all.python <- marginal_effect_for_R(data_part,
            ml_method = ml_method, Y = Y, D = D, X = X, Z = Z,
            trimming_threshold = trimming_threshold, n_folds = n_folds,
            n_estimators = n_estimators,
            solver = solver, max_iter = max_iter, alpha = alpha, hidden_layer_sizes = hidden_layer_sizes, random_state = random_state
        )
        TE.output.all <- data.frame(TE.output.all.python, check.names = FALSE)
        TE.output.all.list[[other.treat.origin[char]]] <- TE.output.all
    }

    final.output <- list(
        diff.info = diff.info,
        treat.info = treat.info,
        est.dml = TE.output.all.list,
        Xlabel = Xlabel,
        Dlabel = Dlabel,
        Ylabel = Ylabel,
        de = de,
        de.tr = treat_den,
        hist.out = hist.out,
        count.tr = treat.hist,
        estimator = "DML"
    )

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
