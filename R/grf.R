interflex.grf <- function(data,
                          Y, # outcome
                          D, # treatment indicator
                          X, # moderator
                          treat.info,
                          diff.info,
                          Z = NULL, # covariates
                          weights = NULL, # weighting variable
                          gate = FALSE,
                          num.trees = 4000,
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
    covariates <- c(X, Z)
    length.covariates <- length(covariates)

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

    TE.output.all.list <- list()
    if (treat.type == "discrete") {
        for (char in other.treat) {
            data_part <- data[data[[D]] %in% c(treat.base, char), ]
            data_part[data_part[[D]] == treat.base, D] <- 0L
            data_part[data_part[[D]] == char, D] <- 1L
            data_part[[D]] <- as.numeric(data_part[[D]])
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

    # GATE estimation
    if (isTRUE(gate)) {
        gate.list <- list()

        for (char in other.treat) {
            data_part <- data[data[[D]] %in% c(treat.base, char), ]
            data_part[data_part[[D]] == treat.base, D] <- 0L
            data_part[data_part[[D]] == char, D] <- 1L
            data_part$D_num <- as.numeric(data_part[[D]])

            # Fit causal forest (reuse same specification as the smooth CATE)
            cf <- causal_forest(
                data_part[covariates],
                data_part[[Y]],
                data_part$D_num,
                num.trees = num.trees
            )

            # Individual CATEs
            cate_i <- predict(cf)$predictions
            var_i <- predict(cf, estimate.variance = TRUE)$variance.estimates

            # Aggregate by X group
            X_vals <- data_part[[X]]
            groups <- sort(unique(X_vals))

            gate_est <- numeric(length(groups))
            gate_se <- numeric(length(groups))
            for (k in seq_along(groups)) {
                idx <- which(X_vals == groups[k])
                gate_est[k] <- mean(cate_i[idx])
                # SE via delta method: SE(mean) = sqrt(mean(var_i) / n_group)
                gate_se[k] <- sqrt(mean(var_i[idx]) / length(idx))
            }

            res <- data.frame(
                X                       = groups,
                ME                      = gate_est,
                sd                      = gate_se,
                `lower CI(95%)`         = gate_est - 1.96 * gate_se,
                `upper CI(95%)`         = gate_est + 1.96 * gate_se,
                check.names = FALSE
            )

            gate.list[[other.treat.origin[char]]] <- res
        }

        final.output$g.est <- gate.list
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
            show.uniform.CI = show.uniform.CI,
            scale = scale,
            height = height,
            width = width
        )
        final.output <- c(final.output, list(figure = figure.output))
    }

    class(final.output) <- "interflex"
    return(final.output)
}
