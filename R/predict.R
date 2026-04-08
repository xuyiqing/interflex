predict.interflex <- function(
    object,
    type = "response", # "response" or "link"
    order = NULL,
    subtitles = NULL,
    show.subtitles = NULL,
    # plot
    Xdistr = "histogram", # can be "histogram","hist","density","none"
    CI = NULL,
    pool = FALSE,
    main = NULL,
    Ylabel = NULL,
    Xlabel = NULL,
    xlab = NULL,
    ylab = NULL,
    xlim = NULL,
    ylim = NULL,
    theme.bw = TRUE,
    show.grid = TRUE,
    cex.main = NULL,
    cex.sub = NULL,
    cex.lab = NULL,
    cex.axis = NULL,
    color = NULL,
    file = NULL,
    interval = NULL,
    legend.title = NULL,
    ncols = NULL,
    show.all = FALSE,
    scale = 1.1,
    height = 7,
    width = 10,
    ...) {
    X <- NULL
    EY <- NULL
    CI_lower <- NULL
    CI_upper <- NULL
    x <- NULL
    y <- NULL
    xmin <- NULL
    xmax <- NULL
    ymin <- NULL
    count1 <- NULL
    ymax <- NULL
    Treatment <- NULL
    end_level <- NULL





    out <- object
    if (!inherits(out, "interflex")) {
        stop("Not an \"interflex\" object.")
    }

    #if (out$use.fe == 1) {
    #    return(0)
    #}

    treat.info <- out$treat.info
    treat.type <- treat.info[["treat.type"]]
    if (treat.type == "discrete") {
        other.treat <- names(treat.info[["other.treat"]])
        all.treat <- names(treat.info[["all.treat"]])
        base <- names(treat.info[["base"]])
    }
    if (treat.type == "continuous") {
        D.sample <- treat.info[["D.sample"]]
        label.name <- names(D.sample)
        all.treat <- label.name
    }

    de <- out$de
    de.tr <- out$de.tr
    hist.out <- out$hist.out
    count.tr <- out$count.tr
    estimator <- out$estimator

    if (!is.null(show.subtitles)) {
        if (!is.logical(show.subtitles) & !is.numeric(show.subtitles)) {
            stop("\"show.subtitles\" is not a logical flag.")
        }
    }

    # CI
    if (!is.null(CI)) {
        if (!is.logical(CI) & !is.numeric(CI)) {
            stop("\"CI\" is not a logical flag.")
        }

        if (estimator == "kernel") {
            if (isTRUE(CI) & isFALSE(out$CI)) {
                stop("Please set CI to FALSE.")
            }
        }
    }

    if (estimator == "kernel") {
        if (is.null(CI)) {
            CI <- out$CI
        }
    }

    if (estimator == "binning" | estimator == "linear") {
        if (is.null(CI)) {
            CI <- TRUE
        }
    }

    # pool
    if (!is.logical(pool) & !is.numeric(pool)) {
        stop("\"pool\" is not a logical flag.")
    }

    # main
    if (!is.null(main)) {
        main <- as.character(main)[1]
    }

    # Ylabel
    if (is.null(Ylabel)) {
        Ylabel <- out$Ylabel
    } else {
        if (!is.character(Ylabel)) {
            stop("\"Ylabel\" is not a string.")
        } else {
            Ylabel <- Ylabel[1]
        }
    }

    # Xlabel
    if (is.null(Xlabel)) {
        Xlabel <- out$Xlabel
    } else {
        if (!is.character(Xlabel)) {
            stop("\"Xlabel\" is not a string.")
        } else {
            Xlabel <- Xlabel[1]
        }
    }

    ## axis labels
    if (is.null(xlab)) {
        xlab <- c(paste("Moderator: ", Xlabel, sep = ""))
    } else {
        if (!is.character(xlab)) {
            stop("\"xlab\" is not a string.")
        }
    }

    if (is.null(ylab)) {
        if (type != "link") {
            ylab <- c(paste("Expected Value of", Ylabel, sep = " "))
        } else {
            ylab <- c(paste("link of", Ylabel, sep = " "))
        }
    } else {
        if (!is.character(ylab)) {
            stop("\"ylab\" is not a string.")
        }
    }

    ## xlim ylim
    if (!is.null(xlim)) {
        if (!is.numeric(xlim)) {
            stop("Some element in \"xlim\" is not numeric.")
        } else {
            if (length(xlim) != 2) {
                stop("\"xlim\" must be of length 2.")
            }
        }
    }

    if (!is.null(ylim)) {
        if (!is.numeric(ylim)) {
            stop("Some element in \"ylim\" is not numeric.")
        } else {
            if (length(ylim) != 2) {
                stop("\"ylim\" must be of length 2.")
            }
        }
    }

    ## theme.bw
    if (!is.logical(theme.bw) & !is.numeric(theme.bw)) {
        stop("\"theme.bw\" is not a logical flag.")
    }

    ## show.grid
    if (!is.logical(show.grid) & !is.numeric(show.grid)) {
        stop("\"show.grid\" is not a logical flag.")
    }

    ## font size
    if (!is.null(cex.main)) {
        if (!is.numeric(cex.main)) {
            stop("\"cex.main\" is not numeric.")
        }
    }
    if (!is.null(cex.sub)) {
        if (!is.numeric(cex.sub)) {
            stop("\"cex.sub\" is not numeric.")
        }
    }
    if (!is.null(cex.lab)) {
        if (!is.numeric(cex.lab)) {
            stop("\"cex.lab\" is not numeric.")
        }
    }
    if (!is.null(cex.axis)) {
        if (!is.numeric(cex.axis)) {
            stop("\"cex.axis\" is not numeric.")
        }
    }

    ## color
    if (!is.null(color)) {
        color <- as.character(color)
        color.in <- c()
        for (char in color) {
            res <- try(col2rgb(char), silent = TRUE)
            if (!"try-error" %in% class(res)) {
                color.in <- c(color.in, char)
            } else {
                stop(paste0(char, " is not one name for a color.\n"))
            }
        }
        color <- color.in
    }

    # file
    if (!is.null(file)) {
        if (!is.character(file)) {
            stop("Wrong file name.")
        }
    }

    # interval
    if (!is.null(interval)) {
        if (!is.numeric(interval)) {
            stop("Some element in \"interval\" is not numeric.")
        }
    }

    ## legend.title
    if (!is.null(legend.title)) {
        legend.title <- as.character(legend.title)[1]
    }

    ## type
    if (!type %in% c("response", "link")) {
        stop("\"type\" should be \"response\" or \"link\".")
    }

    ## ncols
    if (!is.null(ncols)) {
        if (!is.numeric(ncols)) {
            stop("\"ncols\" is not a positive integer.")
        } else {
            ncols <- ncols[1]
            if (ncols %% 1 != 0 | ncols < 1) {
                stop("\"ncols\" is not a positive number.")
            }
        }
    }

    if (treat.type == "discrete") {
        if (!is.null(order)) {
            order <- as.character(order)
        }

        if (!is.null(order)) {
            if (length(order) != length(unique(order))) {
                stop("\"order\" should not contain repeated values z.")
            }

            if (length(order) != length(all.treat)) {
                stop("\"order\" should contain all treatment arms(or treatment arms specified in \"to.show\").")
            }

            if (sum(!is.element(order, all.treat)) != 0 | sum(!is.element(all.treat, order)) != 0) {
                stop("\"order\" should contain all treatment arms(or treatment arms specified in \"to.show\").")
            }
            all.treat <- order
        }

        if (is.null(show.subtitles)) {
            show.subtitles <- TRUE
        }

        if (!is.null(subtitles)) {
            if (length(subtitles) != length(all.treat)) {
                stop("The length of \"subtitles\" should be equal to the number of treatment arms(or treatment arms specified in \"to.show\").")
            }
        }
    }

    if (treat.type == "continuous") {
        if (!is.null(order)) {
            if (!is.numeric(order)) {
                stop("\"order\" should be numeric.")
            }
            if (length(order) != length(unique(order))) {
                stop("\"order\" should not contain repeated values.")
            }
            if (length(order) != length(D.sample)) {
                stop("\"order\" should contain all reference values of D.")
            }
            if (sum(!is.element(order, D.sample)) != 0 | sum(!is.element(D.sample, order)) != 0) {
                stop("\"order\" should contain all reference values of D.")
            }
            all.treat <- label.name
            label.name.order <- c()
            for (a in order) {
                label.name.order <- c(label.name.order, names(D.sample[which(D.sample == a)]))
            }
            all.treat <- label.name.order
        }

        if (is.null(show.subtitles)) {
            show.subtitles <- TRUE
        }

        if (!is.null(subtitles)) {
            if (length(subtitles) != length(all.treat)) {
                stop("The number of elements in \"subtitles\" should be equal to the number of reference values of treatment(D).")
            }
        }
    }

    ## ncols
    ntreat <- length(all.treat)
    if (!is.null(ncols)) {
        if (ncols %% 1 != 0) {
            stop("\"ncols\" is not a positive integer.")
        } else {
            ncols <- ncols[1]
        }
        if (ncols < 1) {
            stop("\"ncols\" is not a positive integer.")
        }
    } else {
        ncols <- length(all.treat)
    }

    if (pool) {
        ncols <- 1
    }

    ## axis labels
    if (is.null(cex.lab)) {
        cex.lab <- 12
    } else {
        cex.lab <- 12 * cex.lab
    }

    if (is.null(cex.axis)) {
        cex.axis <- 12
    } else {
        cex.axis <- 12 * cex.axis
    }


    if (is.null(cex.main)) {
        cex.main <- 18
    } else {
        cex.main <- 18 * cex.main
    }

    if (is.null(cex.sub)) {
        cex.sub <- 10
    } else {
        cex.sub <- 10 * cex.sub
    }


    ## xlim and ylim
    if (!is.null(ylim)) {
        ylim2 <- c(ylim[1] - (ylim[2] - ylim[1]) * 0.25 / 6, ylim[2] + (ylim[2] - ylim[1]) * 0.4 / 6)
    }

    ## color
    requireNamespace("RColorBrewer")
    platte <- rep("#999999", ntreat)
    platte <- c(brewer.pal(n = 8, "Dark2"), platte)
    if (!is.null(color)) {
        platte <- c(color, platte)
    }


    ## yrange
    if (estimator == "kernel") {
        pred <- out$pred.kernel
        if (type == "link") {
            pred <- out$link.kernel
        }
    }
    if (estimator == "binning") {
        pred <- out$pred.bin
        if (type == "link") {
            pred <- out$link.bin
        }
    }
    if (estimator == "linear") {
        pred <- out$pred.lin
        if (type == "link") {
            pred <- out$link.lin
        }
    }



    yrange <- c(0)
    for (char in all.treat) {
        temp.pred <- pred[[char]]
        if (isTRUE(CI)) {
            if(dim(temp.pred)[2]==5){
                colnames(temp.pred) <- c("X", "EY", "sd", "CI_lower", "CI_upper")
                yrange <- c(yrange, na.omit(unlist(temp.pred[, c(4, 5)])))                
            }else{
                colnames(temp.pred) <- c("X", "EY", "sd", "CI_lower", "CI_upper", "CI_uniform_lower", "CI_uniform_upper")
                yrange <- c(yrange, na.omit(unlist(temp.pred[, c(6, 7)])))                 
            }

            pred[[char]] <- as.data.frame(temp.pred)
        }
        if (isFALSE(CI)) {
            temp.pred <- temp.pred[, c(1, 2)]
            colnames(temp.pred) <- c("X", "EY")
            yrange <- c(yrange, na.omit(unlist(temp.pred[, c(2)])))
            pred[[char]] <- as.data.frame(temp.pred)
        }
    }
    if (!is.null(ylim)) {
        yrange <- c(ylim[2], ylim[1] + (ylim[2] - ylim[1]) * 1 / 8)
    }
    maxdiff <- (max(yrange) - min(yrange))
    pos <- max(yrange) - maxdiff / 20

    if (!pool) {
        plot_list <- list()
        k <- 1
        for (char in all.treat) {
            p <- ggplot()
            if (theme.bw) {
                p <- p + theme_bw()
            }
            if (!show.grid) {
                p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
            }
            p <- p + geom_line(data = pred[[char]], aes(x = X, y = EY), linewidth = 1, color = platte[k])
            if (isTRUE(CI)) {
                p <- p + geom_ribbon(
                    data = pred[[char]], aes(x = X, ymin = CI_lower, ymax = CI_upper),
                    alpha = 0.3, fill = platte[k]
                )
                if("CI_uniform_lower" %in% colnames(pred[[char]])){
                    p <- p + geom_ribbon(
                        data = pred[[char]], aes(x = X, ymin = CI_uniform_lower, ymax = CI_uniform_upper),
                        alpha = 0, color = platte[k], linetype = 2
                    )
                }
            }

            if (treat.type == "discrete") {
                if (show.subtitles) {
                    if (is.null(subtitles)) {
                        p <- p + labs(subtitle = paste0("Group:", char)) + theme(plot.subtitle = element_text(hjust = 0.5, size = cex.sub, lineheight = .8))
                    }
                    if (!is.null(subtitles)) {
                        p <- p + labs(subtitle = subtitles[k]) + theme(plot.subtitle = element_text(hjust = 0.5, size = cex.sub, lineheight = .8))
                    }
                }

                if (Xdistr == "density") {
                    deX.ymin <- min(yrange) - maxdiff / 5
                    deX.tr <- data.frame(
                        x = de.tr[[char]]$x,
                        y = de.tr[[char]]$y / max(de.tr[[char]]$y) * maxdiff / 5 + min(yrange) - maxdiff / 5
                    )
                    p <- p + geom_ribbon(
                        data = deX.tr, aes(x = x, ymax = y, ymin = deX.ymin),
                        fill = platte[k], alpha = 0.2
                    )
                }

                if (Xdistr %in% c("histogram", "hist")) {
                    n.hist <- length(hist.out$mids)
                    dist <- hist.out$mids[2] - hist.out$mids[1]
                    hist.max <- max(hist.out$counts)
                    hist.treat <- data.frame(
                        ymin = min(yrange) - maxdiff / 5,
                        # ymax=hist.out$counts/hist.max*maxdiff/5+min(yrange)-maxdiff/5,
                        xmin = hist.out$mids - dist / 2,
                        xmax = hist.out$mids + dist / 2,
                        count1 = count.tr[[char]] / hist.max * maxdiff / 5 + min(yrange) - maxdiff / 5
                    )

                    p <- p + geom_rect(
                        data = hist.treat, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = count1),
                        fill = platte[k], colour = "gray50", alpha = 0.3, linewidth = 0.3
                    )
                }
            }

            if (treat.type == "continuous") {
                if (show.subtitles) {
                    if (is.null(subtitles)) {
                        p <- p + labs(subtitle = char) + theme(plot.subtitle = element_text(hjust = 0.5, size = cex.sub, lineheight = .8))
                    }
                    if (!is.null(subtitles)) {
                        p <- p + labs(subtitle = subtitles[k]) + theme(plot.subtitle = element_text(hjust = 0.5, size = cex.sub, lineheight = .8))
                    }
                }

                if (Xdistr == "density") {
                    deX.ymin <- min(yrange) - maxdiff / 5
                    deX.tr <- data.frame(
                        x = de$x,
                        y = de$y / max(de$y) * maxdiff / 5 + min(yrange) - maxdiff / 5
                    )
                    p <- p + geom_ribbon(
                        data = deX.tr, aes(x = x, ymax = y, ymin = deX.ymin),
                        fill = "gray50", alpha = 0.2
                    )
                }

                if (Xdistr %in% c("histogram", "hist")) {
                    n.hist <- length(hist.out$mids)
                    dist <- hist.out$mids[2] - hist.out$mids[1]
                    hist.max <- max(hist.out$counts)
                    histX <- data.frame(
                        ymin = rep(min(yrange) - maxdiff / 5, n.hist),
                        ymax = hist.out$counts / hist.max * maxdiff / 5 + min(yrange) - maxdiff / 5,
                        xmin = hist.out$mids - dist / 2,
                        xmax = hist.out$mids + dist / 2
                    )

                    p <- p + geom_rect(
                        data = histX, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                        colour = "gray50", fill = "gray50", alpha = 0.3, linewidth = 0.5
                    )
                }
            }
            k <- k + 1
            plot_list[[char]] <- p
        }
    }

    # pool=TRUE
    if (pool) {
        for (char in all.treat) {
            if (char == all.treat[1]) {
                tograph <- pred[[char]]
                tograph[, "Treatment"] <- char
            } else {
                tograph1 <- pred[[char]]
                tograph1[, "Treatment"] <- char
                tograph <- rbind(tograph, tograph1)
            }
        }
        p <- ggplot()
        if (theme.bw) {
            p <- p + theme_bw()
        }
        if (!show.grid) {
            p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        }

        if (treat.type == "discrete") {
            tograph$Treatment <- factor(tograph$Treatment, levels = all.treat)
            p <- p + geom_line(data = tograph, aes(x = X, y = EY, color = Treatment), linewidth = 0.5)
            if (is.null(subtitles)) {
                p <- p + scale_color_manual(values = platte[1:ntreat])
            } else {
                p <- p + scale_color_manual(values = platte[1:ntreat], labels = subtitles)
            }
        }

        if (treat.type == "continuous") {
            tograph$Treatment <- factor(tograph$Treatment, levels = all.treat)
            p <- p + geom_line(data = tograph, aes(x = X, y = EY, color = Treatment), linewidth = 0.5)
            if (is.null(subtitles)) {
                p <- p + scale_color_manual(values = platte[1:ntreat], labels = label.name)
            } else {
                p <- p + scale_color_manual(values = platte[1:ntreat], labels = subtitles)
            }
        }


        if (isTRUE(CI)) {
            p <- p + geom_ribbon(
                data = tograph, aes(x = X, ymin = CI_lower, ymax = CI_upper, fill = Treatment),
                alpha = 0.2, show.legend = T, linewidth = 0
            )
            if("CI_uniform_lower" %in% colnames(tograph)){
                p <- p + geom_ribbon(
                        data = tograph, aes(x = X, ymin = CI_uniform_lower, ymax = CI_uniform_upper, color = Treatment),
                        alpha = 0, linetype = 2
                    )
            }
            if (treat.type == "discrete") {
                if (is.null(subtitles)) {
                    p <- p + scale_fill_manual(values = platte[1:ntreat])
                } else {
                    p <- p + scale_fill_manual(values = platte[1:ntreat], labels = subtitles)
                }
                if("CI_uniform_lower" %in% colnames(tograph)){
                    if (is.null(subtitles)) {
                        p <- p + scale_color_manual(values = platte[1:ntreat])
                    } else {
                        p <- p + scale_color_manual(values = platte[1:ntreat], labels = subtitles)
                    }
                }
            }

            if (treat.type == "continuous") {
                if (is.null(subtitles)) {
                    p <- p + scale_fill_manual(values = platte[1:ntreat], labels = label.name)
                } else {
                    p <- p + scale_fill_manual(values = platte[1:ntreat], labels = subtitles)
                }
                if("CI_uniform_lower" %in% colnames(tograph)){
                    if (is.null(subtitles)) {
                        p <- p + scale_color_manual(values = platte[1:ntreat], labels = label.name)
                    } else {
                        p <- p + scale_color_manual(values = platte[1:ntreat], labels = subtitles)
                    }
                }
            }
        }

        if (!is.null(legend.title)) {
            p <- p + labs(fill = legend.title, color = legend.title)
        }

        if (Xdistr == "density" & treat.type == "discrete") { # density plot
            deX.ymin <- min(yrange) - maxdiff / 5
            k <- 1
            char0 <- all.treat[1]
            start_level <- rep(deX.ymin, length(de.tr[[char0]]$x))
            for (char in all.treat) {
                dex.tr.plot <- data.frame(
                    x = de.tr[[char]]$x,
                    start_level = start_level,
                    end_level = de.tr[[char]]$y / max(de.tr[[char0]]$y) * maxdiff / 10 + start_level
                )

                p <- p + geom_ribbon(
                    data = dex.tr.plot, aes(x = x, ymax = end_level, ymin = start_level), color = platte[k],
                    alpha = 0.0, fill = platte[k], linewidth = 0.3
                )

                k <- k + 1
            }
            p <- p + geom_line(data = dex.tr.plot, aes(x = x, y = min(yrange) - maxdiff / 5), color = "gray50", linewidth = 0.3)
        }

        if (Xdistr == "density" & treat.type == "continuous") {
            deX.ymin <- min(yrange) - maxdiff / 5
            deX.tr <- data.frame(
                x = de$x,
                y = de$y / max(de$y) * maxdiff / 5 + min(yrange) - maxdiff / 5
            )
            p <- p + geom_ribbon(
                data = deX.tr, aes(x = x, ymax = y, ymin = deX.ymin),
                fill = "gray50", alpha = 0.2
            )
        }

        if (Xdistr %in% c("histogram", "hist") & treat.type == "discrete") { # density plot
            deX.ymin <- min(yrange) - maxdiff / 5
            n.hist <- length(hist.out$mids)
            dist <- hist.out$mids[2] - hist.out$mids[1]
            hist.max <- max(hist.out$counts)
            k <- 1
            start_level <- min(yrange) - maxdiff / 5
            for (char in all.treat) {
                hist.treat <- data.frame(
                    ymin = start_level,
                    ymax = count.tr[[char]] / hist.max * maxdiff / 5 + start_level,
                    xmin = hist.out$mids - dist / 2,
                    xmax = hist.out$mids + dist / 2
                )

                start_level <- count.tr[[char]] / hist.max * maxdiff / 5 + start_level

                p <- p + geom_rect(
                    data = hist.treat, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = platte[k], color = "gray50",
                    alpha = 0.3, linewidth = 0.2
                )
                k <- k + 1
            }
        }

        if (Xdistr %in% c("histogram", "hist") & treat.type == "continuous") {
            n.hist <- length(hist.out$mids)
            dist <- hist.out$mids[2] - hist.out$mids[1]
            hist.max <- max(hist.out$counts)
            histX <- data.frame(
                ymin = rep(min(yrange) - maxdiff / 5, n.hist),
                ymax = hist.out$counts / hist.max * maxdiff / 5 + min(yrange) - maxdiff / 5,
                xmin = hist.out$mids - dist / 2,
                xmax = hist.out$mids + dist / 2
            )

            p <- p + geom_rect(
                data = histX, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                colour = "gray50", fill = "gray50", alpha = 0.3, size = 0.5
            )
        }

        plot_list <- p
    }

    if (!pool) {
        for (char in all.treat) {
            p1 <- plot_list[[char]]
            ## mark the original interval (in replicated papers)
            if (!is.null(interval)) {
                p1 <- p1 + geom_vline(xintercept = interval, colour = "steelblue", linetype = 2, linewidth = 1.5)
            }
            ## Other universal options
            p1 <- p1 + xlab(NULL) + ylab(NULL) +
                theme(axis.text = element_text(size = cex.axis))

            plot_list[[char]] <- p1
        }

        ## ylim
        for (char in all.treat) {
            y.limits <- layer_scales(plot_list[[char]])$y$range$range
            if (char == all.treat[1]) {
                ymaxmax <- y.limits[2]
                yminmin <- y.limits[1]
            } else {
                ymaxmax <- max(ymaxmax, y.limits[2])
                yminmin <- min(yminmin, y.limits[1])
            }
        }
        for (char in all.treat) {
            final_ylim <- if (!is.null(ylim)) ylim2 else c(yminmin, ymaxmax)
            final_xlim <- if (!is.null(xlim)) .pad_xlim(xlim) else NULL
            plot_list[[char]] <- plot_list[[char]] +
                coord_cartesian(xlim = final_xlim, ylim = final_ylim)
        }

        requireNamespace("gridExtra")

        suppressMessages(
            graph <- arrangeGrob(
                grobs = plot_list, ncol = ncols, align = "h",
                top = textGrob(main, gp = gpar(fontsize = cex.main, face = "bold")),
                left = textGrob(ylab, rot = 90, vjust = 1, gp = gpar(fontsize = cex.lab)),
                bottom = textGrob(xlab, gp = gpar(fontsize = cex.lab))
            )
        )
    }

    if (pool) {
        p1 <- plot_list
        if (!is.null(interval)) {
            p1 <- p1 + geom_vline(xintercept = interval, colour = "steelblue", linetype = 2, linewidth = 1.5)
        }
        ## Other universal options
        p1 <- p1 + xlab(xlab) + ylab(ylab) +
            theme(axis.text = element_text(size = cex.axis), axis.title = element_text(size = cex.lab))

        if (!is.null(main)) {
            p1 <- p1 + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5, size = cex.main, lineheight = .8, face = "bold"))
        }
        if (!is.null(xlim) & !is.null(ylim)) {
            p1 <- p1 + coord_cartesian(xlim = .pad_xlim(xlim), ylim = ylim2)
        }
        if (is.null(xlim) & !is.null(ylim)) {
            p1 <- p1 + coord_cartesian(ylim = ylim2)
        }
        if (!is.null(xlim) & is.null(ylim)) {
            p1 <- p1 + coord_cartesian(xlim = .pad_xlim(xlim))
        }

        p1 <- p1 + theme(
            legend.title = element_text(colour = "black", size = cex.sub),
            legend.text = element_text(color = "black", size = cex.sub * 0.95)
        )

        graph <- p1
    }


    if (!show.all & !pool) {
        graph <- as.ggplot(graph)
        ## save to file
        if (!is.null(file)) {
            ggsave(file, graph, scale = scale, width = width, height = height)
        }
        return(graph)
    }

    if (pool) {
        ## save to file
        if (!is.null(file)) {
            ggsave(file, graph, scale = scale, width = width, height = height)
        }
        return(graph)
    }
    if (show.all & !pool) {
        cex.axis <- cex.axis
        for (char in all.treat) {
            ## Other universal options
            p1 <- plot_list[[char]] + xlab(xlab) + ylab(ylab) +
                theme(
                    axis.text = element_text(size = cex.axis),
                    axis.title.x = element_text(size = cex.lab),
                    axis.title.y = element_text(size = cex.lab)
                )
            p1 <- p1 + ggtitle(main) +
                theme(plot.title = element_text(hjust = 0.5, size = cex.main, lineheight = .8, face = "bold"))

            plot_list[[char]] <- p1
        }
        return(plot_list)
    }
}
