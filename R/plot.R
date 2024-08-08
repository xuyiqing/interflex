plot.interflex <- function(x,
                           order = NULL,
                           subtitles = NULL,
                           show.subtitles = NULL,
                           CI = NULL,
                           diff.values = NULL,
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
                           bin.labs = TRUE, # bin labels
                           interval = NULL, # interval in replicated papers
                           file = NULL,
                           ncols = NULL,
                           # pool plot
                           pool = FALSE,
                           legend.title = NULL,
                           line.size = 1.2,
                           color = NULL,
                           line.color = "black",
                           CI.color = "black",
                           CI.color.alpha = 0.2,
                           hist.color = c("gray50", "red"),
                           hist.color.alpha = 0.3,
                           density.color = c("gray50", "red"),
                           density.color.alpha = 0.3,
                           show.all = FALSE,
                           scale = 1.1,
                           height = 7,
                           width = 10,
                           ...) {
    y <- NULL
    xmin <- NULL
    xmax <- NULL
    ymin <- NULL
    count1 <- NULL
    ymax <- NULL
    X <- NULL
    TE <- NULL
    CI_lower <- NULL
    CI_upper <- NULL
    ME <- NULL
    x0 <- NULL

    if (pool == TRUE) {
        p <- interflex.plot.pool(
            out = x,
            diff.values = diff.values,
            order = order,
            subtitles = subtitles,
            show.subtitles = show.subtitles,
            legend.title = legend.title,
            CI = CI,
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
            cex.lab = cex.lab,
            cex.axis = cex.axis,
            cex.sub = cex.sub,
            bin.labs = bin.labs, # bin labels
            interval = interval, # interval in replicated papers
            line.size = line.size,
            color = color,
            line.color = line.color,
            CI.color = CI.color,
            CI.color.alpha = CI.color.alpha,
            hist.color = hist.color,
            hist.color.alpha = hist.color.alpha,
            density.color = density.color,
            density.color.alpha = density.color.alpha,
            file = file
        )
        return(p)
    }


    out <- x
    if (!class(out) %in% c("interflex")) {
        stop("Not an \"interflex\" object.")
    }

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
    }

    de <- out$de
    de.tr <- out$de.tr
    hist.out <- out$hist.out
    count.tr <- out$count.tr
    estimator <- out$estimator

    if (is.null(show.subtitles) == FALSE) {
        if (is.logical(show.subtitles) == FALSE & is.numeric(show.subtitles) == FALSE) {
            stop("\"show.subtitles\" is not a logical flag.")
        }
    }

    # CI
    if (is.null(CI) == FALSE) {
        if (is.logical(CI) == FALSE & is.numeric(CI) == FALSE) {
            stop("\"CI\" is not a logical flag.")
        }

        if (estimator == "kernel") {
            if (CI == TRUE & out$CI == FALSE) {
                stop("Please set CI to FALSE.")
            }
        }
    }

    if (estimator == "kernel") {
        if (is.null(CI) == TRUE) {
            CI <- out$CI
        }
    }

    if (estimator == "binning" | estimator == "linear" | estimator == "DML") {
        if (is.null(CI) == TRUE) {
            CI <- TRUE
        }
    }

    # Xdistr
    if (!Xdistr %in% c("hist", "histogram", "density", "none")) {
        stop("\"Xdistr\" must be \"histogram\", \"density\", or \"none\".")
    }

    # main
    if (is.null(main) == FALSE) {
        main <- as.character(main)[1]
    }

    # Ylabel
    if (is.null(Ylabel) == TRUE) {
        Ylabel <- out$Ylabel
    } else {
        if (is.character(Ylabel) == FALSE) {
            stop("\"Ylabel\" is not a string.")
        } else {
            Ylabel <- Ylabel[1]
        }
    }

    # Dlabel
    if (is.null(Dlabel) == TRUE) {
        Dlabel <- out$Dlabel
    } else {
        if (is.character(Dlabel) == FALSE) {
            stop("\"Dlabel\" is not a string.")
        } else {
            Dlabel <- Dlabel[1]
        }
    }

    # Xlabel
    if (is.null(Xlabel) == TRUE) {
        Xlabel <- out$Xlabel
    } else {
        if (is.character(Xlabel) == FALSE) {
            stop("\"Xlabel\" is not a string.")
        } else {
            Xlabel <- Xlabel[1]
        }
    }

    ## axis labels
    if (is.null(xlab) == FALSE) {
        if (is.character(xlab) == FALSE) {
            stop("\"xlab\" is not a string.")
        }
    }
    if (is.null(ylab) == FALSE) {
        if (is.character(ylab) == FALSE) {
            stop("\"ylab\" is not a string.")
        }
    }

    if (is.null(xlab) == TRUE) {
        xlab <- c(paste("Moderator: ", Xlabel, sep = ""))
    } else {
        if (is.character(xlab) == FALSE) {
            stop("\"xlab\" is not a string.")
        }
    }
    if (is.null(ylab) == TRUE) {
        ylab <- c(paste("Marginal Effect of ", Dlabel, " on ", Ylabel, sep = ""))
    } else {
        if (is.character(ylab) == FALSE) {
            stop("\"ylab\" is not a string.")
        }
    }

    ## xlim ylim
    if (is.null(xlim) == FALSE) {
        if (is.numeric(xlim) == FALSE) {
            stop("Some element in \"xlim\" is not numeric.")
        } else {
            if (length(xlim) != 2) {
                stop("\"xlim\" must be of length 2.")
            }
        }
    }

    if (is.null(ylim) == FALSE) {
        if (is.numeric(ylim) == FALSE) {
            stop("Some element in \"ylim\" is not numeric.")
        } else {
            if (length(ylim) != 2) {
                stop("\"ylim\" must be of length 2.")
            }
        }
    }

    ## theme.bw
    if (is.logical(theme.bw) == FALSE & is.numeric(theme.bw) == FALSE) {
        stop("\"theme.bw\" is not a logical flag.")
    }

    ## show.grid
    if (is.logical(show.grid) == FALSE & is.numeric(show.grid) == FALSE) {
        stop("\"show.grid\" is not a logical flag.")
    }

    ## font size
    if (is.null(cex.main) == FALSE) {
        if (is.numeric(cex.main) == FALSE) {
            stop("\"cex.main\" is not numeric.")
        }
    }
    if (is.null(cex.sub) == FALSE) {
        if (is.numeric(cex.sub) == FALSE) {
            stop("\"cex.sub\" is not numeric.")
        }
    }
    if (is.null(cex.lab) == FALSE) {
        if (is.numeric(cex.lab) == FALSE) {
            stop("\"cex.lab\" is not numeric.")
        }
    }
    if (is.null(cex.axis) == FALSE) {
        if (is.numeric(cex.axis) == FALSE) {
            stop("\"cex.axis\" is not numeric.")
        }
    }

    ## bin.labs
    if (is.logical(bin.labs) == FALSE & is.numeric(bin.labs) == FALSE) {
        stop("\"bin.labs\" is not a logical flag.")
    }

    ## interval
    if (is.null(interval) == FALSE) {
        if (is.numeric(interval) == FALSE) {
            stop("Some element in \"interval\" is not numeric.")
        }
    }

    ## file
    if (is.null(file) == FALSE) {
        if (is.character(file) == FALSE) {
            stop("Wrong file name.")
        }
    }

    ## order/subtitles
    if (treat.type == "discrete") {
        other.treat <- sort(all.treat[which(all.treat != base)])
        if (is.null(order) == FALSE) {
            order <- as.character(order)
            if (length(order) != length(unique(order))) {
                stop("\"order\" should not contain repeated values.")
            }

            if (length(order) != length(other.treat)) {
                stop("\"order\" should include all kinds of treatment arms except for the baseline group.")
            }

            if (sum(!is.element(order, other.treat)) != 0 | sum(!is.element(other.treat, order)) != 0) {
                stop("\"order\" should include all kinds of treatment arms except for the baseline group.")
            }
            other.treat <- order
        }

        if (is.null(show.subtitles) == TRUE) {
            if (length(other.treat) == 1) {
                show.subtitles <- FALSE
            } else {
                show.subtitles <- TRUE
            }

            if (is.null(subtitles) == FALSE) {
                show.subtitles <- TRUE
            }
        }

        if (is.null(subtitles) == FALSE) {
            if (length(subtitles) != length(other.treat)) {
                stop("The number of elements in \"subtitles\" should be m-1(m is the number of different treatment arms).")
            }
        }
    }

    if (treat.type == "continuous") {
        if (is.null(order) == FALSE) {
            if (is.numeric(order) == FALSE) {
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
            label.name.order <- c()
            for (a in order) {
                label.name.order <- c(label.name.order, names(D.sample[which(D.sample == a)]))
            }
            label.name <- label.name.order
        }

        if (is.null(show.subtitles) == TRUE) {
            if (length(label.name) == 1) {
                show.subtitles <- FALSE
            } else {
                show.subtitles <- TRUE
            }

            if (is.null(subtitles) == FALSE) {
                show.subtitles <- TRUE
            }
        }

        if (is.null(subtitles) == FALSE) {
            if (length(subtitles) != length(label.name)) {
                stop("The number of elements in \"subtitles\" should equal to the number of values in D.ref.")
            }
        }
    }

    ## ncols
    if (is.null(ncols) == FALSE) {
        if (ncols %% 1 != 0) {
            stop("\"ncols\" is not a positive integer.")
        } else {
            ncols <- ncols[1]
        }
        if (ncols < 1) {
            stop("\"ncols\" is not a positive integer.")
        }
    } else {
        if (treat.type == "discrete") {
            ncols <- length(other.treat)
        }
        if (treat.type == "continuous") {
            ncols <- length(label.name)
        }
    }
    treat_sc <- max(1, ncols - 1)

    ## diff.values
    if (treat.type == "discrete" & (estimator == "linear" | estimator == "binning")) {
        tempxx <- out$est.lin[[other.treat[1]]][, "X"]
    }
    if (treat.type == "discrete" & estimator == "DML") {
        tempxx <- out$est.dml[[other.treat[1]]][, "X"]
    }
    if (treat.type == "discrete" & estimator == "kernel") {
        tempxx <- out$est.kernel[[other.treat[1]]][, "X"]
    }
    if (treat.type == "continuous" & (estimator == "linear" | estimator == "binning")) {
        tempxx <- out$est.lin[[label.name[1]]][, "X"]
    }
    if (treat.type == "continuous" & estimator == "DML") {
        tempxx <- out$est.dml[[label.name[1]]][, "X"]
    }
    if (treat.type == "continuous" & estimator == "kernel") {
        tempxx <- out$est.kernel[[label.name[1]]][, "X"]
    }

    min.XX <- min(tempxx)
    max.XX <- max(tempxx)

    if (is.null(diff.values) == FALSE) {
        if (estimator == "binning") {
            stop("\"diff.values\" can only work after linear or kernel model is applied.")
        }
        if (is.numeric(diff.values) == FALSE) {
            stop("\"diff.values\" is not numeric.")
        }
        if (length(diff.values) < 2) {
            stop("\"diff.values\" must be of length 2 or more.")
        }
        if (treat.type == "discrete" & (estimator == "linear")) {
            tempxx <- out$est.lin[[other.treat[1]]][, "X"]
        }
        if (treat.type == "discrete" & (estimator == "DML")) {
            tempxx <- out$est.dml[[other.treat[1]]][, "X"]
        }
        if (treat.type == "discrete" & estimator == "kernel") {
            tempxx <- out$est.kernel[[other.treat[1]]][, "X"]
        }
        if (treat.type == "continuous" & estimator == "linear") {
            tempxx <- out$est.lin[[label.name[1]]][, "X"]
        }
        if (treat.type == "continuous" & estimator == "DML") {
            tempxx <- out$est.lin[[label.name[1]]][, "X"]
        }
        if (treat.type == "continuous" & estimator == "kernel") {
            tempxx <- out$est.kernel[[label.name[1]]][, "X"]
        }
        min.XX <- min(tempxx)
        max.XX <- max(tempxx)
        for (a in diff.values) {
            if (a < min.XX | a > max.XX) {
                stop("Elements in \"diff.values\" should be within the range of the moderator.")
            }
        }
    } else {
        if (estimator == "binning") {
            diff.values <- NULL
        } else {
            # diff.values <- out$diff.info[["diff.values.plot"]]
            diff.values <- NULL
        }
    }


    # yrange
    if (estimator == "binning") {
        nbins <- out$nbins
        if (treat.type == "discrete") {
            est.lin <- out$est.lin
            est.bin <- out$est.bin
            est.bin2 <- list() ## non missing part
            est.bin3 <- list() ## missing part
            yrange <- c(0)
            for (char in other.treat) {
                est.bin2[[char]] <- as.matrix(est.bin[[char]][which(is.na(est.bin[[char]][, 2]) == FALSE), ])
                est.bin3[[char]] <- as.matrix(est.bin[[char]][which(is.na(est.bin[[char]][, 2]) == TRUE), ])
                if (dim(est.bin2[[char]])[2] == 1) {
                    est.bin2[[char]] <- t(est.bin2[[char]])
                }
                if (dim(est.bin3[[char]])[2] == 1) {
                    est.bin3[[char]] <- t(est.bin3[[char]])
                }

                if (CI == TRUE) {
                    yrange <- c(yrange, na.omit(unlist(c(est.lin[[char]][, c(4, 5)], est.bin[[char]][, c(4, 5)]))))
                } else {
                    yrange <- c(yrange, na.omit(unlist(c(est.lin[[char]][, 2], est.bin[[char]][, 2]))))
                }
            }

            if (is.null(ylim) == FALSE) {
                yrange <- c(ylim[2], ylim[1] + (ylim[2] - ylim[1]) * 1 / 8)
            }
            X.lvls <- est.lin[[other.treat[1]]][, 1]
            errorbar.width <- (max(X.lvls) - min(X.lvls)) / 20
            maxdiff <- (max(yrange) - min(yrange))
            pos <- max(yrange) - maxdiff / 20
        }

        if (treat.type == "continuous") {
            est.lin <- out$est.lin
            est.bin <- out$est.bin
            est.bin2 <- list() ## non missing part
            est.bin3 <- list() ## missing part
            yrange <- c(0)
            for (label in label.name) {
                est.bin2[[label]] <- as.matrix(est.bin[[label]][which(is.na(est.bin[[label]][, 2]) == FALSE), ])
                est.bin3[[label]] <- as.matrix(est.bin[[label]][which(is.na(est.bin[[label]][, 2]) == TRUE), ])
                if (dim(est.bin2[[label]])[2] == 1) {
                    est.bin2[[label]] <- t(est.bin2[[label]])
                }
                if (dim(est.bin3[[label]])[2] == 1) {
                    est.bin3[[label]] <- t(est.bin3[[label]])
                }
                if (CI == TRUE) {
                    yrange <- c(yrange, na.omit(unlist(c(est.lin[[label]][, c(4, 5)], est.bin[[label]][, c(4, 5)]))))
                } else {
                    yrange <- c(yrange, na.omit(unlist(c(est.lin[[label]][, 2], est.bin[[label]][, 2]))))
                }
            }

            X.lvls <- est.lin[[label.name[1]]][, 1]
            errorbar.width <- (max(X.lvls) - min(X.lvls)) / 20
            if (is.null(ylim) == FALSE) {
                yrange <- c(ylim[2], ylim[1] + (ylim[2] - ylim[1]) * 1 / 8)
            }
            maxdiff <- (max(yrange) - min(yrange))
            pos <- max(yrange) - maxdiff / 20
        }
    }

    if (estimator == "DML") {
        if (treat.type == "discrete") {
            est.dml <- out$est.dml
            yrange <- c(0)
            for (char in other.treat) {
                if (CI == TRUE) {
                    yrange <- c(yrange, na.omit(unlist(c(est.dml[[char]][, c(4, 5)]))))
                    if(ncol(est.dml[[char]])>5){
                        yrange <- c(yrange, na.omit(unlist(c(est.dml[[char]][, c(6, 7)]))))
                    }
                } else {
                    yrange <- c(yrange, na.omit(unlist(c(est.dml[[char]][, 2]))))
                }
            }
            X.lvls <- est.dml[[other.treat[1]]][, 1]
        }
        if (treat.type == "continuous") {
            est.dml <- out$est.dml
            yrange <- c(0)
            for (label in label.name) {
                if (CI == TRUE) {
                    yrange <- c(yrange, na.omit(unlist(c(est.dml[[label]][, c(4, 5)]))))
                    if(ncol(est.dml[[label]])>5){
                        yrange <- c(yrange, na.omit(unlist(c(est.dml[[label]][, c(6, 7)]))))
                    }
                } else {
                    yrange <- c(yrange, na.omit(unlist(c(est.dml[[label]][, 2]))))
                }
            }
            X.lvls <- est.dml[[label.name[1]]][, 1]
        }
        errorbar.width <- (max(X.lvls) - min(X.lvls)) / 20
        if (is.null(ylim) == FALSE) {
            yrange <- c(ylim[2], ylim[1] + (ylim[2] - ylim[1]) * 1 / 8)
        }
        maxdiff <- (max(yrange) - min(yrange))
        pos <- max(yrange) - maxdiff / 20
    }

    if (estimator == "linear") {
        if (treat.type == "discrete") {
            est.lin <- out$est.lin
            yrange <- c(0)
            for (char in other.treat) {
                if (CI == TRUE) {
                    yrange <- c(yrange, na.omit(unlist(c(est.lin[[char]][, c(4, 5)]))))
                    if(ncol(est.lin[[char]])>5){
                        yrange <- c(yrange, na.omit(unlist(c(est.lin[[char]][, c(6, 7)]))))
                    }
                } else {
                    yrange <- c(yrange, na.omit(unlist(c(est.lin[[char]][, 2]))))
                }
            }
            X.lvls <- est.lin[[other.treat[1]]][, 1]
        }

        if (treat.type == "continuous") {
            est.lin <- out$est.lin
            yrange <- c(0)
            for (label in label.name) {
                if (CI == TRUE) {
                    yrange <- c(yrange, na.omit(unlist(c(est.lin[[label]][, c(4, 5)]))))
                    if(ncol(est.lin[[label]])>5){
                        yrange <- c(yrange, na.omit(unlist(c(est.lin[[label]][, c(6, 7)]))))
                    }
                } else {
                    yrange <- c(yrange, na.omit(unlist(c(est.lin[[label]][, 2]))))
                }
            }
            X.lvls <- est.lin[[label.name[1]]][, 1]
        }

        errorbar.width <- (max(X.lvls) - min(X.lvls)) / 20
        if (is.null(ylim) == FALSE) {
            yrange <- c(ylim[2], ylim[1] + (ylim[2] - ylim[1]) * 1 / 8)
        }
        maxdiff <- (max(yrange) - min(yrange))
        pos <- max(yrange) - maxdiff / 20
    }

    if (estimator == "kernel") {
        est.kernel <- out$est.kernel
        yrange <- c(0)
        if (CI == FALSE) {
            if (treat.type == "discrete") {
                for (char in other.treat) {
                    yrange <- c(yrange, na.omit(unlist(c(est.kernel[[char]][, 2]))))
                }
                X.lvls <- est.kernel[[other.treat[1]]][, 1]
            }

            if (treat.type == "continuous") {
                for (label in label.name) {
                    yrange <- c(yrange, na.omit(unlist(c(est.kernel[[label]][, 2]))))
                }
                X.lvls <- est.kernel[[label.name[1]]][, 1]
            }
        }

        if (CI == TRUE) {
            if (treat.type == "discrete") {
                for (char in other.treat) {
                    yrange <- c(yrange, na.omit(unlist(c(est.kernel[[char]][, c(4, 5)]))))
                    if(ncol(est.kernel[[char]])>5){
                        yrange <- c(yrange, na.omit(unlist(c(est.kernel[[char]][, c(6, 7)]))))
                    }
                }
                X.lvls <- est.kernel[[other.treat[1]]][, 1]
            }

            if (treat.type == "continuous") {
                for (label in label.name) {
                    yrange <- c(yrange, na.omit(unlist(c(est.kernel[[label]][, c(4, 5)]))))
                    if(ncol(est.kernel[[label]])>5){
                        yrange <- c(yrange, na.omit(unlist(c(est.kernel[[label]][, c(6, 7)]))))
                    }
                }
                X.lvls <- est.kernel[[label.name[1]]][, 1]
            }
        }

        if (is.null(ylim) == FALSE) {
            yrange <- c(ylim[2], ylim[1] + (ylim[2] - ylim[1]) * 1 / 8)
        }
        errorbar.width <- (max(X.lvls) - min(X.lvls)) / 20
        maxdiff <- (max(yrange) - min(yrange))
        pos <- max(yrange) - maxdiff / 20
    }


    # plot initialization
    p.group <- list()
    if (treat.type == "discrete") {
        for (char in other.treat) {
            p1 <- ggplot()
            ## black white theme and mark zero
            if (theme.bw == FALSE) {
                p1 <- p1 + geom_hline(yintercept = 0, colour = "white", size = 2)
            } else {
                p1 <- p1 + theme_bw() + geom_hline(yintercept = 0, colour = "#AAAAAA50", size = 2)
            }
            if (show.grid == FALSE) {
                p1 <- p1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
            }
            p.group[[char]] <- p1
        }
    }

    if (treat.type == "continuous") {
        for (label in label.name) {
            p1 <- ggplot()
            ## black white theme and mark zero
            if (theme.bw == FALSE) {
                p1 <- p1 + geom_hline(yintercept = 0, colour = "white", size = 2)
            } else {
                p1 <- p1 + theme_bw() + geom_hline(yintercept = 0, colour = "#AAAAAA50", size = 2)
            }
            if (show.grid == FALSE) {
                p1 <- p1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
            }
            p.group[[label]] <- p1
        }
    }

    # density
    if (Xdistr == "density") {
        if (treat.type == "discrete") {
            deX.ymin <- min(yrange) - maxdiff / 5
            deX.co <- data.frame(
                x = de.tr[[base]]$x,
                y = de.tr[[base]]$y / max(de.tr[[base]]$y) * maxdiff / 5 + min(yrange) - maxdiff / 5
            )

            ## color
            col.co <- density.color[1]
            col.tr <- density.color[1]
            if (length(density.color) == 2) {
                col.tr <- density.color[2]
            }

            for (char in other.treat) {
                deX.tr <- data.frame(
                    x = de.tr[[char]]$x,
                    y = de.tr[[char]]$y / max(de.tr[[char]]$y) * maxdiff / 5 + min(yrange) - maxdiff / 5
                )

                p1 <- p.group[[char]] + geom_ribbon(
                    data = deX.co, aes(x = x, ymax = y, ymin = deX.ymin),
                    fill = col.co, alpha = density.color.alpha
                ) +
                    geom_ribbon(
                        data = deX.tr, aes(x = x, ymax = y, ymin = deX.ymin),
                        fill = col.tr, alpha = density.color.alpha
                    )
                p.group[[char]] <- p1
            }
        }

        if (treat.type == "continuous") {
            deX.ymin <- min(yrange) - maxdiff / 5
            deX <- data.frame(
                x = de$x,
                y = de$y / max(de$y) * maxdiff / 5 + min(yrange) - maxdiff / 5
            )

            for (label in label.name) {
                ## plotting
                p1 <- p.group[[label]]
                p1 <- p1 + geom_ribbon(
                    data = deX, aes(x = x, ymax = y, ymin = deX.ymin),
                    fill = density.color[1], alpha = density.color.alpha
                )
                p.group[[label]] <- p1
            }
        }
    }

    # hist
    if (Xdistr %in% c("histogram", "hist")) {
        if (treat.type == "discrete") {
            n.hist <- length(hist.out$mids)
            dist <- hist.out$mids[2] - hist.out$mids[1]
            hist.max <- max(hist.out$counts)
            hist.col <- data.frame(
                ymin = rep(min(yrange) - maxdiff / 5, n.hist),
                # ymax=hist.out$counts/hist.max*maxdiff/5+min(yrange)-maxdiff/5,
                xmin = hist.out$mids - dist / 2,
                xmax = hist.out$mids + dist / 2,
                count1 = count.tr[[base]] / hist.max * maxdiff / 5 + min(yrange) - maxdiff / 5
            )

            for (char in other.treat) {
                hist.treat <- data.frame(
                    ymin = hist.col[, "count1"],
                    # ymax=hist.out$counts/hist.max*maxdiff/5+min(yrange)-maxdiff/5,
                    xmin = hist.out$mids - dist / 2,
                    xmax = hist.out$mids + dist / 2,
                    count1 = count.tr[[char]] / hist.max * maxdiff / 5 + hist.col[, "count1"]
                )

                fill1 <- hist.color[1]
                fill2 <- hist.color[1]
                if (length(hist.color) == 2) {
                    fill2 <- hist.color[2]
                }
                p1 <- p.group[[char]] + geom_rect(
                    data = hist.col, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = count1),
                    fill = fill1, colour = "gray50", alpha = hist.color.alpha, size = 0.3
                ) + # control
                    geom_rect(
                        data = hist.treat, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = count1),
                        fill = fill2, colour = "gray50", alpha = hist.color.alpha, size = 0.3
                    )
                p.group[[char]] <- p1
            }
        }

        if (treat.type == "continuous") {
            n.hist <- length(hist.out$mids)
            dist <- hist.out$mids[2] - hist.out$mids[1]
            hist.max <- max(hist.out$counts)
            histX <- data.frame(
                ymin = rep(min(yrange) - maxdiff / 5, n.hist),
                ymax = hist.out$counts / hist.max * maxdiff / 5 + min(yrange) - maxdiff / 5,
                xmin = hist.out$mids - dist / 2,
                xmax = hist.out$mids + dist / 2
            )
            for (label in label.name) {
                p1 <- p.group[[label]]
                p1 <- p1 + geom_rect(
                    data = histX, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                    fill = hist.color[1], colour = "gray50", alpha = hist.color.alpha, size = 0.5
                )
                p.group[[label]] <- p1
            }
        }
    }

    # ME/TE in kernel/linear
    if (estimator == "kernel" | estimator == "linear" | estimator == "DML") {
        if (estimator == "kernel") {
            est <- est.kernel
        } else if (estimator == "DML") {
            est <- est.dml
        } else {
            est <- est.lin
        }

        if (treat.type == "discrete") {
            for (char in other.treat) {
                p1 <- p.group[[char]]
                tempest <- est[[char]]
                if (CI == TRUE) {
                    if(ncol(tempest)==5) {
                        colnames(tempest) <- c("X", "TE", "sd", "CI_lower", "CI_upper")
                    } else{
                        colnames(tempest) <- c("X", "TE", "sd", "CI_lower", "CI_upper","CI_uniform_lower","CI_uniform_upper")
                    }
                }
                if (CI == FALSE) {
                    tempest <- tempest[, c(1, 2)]
                    colnames(tempest) <- c("X", "TE")
                }
                tempest <- as.data.frame(tempest)
                p1 <- p1 + geom_line(data = tempest, aes(X, TE), color = line.color, size = line.size)
                if (CI == TRUE) {
                    if(estimator == "kernel" | estimator == "linear" | estimator == "DML"){
                        p1 <- p1 + geom_ribbon(
                            data = tempest, aes(x = X, ymin = CI_lower, ymax = CI_upper),
                            fill = CI.color, alpha = CI.color.alpha
                        )                        
                    }
                    #if(estimator == "DML"){
                    #    p1 <- p1 + geom_line(data = tempest,aes(x=X, y = CI_lower),linetype = 'dashed',color = 'gray50') + geom_line(data = tempest,aes(x=X, y = CI_upper),linetype = 'dashed',color = 'gray50')
                    #}

                    if("CI_uniform_lower" %in% colnames(tempest)){
                        p1 <- p1 + geom_line(data = tempest,aes(x=X, y = CI_uniform_lower),linetype = 'dashed',color = 'gray50') + geom_line(data = tempest,aes(x=X, y = CI_uniform_upper),linetype = 'dashed',color = 'gray50')
                    }
                }
                # ymin=min(yrange)-maxdiff/5

                if (is.null(diff.values) == FALSE) {
                    for (target.value in diff.values) {
                        Xnew <- abs(tempest[, "X"] - target.value)
                        d1 <- min(Xnew)
                        label1 <- which.min(Xnew)
                        Xnew[label1] <- Inf
                        d2 <- min(Xnew)
                        label2 <- which.min(Xnew)
                        if (d1 == 0) {
                            est.mark <- tempest[label1, "TE"]
                            if (CI == TRUE) {
                                lb.mark <- tempest[label1, "CI_lower"]
                                ub.mark <- tempest[label1, "CI_upper"]
                            }
                        } else if (d2 == 0) {
                            est.mark <- tempest[label2, "TE"]
                            if (CI == TRUE) {
                                lb.mark <- tempest[label2, "CI_lower"]
                                ub.mark <- tempest[label2, "CI_upper"]
                            }
                        } else { ## weighted average
                            est.mark1 <- tempest[label1, "TE"]
                            est.mark2 <- tempest[label2, "TE"]
                            est.mark <- ((est.mark1 * d2 + est.mark2 * d1) / (d1 + d2))
                            if (CI == TRUE) {
                                lb.mark1 <- tempest[label1, "CI_lower"]
                                ub.mark1 <- tempest[label1, "CI_upper"]
                                lb.mark2 <- tempest[label2, "CI_lower"]
                                ub.mark2 <- tempest[label2, "CI_upper"]
                                lb.mark <- ((lb.mark1 * d2 + lb.mark2 * d1) / (d1 + d2))
                                ub.mark <- ((ub.mark1 * d2 + ub.mark2 * d1) / (d1 + d2))
                            }
                        }

                        p1 <- p1 + annotate("point", x = target.value, y = est.mark, size = 1, colour = "red")
                        if (CI == TRUE) {
                            p1 <- p1 + annotate("errorbar", x = target.value, ymin = lb.mark, ymax = ub.mark, colour = "red", size = 0.5, width = (max(tempxx) - min(tempxx)) / 30)
                        }
                    }
                }
                p.group[[char]] <- p1
            }
        }

        if (treat.type == "continuous") {
            for (label in label.name) {
                p1 <- p.group[[label]]
                tempest <- est[[label]]
                if (CI == TRUE) {
                    if(ncol(tempest)==5){
                        colnames(tempest) <- c("X", "ME", "sd", "CI_lower", "CI_upper")
                    }else{
                        colnames(tempest) <- c("X", "ME", "sd", "CI_lower", "CI_upper","CI_uniform_lower","CI_uniform_upper")
                    }
                    
                    tempest <- as.data.frame(tempest)
                }
                if (CI == FALSE) {
                    tempest <- tempest[, c(1, 2)]
                    colnames(tempest) <- c("X", "ME")
                    tempest <- as.data.frame(tempest)
                }
                p1 <- p1 + geom_line(data = tempest, aes(X, ME), color = line.color, size = line.size)
                if (CI == TRUE) {
                    if(estimator == "kernel" | estimator == "linear" | estimator == "DML"){
                        p1 <- p1 + geom_ribbon(
                            data = tempest, aes(x = X, ymin = CI_lower, ymax = CI_upper),
                            fill = CI.color, alpha = CI.color.alpha
                        )                        
                    }
                    #if(estimator == "DML"){
                    #    p1 <- p1 + geom_line(data = tempest,aes(x=X, y = CI_lower),linetype = 'dashed',color = 'gray50') + geom_line(data = tempest,aes(x=X, y = CI_upper),linetype = 'dashed',color = 'gray50')
                    #}
                    if("CI_uniform_lower" %in% colnames(tempest)){
                        p1 <- p1 + geom_line(data = tempest,aes(x=X, y = CI_uniform_lower),linetype = 'dashed',color = 'gray50') + geom_line(data = tempest,aes(x=X, y = CI_uniform_upper),linetype = 'dashed',color = 'gray50')
                    }
                }
                # ymin=min(yrange)-maxdiff/5

                if (is.null(diff.values) == FALSE) {
                    for (target.value in diff.values) {
                        Xnew <- abs(tempest[, "X"] - target.value)
                        d1 <- min(Xnew)
                        label1 <- which.min(Xnew)
                        Xnew[label1] <- Inf
                        d2 <- min(Xnew)
                        label2 <- which.min(Xnew)
                        if (d1 == 0) {
                            est.mark <- tempest[label1, "ME"]
                            if (CI == TRUE) {
                                lb.mark <- tempest[label1, "CI_lower"]
                                ub.mark <- tempest[label1, "CI_upper"]
                            }
                        } else if (d2 == 0) {
                            est.mark <- tempest[label2, "ME"]
                            if (CI == TRUE) {
                                lb.mark <- tempest[label2, "CI_lower"]
                                ub.mark <- tempest[label2, "CI_upper"]
                            }
                        } else { ## weighted average
                            est.mark1 <- tempest[label1, "ME"]
                            est.mark2 <- tempest[label2, "ME"]
                            est.mark <- ((est.mark1 * d2 + est.mark2 * d1) / (d1 + d2))
                            if (CI == TRUE) {
                                lb.mark1 <- tempest[label1, "CI_lower"]
                                ub.mark1 <- tempest[label1, "CI_upper"]
                                lb.mark2 <- tempest[label2, "CI_lower"]
                                ub.mark2 <- tempest[label2, "CI_upper"]
                                lb.mark <- ((lb.mark1 * d2 + lb.mark2 * d1) / (d1 + d2))
                                ub.mark <- ((ub.mark1 * d2 + ub.mark2 * d1) / (d1 + d2))
                            }
                        }

                        p1 <- p1 + annotate("point", x = target.value, y = est.mark, size = 1, colour = "red")
                        if (CI == TRUE) {
                            p1 <- p1 + annotate("errorbar", x = target.value, ymin = lb.mark, ymax = ub.mark, colour = "red", size = 0.5, width = (max(tempxx) - min(tempxx)) / 30)
                        }
                    }
                }

                p.group[[label]] <- p1
            }
        }
    }

    if (estimator == "binning") {
        if (treat.type == "discrete") {
            for (char in other.treat) {
                p1 <- p.group[[char]]
                tempest <- est.lin[[char]]
                tempest.bin2 <- est.bin2[[char]]
                if (CI == TRUE) {
                    colnames(tempest) <- c("X", "TE", "sd", "CI_lower", "CI_upper")
                    colnames(tempest.bin2) <- c("x0", "TE", "sd", "CI_lower", "CI_upper")
                    tempest <- as.data.frame(tempest)
                    tempest.bin2 <- as.data.frame(tempest.bin2)
                }
                if (CI == FALSE) {
                    tempest <- tempest[, c(1, 2)]
                    colnames(tempest) <- c("X", "TE")
                    tempest.bin2 <- tempest.bin2[, c(1, 2)]
                    colnames(tempest.bin2) <- c("x0", "TE")
                    tempest <- as.data.frame(tempest)
                    tempest.bin2 <- as.data.frame(tempest.bin2)
                }

                p1 <- p1 + geom_line(data = tempest, aes(X, TE), color = line.color, size = line.size)
                if (CI == TRUE) {
                    p1 <- p1 + geom_ribbon(
                        data = tempest, aes(x = X, ymin = CI_lower, ymax = CI_upper),
                        fill = CI.color, alpha = CI.color.alpha
                    )
                }
                ## bin estimates
                p1 <- p1 + geom_point(data = tempest.bin2, aes(x0, TE), size = 4 / treat_sc, shape = 21, fill = "white", colour = "red")
                if (CI == TRUE) {
                    p1 <- p1 + geom_errorbar(
                        data = tempest.bin2, aes(x = x0, ymin = CI_lower, ymax = CI_upper), colour = "red", size = 1,
                        width = errorbar.width
                    )
                }

                ## in case there's non-overlap


                p1 <- p1 + annotate(
                    geom = "text", x = est.bin3[[char]][, 1], y = rep(0, dim(est.bin3[[char]])[1]),
                    label = "NaN", colour = "red"
                )

                ## labels: L, M, H and so on
                if (bin.labs == TRUE) {
                    if (nbins == 3) {
                        p1 <- p1 + annotate(
                            geom = "text", x = est.bin[[char]][1, 1], y = pos,
                            label = "L", colour = "gray50", size = 10 / treat_sc
                        ) +
                            annotate(
                                geom = "text", x = est.bin[[char]][2, 1], y = pos,
                                label = "M", colour = "gray50", size = 10 / treat_sc
                            ) +
                            annotate(
                                geom = "text", x = est.bin[[char]][3, 1], y = pos,
                                label = "H", colour = "gray50", size = 10 / treat_sc
                            )
                    } else if (nbins == 4) {
                        p1 <- p1 + annotate(
                            geom = "text", x = est.bin[[char]][1, 1], y = pos,
                            label = "L", colour = "gray50", size = 10 / treat_sc
                        ) +
                            annotate(
                                geom = "text", x = est.bin[[char]][2, 1], y = pos,
                                label = "M1", colour = "gray50", size = 10 / treat_sc
                            ) +
                            annotate(
                                geom = "text", x = est.bin[[char]][3, 1], y = pos,
                                label = "M2", colour = "gray50", size = 10 / treat_sc
                            ) +
                            annotate(
                                geom = "text", x = est.bin[[char]][4, 1], y = pos,
                                label = "H", colour = "gray50", size = 10 / treat_sc
                            )
                    } else if (nbins == 2) {
                        p1 <- p1 + annotate(
                            geom = "text", x = est.bin[[char]][1, 1], y = pos,
                            label = "L", colour = "gray50", size = 10 / treat_sc
                        ) +
                            annotate(
                                geom = "text", x = est.bin[[char]][2, 1], y = pos,
                                label = "H", colour = "gray50", size = 10 / treat_sc
                            )
                    }
                }
                p.group[[char]] <- p1
            }
        }

        if (treat.type == "continuous") {
            for (label in label.name) {
                p1 <- p.group[[label]]
                tempest <- est.lin[[label]]
                tempest.bin2 <- est.bin2[[label]]
                if (CI == TRUE) {
                    colnames(tempest) <- c("X", "ME", "sd", "CI_lower", "CI_upper")
                    colnames(tempest.bin2) <- c("x0", "ME", "sd", "CI_lower", "CI_upper")
                    tempest <- as.data.frame(tempest)
                    tempest.bin2 <- as.data.frame(tempest.bin2)
                }
                if (CI == FALSE) {
                    tempest <- tempest[, c(1, 2)]
                    colnames(tempest) <- c("X", "ME")
                    tempest.bin2 <- tempest.bin2[, c(1, 2)]
                    colnames(tempest.bin2) <- c("x0", "ME")
                    tempest <- as.data.frame(tempest)
                    tempest.bin2 <- as.data.frame(tempest.bin2)
                }

                p1 <- p1 + geom_line(data = tempest, aes(X, ME), color = line.color, size = line.size)
                if (CI == TRUE) {
                    p1 <- p1 + geom_ribbon(
                        data = tempest, aes(x = X, ymin = CI_lower, ymax = CI_upper),
                        fill = CI.color, alpha = CI.color.alpha
                    )
                }
                ## bin estimates
                p1 <- p1 + geom_point(data = tempest.bin2, aes(x0, ME), size = 4 / treat_sc, shape = 21, fill = "white", colour = "red")
                if (CI == TRUE) {
                    p1 <- p1 + geom_errorbar(
                        data = tempest.bin2, aes(x = x0, ymin = CI_lower, ymax = CI_upper), colour = "red", size = 1,
                        width = errorbar.width
                    )
                }
                ## in case there's non-overlap
                p1 <- p1 + annotate(
                    geom = "text", x = est.bin3[[label]][, 1], y = rep(0, dim(est.bin3[[label]])[1]),
                    label = "NaN", colour = "red"
                )

                ## labels: L, M, H and so on
                if (bin.labs == TRUE) {
                    if (nbins == 3) {
                        p1 <- p1 + annotate(
                            geom = "text", x = est.bin[[label]][1, 1], y = pos,
                            label = "L", colour = "gray50", size = 10 / treat_sc
                        ) +
                            annotate(
                                geom = "text", x = est.bin[[label]][2, 1], y = pos,
                                label = "M", colour = "gray50", size = 10 / treat_sc
                            ) +
                            annotate(
                                geom = "text", x = est.bin[[label]][3, 1], y = pos,
                                label = "H", colour = "gray50", size = 10 / treat_sc
                            )
                    } else if (nbins == 4) {
                        p1 <- p1 + annotate(
                            geom = "text", x = est.bin[[label]][1, 1], y = pos,
                            label = "L", colour = "gray50", size = 10 / treat_sc
                        ) +
                            annotate(
                                geom = "text", x = est.bin[[label]][2, 1], y = pos,
                                label = "M1", colour = "gray50", size = 10 / treat_sc
                            ) +
                            annotate(
                                geom = "text", x = est.bin[[label]][3, 1], y = pos,
                                label = "M2", colour = "gray50", size = 10 / treat_sc
                            ) +
                            annotate(
                                geom = "text", x = est.bin[[label]][4, 1], y = pos,
                                label = "H", colour = "gray50", size = 10 / treat_sc
                            )
                    } else if (nbins == 2) {
                        p1 <- p1 + annotate(
                            geom = "text", x = est.bin[[label]][1, 1], y = pos,
                            label = "L", colour = "gray50", size = 10 / treat_sc
                        ) +
                            annotate(
                                geom = "text", x = est.bin[[label]][2, 1], y = pos,
                                label = "H", colour = "gray50", size = 10 / treat_sc
                            )
                    }
                }
                p.group[[label]] <- p1
            }
        }
    }

    # cex/title...
    if (is.null(cex.lab) == TRUE) {
        cex.lab <- 15
    } else {
        cex.lab <- 15 * cex.lab
    }
    if (is.null(cex.axis) == TRUE) {
        cex.axis <- 15 / treat_sc
    } else {
        cex.axis <- 15 * cex.axis / treat_sc
    }
    ## title
    if (is.null(cex.main) == TRUE) {
        cex.main <- 18
    } else {
        cex.main <- 18 * cex.main
    }

    if (is.null(cex.sub) == TRUE) {
        cex.sub <- 12
    } else {
        cex.sub <- 12 * cex.sub
    }

    ## xlim and ylim
    if (is.null(ylim) == FALSE) {
        ylim2 <- c(ylim[1] - (ylim[2] - ylim[1]) * 0.25 / 6, ylim[2] + (ylim[2] - ylim[1]) * 0.4 / 6)
    }

    if (treat.type == "discrete") {
        k <- 1
        for (char in other.treat) {
            p1 <- p.group[[char]]
            ## mark the original interval (in replicated papers)
            if (is.null(interval) == FALSE) {
                p1 <- p1 + geom_vline(xintercept = interval, colour = "steelblue", linetype = 2, size = 1.5)
            }

            ## Other universal options
            p1 <- p1 + xlab(NULL) + ylab(NULL) +
                theme(axis.text = element_text(size = cex.axis))

            if (show.subtitles == TRUE) {
                if (is.null(subtitles) == TRUE) {
                    subtitle.temp <- paste0("Treated = ", char, ", Baseline = ", base)
                    p1 <- p1 + labs(subtitle = subtitle.temp) + theme(plot.subtitle = element_text(hjust = 0.5, size = cex.sub, lineheight = .8))
                }

                if (is.null(subtitles) == FALSE) {
                    subtitle.temp <- subtitles[k]
                    p1 <- p1 + labs(subtitle = subtitle.temp) + theme(plot.subtitle = element_text(hjust = 0.5, size = cex.sub, lineheight = .8))
                }
            }

            if (is.null(xlim) == FALSE & is.null(ylim) == FALSE) {
                p1 <- p1 + coord_cartesian(xlim = xlim, ylim = ylim2)
            }
            if (is.null(xlim) == TRUE & is.null(ylim) == FALSE) {
                p1 <- p1 + coord_cartesian(ylim = ylim2)
            }
            if (is.null(xlim) == FALSE & is.null(ylim) == TRUE) {
                p1 <- p1 + coord_cartesian(xlim = xlim)
            }

            p.group[[char]] <- p1
            k <- k + 1
        }

        ## ylim
        for (char in other.treat) {
            y.limits <- layer_scales(p.group[[char]])$y$range$range
            if (char == other.treat[1]) {
                ymaxmax <- y.limits[2]
                yminmin <- y.limits[1]
            } else {
                ymaxmax <- max(ymaxmax, y.limits[2])
                yminmin <- min(yminmin, y.limits[1])
            }
        }
        for (char in other.treat) {
            p.group[[char]] <- p.group[[char]] + ylim(c(yminmin, ymaxmax))
        }

        requireNamespace("gridExtra")
        requireNamespace("grid")
        requireNamespace("ggplotify")
        suppressMessages(
            graph <- arrangeGrob(
                grobs = p.group, ncol = ncols, align = "h",
                top = textGrob(main, gp = gpar(fontsize = cex.main, face = "bold")),
                left = textGrob(ylab, rot = 90, vjust = 1, gp = gpar(fontsize = cex.lab)),
                bottom = textGrob(xlab, gp = gpar(fontsize = cex.lab))
            )
        )
    }

    if (treat.type == "continuous") {
        k <- 1
        for (label in label.name) {
            p1 <- p.group[[label]]
            ## mark the original interval (in replicated papers)
            if (is.null(interval) == FALSE) {
                p1 <- p1 + geom_vline(xintercept = interval, colour = "steelblue", linetype = 2, size = 1.5)
            }

            ## Other universal options
            p1 <- p1 + xlab(NULL) + ylab(NULL) +
                theme(axis.text = element_text(size = cex.axis))

            if (show.subtitles == TRUE) {
                if (is.null(subtitles) == TRUE) {
                    subtitle.temp <- label
                    p1 <- p1 + labs(subtitle = subtitle.temp) + theme(plot.subtitle = element_text(hjust = 0.5, size = cex.sub, lineheight = .8))
                }

                if (is.null(subtitles) == FALSE) {
                    subtitle.temp <- subtitles[k]
                    p1 <- p1 + labs(subtitle = subtitle.temp) + theme(plot.subtitle = element_text(hjust = 0.5, size = cex.sub, lineheight = .8))
                }
            }

            if (is.null(xlim) == FALSE & is.null(ylim) == FALSE) {
                p1 <- p1 + coord_cartesian(xlim = xlim, ylim = ylim2)
            }
            if (is.null(xlim) == TRUE & is.null(ylim) == FALSE) {
                p1 <- p1 + coord_cartesian(ylim = ylim2)
            }
            if (is.null(xlim) == FALSE & is.null(ylim) == TRUE) {
                p1 <- p1 + coord_cartesian(xlim = xlim)
            }

            p.group[[label]] <- p1
            k <- k + 1
        }

        ## ylim
        for (label in label.name) {
            y.limits <- layer_scales(p.group[[label]])$y$range$range
            if (label == label.name[1]) {
                ymaxmax <- y.limits[2]
                yminmin <- y.limits[1]
            } else {
                ymaxmax <- max(ymaxmax, y.limits[2])
                yminmin <- min(yminmin, y.limits[1])
            }
        }
        for (label in label.name) {
            p.group[[label]] <- p.group[[label]] + ylim(c(yminmin, ymaxmax))
        }

        requireNamespace("gridExtra")
        requireNamespace("grid")
        requireNamespace("ggplotify")
        suppressMessages(
            graph <- arrangeGrob(
                grobs = p.group, ncol = ncols, align = "h",
                top = textGrob(main, gp = gpar(fontsize = cex.main, face = "bold")),
                left = textGrob(ylab, rot = 90, vjust = 1, gp = gpar(fontsize = cex.lab)),
                bottom = textGrob(xlab, gp = gpar(fontsize = cex.lab))
            )
        )
    }



    ## save to file
    if (is.null(file) == FALSE) {
        graph <- as.ggplot(graph)
        ggsave(file, graph, scale = scale, width = width, height = height)
    }

    if (show.all == FALSE) {
        graph <- as.ggplot(graph)
        return(graph)
    } else {
        cex.axis <- cex.axis
        if (treat.type == "discrete") {
            for (char in other.treat) {
                ## Other universal options
                p1 <- p.group[[char]] + xlab(xlab) + ylab(ylab) +
                    theme(
                        axis.text = element_text(size = cex.axis),
                        axis.title.x = element_text(size = cex.lab),
                        axis.title.y = element_text(size = cex.lab)
                    )
                p1 <- p1 + ggtitle(main) +
                    theme(plot.title = element_text(hjust = 0.5, size = cex.main, lineheight = .8, face = "bold"))
                p.group[[char]] <- p1
            }
        }

        if (treat.type == "continuous") {
            for (label in label.name) {
                ## Other universal options
                p1 <- p.group[[label]] + xlab(xlab) + ylab(ylab) +
                    theme(
                        axis.text = element_text(size = cex.axis),
                        axis.title.x = element_text(size = cex.lab),
                        axis.title.y = element_text(size = cex.lab)
                    )
                p1 <- p1 + ggtitle(main) +
                    theme(plot.title = element_text(hjust = 0.5, size = cex.main, lineheight = .8, face = "bold"))
                p.group[[label]] <- p1
            }
        }
        return(p.group)
    }
}
