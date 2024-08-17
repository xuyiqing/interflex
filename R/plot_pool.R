interflex.plot.pool <- function(
    # only for discrete treatments
    out,
    diff.values = NULL,
    order = NULL,
    subtitles = NULL,
    show.subtitles = NULL,
    legend.title = NULL,
    CI = TRUE,
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
    cex.lab = NULL,
    cex.axis = NULL,
    cex.sub = NULL,
    bin.labs = TRUE, # bin labels
    interval = NULL, # interval in replicated papers
    line.size = NULL,
    color = NULL,
    line.color = NULL,
    CI.color = NULL,
    CI.color.alpha = NULL,
    hist.color = NULL,
    hist.color.alpha = NULL,
    density.color = NULL,
    density.color.alpha = NULL,
    file = NULL,
    scale = 1.1,
    height = 7,
    width = 10) {
    X <- NULL
    TE <- NULL
    Treatment <- NULL
    CI_lower <- NULL
    CI_upper <- NULL
    ME <- NULL
    x0 <- NULL
    CI.lower <- NULL
    CI.upper <- NULL
    x <- NULL
    y <- NULL
    end_level <- NULL
    xmin <- NULL
    xmax <- NULL
    count1 <- NULL
    ymax <- NULL
    D <- NULL
    name <- NULL
    r <- NULL


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
    } else {
        show.subtitles <- TRUE
    }

    # CI
    if (is.null(CI) == FALSE) {
        if (is.logical(CI) == FALSE & is.numeric(CI) == FALSE) {
            stop("\"CI\" is not a logical flag.")
        }
        if (estimator == "kernel") {
            if (CI == TRUE & out$CI == FALSE) {
                stop("Confidence intervals are not estimated, please set CI to FALSE.")
            }
        }
    }

    if (estimator == "kernel") {
        if (is.null(CI) == TRUE) {
            CI <- out$CI
        }
    }

    if (estimator == "binning" | estimator == "linear" | estimator == "DML" | estimator == "grf") {
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

    ## color
    if (is.null(color) == FALSE) {
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

    if (is.null(legend.title) == FALSE) {
        legend.title <- as.character(legend.title)[1]
    }

    if (treat.type == "discrete" & (estimator == "linear" | estimator == "binning")) {
        tempxx <- out$est.lin[[other.treat[1]]][, "X"]
    }
    if (treat.type == "discrete" & estimator == "DML") {
        tempxx <- out$est.dml[[other.treat[1]]][, "X"]
    }
    if (treat.type == "discrete" & estimator == "grf") {
        tempxx <- out$est.grf[[other.treat[1]]][, "X"]
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
            show.subtitles <- TRUE
        }

        if (is.null(subtitles) == FALSE) {
            if (length(subtitles) != length(all.treat)) {
                stop("The number of elements in \"subtitles\" should be m(m is the number of different treatment arms including the baseline group).")
            }
        }

        if (is.null(subtitles) == TRUE) {
            base.name <- paste0("Base Group (", base, ")")
            subtitles <- c(base.name, other.treat)
        } else {
            base.name <- subtitles[1]
        }

        subtitles.all <- as.character(subtitles)
        subtitles <- subtitles.all[2:length(subtitles.all)]
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
            show.subtitles <- TRUE
        }

        if (is.null(subtitles) == FALSE) {
            if (length(subtitles) != length(label.name)) {
                stop("The number of elements in \"subtitles\" should equal to the number of values in D.ref.")
            }
        }

        if (is.null(subtitles) == TRUE) {
            subtitles <- label.name
        }

        subtitles.all <- subtitles
    }

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
        if (treat.type == "discrete" & estimator == "linear") {
            tempxx <- out$est.lin[[other.treat[1]]][, "X"]
        }
        if (treat.type == "discrete" & (estimator == "DML")) {
            tempxx <- out$est.dml[[other.treat[1]]][, "X"]
        }
        if (treat.type == "discrete" & (estimator == "grf")) {
            tempxx <- out$est.grf[[other.treat[1]]][, "X"]
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
        ymin <- min(yrange) - maxdiff / 5
    }

    if (estimator == "DML") {
        if (treat.type == "discrete") {
            est.dml <- out$est.dml
            yrange <- c(0)
            for (char in other.treat) {
                if (CI == TRUE) {
                    yrange <- c(yrange, na.omit(unlist(c(est.dml[[char]][, c(4, 5)]))))
                    if (ncol(est.dml[[char]]) > 5) {
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
                    if (ncol(est.dml[[label]]) > 5) {
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

    if (estimator == "grf") {
        if (treat.type == "discrete") {
            est.grf <- out$est.grf
            yrange <- c(0)
            for (char in other.treat) {
                if (CI == TRUE) {
                    yrange <- c(yrange, na.omit(unlist(c(est.grf[[char]][, c(4, 5)]))))
                    if (ncol(est.grf[[char]]) > 5) {
                        yrange <- c(yrange, na.omit(unlist(c(est.grf[[char]][, c(6, 7)]))))
                    }
                } else {
                    yrange <- c(yrange, na.omit(unlist(c(est.grf[[char]][, 2]))))
                }
            }
            X.lvls <- est.grf[[other.treat[1]]][, 1]
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
        ymin <- min(yrange) - maxdiff / 5
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
                }
                X.lvls <- est.kernel[[other.treat[1]]][, 1]
            }

            if (treat.type == "continuous") {
                for (label in label.name) {
                    yrange <- c(yrange, na.omit(unlist(c(est.kernel[[label]][, c(4, 5)]))))
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
        ymin <- min(yrange) - maxdiff / 5
    }

    # color
    # get base.color & platte(for discrete data)
    requireNamespace("RColorBrewer")
    platte <- brewer.pal(n = 8, "Set2")

    if (is.null(color) == TRUE) {
        base.color <- "gray50"
    }

    if (is.null(color) == FALSE) {
        if (treat.type == "discrete") {
            base.color <- color[1]
            if (length(color) == 1) {
                platte <- platte
            } else {
                platte <- c(color[2:length(color)], platte)
            }
        }

        if (treat.type == "continuous") {
            platte <- c(color, platte)
        }
    }

    if (treat.type == "discrete") {
        num.treat <- length(other.treat)
    }

    if (treat.type == "continuous") {
        num.treat <- length(label.name)
    }

    platte <- platte[1:num.treat]
    # initialize
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

    if (estimator == "kernel" | estimator == "linear" | estimator == "DML" | estimator == "grf") {
        if (estimator == "kernel") {
            est <- est.kernel
        } else if (estimator == "DML") {
            est <- est.dml
        } else if (estimator == "grf") {
            est <- est.grf
        } else {
            est <- est.lin
        }

        if (treat.type == "discrete") {
            for (char in other.treat) {
                est.touse <- est[[char]]
                if (CI == TRUE) {
                    if (estimator == "DML") {
                        colnames(est.touse) <- c("X", "TE", "sd", "CI_lower", "CI_upper")
                    } else {
                        if (dim(est.touse)[2] == 5) {
                            colnames(est.touse) <- c("X", "TE", "sd", "CI_lower", "CI_upper")
                        } else {
                            colnames(est.touse) <- c("X", "TE", "sd", "CI_lower", "CI_upper", "CI_uniform_lower", "CI_uniform_upper")
                        }
                    }
                }
                if (CI == FALSE) {
                    est.touse <- est.touse[, c(1, 2)]
                    colnames(est.touse) <- c("X", "TE")
                }
                est.touse <- as.data.frame(est.touse)
                est.touse[["Treatment"]] <- rep(char, dim(est.touse)[1])

                if (char == other.treat[1]) {
                    tempest <- est.touse
                } else {
                    tempest <- rbind(tempest, est.touse)
                }
            }
            tempest <- as.data.frame(tempest)
            tempest$Treatment <- factor(tempest$Treatment, levels = other.treat)
            p1 <- p1 + geom_line(data = tempest, aes(x = X, y = TE, color = Treatment), show.legend = FALSE)
            p1 <- p1 + scale_color_manual(values = platte, labels = subtitles)
            if (CI == TRUE) {
                # if(estimator!='DML'){
                p1 <- p1 + geom_ribbon(data = tempest, aes(x = X, ymin = CI_lower, ymax = CI_upper, fill = Treatment), alpha = 0.2, show.legend = F)
                p1 <- p1 + scale_fill_manual(values = platte, labels = subtitles)
                if ("CI_uniform_lower" %in% colnames(tempest)) {
                    p1 <- p1 + geom_ribbon(data = tempest, aes(x = X, ymin = CI_uniform_lower, ymax = CI_uniform_upper, color = Treatment), linetype = 2, alpha = 0, show.legend = F)
                    p1 <- p1 + scale_color_manual(values = platte, labels = subtitles)
                }
                # }
                # else{
                #    p1 <- p1 + geom_ribbon(data = tempest, aes(x = X, ymin = CI_uniform_lower, ymax = CI_uniform_upper, color = Treatment), linetype = 2,alpha = 0, show.legend = F)
                #    p1 <- p1 + scale_color_manual(values = platte, labels = subtitles)
                # }
            }

            if (is.null(diff.values) == FALSE) {
                k <- 1
                for (char in other.treat) {
                    tempest <- est[[char]]
                    for (target.value in diff.values) {
                        Xnew <- abs(tempest[, "X"] - target.value)
                        d1 <- min(Xnew)
                        label1 <- which.min(Xnew)
                        Xnew[label1] <- Inf
                        d2 <- min(Xnew)
                        label2 <- which.min(Xnew)
                        if (d1 == 0) {
                            est.mark <- tempest[label1, 2]
                            if (CI == TRUE) {
                                lb.mark <- tempest[label1, 4]
                                ub.mark <- tempest[label1, 5]
                            }
                        } else if (d2 == 0) {
                            est.mark <- tempest[label2, 2]
                            if (CI == TRUE) {
                                lb.mark <- tempest[label2, 4]
                                ub.mark <- tempest[label2, 5]
                            }
                        } else { ## weighted average
                            est.mark1 <- tempest[label1, 2]
                            est.mark2 <- tempest[label2, 2]
                            est.mark <- ((est.mark1 * d2 + est.mark2 * d1) / (d1 + d2))
                            if (CI == TRUE) {
                                lb.mark1 <- tempest[label1, 4]
                                ub.mark1 <- tempest[label1, 5]
                                lb.mark2 <- tempest[label2, 4]
                                ub.mark2 <- tempest[label2, 5]
                                lb.mark <- ((lb.mark1 * d2 + lb.mark2 * d1) / (d1 + d2))
                                ub.mark <- ((ub.mark1 * d2 + ub.mark2 * d1) / (d1 + d2))
                            }
                        }

                        p1 <- p1 + annotate("point", x = target.value, y = est.mark, size = 1, colour = platte[k])
                        if (CI == TRUE) {
                            p1 <- p1 + annotate("errorbar", x = target.value, ymin = lb.mark, ymax = ub.mark, colour = platte[k], size = 0.5, width = (max(tempxx) - min(tempxx)) / 20)
                        }
                    }
                    k <- k + 1
                }
            }
        }

        if (treat.type == "continuous") {
            for (label in label.name) {
                est.touse <- est[[label]]
                if (CI == TRUE) {
                    if (dim(est.touse)[2] == 5) {
                        if (estimator == "DML") {
                            colnames(est.touse) <- c("X", "TE", "sd", "CI_lower", "CI_upper")
                        } else {
                            colnames(est.touse) <- c("X", "TE", "sd", "CI_lower", "CI_upper")
                        }
                    } else {
                        colnames(est.touse) <- c("X", "TE", "sd", "CI_lower", "CI_upper", "CI_uniform_lower", "CI_uniform_upper")
                    }
                }
                if (CI == FALSE) {
                    est.touse <- est.touse[, c(1, 2)]
                    colnames(est.touse) <- c("X", "ME")
                }

                est.touse <- as.data.frame(est.touse)
                est.touse[["Treatment"]] <- rep(label, dim(est.touse)[1])

                if (label == label.name[1]) {
                    tempest <- est.touse
                } else {
                    tempest <- rbind(tempest, est.touse)
                }
            }
            tempest <- as.data.frame(tempest)
            tempest$Treatment <- factor(tempest$Treatment, levels = label.name)
            p1 <- p1 + geom_line(data = tempest, aes(x = X, y = ME, color = Treatment), show.legend = FALSE)
            p1 <- p1 + scale_color_manual(values = platte, labels = subtitles)
            if (CI == TRUE) {
                # if(estimator!='DML'){
                p1 <- p1 + geom_ribbon(data = tempest, aes(x = X, ymin = CI_lower, ymax = CI_upper, fill = Treatment), alpha = 0.2, show.legend = F)
                p1 <- p1 + scale_fill_manual(values = platte, labels = subtitles)
                if ("CI_uniform_lower" %in% colnames(tempest)) {
                    p1 <- p1 + geom_ribbon(data = tempest, aes(x = X, ymin = CI_uniform_lower, ymax = CI_uniform_upper, color = Treatment), linetype = 2, alpha = 0, show.legend = F)
                    p1 <- p1 + scale_color_manual(values = platte, labels = subtitles)
                }
                # }
                # else{
                #    p1 <- p1 + geom_ribbon(data = tempest, aes(x = X, ymin = CI_uniform_lower, ymax = CI_uniform_upper, color = Treatment), linetype = 2,alpha = 0, show.legend = F)
                #    p1 <- p1 + scale_color_manual(values = platte, labels = subtitles)
                # }
            }

            if (is.null(diff.values) == FALSE) {
                k <- 1
                for (label in label.name) {
                    tempest <- est[[label]]
                    for (target.value in diff.values) {
                        Xnew <- abs(tempest[, "X"] - target.value)
                        d1 <- min(Xnew)
                        label1 <- which.min(Xnew)
                        Xnew[label1] <- Inf
                        d2 <- min(Xnew)
                        label2 <- which.min(Xnew)
                        if (d1 == 0) {
                            est.mark <- tempest[label1, 2]
                            if (CI == TRUE) {
                                lb.mark <- tempest[label1, 4]
                                ub.mark <- tempest[label1, 5]
                            }
                        } else if (d2 == 0) {
                            est.mark <- tempest[label2, 2]
                            if (CI == TRUE) {
                                lb.mark <- tempest[label2, 4]
                                ub.mark <- tempest[label2, 5]
                            }
                        } else { ## weighted average
                            est.mark1 <- tempest[label1, 2]
                            est.mark2 <- tempest[label2, 2]
                            est.mark <- ((est.mark1 * d2 + est.mark2 * d1) / (d1 + d2))
                            if (CI == TRUE) {
                                lb.mark1 <- tempest[label1, 4]
                                ub.mark1 <- tempest[label1, 5]
                                lb.mark2 <- tempest[label2, 4]
                                ub.mark2 <- tempest[label2, 5]
                                lb.mark <- ((lb.mark1 * d2 + lb.mark2 * d1) / (d1 + d2))
                                ub.mark <- ((ub.mark1 * d2 + ub.mark2 * d1) / (d1 + d2))
                            }
                        }

                        p1 <- p1 + annotate("point", x = target.value, y = est.mark, size = 1, colour = platte[k])
                        if (CI == TRUE) {
                            p1 <- p1 + annotate("errorbar", x = target.value, ymin = lb.mark, ymax = ub.mark, colour = platte[k], size = 0.5, width = (max(tempxx) - min(tempxx)) / 20)
                        }
                    }
                    k <- k + 1
                }
            }
        }
    }

    if (estimator == "binning") {
        # est <- est.lin
        if (treat.type == "discrete") {
            for (char in other.treat) {
                est.touse <- est.lin[[char]]
                if (CI == TRUE) {
                    colnames(est.touse) <- c("X", "TE", "sd", "CI_lower", "CI_upper")
                }
                if (CI == FALSE) {
                    est.touse <- est.touse[, c(1, 2)]
                    colnames(est.touse) <- c("X", "TE")
                }
                est.touse <- as.data.frame(est.touse)
                est.touse[["Treatment"]] <- rep(char, dim(est.touse)[1])

                if (char == other.treat[1]) {
                    tempest <- est.touse
                } else {
                    tempest <- rbind(tempest, est.touse)
                }
            }
            tempest$Treatment <- factor(tempest$Treatment, levels = other.treat)
            p1 <- p1 + geom_line(data = tempest, aes(x = X, y = TE, color = Treatment), show.legend = FALSE)
            p1 <- p1 + scale_color_manual(values = platte, labels = subtitles)
            if (CI == TRUE) {
                p1 <- p1 + geom_ribbon(data = tempest, aes(x = X, ymin = CI_lower, ymax = CI_upper, fill = Treatment), alpha = 0.2, show.legend = FALSE)
                p1 <- p1 + scale_fill_manual(values = platte, labels = subtitles)
            }

            k <- 1
            for (char in other.treat) {
                tempest2 <- as.data.frame(est.bin2[[char]])
                tempest3 <- as.data.frame(est.bin3[[char]])
                p1 <- p1 + geom_errorbar(
                    data = tempest2, aes(x = x0, ymin = CI.lower, ymax = CI.upper), color = platte[k],
                    width = errorbar.width / 3
                ) +
                    geom_point(data = tempest2, aes(x = x0, y = coef), size = 3, shape = 21, fill = platte[k], color = platte[k])

                if (dim(tempest3)[1] != 0) {
                    p1 <- p1 + geom_text(data = tempest3, aes(x = x0, y = 0), label = "NaN", color = platte[k])
                }
                k <- k + 1
            }
        }

        if (treat.type == "continuous") {
            for (label in label.name) {
                est.touse <- est.lin[[label]]
                if (CI == TRUE) {
                    colnames(est.touse) <- c("X", "ME", "sd", "CI_lower", "CI_upper")
                }
                if (CI == FALSE) {
                    est.touse <- est.touse[, c(1, 2)]
                    colnames(est.touse) <- c("X", "ME")
                }
                est.touse <- as.data.frame(est.touse)
                est.touse[["Treatment"]] <- rep(label, dim(est.touse)[1])

                if (label == label.name[1]) {
                    tempest <- est.touse
                } else {
                    tempest <- rbind(tempest, est.touse)
                }
            }
            tempest$Treatment <- factor(tempest$Treatment, levels = label.name)
            p1 <- p1 + geom_line(data = tempest, aes(x = X, y = ME, color = Treatment), show.legend = FALSE)
            p1 <- p1 + scale_color_manual(values = platte, labels = subtitles)
            if (CI == TRUE) {
                p1 <- p1 + geom_ribbon(data = tempest, aes(x = X, ymin = CI_lower, ymax = CI_upper, fill = Treatment), alpha = 0.2, show.legend = FALSE)
                p1 <- p1 + scale_fill_manual(values = platte, labels = subtitles)
            }

            k <- 1
            for (label in label.name) {
                tempest2 <- as.data.frame(est.bin2[[label]])
                tempest3 <- as.data.frame(est.bin3[[label]])

                p1 <- p1 + geom_errorbar(
                    data = tempest2, aes(x = x0, ymin = CI.lower, ymax = CI.upper), color = platte[k],
                    width = errorbar.width / 3
                ) +
                    geom_point(data = tempest2, aes(x = x0, y = coef), size = 3, shape = 21, fill = platte[k], color = platte[k])

                if (dim(tempest3)[1] != 0) {
                    p1 <- p1 + geom_text(data = tempest3, aes(x = x0, y = 0), label = "NaN", color = platte[k])
                }
                k <- k + 1
            }
        }

        if (bin.labs == TRUE) {
            if (treat.type == "discrete") {
                char0 <- other.treat[1]
            }
            if (treat.type == "continuous") {
                char0 <- label.name[1]
            }
            if (nbins == 3) {
                p1 <- p1 + annotate(
                    geom = "text", x = est.bin[[char0]][1, 1], y = pos,
                    label = "L", colour = "gray50", size = 10
                ) +
                    annotate(
                        geom = "text", x = est.bin[[char0]][2, 1], y = pos,
                        label = "M", colour = "gray50", size = 10
                    ) +
                    annotate(
                        geom = "text", x = est.bin[[char0]][3, 1], y = pos,
                        label = "H", colour = "gray50", size = 10
                    )
            } else if (nbins == 4) {
                p1 <- p1 + annotate(
                    geom = "text", x = est.bin[[char0]][1, 1], y = pos,
                    label = "L", colour = "gray50", size = 10
                ) +
                    annotate(
                        geom = "text", x = est.bin[[char0]][2, 1], y = pos,
                        label = "M1", colour = "gray50", size = 10
                    ) +
                    annotate(
                        geom = "text", x = est.bin[[char0]][3, 1], y = pos,
                        label = "M2", colour = "gray50", size = 10
                    ) +
                    annotate(
                        geom = "text", x = est.bin[[char0]][4, 1], y = pos,
                        label = "H", colour = "gray50", size = 10
                    )
            } else if (nbins == 2) {
                p1 <- p1 + annotate(
                    geom = "text", x = est.bin[[char0]][1, 1], y = pos,
                    label = "L", colour = "gray50", size = 10
                ) +
                    annotate(
                        geom = "text", x = est.bin[[char0]][2, 1], y = pos,
                        label = "H", colour = "gray50", size = 10
                    )
            }
        }
    }

    if (Xdistr == "density") { # density plot
        if (treat.type == "discrete") {
            ## put in data frames
            dist <- hist.out$mids[2] - hist.out$mids[1]
            deX.ymin <- min(yrange) - maxdiff / 5
            deX.co <- data.frame(
                x = de.tr[[base]]$x,
                y = de.tr[[base]]$y / max(de.tr[[base]]$y) * maxdiff / 10 + min(yrange) - maxdiff / 5
            )
            ## color
            p1 <- p1 + geom_ribbon(
                data = deX.co, aes(x = x, ymax = y, ymin = deX.ymin), color = base.color,
                fill = base.color, alpha = 0.0, size = 0.3
            )
            k <- 1
            char0 <- other.treat[1]
            start_level <- rep(deX.ymin, length(de.tr[[char0]]$x))
            for (char in other.treat) {
                dex.tr.plot <- data.frame(
                    x = de.tr[[char]]$x,
                    start_level = start_level,
                    end_level = de.tr[[char]]$y / max(de.tr[[char0]]$y) * maxdiff / 10 + start_level
                )

                p1 <- p1 + geom_ribbon(
                    data = dex.tr.plot, aes(x = x, ymax = end_level, ymin = start_level), color = platte[k],
                    alpha = 0.0, fill = platte[k], size = 0.3
                )

                k <- k + 1
            }
            p1 <- p1 + geom_line(data = dex.tr.plot, aes(x = x, y = ymin), color = "gray50", size = 0.3)
        }

        if (treat.type == "continuous") {
            deX.ymin <- min(yrange) - maxdiff / 5
            deX <- data.frame(
                x = de$x,
                y = de$y / max(de$y) * maxdiff / 5 + min(yrange) - maxdiff / 5
            )

            ## color
            feed.col <- col2rgb("gray50")
            col <- rgb(feed.col[1] / 1000, feed.col[2] / 1000, feed.col[3] / 1000)
            p1 <- p1 + geom_ribbon(
                data = deX, aes(x = x, ymax = y, ymin = deX.ymin),
                fill = col, alpha = 0.2
            )
        }
    }

    if (Xdistr %in% c("histogram", "hist")) { # histogram plot
        if (treat.type == "discrete") {
            n.hist <- length(hist.out$mids)
            dist <- hist.out$mids[2] - hist.out$mids[1]
            hist.max <- max(hist.out$counts)
            hist.col <- data.frame(
                ymin = rep(min(yrange) - maxdiff / 5, n.hist),
                ymax = hist.out$counts / hist.max * maxdiff / 5 + min(yrange) - maxdiff / 5,
                xmin = hist.out$mids - dist / 2,
                xmax = hist.out$mids + dist / 2,
                count1 = count.tr[[base]] / hist.max * maxdiff / 5 + min(yrange) - maxdiff / 5
            )

            p1 <- p1 + geom_rect(
                data = hist.col, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = count1), fill = base.color, color = "gray50",
                alpha = 0.3, size = 0.5
            )

            k <- 1
            start_level <- count.tr[[base]] / hist.max * maxdiff / 5 + min(yrange) - maxdiff / 5
            for (char in other.treat) {
                hist.treat <- data.frame(
                    ymin = start_level,
                    ymax = count.tr[[char]] / hist.max * maxdiff / 5 + start_level,
                    xmin = hist.out$mids - dist / 2,
                    xmax = hist.out$mids + dist / 2
                )

                start_level <- count.tr[[char]] / hist.max * maxdiff / 5 + start_level

                p1 <- p1 + geom_rect(
                    data = hist.treat, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = platte[k], color = "gray50",
                    alpha = 0.5, size = 0.5
                )
                k <- k + 1
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
            p1 <- p1 + geom_rect(
                data = histX, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                colour = "gray50", alpha = 0.3, size = 0.5
            )
        }
    }


    ## other properties
    if (is.null(legend.title) == FALSE) {
        p1 <- p1 + labs(fill = legend.title, color = legend.title)
    }

    ## mark the original interval (in replicated papers)
    if (is.null(interval) == FALSE) {
        p1 <- p1 + geom_vline(xintercept = interval, colour = "steelblue", linetype = 2, size = 1.5)
    }

    ## Other universal options
    ## axis labels
    if (is.null(cex.lab) == TRUE) {
        cex.lab <- 15
    } else {
        cex.lab <- 15 * cex.lab
    }
    if (is.null(cex.axis) == TRUE) {
        cex.axis <- 15
    } else {
        cex.axis <- 15 * cex.axis
    }
    p1 <- p1 + xlab(xlab) + ylab(ylab) +
        theme(axis.text = element_text(size = cex.axis), axis.title = element_text(size = cex.lab))

    ## title
    if (is.null(cex.main) == TRUE) {
        cex.main <- 18
    } else {
        cex.main <- 18 * cex.main
    }

    if (is.null(cex.sub) == TRUE) {
        cex.sub <- 10
    } else {
        cex.sub <- 10 * cex.sub
    }

    if (is.null(main) == FALSE) {
        p1 <- p1 + ggtitle(main) +
            theme(plot.title = element_text(
                hjust = 0.5, size = cex.main,
                lineheight = .8, face = "bold"
            ))
    }

    ## xlim and ylim
    if (is.null(ylim) == FALSE) {
        ylim2 <- c(ylim[1] - (ylim[2] - ylim[1]) * 0.25 / 6, ylim[2] + (ylim[2] - ylim[1]) * 0.4 / 6)
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

    # legend
    if (show.subtitles == TRUE) {
        if (treat.type == "discrete") {
            p1_table <- ggplot_gtable(ggplot_build(p1))
            data.touse3 <- data.frame(X = rep(1, length(all.treat)), ymin = -1, ymax = 1, D = all.treat)
            data.touse3$D <- factor(data.touse3$D, levels = all.treat)
            p0 <- ggplot() +
                geom_ribbon(
                    data = data.touse3, aes(x = X, ymin = ymin, ymax = ymax, fill = D),
                    alpha = 0.3
                )
            p0 <- p0 + scale_fill_manual(
                values = c(base.color, platte[1:length(other.treat)]),
                labels = as.character(subtitles.all)
            )
            if (is.null(legend.title) == FALSE) {
                p0 <- p0 + labs(fill = legend.title, color = legend.title)
            } else {
                p0 <- p0 + labs(fill = "Treatment", color = "Treatment")
            }

            p0 <- p0 + theme(
                legend.title = element_text(colour = "black", size = cex.sub),
                legend.text = element_text(color = "black", size = cex.sub * 0.95)
            )

            p0 <- p0 + xlab(xlab) + ylab(ylab) + theme(axis.text = element_text(size = cex.axis), axis.title = element_text(size = cex.lab))

            if (is.null(main) == FALSE) {
                p0 <- p0 + ggtitle(main) + theme(plot.title = element_text(
                    hjust = 0.5, size = cex.main,
                    lineheight = .8, face = "bold"
                ))
            }

            if (is.null(ylim) == FALSE) {
                ylim2 <- c(ylim[1] - (ylim[2] - ylim[1]) * 0.25 / 6, ylim[2] + (ylim[2] - ylim[1]) * 0.4 / 6)
            }
            if (is.null(xlim) == FALSE & is.null(ylim) == FALSE) {
                p0 <- p0 + coord_cartesian(xlim = xlim, ylim = ylim2)
            }
            if (is.null(xlim) == TRUE & is.null(ylim) == FALSE) {
                p0 <- p0 + coord_cartesian(ylim = ylim2)
            }
            if (is.null(xlim) == FALSE & is.null(ylim) == TRUE) {
                p0 <- p0 + coord_cartesian(xlim = xlim)
            }
            y.limits <- layer_scales(p1)$y$range$range
            x.limits <- layer_scales(p1)$x$range$range
            ymaxmax <- y.limits[2]
            yminmin <- y.limits[1]
            xmaxmax <- x.limits[2]
            xminmin <- x.limits[1]

            suppressWarnings(
                p0 <- p0 + ylim(c(yminmin, ymaxmax)) + xlim(c(xminmin, xmaxmax))
            )
            suppressWarnings(
                p0 <- ggplot_gtable(ggplot_build(p0))
            )
            suppressWarnings(
                pp <- c(subset(p0$layout, name == "panel", se = t:r))
            )

            gt <- gtable_add_grob(
                p0,
                p1_table$grobs[[which(p1_table$layout$name == "panel")]],
                pp$t, pp$l, pp$b, pp$l
            )

            gt <- as.ggplot(gt)
            p1 <- gt
        }

        if (treat.type == "continuous") {
            p1_table <- ggplot_gtable(ggplot_build(p1))
            data.touse3 <- data.frame(X = rep(1, length(label.name)), ymin = -1, ymax = 1, D = label.name)
            data.touse3$D <- factor(data.touse3$D, levels = label.name)
            p0 <- ggplot() +
                geom_ribbon(
                    data = data.touse3, aes(x = X, ymin = ymin, ymax = ymax, fill = D),
                    alpha = 0.3
                )
            p0 <- p0 + scale_fill_manual(
                values = platte[1:length(label.name)],
                labels = as.character(subtitles.all)
            )
            if (is.null(legend.title) == FALSE) {
                p0 <- p0 + labs(fill = legend.title, color = legend.title)
            } else {
                p0 <- p0 + labs(fill = "Treatment", color = "Treatment")
            }

            p0 <- p0 + theme(
                legend.title = element_text(colour = "black", size = cex.sub),
                legend.text = element_text(color = "black", size = cex.sub * 0.95)
            )

            p0 <- p0 + xlab(xlab) + ylab(ylab) + theme(axis.text = element_text(size = cex.axis), axis.title = element_text(size = cex.lab))

            if (is.null(main) == FALSE) {
                p0 <- p0 + ggtitle(main) + theme(plot.title = element_text(
                    hjust = 0.5, size = cex.main,
                    lineheight = .8, face = "bold"
                ))
            }

            ## xlim and ylim
            if (is.null(ylim) == FALSE) {
                ylim2 <- c(ylim[1] - (ylim[2] - ylim[1]) * 0.25 / 6, ylim[2] + (ylim[2] - ylim[1]) * 0.4 / 6)
            }
            if (is.null(xlim) == FALSE & is.null(ylim) == FALSE) {
                p0 <- p0 + coord_cartesian(xlim = xlim, ylim = ylim2)
            }
            if (is.null(xlim) == TRUE & is.null(ylim) == FALSE) {
                p0 <- p0 + coord_cartesian(ylim = ylim2)
            }
            if (is.null(xlim) == FALSE & is.null(ylim) == TRUE) {
                p0 <- p0 + coord_cartesian(xlim = xlim)
            }
            y.limits <- layer_scales(p1)$y$range$range
            x.limits <- layer_scales(p1)$x$range$range
            ymaxmax <- y.limits[2]
            yminmin <- y.limits[1]
            xmaxmax <- x.limits[2]
            xminmin <- x.limits[1]
            suppressWarnings(
                p0 <- p0 + ylim(c(yminmin, ymaxmax)) + xlim(c(xminmin, xmaxmax))
            )
            suppressWarnings(
                p0 <- ggplot_gtable(ggplot_build(p0))
            )
            suppressWarnings(
                pp <- c(subset(p0$layout, name == "panel", se = t:r))
            )

            gt <- gtable_add_grob(
                p0,
                p1_table$grobs[[which(p1_table$layout$name == "panel")]],
                pp$t, pp$l, pp$b, pp$l
            )

            gt <- as.ggplot(gt)
            p1 <- gt
        }
    }

    ## save to file
    if (is.null(file) == FALSE) {
        ggsave(file, p1, scale = scale, width = width, height = height)
    }

    return(p1)
}
