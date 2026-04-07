# raw
interflex.raw <- function(data,
                          Y,
                          D,
                          X,
                          treat.info,
                          order = NULL,
                          subtitles = NULL,
                          weights = NULL,
                          nbins = 3,
                          cutoffs = NULL,
                          span = NULL,
                          main = NULL,
                          Ylabel = NULL,
                          Dlabel = NULL,
                          Xlabel = NULL,
                          theme.bw = TRUE,
                          show.grid = TRUE,
                          cex.main = NULL,
                          cex.lab = NULL,
                          cex.axis = NULL,
                          ncols = NULL,
                          file = NULL,
                          scale = 1.1,
                          height = 7,
                          width = 10,
                          box.pos = "down",
                          xlim = NULL,
                          ylim = NULL) {
    n <- dim(data)[1]
    ti <- .extract_treat_info(treat.info)
    treat.type <- ti$treat.type
    if (treat.type == "discrete") {
        other.treat <- ti$other.treat
        other.treat.origin <- ti$other.treat.origin
        all.treat <- ti$all.treat
        all.treat.origin <- ti$all.treat.origin
        base <- ti$base
    }
    if (treat.type == "continuous") {
        D.sample <- ti$D.sample
        label.name <- ti$label.name
    }
    ncols <- ti$ncols

    if (!is.null(span)) {
        if (!is.numeric(span)) {
            stop("\"span\" is not numeric.\n")
        } else {
            span <- span[1]
            if (span <= 0) {
                stop("\"span\" is not a positive number.\n")
            }
        }
    }

    data <- na.omit(data[, c(Y, D, X)])

    if (treat.type == "discrete") {
        for (i in 1:length(all.treat)) {
            data[which(data[, D] == all.treat[i]), D] <- all.treat.origin[i]
        }
        all.treat <- unique(data[, D])
        if (!is.null(order)) {
            order <- as.character(order)
            if (length(order) != length(unique(order))) {
                stop("\"order\" should not contain repeated values.")
            }
            if (length(order) != length(all.treat)) {
                stop("\"order\" should contain all kinds of treatment arms.")
            }
            if (sum(!is.element(order, all.treat)) != 0 | sum(!is.element(all.treat, order)) != 0) {
                stop("\"order\" should contain all kinds of treatment arms.")
            }
        } else {
            order <- sort(all.treat)
        }
    }

    ## plotting
    if (treat.type == "discrete") { ## discrte
        all_treatment <- as.character(unique(data[, D]))
        data.aug <- data
        data.aug[, "treat"] <- paste("Treatment =", data[, D])
        order.lab <- paste("Treatment =", order)

        if (is.null(subtitles)) {
            subtitles <- order.lab
        }

        if (box.pos %in% c("down", "1", 1)) {
            box.pos <- min(data[, Y])
        } else if (box.pos %in% c("up", "2", 2)) {
            box.pos <- max(data[, Y])
        } else {
            stop("\"box.pos\" should only be \"down/1\" or \"up/2\".")
        }

        data.aug$qt90 <- NA
        data.aug$qt50 <- NA
        data.aug$med <- NA
        for (char in all_treatment) {
            treat <- paste0("Treat", char)
            treat_90 <- paste0(treat, "90")
            treat_50 <- paste0(treat, "50")
            treat_med <- paste0(treat, "med")
            assign(treat, which(data.aug[, D] == char))
            assign(treat_90, quantile(data[get(treat), X], c(0.05, 0.95)))
            assign(treat_50, quantile(data[get(treat), X], c(0.25, 0.75)))
            assign(treat_med, median(data[get(treat), X]))
            data.aug$qt90[which(data.aug[, D] == char &
                data.aug[, X] >= get(treat_90)[1] &
                data.aug[, X] <= get(treat_90)[2])] <- box.pos
            data.aug$qt50[which(data.aug[, D] == char &
                data.aug[, X] >= get(treat_50)[1] &
                data.aug[, X] <= get(treat_50)[2])] <- box.pos
            data.aug$med[get(treat)[which.min(abs(data.aug[get(treat), X] - get(treat_med)))]] <- box.pos
        }

        ## plotting
        if (is.null(weights)) {
            p1 <- ggplot(transform(data.aug, treat = factor(treat, levels = order.lab, labels = subtitles)), aes(.data[[X]], .data[[Y]]))
        } else {
            p1 <- ggplot(transform(data.aug, treat = factor(treat, levels = order.lab, labels = subtitles)), aes(.data[[X]], .data[[Y]], weight = .data[[weights]]))
        }
        if (theme.bw) {
            p1 <- p1 + theme_bw()
        }
        if (!show.grid) {
            p1 <- p1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        }
        p1 <- p1 + geom_point(na.rm = TRUE) + geom_smooth(method = "lm", formula = y ~ x, se = F, fullrange = T, colour = "steelblue", linewidth = 1, na.rm = TRUE)

        if (is.null(span)) {
            p1 <- p1 + geom_smooth(method = "loess", formula = y ~ x, se = F, colour = "red", na.rm = TRUE)
        } else {
            p1 <- p1 + geom_smooth(method = "loess", formula = y ~ x, se = F, colour = "red", span = span, na.rm = TRUE)
        }
        p1 <- p1 + xlab(Xlabel) + ylab(Ylabel)

        p1 <- p1 + geom_line(aes(.data[[X]], .data[["qt90"]]), linewidth = 1, colour = "grey50", na.rm = TRUE)
        p1 <- p1 + geom_line(aes(.data[[X]], .data[["qt50"]]), linewidth = 3, colour = "grey50", na.rm = TRUE)
        p1 <- p1 + geom_point(aes(.data[[X]], .data[["med"]]), size = 3, shape = 21, fill = "white", colour = "red", na.rm = TRUE)
        p1 <- p1 + facet_wrap(~treat, ncol = ncols)
    } else { # continuous case
        if (!is.numeric(data[, D])) {
            stop("\"D\" is not a numeric variable")
        }

        ## grouping by X
        if (is.null(cutoffs)) {
            cutoff <- quantile(data[, X], probs = seq(0, 1, 1 / nbins))
            while (length(unique(cutoff)) != nbins + 1) {
                nbins <- nbins - 1
                cutoff <- quantile(data[, X], probs = seq(0, 1, 1 / nbins))
            }
        } else {
            cutoffs <- cutoffs[which(cutoffs > min(data[, X]) & cutoffs < max(data[, X]))]
            cutoff <- sort(unique(c(min(data[, X]), cutoffs, max(data[, X]))))
        }
        groupID <- cut(data[, X], breaks = cutoff, labels = FALSE)
        groupID[which(data[, X] == min(data[, X]))] <- 1

        ## X labels
        groupID2 <- cut(data[, X], breaks = cutoff)
        gp.lab <- paste(Xlabel, ": ", levels(groupID2), sep = "")
        gp.lab[1] <- paste(Xlabel, ": [", substring(levels(groupID2)[1], 2), sep = "")
        nbins <- length(unique(groupID))
        groupID <- factor(groupID, labels = gp.lab)
        data.aug <- data
        data.aug$groupID <- groupID

        ## plotting
        if (is.null(weights)) {
            p1 <- ggplot(data.aug, aes(.data[[D]], .data[[Y]]))
        } else {
            p1 <- ggplot(data.aug, aes(.data[[D]], .data[[Y]], weight = .data[[weights]]))
        }
        if (theme.bw) {
            p1 <- p1 + theme_bw()
        }
        if (!show.grid) {
            p1 <- p1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        }
        p1 <- p1 + geom_point(na.rm = TRUE) + geom_smooth(formula = y ~ x, method = "lm", se = F, fullrange = T, colour = "steelblue", linewidth = 1, na.rm = TRUE)

        if (is.null(span)) {
            p1 <- p1 + geom_smooth(method = "loess", formula = y ~ x, se = F, colour = "red", na.rm = TRUE)
        } else {
            p1 <- p1 + geom_smooth(method = "loess", formula = y ~ x, se = F, colour = "red", span = span, na.rm = TRUE)
        }
        p1 <- p1 + xlab(Dlabel) + ylab(Ylabel) + facet_wrap(. ~ groupID, ncol = ncols)
    }

    ## axis labels
    if (is.null(cex.lab)) {
        cex.lab <- 15
    } else {
        cex.lab <- 15 * cex.lab
    }
    if (is.null(cex.axis)) {
        cex.axis <- 15
    } else {
        cex.axis <- 15 * cex.axis
    }
    p1 <- p1 + theme(axis.text = element_text(size = cex.axis), axis.title = element_text(size = cex.lab))

    ## title
    if (is.null(cex.main)) {
        cex.main <- 18
    } else {
        cex.main <- 18 * cex.main
    }
    if (!is.null(main)) {
        p1 <- p1 + ggtitle(main) +
            theme(plot.title = element_text(hjust = 0.5, size = cex.main, lineheight = .8, face = "bold"))
    }

    if (!is.null(xlim)) {
        p1 <- p1 + xlim(xlim[1], xlim[2])
    }

    if (!is.null(ylim)) {
        p1 <- p1 + ylim(ylim[1], ylim[2])
    }

    ## save to file
    if (!is.null(file)) {
        ggsave(file, p1, scale = scale, width = width, height = height)
    }

    return(p1)
}
