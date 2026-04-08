# Internal: add a small symmetric pad to a user-supplied xlim so curves do not
# clip exactly at the panel edge. Returns `xlim` unchanged when NULL or when
# the range is degenerate (non-finite or zero-width).
.pad_xlim <- function(xlim, mult = 0.04) {
    if (is.null(xlim)) return(NULL)
    if (length(xlim) != 2L) return(xlim)
    if (any(!is.finite(xlim))) return(xlim)
    span <- xlim[2] - xlim[1]
    if (!is.finite(span) || span <= 0) return(xlim)
    pad <- span * mult
    c(xlim[1] - pad, xlim[2] + pad)
}

# Internal: defensively append CI columns from an estimator output table to a
# running yrange vector. Mirrors the existing `ncol(...) > 5` guard pattern but
# also guards the pointwise CI cols (4, 5). When the table is too narrow to
# carry CIs, falls back to the point-estimate column (col 2), matching what the
# CI=FALSE branch already does upstream.
.append_yrange_ci <- function(yrange, tab) {
    if (is.null(tab)) return(yrange)
    nc <- ncol(tab)
    if (is.null(nc) || nc < 2L) return(yrange)
    if (nc >= 5L) {
        yrange <- c(yrange, na.omit(unlist(c(tab[, c(4, 5)]))))
        if (nc >= 7L) {
            yrange <- c(yrange, na.omit(unlist(c(tab[, c(6, 7)]))))
        }
    } else {
        yrange <- c(yrange, na.omit(unlist(c(tab[, 2]))))
    }
    yrange
}

# Internal: defensively rename the columns of an estimator output table to the
# canonical (X, <point>, sd, CI_lower, CI_upper, [CI_uniform_lower, CI_uniform_upper])
# schema. Tolerates narrow tables (2-4 cols, e.g. a DML output that carries only
# X / point / sd with no CI columns) by assigning only the names that exist.
# Mirrors the .append_yrange_ci guard so that the colnames-rename site does not
# crash on the same narrow tables that helper was added to handle.
# `point` is the label for column 2 ("TE" or "ME").
.rename_est_ci <- function(tab, point = "TE") {
    if (is.null(tab)) return(tab)
    nc <- ncol(tab)
    if (is.null(nc) || nc < 1L) return(tab)
    full <- c("X", point, "sd", "CI_lower", "CI_upper", "CI_uniform_lower", "CI_uniform_upper")
    if (nc >= 7L) {
        colnames(tab) <- full[1:7]
    } else if (nc >= 5L) {
        colnames(tab) <- full[1:5]
    } else {
        # Narrow table (2-4 cols): assign only the names that exist. No CI ribbon
        # columns are created -- downstream CI/ribbon code is gated on
        # `"CI_lower" %in% colnames(tab)` so a narrow table degrades gracefully
        # to a curve-only plot.
        colnames(tab) <- full[seq_len(nc)]
    }
    tab
}

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
                           theme.bw = TRUE,
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
                           by.group = FALSE, # if g.est is in the output
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
                           show.uniform.CI = TRUE,
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



    if (pool) {
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
    if (!inherits(out, "interflex")) {
        stop("Not an \"interflex\" object.")
    }

    ## Snapshot user-supplied limits before any downstream code mutates them
    ## (e.g. defaults derived from the data range). These snapshots are stamped
    ## onto the returned graph as attributes so a downstream re-plot can recover
    ## the original window the user requested.
    ## Distinguish user-supplied limits from auto-trim defaults applied
    ## upstream by interflex(). When the option flags are unset (e.g. when
    ## plot.interflex is invoked directly on an existing object), treat any
    ## non-NULL limit as user-supplied (backward compatible).
    .xlim_explicit <- getOption("interflex.user_xlim_explicit", default = NA)
    .ylim_explicit <- getOption("interflex.user_ylim_explicit", default = NA)
    .user_xlim_in <- if (isTRUE(.xlim_explicit) || is.na(.xlim_explicit)) xlim else NULL
    .user_ylim_in <- if (isTRUE(.ylim_explicit) || is.na(.ylim_explicit)) ylim else NULL

    ## Recover user-supplied limits from a previously-built figure when the
    ## caller re-plots an existing interflex object without re-passing them.
    ## plot.interflex stamps the original xlim/ylim onto the returned graph as
    ## attributes (see end of function), so a downstream re-plot via
    ## plot.interflex(out, show.all = TRUE) honors the same window the user
    ## originally requested.
    if (is.null(xlim) && !is.null(out$figure)) {
        prior_xlim <- attr(out$figure, "interflex_xlim")
        if (!is.null(prior_xlim)) {
            xlim <- prior_xlim
            .user_xlim_in <- prior_xlim
        }
    }
    if (is.null(ylim) && !is.null(out$figure)) {
        prior_ylim <- attr(out$figure, "interflex_ylim")
        if (!is.null(prior_ylim)) {
            ylim <- prior_ylim
            .user_ylim_in <- prior_ylim
        }
    }

    ## --- PAD-001 PASS 2: plot-time row filter for Path 2 ---
    ## When the user supplies xlim at plot time on a fit built without xlim,
    ## the stored est.<estimator> tables span the full data range. We drop
    ## rows whose moderator value lies outside [lo, hi] BEFORE the rest of
    ## plot.interflex consumes them, so geom_line / geom_ribbon /
    ## geom_path / uniform-CI layers all see the restricted window without
    ## any post-build layer manipulation. We never re-fit; we never modify
    ## `out` in place. PASS 1 (fit-time grid restriction) is idempotent
    ## with PASS 2 -- when both fired the filter is a no-op.
    .filter_xlim_rows <- function(tbl, lo, hi) {
        if (is.null(tbl)) return(tbl)
        if (is.null(dim(tbl))) return(tbl)
        if (!"X" %in% colnames(tbl)) return(tbl)
        eps <- 1e-9 * max(1, abs(hi - lo))
        keep <- which(tbl[, "X"] >= lo - eps & tbl[, "X"] <= hi + eps)
        if (length(keep) == 0L) return(tbl)        # never empty out the table
        tbl[keep, , drop = FALSE]
    }
    .filter_xlim_field <- function(field, lo, hi) {
        if (is.null(field)) return(field)
        if (is.list(field) && is.null(dim(field))) {
            return(lapply(field, .filter_xlim_rows, lo = lo, hi = hi))
        }
        .filter_xlim_rows(field, lo, hi)
    }
    ## PAD-001 PASS 2b: filter density objects (with $x and $y vectors,
    ## as produced by stats::density) so continuous density ribbons do
    ## not extend past the user xlim window.
    .filter_xlim_density <- function(dens, lo, hi) {
        if (is.null(dens)) return(dens)
        if (!is.list(dens)) return(dens)
        if (is.null(dens$x) || is.null(dens$y)) return(dens)
        eps <- 1e-9 * max(1, abs(hi - lo))
        keep <- which(dens$x >= lo - eps & dens$x <= hi + eps)
        if (length(keep) == 0L) return(dens)  # empty guard
        dens$x <- dens$x[keep]
        dens$y <- dens$y[keep]
        dens
    }
    ## PAD-001 PASS 2b: filter histogram object (hist.out has $mids,
    ## $counts, possibly $breaks/$density). count.tr (per-treatment
    ## count vectors) is index-aligned with mids and must be filtered
    ## with the same mask -- but we only filter it here; the continuous
    ## branch does not use count.tr (discrete branch does), so touching
    ## count.tr is a safety measure consistent with the mask.
    .filter_xlim_histlike <- function(h, lo, hi) {
        if (is.null(h)) return(h)
        if (!is.list(h) || is.null(h$mids)) return(h)
        eps <- 1e-9 * max(1, abs(hi - lo))
        keep <- which(h$mids >= lo - eps & h$mids <= hi + eps)
        if (length(keep) == 0L) return(h)  # empty guard
        if (!is.null(h$counts) && length(h$counts) == length(h$mids)) {
            h$counts <- h$counts[keep]
        }
        if (!is.null(h$density) && length(h$density) == length(h$mids)) {
            h$density <- h$density[keep]
        }
        h$mids <- h$mids[keep]
        h
    }
    .pad_treat_type <- out$treat.info[["treat.type"]]
    .pad_xlim_gate <- (!is.null(.pad_treat_type)
        && .pad_treat_type == "continuous"
        && !is.null(.user_xlim_in)
        && is.numeric(.user_xlim_in)
        && length(.user_xlim_in) == 2L
        && all(is.finite(.user_xlim_in))
        && .user_xlim_in[1] < .user_xlim_in[2])
    .out_filt <- out
    if (isTRUE(.pad_xlim_gate)) {
        .lo <- .user_xlim_in[1]
        .hi <- .user_xlim_in[2]
        for (.fld in c("est.lin", "est.bin", "est.kernel", "est.dml",
                       "est.lasso", "est.grf")) {
            if (.fld %in% names(.out_filt)) {
                .out_filt[[.fld]] <- .filter_xlim_field(.out_filt[[.fld]], .lo, .hi)
            }
        }
        ## PAD-001 PASS 2b: also filter bar/density/histogram source
        ## fields so continuous Xdistr layers (geom_ribbon on de, geom_rect
        ## on hist.out) are bounded to the user xlim window. Discrete
        ## branches use out$de.tr / out$count.tr directly and are untouched.
        if ("de" %in% names(.out_filt)) {
            .out_filt$de <- .filter_xlim_density(.out_filt$de, .lo, .hi)
        }
        if ("hist.out" %in% names(.out_filt)) {
            .out_filt$hist.out <- .filter_xlim_histlike(.out_filt$hist.out, .lo, .hi)
        }
    }
    ## --- end PAD-001 PASS 2 hoisted filter ---

    ## Auto by.group when gate estimates are available and non-empty
    .has_gate <- function(obj, field) {
        if (!field %in% names(obj)) return(FALSE)
        g <- obj[[field]]
        is.list(g) && length(g) > 0 && any(vapply(g, NROW, 0L) > 0L)
    }
    if (!by.group && (.has_gate(out, "g.est") || .has_gate(out, "g.est.dml"))) {
        by.group <- TRUE
    }
    if (by.group) {
        if (!"g.est" %in% names(out) && !"g.est.dml" %in% names(out)) {
            stop("Group-specific Average Treatment Effects have to be estimated first. Use gate = TRUE.\n")
        }
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

    if (estimator == "binning" | estimator == "linear" | estimator == "dml" | estimator == "grf" | estimator == "lasso") {
        if (is.null(CI)) {
            CI <- TRUE
        }
    }

    # Xdistr
    if (!Xdistr %in% c("hist", "histogram", "density", "none")) {
        stop("\"Xdistr\" must be \"histogram\", \"density\", or \"none\".")
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

    # Dlabel
    if (is.null(Dlabel)) {
        Dlabel <- out$Dlabel
    } else {
        if (!is.character(Dlabel)) {
            stop("\"Dlabel\" is not a string.")
        } else {
            Dlabel <- Dlabel[1]
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
    if (!is.null(xlab)) {
        if (!is.character(xlab)) {
            stop("\"xlab\" is not a string.")
        }
    }
    if (!is.null(ylab)) {
        if (!is.character(ylab)) {
            stop("\"ylab\" is not a string.")
        }
    }

    if (is.null(xlab)) {
        xlab <- c(paste("Moderator: ", Xlabel, sep = ""))
    } else {
        if (!is.character(xlab)) {
            stop("\"xlab\" is not a string.")
        }
    }
    if (is.null(ylab)) {
        if (by.group && (!is.null(out$g.est) || !is.null(out$g.est.dml))) {
            ylab <- paste("GATE of ", Dlabel, " on ", Ylabel, sep = "")
        } else {
            ylab <- paste("CME of ", Dlabel, " on ", Ylabel, sep = "")
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

    ## bin.labs
    if (!is.logical(bin.labs) & !is.numeric(bin.labs)) {
        stop("\"bin.labs\" is not a logical flag.")
    }

    ## interval
    if (!is.null(interval)) {
        if (!is.numeric(interval)) {
            stop("Some element in \"interval\" is not numeric.")
        }
    }

    ## file
    if (!is.null(file)) {
        if (!is.character(file)) {
            stop("Wrong file name.")
        }
    }

    ## order/subtitles
    if (treat.type == "discrete") {
        other.treat <- sort(all.treat[which(all.treat != base)])
        if (!is.null(order)) {
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

        if (is.null(show.subtitles)) {
            if (length(other.treat) == 1) {
                show.subtitles <- FALSE
            } else {
                show.subtitles <- TRUE
            }

            if (!is.null(subtitles)) {
                show.subtitles <- TRUE
            }
        }

        if (!is.null(subtitles)) {
            if (length(subtitles) != length(other.treat)) {
                stop("The number of elements in \"subtitles\" should be m-1(m is the number of different treatment arms).")
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
            label.name.order <- c()
            for (a in order) {
                label.name.order <- c(label.name.order, names(D.sample[which(D.sample == a)]))
            }
            label.name <- label.name.order
        }

        if (is.null(show.subtitles)) {
            if (length(label.name) == 1) {
                show.subtitles <- FALSE
            } else {
                show.subtitles <- TRUE
            }

            if (!is.null(subtitles)) {
                show.subtitles <- TRUE
            }
        }

        if (!is.null(subtitles)) {
            if (length(subtitles) != length(label.name)) {
                stop("The number of elements in \"subtitles\" should equal to the number of values in D.ref.")
            }
        }
    }

    ## ncols
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
    if (treat.type == "discrete" & estimator == "dml") {
        tempxx <- out$est.dml[[other.treat[1]]][, "X"]
    }
    if (treat.type == "discrete" & estimator == "grf") {
        tempxx <- out$est.grf[[other.treat[1]]][, "X"]
    }
    
    if (treat.type == "discrete" & estimator == "lasso") {
        tempxx <- out$est.lasso[[other.treat[1]]][, "X"]
    }
    if (treat.type == "discrete" & estimator == "kernel") {
        tempxx <- out$est.kernel[[other.treat[1]]][, "X"]
    }
    if (treat.type == "continuous" & (estimator == "linear" | estimator == "binning")) {
        tempxx <- out$est.lin[[label.name[1]]][, "X"]
    }
    if (treat.type == "continuous" & estimator == "dml") {
        tempxx <- out$est.dml[[label.name[1]]][, "X"]
    }
    if (treat.type == "continuous" & estimator == "lasso") {
        tempxx <- out$est.lasso[[label.name[1]]][, "X"]
    }
    if (treat.type == "continuous" & estimator == "kernel") {
        tempxx <- out$est.kernel[[label.name[1]]][, "X"]
    }

    min.XX <- min(tempxx)
    max.XX <- max(tempxx)

    if (!is.null(diff.values)) {
        if (estimator == "binning") {
            stop("\"diff.values\" can only work after linear or kernel model is applied.")
        }
        if (!is.numeric(diff.values)) {
            stop("\"diff.values\" is not numeric.")
        }
        if (length(diff.values) < 2) {
            stop("\"diff.values\" must be of length 2 or more.")
        }
        if (treat.type == "discrete" & (estimator == "linear")) {
            tempxx <- out$est.lin[[other.treat[1]]][, "X"]
        }
        if (treat.type == "discrete" & (estimator == "dml")) {
            tempxx <- out$est.dml[[other.treat[1]]][, "X"]
        }
        if (treat.type == "discrete" & (estimator == "grf")) {
            tempxx <- out$est.grf[[other.treat[1]]][, "X"]
        }
        if (treat.type == "discrete" & (estimator == "lasso")) {
            tempxx <- out$est.lasso[[other.treat[1]]][, "X"]
        }
        if (treat.type == "discrete" & estimator == "kernel") {
            tempxx <- out$est.kernel[[other.treat[1]]][, "X"]
        }
        if (treat.type == "continuous" & estimator == "linear") {
            tempxx <- out$est.lin[[label.name[1]]][, "X"]
        }
        if (treat.type == "continuous" & estimator == "dml") {
            tempxx <- out$est.dml[[label.name[1]]][, "X"]
        }
        if (treat.type == "continuous" & estimator == "lasso") {
            tempxx <- out$est.lasso[[label.name[1]]][, "X"]
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
                est.bin2[[char]] <- as.matrix(est.bin[[char]][which(!is.na(est.bin[[char]][, 2])), ])
                est.bin3[[char]] <- as.matrix(est.bin[[char]][which(is.na(est.bin[[char]][, 2])), ])
                if (dim(est.bin2[[char]])[2] == 1) {
                    est.bin2[[char]] <- t(est.bin2[[char]])
                }
                if (dim(est.bin3[[char]])[2] == 1) {
                    est.bin3[[char]] <- t(est.bin3[[char]])
                }

                if (isTRUE(CI)) {
                    yrange <- .append_yrange_ci(yrange, est.lin[[char]])
                    yrange <- .append_yrange_ci(yrange, est.bin[[char]])
                } else {
                    yrange <- c(yrange, na.omit(unlist(c(est.lin[[char]][, 2], est.bin[[char]][, 2]))))
                }
            }

            if (!is.null(ylim)) {
                yrange <- c(ylim[2], ylim[1] + (ylim[2] - ylim[1]) * 1 / 8)
            }
            X.lvls <- est.lin[[other.treat[1]]][, 1]
            errorbar.width <- (max(X.lvls) - min(X.lvls)) / 20
            maxdiff <- (max(yrange) - min(yrange))
            pos <- max(yrange) - maxdiff / 20
        }

        if (treat.type == "continuous") {
            est.lin <- .out_filt$est.lin   # PAD-001 PASS 2: filtered Path-2 view
            est.bin <- .out_filt$est.bin   # PAD-001 PASS 2: filtered Path-2 view
            est.bin2 <- list() ## non missing part
            est.bin3 <- list() ## missing part
            yrange <- c(0)
            for (label in label.name) {
                est.bin2[[label]] <- as.matrix(est.bin[[label]][which(!is.na(est.bin[[label]][, 2])), ])
                est.bin3[[label]] <- as.matrix(est.bin[[label]][which(is.na(est.bin[[label]][, 2])), ])
                if (dim(est.bin2[[label]])[2] == 1) {
                    est.bin2[[label]] <- t(est.bin2[[label]])
                }
                if (dim(est.bin3[[label]])[2] == 1) {
                    est.bin3[[label]] <- t(est.bin3[[label]])
                }
                if (isTRUE(CI)) {
                    yrange <- .append_yrange_ci(yrange, est.lin[[label]])
                    yrange <- .append_yrange_ci(yrange, est.bin[[label]])
                } else {
                    yrange <- c(yrange, na.omit(unlist(c(est.lin[[label]][, 2], est.bin[[label]][, 2]))))
                }
            }

            X.lvls <- est.lin[[label.name[1]]][, 1]
            errorbar.width <- (max(X.lvls) - min(X.lvls)) / 20
            if (!is.null(ylim)) {
                yrange <- c(ylim[2], ylim[1] + (ylim[2] - ylim[1]) * 1 / 8)
            }
            maxdiff <- (max(yrange) - min(yrange))
            pos <- max(yrange) - maxdiff / 20
        }
    }

    if (estimator == "dml") {
        if (treat.type == "discrete") {
            est.dml <- out$est.dml
            if (by.group) {
                est.dml <- if (!is.null(out$g.est)) out$g.est else out$g.est.dml
            }

            yrange <- c(0)
            for (char in other.treat) {
                if (isTRUE(CI)) {
                    yrange <- .append_yrange_ci(yrange, est.dml[[char]])
                } else {
                    yrange <- c(yrange, na.omit(unlist(c(est.dml[[char]][, 2]))))
                }
            }
            X.lvls <- est.dml[[other.treat[1]]][, 1]
        }
        if (treat.type == "continuous") {
            est.dml <- .out_filt$est.dml   # PAD-001 PASS 2: filtered Path-2 view
            if (by.group) {
                est.dml <- if (!is.null(out$g.est)) out$g.est else out$g.est.dml
            }
            yrange <- c(0)
            for (label in label.name) {
                if (isTRUE(CI)) {
                    yrange <- .append_yrange_ci(yrange, est.dml[[label]])
                } else {
                    yrange <- c(yrange, na.omit(unlist(c(est.dml[[label]][, 2]))))
                }
            }
            X.lvls <- est.dml[[label.name[1]]][, 1]
        }
        errorbar.width <- (max(X.lvls) - min(X.lvls)) / 20
        if (!is.null(ylim)) {
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
                if (isTRUE(CI)) {
                    yrange <- .append_yrange_ci(yrange, est.lin[[char]])
                } else {
                    yrange <- c(yrange, na.omit(unlist(c(est.lin[[char]][, 2]))))
                }
            }
            X.lvls <- est.lin[[other.treat[1]]][, 1]
        }

        if (treat.type == "continuous") {
            est.lin <- .out_filt$est.lin   # PAD-001 PASS 2: filtered Path-2 view
            yrange <- c(0)
            for (label in label.name) {
                if (isTRUE(CI)) {
                    yrange <- .append_yrange_ci(yrange, est.lin[[label]])
                } else {
                    yrange <- c(yrange, na.omit(unlist(c(est.lin[[label]][, 2]))))
                }
            }
            X.lvls <- est.lin[[label.name[1]]][, 1]
        }

        errorbar.width <- (max(X.lvls) - min(X.lvls)) / 20
        if (!is.null(ylim)) {
            yrange <- c(ylim[2], ylim[1] + (ylim[2] - ylim[1]) * 1 / 8)
        }
        maxdiff <- (max(yrange) - min(yrange))
        pos <- max(yrange) - maxdiff / 20
    }

    if (estimator == "kernel") {
        est.kernel <- .out_filt$est.kernel   # PAD-001 PASS 2: filtered Path-2 view (no-op for discrete)
        yrange <- c(0)
        if (isFALSE(CI)) {
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

        if (isTRUE(CI)) {
            if (treat.type == "discrete") {
                for (char in other.treat) {
                    yrange <- .append_yrange_ci(yrange, est.kernel[[char]])
                }
                X.lvls <- est.kernel[[other.treat[1]]][, 1]
            }

            if (treat.type == "continuous") {
                for (label in label.name) {
                    yrange <- .append_yrange_ci(yrange, est.kernel[[label]])
                }
                X.lvls <- est.kernel[[label.name[1]]][, 1]
            }
        }

        if (!is.null(ylim)) {
            yrange <- c(ylim[2], ylim[1] + (ylim[2] - ylim[1]) * 1 / 8)
        }
        errorbar.width <- (max(X.lvls) - min(X.lvls)) / 20
        maxdiff <- (max(yrange) - min(yrange))
        pos <- max(yrange) - maxdiff / 20
    }

    if (estimator == "grf") {
        if (treat.type == "discrete") {
            est.grf <- out$est.grf
            if (by.group) {
                est.grf <- out$est.grf
            }
            yrange <- c(0)
            for (char in other.treat) {
                if (isTRUE(CI)) {
                    yrange <- .append_yrange_ci(yrange, est.grf[[char]])
                } else {
                    yrange <- c(yrange, na.omit(unlist(c(est.grf[[char]][, 2]))))
                }
            }
            X.lvls <- est.grf[[other.treat[1]]][, 1]
        }
        errorbar.width <- (max(X.lvls) - min(X.lvls)) / 20
        if (!is.null(ylim)) {
            yrange <- c(ylim[2], ylim[1] + (ylim[2] - ylim[1]) * 1 / 8)
        }
        maxdiff <- (max(yrange) - min(yrange))
        pos <- max(yrange) - maxdiff / 20
    }

    if (estimator == "lasso") {
        if (treat.type == "discrete") {
            est.lasso <- out$est.lasso

            yrange <- c(0)
            for (char in other.treat) {
                if (isTRUE(CI)) {
                    yrange <- .append_yrange_ci(yrange, est.lasso[[char]])
                } else {
                    yrange <- c(yrange, na.omit(unlist(c(est.lasso[[char]][, 2]))))
                }
            }
            X.lvls <- est.lasso[[other.treat[1]]][, 1]
        }
        if (treat.type == "continuous") {
            est.lasso <- .out_filt$est.lasso   # PAD-001 PASS 2: filtered Path-2 view
            if (by.group) {
                est.lasso <- out$est.lasso
            }
            yrange <- c(0)
            for (label in label.name) {
                if (isTRUE(CI)) {
                    yrange <- .append_yrange_ci(yrange, est.lasso[[label]])
                } else {
                    yrange <- c(yrange, na.omit(unlist(c(est.lasso[[label]][, 2]))))
                }
            }
            X.lvls <- est.lasso[[label.name[1]]][, 1]
        }
        errorbar.width <- (max(X.lvls) - min(X.lvls)) / 20
        if (!is.null(ylim)) {
            yrange <- c(ylim[2], ylim[1] + (ylim[2] - ylim[1]) * 1 / 8)
        }
        maxdiff <- (max(yrange) - min(yrange))
        pos <- max(yrange) - maxdiff / 20
    }

    # Histogram height fraction: smaller for GATE plots to avoid overlap
    hist_frac <- if (by.group) 8 else 5

    # plot initialization
    p.group <- list()
    if (treat.type == "discrete") {
        for (char in other.treat) {
            p1 <- ggplot()
            ## black white theme and mark zero
            if (!theme.bw) {
                p1 <- p1 + geom_hline(yintercept = 0, colour = "white", linewidth = 2)
            } else {
                p1 <- p1 + theme_bw() + geom_hline(yintercept = 0, colour = "#AAAAAA50", linewidth = 2)
            }
            if (!show.grid) {
                p1 <- p1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
            }
            p.group[[char]] <- p1
        }
    }

    if (treat.type == "continuous") {
        for (label in label.name) {
            p1 <- ggplot()
            ## black white theme and mark zero
            if (!theme.bw) {
                p1 <- p1 + geom_hline(yintercept = 0, colour = "white", linewidth = 2)
            } else {
                p1 <- p1 + theme_bw() + geom_hline(yintercept = 0, colour = "#AAAAAA50", linewidth = 2)
            }
            if (!show.grid) {
                p1 <- p1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
            }
            p.group[[label]] <- p1
        }
    }

    # density
    if (Xdistr == "density") {
        if (treat.type == "discrete") {
            deX.ymin <- min(yrange) - maxdiff / hist_frac
            deX.co <- data.frame(
                x = de.tr[[base]]$x,
                y = de.tr[[base]]$y / max(de.tr[[base]]$y) * maxdiff / hist_frac + min(yrange) - maxdiff / hist_frac
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
                    y = de.tr[[char]]$y / max(de.tr[[char]]$y) * maxdiff / hist_frac + min(yrange) - maxdiff / hist_frac
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
            de <- .out_filt$de   # PAD-001 PASS 2b: filtered density view
            deX.ymin <- min(yrange) - maxdiff / hist_frac
            deX <- data.frame(
                x = de$x,
                y = de$y / max(de$y) * maxdiff / hist_frac + min(yrange) - maxdiff / hist_frac
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
                ymin = rep(min(yrange) - maxdiff / hist_frac, n.hist),
                # ymax=hist.out$counts/hist.max*maxdiff/5+min(yrange)-maxdiff/5,
                xmin = hist.out$mids - dist / 2,
                xmax = hist.out$mids + dist / 2,
                count1 = count.tr[[base]] / hist.max * maxdiff / hist_frac + min(yrange) - maxdiff / hist_frac
            )

            for (char in other.treat) {
                hist.treat <- data.frame(
                    ymin = hist.col[, "count1"],
                    # ymax=hist.out$counts/hist.max*maxdiff/5+min(yrange)-maxdiff/5,
                    xmin = hist.out$mids - dist / 2,
                    xmax = hist.out$mids + dist / 2,
                    count1 = count.tr[[char]] / hist.max * maxdiff / hist_frac + hist.col[, "count1"]
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
            .hist.out.full <- hist.out  # PAD-001 PASS 2b: keep full ref for bin width
            hist.out <- .out_filt$hist.out   # PAD-001 PASS 2b: filtered hist view
            n.hist <- length(hist.out$mids)
            dist <- if (length(hist.out$mids) >= 2L) {
                hist.out$mids[2] - hist.out$mids[1]
            } else if (length(.hist.out.full$mids) >= 2L) {
                .hist.out.full$mids[2] - .hist.out.full$mids[1]
            } else {
                1
            }
            hist.max <- max(hist.out$counts)
            histX <- data.frame(
                ymin = rep(min(yrange) - maxdiff / hist_frac, n.hist),
                ymax = hist.out$counts / hist.max * maxdiff / hist_frac + min(yrange) - maxdiff / hist_frac,
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
    if (estimator == "kernel" | estimator == "linear" | estimator == "dml" | estimator == "grf" | estimator == "lasso") {
        if (estimator == "kernel") {
            est <- est.kernel
        } else if (estimator == "dml") {
            est <- est.dml
        } else if (estimator == "grf") {
            est <- est.grf
        } else if (estimator == "lasso") {
            est <- est.lasso
        } else {
            est <- est.lin
        }

        # When by.group = TRUE and g.est is available, use GATE values
        # instead of the smooth curve -- GATE only has actual moderator values
        if (by.group && !is.null(out$g.est)) {
            est <- out$g.est
        } else if (by.group && !is.null(out$g.est.dml)) {
            est <- out$g.est.dml
        }

        if (treat.type == "discrete") {
            for (char in other.treat) {
                p1 <- p.group[[char]]
                tempest <- est[[char]]
                if (isTRUE(CI)) {
                    tempest <- .rename_est_ci(tempest, point = "TE")
                }
                if (isFALSE(CI)) {
                    tempest <- tempest[, c(1, 2)]
                    colnames(tempest) <- c("X", "TE")
                }
                tempest <- as.data.frame(tempest)
                if (!by.group) {
                    p1 <- p1 + geom_line(data = tempest, aes(X, TE), color = line.color, linewidth = line.size)
                }

                if (isTRUE(CI) && "CI_lower" %in% colnames(tempest)) {
                    if (estimator == "kernel" | estimator == "linear" | estimator == "dml" | estimator == "grf" | estimator == "lasso") {
                        if (!by.group) {
                            p1 <- p1 + geom_ribbon(
                                data = tempest, aes(x = X, ymin = CI_lower, ymax = CI_upper),
                                fill = CI.color, alpha = CI.color.alpha
                            )
                        } else {
                            p1 <- p1 + geom_errorbar(
                                data = tempest, aes(x = X, ymin = CI_lower, ymax = CI_upper), linewidth = line.size, width = 0.2
                            )
                        }
                    }

                    if ("CI_uniform_lower" %in% colnames(tempest) & show.uniform.CI) {
                        if (!by.group) {
                            p1 <- p1 + geom_line(data = tempest, aes(x = X, y = CI_uniform_lower), linetype = "dashed", color = "gray50") + geom_line(data = tempest, aes(x = X, y = CI_uniform_upper), linetype = "dashed", color = "gray50")
                        } else {
                            p1 <- p1 + geom_errorbar(
                                data = tempest, aes(x = X, ymin = CI_uniform_lower, ymax = CI_uniform_upper), linewidth = line.size, width = 0.2, color = "gray75"
                            )
                            p1 <- p1 + geom_errorbar(
                                data = tempest, aes(x = X, ymin = CI_lower, ymax = CI_upper), linewidth = line.size, width = 0.2
                            )
                        }
                    }
                }

                if (by.group) {
                    p1 <- p1 + geom_point(
                        data = tempest, aes(x = X, y = TE),
                        shape = 21, fill = line.color, colour = "white",
                        size = 3 * line.size, stroke = 0.8
                    )
                }
                # ymin=min(yrange)-maxdiff/5

                if (!is.null(diff.values)) {
                    for (target.value in diff.values) {
                        Xnew <- abs(tempest[, "X"] - target.value)
                        d1 <- min(Xnew)
                        label1 <- which.min(Xnew)
                        Xnew[label1] <- Inf
                        d2 <- min(Xnew)
                        label2 <- which.min(Xnew)
                        if (d1 == 0) {
                            est.mark <- tempest[label1, "TE"]
                            if (isTRUE(CI)) {
                                lb.mark <- tempest[label1, "CI_lower"]
                                ub.mark <- tempest[label1, "CI_upper"]
                            }
                        } else if (d2 == 0) {
                            est.mark <- tempest[label2, "TE"]
                            if (isTRUE(CI)) {
                                lb.mark <- tempest[label2, "CI_lower"]
                                ub.mark <- tempest[label2, "CI_upper"]
                            }
                        } else { ## weighted average
                            est.mark1 <- tempest[label1, "TE"]
                            est.mark2 <- tempest[label2, "TE"]
                            est.mark <- ((est.mark1 * d2 + est.mark2 * d1) / (d1 + d2))
                            if (isTRUE(CI)) {
                                lb.mark1 <- tempest[label1, "CI_lower"]
                                ub.mark1 <- tempest[label1, "CI_upper"]
                                lb.mark2 <- tempest[label2, "CI_lower"]
                                ub.mark2 <- tempest[label2, "CI_upper"]
                                lb.mark <- ((lb.mark1 * d2 + lb.mark2 * d1) / (d1 + d2))
                                ub.mark <- ((ub.mark1 * d2 + ub.mark2 * d1) / (d1 + d2))
                            }
                        }

                        p1 <- p1 + annotate("point", x = target.value, y = est.mark, size = 1, colour = "red")
                        if (isTRUE(CI)) {
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
                if (isTRUE(CI)) {
                    tempest <- .rename_est_ci(tempest, point = "ME")
                    tempest <- as.data.frame(tempest)
                }
                if (isFALSE(CI)) {
                    tempest <- tempest[, c(1, 2)]
                    colnames(tempest) <- c("X", "ME")
                    tempest <- as.data.frame(tempest)
                }
                if (!by.group) {
                    p1 <- p1 + geom_line(data = tempest, aes(X, ME), color = line.color, linewidth = line.size)
                }

                if (isTRUE(CI) && "CI_lower" %in% colnames(tempest)) {
                    if (estimator == "kernel" | estimator == "linear" | estimator == "dml" | estimator == "grf" | estimator == "lasso") {
                        if (!by.group) {
                            p1 <- p1 + geom_ribbon(
                                data = tempest, aes(x = X, ymin = CI_lower, ymax = CI_upper),
                                fill = CI.color, alpha = CI.color.alpha
                            )
                        } else {
                            p1 <- p1 + geom_errorbar(
                                data = tempest, aes(x = X, ymin = CI_lower, ymax = CI_upper), linewidth = line.size, width = 0.2
                            )
                        }
                    }
                    if ("CI_uniform_lower" %in% colnames(tempest) & show.uniform.CI) {
                        if (!by.group) {
                            p1 <- p1 + geom_line(data = tempest, aes(x = X, y = CI_uniform_lower), linetype = "dashed", color = "gray50") + geom_line(data = tempest, aes(x = X, y = CI_uniform_upper), linetype = "dashed", color = "gray50")
                        } else {
                            p1 <- p1 + geom_errorbar(
                                data = tempest, aes(x = X, ymin = CI_uniform_lower, ymax = CI_uniform_upper), linewidth = line.size, width = 0.2, color = "gray75"
                            )
                            p1 <- p1 + geom_errorbar(
                                data = tempest, aes(x = X, ymin = CI_lower, ymax = CI_upper), linewidth = line.size, width = 0.2
                            )
                        }
                    }
                }

                if (by.group) {
                    p1 <- p1 + geom_point(
                        data = tempest, aes(x = X, y = ME),
                        shape = 21, fill = line.color, colour = "white",
                        size = 3 * line.size, stroke = 0.8
                    )
                }
                # ymin=min(yrange)-maxdiff/5

                if (!is.null(diff.values)) {
                    for (target.value in diff.values) {
                        Xnew <- abs(tempest[, "X"] - target.value)
                        d1 <- min(Xnew)
                        label1 <- which.min(Xnew)
                        Xnew[label1] <- Inf
                        d2 <- min(Xnew)
                        label2 <- which.min(Xnew)
                        if (d1 == 0) {
                            est.mark <- tempest[label1, "ME"]
                            if (isTRUE(CI)) {
                                lb.mark <- tempest[label1, "CI_lower"]
                                ub.mark <- tempest[label1, "CI_upper"]
                            }
                        } else if (d2 == 0) {
                            est.mark <- tempest[label2, "ME"]
                            if (isTRUE(CI)) {
                                lb.mark <- tempest[label2, "CI_lower"]
                                ub.mark <- tempest[label2, "CI_upper"]
                            }
                        } else { ## weighted average
                            est.mark1 <- tempest[label1, "ME"]
                            est.mark2 <- tempest[label2, "ME"]
                            est.mark <- ((est.mark1 * d2 + est.mark2 * d1) / (d1 + d2))
                            if (isTRUE(CI)) {
                                lb.mark1 <- tempest[label1, "CI_lower"]
                                ub.mark1 <- tempest[label1, "CI_upper"]
                                lb.mark2 <- tempest[label2, "CI_lower"]
                                ub.mark2 <- tempest[label2, "CI_upper"]
                                lb.mark <- ((lb.mark1 * d2 + lb.mark2 * d1) / (d1 + d2))
                                ub.mark <- ((ub.mark1 * d2 + ub.mark2 * d1) / (d1 + d2))
                            }
                        }

                        p1 <- p1 + annotate("point", x = target.value, y = est.mark, size = 1, colour = "red")
                        if (isTRUE(CI)) {
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
                if (isTRUE(CI)) {
                    colnames(tempest) <- c("X", "TE", "sd", "CI_lower", "CI_upper")
                    colnames(tempest.bin2) <- c("x0", "TE", "sd", "CI_lower", "CI_upper")
                    tempest <- as.data.frame(tempest)
                    tempest.bin2 <- as.data.frame(tempest.bin2)
                }
                if (isFALSE(CI)) {
                    tempest <- tempest[, c(1, 2)]
                    colnames(tempest) <- c("X", "TE")
                    tempest.bin2 <- tempest.bin2[, c(1, 2)]
                    colnames(tempest.bin2) <- c("x0", "TE")
                    tempest <- as.data.frame(tempest)
                    tempest.bin2 <- as.data.frame(tempest.bin2)
                }

                p1 <- p1 + geom_line(data = tempest, aes(X, TE), color = line.color, linewidth = line.size)
                if (isTRUE(CI)) {
                    p1 <- p1 + geom_ribbon(
                        data = tempest, aes(x = X, ymin = CI_lower, ymax = CI_upper),
                        fill = CI.color, alpha = CI.color.alpha
                    )
                }
                ## bin estimates
                p1 <- p1 + geom_point(data = tempest.bin2, aes(x0, TE), size = 4 / treat_sc, shape = 21, fill = "white", colour = "red")
                if (isTRUE(CI)) {
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
                if (bin.labs) {
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
                if (isTRUE(CI)) {
                    colnames(tempest) <- c("X", "ME", "sd", "CI_lower", "CI_upper")
                    colnames(tempest.bin2) <- c("x0", "ME", "sd", "CI_lower", "CI_upper")
                    tempest <- as.data.frame(tempest)
                    tempest.bin2 <- as.data.frame(tempest.bin2)
                }
                if (isFALSE(CI)) {
                    tempest <- tempest[, c(1, 2)]
                    colnames(tempest) <- c("X", "ME")
                    tempest.bin2 <- tempest.bin2[, c(1, 2)]
                    colnames(tempest.bin2) <- c("x0", "ME")
                    tempest <- as.data.frame(tempest)
                    tempest.bin2 <- as.data.frame(tempest.bin2)
                }

                p1 <- p1 + geom_line(data = tempest, aes(X, ME), color = line.color, linewidth = line.size)
                if (isTRUE(CI)) {
                    p1 <- p1 + geom_ribbon(
                        data = tempest, aes(x = X, ymin = CI_lower, ymax = CI_upper),
                        fill = CI.color, alpha = CI.color.alpha
                    )
                }
                ## bin estimates
                p1 <- p1 + geom_point(data = tempest.bin2, aes(x0, ME), size = 4 / treat_sc, shape = 21, fill = "white", colour = "red")
                if (isTRUE(CI)) {
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
                if (bin.labs) {
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
    if (is.null(cex.lab)) {
        cex.lab <- 15
    } else {
        cex.lab <- 15 * cex.lab
    }
    if (is.null(cex.axis)) {
        cex.axis <- 15 / treat_sc
    } else {
        cex.axis <- 15 * cex.axis / treat_sc
    }
    ## title
    if (is.null(cex.main)) {
        cex.main <- 18
    } else {
        cex.main <- 18 * cex.main
    }

    if (is.null(cex.sub)) {
        cex.sub <- 12
    } else {
        cex.sub <- 12 * cex.sub
    }

    ## xlim and ylim
    if (!is.null(ylim)) {
        ylim2 <- c(ylim[1] - (ylim[2] - ylim[1]) * 0.25 / 6, ylim[2] + (ylim[2] - ylim[1]) * 0.4 / 6)
    }

    if (treat.type == "discrete") {
        k <- 1
        for (char in other.treat) {
            p1 <- p.group[[char]]
            ## mark the original interval (in replicated papers)
            if (!is.null(interval)) {
                p1 <- p1 + geom_vline(xintercept = interval, colour = "steelblue", linetype = 2, linewidth = 1.5)
            }

            ## Other universal options
            p1 <- p1 + xlab(NULL) + ylab(NULL) +
                theme(axis.text = element_text(size = cex.axis))

            if (show.subtitles) {
                if (is.null(subtitles)) {
                    subtitle.temp <- paste0("Treated = ", char, ", Baseline = ", base)
                    p1 <- p1 + labs(subtitle = subtitle.temp) + theme(plot.subtitle = element_text(hjust = 0.5, size = cex.sub, lineheight = .8))
                }

                if (!is.null(subtitles)) {
                    subtitle.temp <- subtitles[k]
                    p1 <- p1 + labs(subtitle = subtitle.temp) + theme(plot.subtitle = element_text(hjust = 0.5, size = cex.sub, lineheight = .8))
                }
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
            final_ylim <- if (!is.null(.user_ylim_in)) .user_ylim_in else c(yminmin, ymaxmax)
            final_xlim <- if (!is.null(.user_xlim_in)) .pad_xlim(.user_xlim_in) else NULL
            p.group[[char]] <- p.group[[char]] +
                coord_cartesian(xlim = final_xlim, ylim = final_ylim)
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
            if (!is.null(interval)) {
                p1 <- p1 + geom_vline(xintercept = interval, colour = "steelblue", linetype = 2, linewidth = 1.5)
            }

            ## Other universal options
            p1 <- p1 + xlab(NULL) + ylab(NULL) +
                theme(axis.text = element_text(size = cex.axis))

            if (show.subtitles) {
                if (is.null(subtitles)) {
                    subtitle.temp <- label
                    p1 <- p1 + labs(subtitle = subtitle.temp) + theme(plot.subtitle = element_text(hjust = 0.5, size = cex.sub, lineheight = .8))
                }

                if (!is.null(subtitles)) {
                    subtitle.temp <- subtitles[k]
                    p1 <- p1 + labs(subtitle = subtitle.temp) + theme(plot.subtitle = element_text(hjust = 0.5, size = cex.sub, lineheight = .8))
                }
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
            final_ylim <- if (!is.null(.user_ylim_in)) .user_ylim_in else c(yminmin, ymaxmax)
            final_xlim <- if (!is.null(.user_xlim_in)) .pad_xlim(.user_xlim_in) else NULL
            p.group[[label]] <- p.group[[label]] +
                coord_cartesian(xlim = final_xlim, ylim = final_ylim)
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
    if (!is.null(file)) {
        graph <- as.ggplot(graph)
        ggsave(file, graph, scale = scale, width = width, height = height)
    }

    if (!show.all) {
        graph <- as.ggplot(graph)
        attr(graph, "interflex_xlim") <- .user_xlim_in
        attr(graph, "interflex_ylim") <- .user_ylim_in
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
