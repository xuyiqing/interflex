#' @rdname ttest.interflex
#' @export
inter.test <- function(out,
                       diff.values,
                       percentile = FALSE,
                       k = 16) {
    if (!class(out) %in% c("interflex")) {
        stop("Not an \"interflex\" object.")
    }

    requireNamespace("mgcv")

    if (is.numeric(k) == FALSE) {
        stop("\"k\" is not numeric.")
    }

    if (is.numeric(diff.values) == FALSE) {
        stop("\"diff.values\" is not numeric.")
    }
    if (length(diff.values) != 3 & length(diff.values) != 2) {
        stop("\"diff.values\" must be of length 2 or 3.")
    }

    if (is.logical(percentile) == FALSE & is.numeric(percentile) == FALSE) {
        stop("\"percentile\" is not a logical flag.")
    }

    if (k %% 1 != 0) {
        stop("\"k\" is not a positive integer.")
    } else {
        k <- k[1]
    }
    if (k < 1) {
        stop("\"k\" should be a positive integer larger than 1.")
    }

    if (percentile == TRUE) {
        for (a in diff.values) {
            if (a < 0 | a > 1) {
                stop("Elements in \"diff.values\" should be between 0 and 1 when percentile==TRUE.")
            }
        }
    }

    treat.info <- out$treat.info
    treat.type <- treat.info[["treat.type"]]
    if (treat.type == "discrete") {
        other.treat <- treat.info[["other.treat"]]
        other.treat.origin <- names(other.treat)
        other.treat <- other.treat.origin
        all.treat <- treat.info[["all.treat"]]
        all.treat.origin <- names(all.treat)
        all.treat <- all.treat.origin
        base <- treat.info[["base"]]
    }
    if (treat.type == "continuous") {
        D.sample <- treat.info[["D.sample"]]
        label.name <- names(D.sample)
        # names(label.name) <- D.sample
    }

    estimator <- out$estimator

    if (estimator == "kernel") {
        if (out$CI == FALSE | length(out$vcov.matrix) == 0) {
            stop("ttest() can't work without vcov matrix. You may try set vartype as bootstrap for kernel method.")
        }
    }

    if (estimator == "binning" | estimator == "dml") {
        stop("ttest() can only work after linear or kernel estimations.")
    }


    if (treat.type == "discrete" & estimator == "linear") {
        data.x <- out$est.lin[[other.treat[1]]][, "X"]
    }
    if (treat.type == "discrete" & estimator == "kernel") {
        data.x <- out$est.kernel[[other.treat[1]]][, "X"]
    }
    if (treat.type == "continuous" & estimator == "linear") {
        data.x <- out$est.lin[[label.name[1]]][, "X"]
    }
    if (treat.type == "continuous" & estimator == "kernel") {
        data.x <- out$est.kernel[[label.name[1]]][, "X"]
    }
    min.XX <- min(data.x)
    max.XX <- max(data.x)

    if (percentile == TRUE) {
        diff.pc <- diff.values
        diff.values <- quantile(data.x, probs = diff.values)
    }

    for (a in diff.values) {
        if (a < min.XX | a > max.XX) {
            stop("Elements in \"diff.values\" should be greater than the minimum and less than the maximum of the moderator.")
        }
    }

    if (estimator == "linear") {
        est <- out$est.lin
    }
    if (estimator == "kernel") {
        est <- out$est.kernel
    }
    vcov.matrix <- out$vcov.matrix

    if (length(diff.values) == 3) {
        if (percentile == FALSE) {
            diff.name <- c(
                paste0(round(diff.values[2], 3), " vs ", round(diff.values[1], 3)),
                paste0(round(diff.values[3], 3), " vs ", round(diff.values[2], 3)),
                paste0(round(diff.values[3], 3), " vs ", round(diff.values[1], 3))
            )
        }
        if (percentile == TRUE) {
            diff.name <- c(
                paste0(round(100 * diff.pc[2], 3), "%", " vs ", round(100 * diff.pc[1], 3), "%"),
                paste0(round(100 * diff.pc[3], 3), "%", " vs ", round(100 * diff.pc[2], 3), "%"),
                paste0(round(100 * diff.pc[3], 3), "%", " vs ", round(100 * diff.pc[1], 3), "%")
            )
        }
    }
    if (length(diff.values) == 2) {
        if (percentile == FALSE) {
            diff.name <- c(paste0(round(diff.values[2], 3), " vs ", round(diff.values[1], 3)))
        }
        if (percentile == TRUE) {
            diff.name <- c(paste0(round(100 * diff.pc[2], 3), "%", " vs ", round(100 * diff.pc[1], 3), "%"))
        }
    }

    get_stat <- function(point1, point2) {
        var1 <- predict(model_use, newdata = cbind.data.frame(x1 = point1, x2 = point1))
        var2 <- predict(model_use, newdata = cbind.data.frame(x1 = point2, x2 = point2))
        cov12 <- predict(model_use, newdata = cbind.data.frame(x1 = point1, x2 = point2))
        marg1 <- predict(model_use2, newdata = cbind.data.frame(X = point1))
        marg2 <- predict(model_use2, newdata = cbind.data.frame(X = point2))
        diff <- marg2 - marg1
        var.diff <- var2 + var1 - 2 * cov12
        se.diff <- sqrt(var.diff)
        z.value <- diff / se.diff
        pvalue2sided <- 2 * pnorm(-abs(z.value))
        lbound <- diff - 1.96 * se.diff
        ubound <- diff + 1.96 * se.diff
        diff.table <- c(diff, se.diff, z.value, pvalue2sided, lbound, ubound)
        diff.table <- sprintf("%.3f", diff.table)
        return(diff.table)
    }


    if (treat.type == "discrete") {
        x1 <- rep(data.x, length(data.x))
        x2 <- rep(data.x, each = length(data.x))
        diff.table.list <- list()

        for (char in other.treat) {
            cov_base <- c(vcov.matrix[[char]])

            data_touse <- cbind.data.frame(x1 = x1, x2 = x2, y = cov_base)
            data_touse2 <- cbind.data.frame(X = est[[char]][, 1], ME = est[[char]][, 2])

            model_use <- gam(y ~ te(x1, x2, k = k), data = data_touse)
            model_use2 <- gam(ME ~ s(X, k = k), data = data_touse2)

            diff.table <- matrix(0, nrow = 0, ncol = 6)
            colnames(diff.table) <- c("diff.estimate", "sd", "z-value", "p-value", "lower CI(95%)", "upper CI(95%)")

            if (length(diff.values) == 3) {
                diff.table <- rbind(diff.table, get_stat(diff.values[1], diff.values[2]))
                diff.table <- rbind(diff.table, get_stat(diff.values[2], diff.values[3]))
                diff.table <- rbind(diff.table, get_stat(diff.values[1], diff.values[3]))
                rownames(diff.table) <- diff.name
            }

            if (length(diff.values) == 2) {
                diff.table <- rbind(diff.table, get_stat(diff.values[1], diff.values[2]))
                rownames(diff.table) <- diff.name
            }
            diff.table <- as.data.frame(diff.table)
            diff.table.list[[char]] <- diff.table
        }
        diff.table <- diff.table.list
    }

    if (treat.type == "continuous") {
        x1 <- rep(data.x, length(data.x))
        x2 <- rep(data.x, each = length(data.x))
        diff.table.list <- list()

        for (label in label.name) {
            cov_base <- c(vcov.matrix[[label]])

            data_touse <- cbind.data.frame(x1 = x1, x2 = x2, y = cov_base)
            data_touse2 <- cbind.data.frame(X = est[[label]][, 1], ME = est[[label]][, 2])

            model_use <- gam(y ~ te(x1, x2, k = k), data = data_touse)
            model_use2 <- gam(ME ~ s(X, k = k), data = data_touse2)

            diff.table <- matrix(0, nrow = 0, ncol = 6)
            colnames(diff.table) <- c("diff.estimate", "sd", "z-value", "p-value", "lower CI(95%)", "upper CI(95%)")

            if (length(diff.values) == 3) {
                diff.table <- rbind(diff.table, get_stat(diff.values[1], diff.values[2]))
                diff.table <- rbind(diff.table, get_stat(diff.values[2], diff.values[3]))
                diff.table <- rbind(diff.table, get_stat(diff.values[1], diff.values[3]))
                rownames(diff.table) <- diff.name
            }

            if (length(diff.values) == 2) {
                diff.table <- rbind(diff.table, get_stat(diff.values[1], diff.values[2]))
                rownames(diff.table) <- diff.name
            }
            diff.table <- as.data.frame(diff.table)
            diff.table.list[[label]] <- diff.table
        }
        diff.table <- diff.table.list
    }


    return(diff.table)
}
