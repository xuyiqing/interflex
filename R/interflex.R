## Input Check
interflex <- function(estimator, # "linear", "kernel", "binning" , "gam", "raw", "DML"
                      data,
                      Y, # outcome
                      D, # treatment indicator
                      X, # moderator
                      treat.type = NULL, # discrete or continuous
                      base = NULL, # base group when treatments are discrete
                      Z = NULL, # covariates
                      IV = NULL, # instrument variables
                      FE = NULL, # fixed effects for linear(method) model
                      full.moderate = FALSE, # fully moderated model
                      weights = NULL, # weighting variable
                      na.rm = FALSE,
                      Xunif = FALSE,
                      CI = TRUE,
                      neval = 50,
                      X.eval = NULL,
                      method = "linear", ## "probit"; "logit";  "poisson"; "nbinom"; "linear"
                      vartype = "delta", ## variance type "simu"; "bootstrap"; "delta"
                      vcov.type = "robust", ## "homoscedastic"; "robust"; "cluster"; "pcse"
                      time = NULL,
                      pairwise = TRUE,
                      nboots = 200,
                      nsimu = 1000,
                      parallel = FALSE,
                      cores = 4,
                      cl = NULL, # variable to be clustered on
                      Z.ref = NULL, # same length as Z, set the value of Z when estimating marginal effects/predicted value
                      D.ref = NULL, # default to the mean of D when plotting marginal effects
                      nbins = 3,
                      cutoffs = NULL,
                      wald = TRUE,
                      bw = NULL,
                      kfold = 10,
                      grid = 30,
                      metric = NULL,
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
                      gate = FALSE,
                      grf.num.trees = 4000,
                      figure = TRUE,
                      bin.labs = TRUE,
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
                      show.all = FALSE,
                      legend.title = NULL,
                      diff.values = NULL,
                      percentile = FALSE,
                      gam.k = 10,
                      angle = c(30, 100, -30, -120),
                      span = NULL,
                      scale = 1.1,
                      height = 7,
                      width = 10,
                      box.pos = "down") {
    ##################################################### CHECK ##############################################################
    ## in case data is in tibble format
    if (is.data.frame(data) == FALSE) {
        data <- as.data.frame(data)
    }
    n <- dim(data)[1]

    ## estimator
    if (!estimator %in% c("linear", "binning", "kernel", "gam", "raw", "grf", "dml", "DML")) {
        stop("estimator must be one of the following: raw, linear, binning, kernel, gam, grf or DML.\n")
    }

    ## Y
    if (is.character(Y) == FALSE) {
        stop("\"Y\" is not a string.")
    } else {
        Y <- Y[1]
    }
    if (!Y %in% colnames(data)) {
        f <- sprintf("%s is not in the data. \n", Y)
        stop(f)
    }


    ## D
    if (is.character(D) == FALSE) {
        stop("\"D\" is not a string.")
    } else {
        D <- D[1]
    }
    if (!D %in% colnames(data)) {
        f <- sprintf("%s is not in the data. \n", D)
        stop(f)
    }

    ## X
    if (is.character(X) == FALSE) {
        stop("\"X\" is not a string.")
    } else {
        X <- X[1]
    }
    if (!X %in% colnames(data)) {
        f <- sprintf("%s is not in the data. \n", X)
        stop(f)
    }

    ## Z
    if (is.null(Z) == FALSE) {
        for (i in 1:length(Z)) {
            if (is.character(Z[i]) == FALSE) {
                stop("Some element in \"Z\" is not a string.")
            }
            if (!Z[i] %in% colnames(data)) {
                f <- sprintf("%s is not in the data. \n", Z[i])
                stop(f)
            }
        }
    }

    ## FE
    if (is.null(FE) == FALSE) {
        if (method != "linear") {
            stop("glm models with fixed effects are not supported in this version.")
        }
        requireNamespace("lfe")
        for (i in 1:length(FE)) {
            if (is.character(FE[i]) == FALSE) {
                stop("Some element in \"FE\" is not a string.")
            }
            if (!FE[i] %in% colnames(data)) {
                f <- sprintf("%s is not in the data. \n", FE[i])
                stop(f)
            }
        }
    }

    ## IV
    if (is.null(IV) == FALSE) {
        if (method != "linear") {
            stop("glm models using instrumental variables are not supported in this version.")
        }
        if (vcov.type == "pcse") {
            stop("instrumental regression doesn't stop pcse variances.\n")
        }
        requireNamespace("AER")
        for (i in 1:length(IV)) {
            if (is.character(IV[i]) == FALSE) {
                stop("Some element in \"IV\" is not a string.")
            }
            if (!IV[i] %in% colnames(data)) {
                f <- sprintf("%s is not in the data. \n", IV[i])
                stop(f)
            }
        }
    }

    ## fully moderated model
    if (is.logical(full.moderate) == FALSE & is.numeric(full.moderate) == FALSE) {
        stop("\"CI\" is not a logical flag.")
    }

    ## Weights
    if (is.null(weights) == FALSE) {
        if (is.character(weights) == FALSE) {
            stop("\"weights\" is not a string.")
        }
        if (!weights %in% colnames(data)) {
            f <- sprintf("%s is not in the data. \n", weights)
            stop(f)
        }
    }

    # cl
    if (is.null(cl) == FALSE) {
        if (is.character(cl) == FALSE) {
            stop("\"cl\" is not a string.")
        } else {
            cl <- cl[1]
        }
        if (!cl %in% colnames(data)) {
            f <- sprintf("%s is not in the data. \n", cl)
            stop(f)
        }
    }

    ## na.rm
    if (is.logical(na.rm) == FALSE && is.numeric(na.rm) == FALSE) {
        stop("\"na.rm\" is not a logical flag.")
    }

    ## Xunif
    if (is.logical(Xunif) == FALSE && is.numeric(Xunif) == FALSE) {
        stop("\"Xunif\" is not a logical flag.")
    }

    ## CI
    if (is.logical(CI) == FALSE && is.numeric(CI) == FALSE) {
        stop("\"CI\" is not a logical flag.")
    }

    ## method
    if (is.null(method) == TRUE) {
        method <- "linear"
    }

    if (!method %in% c("linear", "logit", "probit", "poisson", "nbinom")) {
        stop("\"method\" must be one of the following: \"linear\",\"logit\",\"probit\",\"poisson\",\"nbinom\".")
    }

    ## vartype
    if (is.null(vartype) == TRUE) {
        vartype <- "simu"
    }
    if (!vartype %in% c("simu", "delta", "bootstrap")) {
        stop("\"vartype\" must be one of the following: \"simu\",\"delta\".")
    }

    ## vcov.type
    if (is.null(vcov.type) == TRUE) {
        vcov.type <- "robust"
    }
    if (!vcov.type %in% c("robust", "homoscedastic", "cluster", "pcse")) {
        stop("\"vcov.type\" must be one of the following: \"robust\",\"homoscedastic\",\"cluster\",\"pcse\".")
    }

    if (vcov.type == "pcse") {
        if (is.null(FE) == FALSE) {
            vcov.type <- "cluster"
            warning("Fixed-effect models do not allow panel corrected standard errors; changed to clustered standard errors.")
        }
        if (method != "linear") {
            stop("panel corrected standard error can only be used in linear model.")
        }
        if (is.null(cl) == TRUE | is.null(time) == TRUE) {
            stop("\"cl\" and \"time\" should be specified when using panel corrected standard error.")
        }
        if (is.logical(pairwise) == FALSE & is.numeric(pairwise) == FALSE) {
            stop("\"pairwise\" should be a logical flag.")
        }
    }

    if (vcov.type == "cluster" & is.null(cl) == TRUE) {
        stop("\"cl\" should be specified when using clustered standard error.")
    }

    # nboots
    if (is.null(nboots) == FALSE) {
        if (is.numeric(nboots) == FALSE) {
            stop("\"nboots\" is not a positive integer.")
        } else {
            nboots <- nboots[1]
            if (nboots %% 1 != 0 | nboots < 1) {
                stop("\"nboots\" is not a positive number.")
            }
        }
    } else {
        nboots <- 200
    }

    # nsimu
    if (is.null(nsimu) == FALSE) {
        if (is.numeric(nsimu) == FALSE) {
            stop("\"nsimu\" is not a positive integer.")
        } else {
            nsimu <- nsimu[1]
            if (nsimu %% 1 != 0 | nsimu < 1) {
                stop("\"nsimu\" is not a positive number.")
            }
        }
    } else {
        nsimu <- 200
    }

    # Parallel
    if (is.logical(parallel) == FALSE & is.numeric(parallel) == FALSE) {
        stop("\"paralell\" is not a logical flag.")
    }

    # Cores
    if (is.numeric(cores) == FALSE) {
        stop("\"cores\" is not a positive integer.")
    } else {
        cores <- cores[1]
        if (cores %% 1 != 0 | cores <= 0) {
            stop("\"cores\" is not a positive integer.")
        }
    }


    ## check missing values
    vars <- c(Y, D, X, Z, cl, time, weights, FE, IV)
    if (na.rm == TRUE) {
        data <- na.omit(data[, vars])
    } else {
        if (sum(is.na(data[, vars])) > 0) {
            stop("Have Missing values. Try option na.rm = TRUE\n")
        }
    }

    # Z.ref
    if (is.null(Z.ref) == FALSE) {
        if (length(Z.ref) != length(Z)) {
            stop("\"Z.ref\" should have the same length as \"Z\".")
        }
    }

    # D.ref
    if (is.null(D.ref) == FALSE) {
        if (is.numeric(D.ref) == FALSE) {
            stop("\"D.ref\" is not a numeric variable.")
        }
        if (length(D.ref) > 9) {
            stop("Too many values in \"D.ref\".")
        }
    }

    # figure
    if (is.logical(figure) == FALSE & is.numeric(figure) == FALSE) {
        stop("\"figure\" is not a logical flag.")
    }

    # show.subtitles
    if (is.null(show.subtitles) == FALSE) {
        if (is.logical(show.subtitles) == FALSE & is.numeric(show.subtitles) == FALSE) {
            stop("\"show.subtitles\" is not a logical flag.")
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
        Ylabel <- Y
    } else {
        if (is.character(Ylabel) == FALSE) {
            stop("\"Ylabel\" is not a string.")
        } else {
            Ylabel <- Ylabel[1]
        }
    }

    # Dlabel
    if (is.null(Dlabel) == TRUE) {
        Dlabel <- D
    } else {
        if (is.character(Dlabel) == FALSE) {
            stop("\"Dlabel\" is not a string.")
        } else {
            Dlabel <- Dlabel[1]
        }
    }

    # Xlabel
    if (is.null(Xlabel) == TRUE) {
        Xlabel <- X
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

    # interval
    if (is.null(interval) == FALSE) {
        if (is.numeric(interval) == FALSE) {
            stop("Some element in \"interval\" is not numeric.")
        }
    }

    # file
    if (is.null(file) == FALSE) {
        if (is.character(file) == FALSE) {
            stop("Wrong file name.")
        }
    }

    # ncols
    if (is.null(ncols) == FALSE) {
        if (is.numeric(ncols) == FALSE) {
            stop("\"ncols\" is not a positive integer.")
        } else {
            ncols <- ncols[1]
            if (ncols %% 1 != 0 | ncols < 1) {
                stop("\"ncols\" is not a positive number.")
            }
        }
    }

    ## pool
    if (is.logical(pool) == FALSE & is.numeric(pool) == FALSE) {
        stop("\"pool\" is not a logical flag.")
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

    ## legend.title
    if (is.null(legend.title) == FALSE) {
        legend.title <- as.character(legend.title)[1]
    }

    ## diff.values
    if (is.null(diff.values) == FALSE) {
        if (is.numeric(diff.values) == FALSE) {
            stop("\"diff.values\" is not numeric.")
        }
        if (length(diff.values) != 3 & length(diff.values) != 2) {
            stop("\"diff.values\" must be of length 3.")
        }
        min.XX <- min(data[, X])
        max.XX <- max(data[, X])
        for (a in diff.values) {
            if (a < min.XX | a > max.XX) {
                warning("Elements in \"diff.values\" should be within the range of the moderator.")
            }
        }
    }

    ## percentile
    if (is.logical(percentile) == FALSE & is.numeric(percentile) == FALSE) {
        stop("\"percentile\" is not a logical flag.")
    }

    if (percentile == TRUE) {
        for (a in diff.values) {
            if (a < 0 | a > 1) {
                stop("Elements in \"diff.values\" should be between 0 and 1 when percentile==TRUE.")
            }
        }
    }

    ##################################################### Preprocess ##############################################################

    ## method
    unique.y <- sort(unique(data[, Y]))
    y <- data[, Y]
    if (length(unique.y) == 2) {
        if (unique.y[1] != 0 | unique.y[2] != 1) {
            stop("\"Y\" should only contain 0 and 1 in binary model.")
        }
    }

    if (method == "probit" | method == "logit") {
        if (length(unique.y) != 2) {
            stop("\"Y\" should only contain 0 and 1 when using logit/probit methods.")
        }
        if (unique.y[1] != 0 | unique.y[2] != 1) {
            stop("\"Y\" should only contain 0 and 1 when using logit/probit methods.")
        }
    }

    if (method == "poisson" | method == "nbinom") {
        if (FALSE %in% (y == round(y, 0))) {
            stop("\"Y\" should only contain non-negative integers when using poisson/nbinomial methods.")
        }
        if (FALSE %in% (y >= 0)) {
            stop("\"Y\" should only contain non-negative integers when using poisson/nbinomial methods.")
        }
    }

    # metric
    if (is.null(metric) == TRUE) {
        if (length(unique(y)) == 2) {
            metric <- c("Cross Entropy")
        } else {
            metric <- c("MSE")
        }
    } else {
        if (!metric %in% c("MSE", "MAE", "Cross Entropy", "AUC")) {
            stop("metric must be one of the following: MSE, MAE, Cross Entropy and AUC.")
        }
        if (metric %in% c("Cross Entropy", "AUC") & length(unique(y)) > 2) {
            stop("Not binary outcome model, can't use Cross Entropy or AUC criteria.")
        }
    }


    # Xunif
    if (Xunif == TRUE) {
        x <- data[, X]
        data[, X] <- rank(x, ties.method = "average") / length(x) * 100
        Xlabel <- paste(Xlabel, "(Percentile)")
    }

    # Covariates & Z.ref
    to_dummy_var <- c()
    k <- 1
    Z.origin <- Z

    for (a in Z) {
        if (is.factor(data[, a]) == TRUE) {
            to_dummy_var <- c(to_dummy_var, a)
            unique.a <- unique(as.character(data[, a]))
            if (is.null(Z.ref) == FALSE) {
                if (!(Z.ref[k] %in% unique.a)) {
                    stop("A value in \"Z.ref\" doesn't exist in the data...")
                }
            }
        } else if (is.character(data[, a]) == TRUE) {
            stop("\"Z\" should be numeric or factorial.")
        } else {
            if (is.null(Z.ref) == FALSE) {
                suppressWarnings(
                    ww <- as.numeric(Z.ref[k])
                )

                if (is.na(ww) == TRUE) {
                    stop("A value in \"Z.ref\" is not numeric... ")
                }
            }
        }
        k <- k + 1
    }

    if (length(to_dummy_var) > 0) {
        fnames <- paste("factor(", to_dummy_var, ")", sep = "")
        # contr.list <- list(contr.sum, contr.sum,contr.sum,contr.sum)
        contr.list <- rep("contr.sum", length(to_dummy_var))
        names(contr.list) <- fnames
        to_dummy_form <- as.formula(paste("~", paste(fnames, collapse = " + ")))

        if (is.null(Z.ref) == FALSE) {
            first_row <- as.data.frame(data[1, ])
            k <- 1
            for (a in Z) {
                if (a %in% to_dummy_var) {
                    first_row[1, a] <- Z.ref[k]
                } else {
                    first_row[1, a] <- as.numeric(Z.ref[k])
                }
                k <- k + 1
            }
            data <- rbind(first_row, data)
        }

        suppressWarnings(
            to_dummy_mat <- model.matrix(to_dummy_form,
                data = data,
                contrasts.arg = contr.list
            )[, -1]
        )
        to_dummy_mat <- as.matrix(to_dummy_mat)
        dummy_colnames <- c()
        for (i in 1:dim(to_dummy_mat)[2]) {
            dummy_colnames <- c(dummy_colnames, paste0("Dummy.Covariate.", i))
        }
        colnames(to_dummy_mat) <- dummy_colnames
        data <- cbind(data, to_dummy_mat)
        Z <- Z[!Z %in% to_dummy_var]
        Z <- c(Z, dummy_colnames)
        if (is.null(Z.ref) == FALSE) {
            Z.ref <- data[1, Z]
            data <- data[2:dim(data)[1], ]
        }
    }

    if (is.null(Z.ref) == TRUE) {
        Z.ref <- c()
        for (a in Z.origin) {
            if (!(a %in% to_dummy_var)) {
                Z.ref <- c(Z.ref, mean(data[, a]))
            }
        }
        if (length(to_dummy_var) > 0) {
            Z.ref <- c(Z.ref, rep(0, length(dummy_colnames)))
        }
    }

    Z.ref <- as.numeric(Z.ref)
    names(Z.ref) <- Z
    rownames(data) <- NULL

    Z.X <- c()
    if (full.moderate == TRUE & estimator != "kernel") {
        for (sub.Z in Z) {
            Z.X <- c(Z.X, paste0(sub.Z, ".X"))
            data[, paste0(sub.Z, ".X")] <- data[, sub.Z] * data[, X]
        }
    }
    if (full.moderate == TRUE) {
        cat("Use fully moderated model.\n")
    }

    ## treat
    if (is.null(treat.type) == TRUE) {
        if (is.numeric(data[, D]) == TRUE) {
            if (length(unique(data[, D])) > 5) {
                treat.type <- "continuous"
            } else {
                treat.type <- "discrete"
            }
        } else {
            treat.type <- "discrete"
        }
    }

    if (!treat.type %in% c("discrete", "continuous")) {
        stop("\"treat.type\" must be one of the following: \"discrete\",\"continuous\".")
    }

    if (treat.type == "continuous") {
        if (estimator == "grf") {
            stop("\"treat.type\" must be \"discrete\" when \"estimator\" is \"grf\".")
        }
    }

    # if treat is discrete
    if (treat.type == "discrete") {
        D.sample <- NULL

        if (length(unique(data[, D])) > 9) {
            stop("Too many kinds of treatment arms")
        }
        data[, D] <- as.character(data[, D])
        if (is.null(base) == TRUE) {
            base <- sort(unique(data[, D]))[1]
            f <- sprintf("Baseline group not specified; choose treat = %s as the baseline group. \n", base)
            cat(f)
        } else {
            base <- as.character(base)
            if (!base %in% unique(data[, D])) {
                stop("\"base\" must be one kind of treatments.")
            }
            f <- sprintf("Baseline group: treat = %s \n", base)
            cat(f)
        }

        ## in case there are special characters in D
        all.treat <- sort(unique(data[, D]))
        names(all.treat) <- paste("Group", c(1:length(all.treat)), sep = ".")
        other.treat <- all.treat[which(all.treat != base)]
        other.treat <- sort(other.treat)

        if (is.null(order) == FALSE) {
            order <- as.character(order)
            if (length(order) != length(unique(order))) {
                stop("\"order\" should not contain repeated values.")
            }

            if (length(order) != length(other.treat) & estimator != "raw") {
                stop("\"order\" should include all kinds of treatment arms except for the baseline group.")
            }

            if (length(order) != length(all.treat) & estimator == "raw") {
                stop("\"order\" should include all kinds of treatment arms.")
            }

            if (estimator != "raw") {
                if (sum(!is.element(order, other.treat)) != 0 | sum(!is.element(other.treat, order)) != 0) {
                    stop("\"order\" should include all kinds of treatment arms except for the baseline group.")
                }
            }

            if (estimator == "raw") {
                if (sum(!is.element(order, all.treat)) != 0 | sum(!is.element(all.treat, order)) != 0) {
                    stop("\"order\" should include all kinds of treatment arms.")
                }
            }

            other.treat <- order
            colnames.p <- c()
            for (char in other.treat) {
                colnames.p <- c(colnames.p, names(all.treat[which(all.treat == char)]))
            }
            names(other.treat) <- colnames.p
        }

        all.treat.origin <- all.treat
        other.treat.origin <- other.treat
        base.origin <- base
        for (char in names(all.treat)) {
            data[which(data[, D] == all.treat[char]), D] <- char
        }

        all.treat <- names(all.treat.origin)
        other.treat <- names(other.treat.origin)
        base <- names(all.treat.origin[which(all.treat.origin == base.origin)])
        names(all.treat) <- all.treat.origin
        names(other.treat) <- other.treat.origin
        names(base) <- base.origin

        if (is.null(subtitles) == FALSE) {
            if (length(subtitles) != length(other.treat) & estimator != "raw") {
                stop("The number of elements in \"subtitles\" should be m-1(m is the number of different treatment arms).\n")
            }
            if (length(subtitles) != length(all.treat) & estimator == "raw") {
                stop("The number of elements in \"subtitles\" should be the number of different treatment arms.\n")
            }
        }

        if (is.logical(show.subtitles) == FALSE & is.numeric(show.subtitles) == FALSE & is.null(show.subtitles) == FALSE) {
            stop("\"show.subtitles\" is not a logical flag.")
        }

        if (is.null(IV) == FALSE) {
            if (length(IV) < length(other.treat)) {
                stop("Insufficient number of IVs.\n")
            }
        }
    }

    if (treat.type == "continuous") {
        base <- NULL
        other.treat <- NULL
        if (is.numeric(data[, D]) == FALSE) {
            stop("\"D\" is not a numeric variable")
        }
        if (is.null(D.ref) == TRUE) {
            D.sample <- quantile(data[, D], probs = c(0.5), na.rm = T)
            all.treat <- names(D.sample)
            ntreat <- length(D.sample)
            labelname <- c()
            for (targetD in D.sample) {
                labelname <- c(labelname, paste0(D, "=", round(targetD, 3)))
            }
            labelname <- paste0(labelname, " (", all.treat, ")")
            names(D.sample) <- labelname
        } else {
            if (is.numeric(D.ref) == FALSE) {
                stop("\"D.ref\" is not a numeric variable")
            }
            D.sample <- D.ref
            labelname <- c()
            for (targetD in D.sample) {
                labelname <- c(labelname, paste0(D, "=", round(targetD, 3)))
            }
            names(D.sample) <- labelname
            all.treat <- labelname
            ntreat <- length(D.sample)
        }
        if (is.null(IV) == FALSE) {
            if (length(IV) < 1) {
                stop("Insufficient number of IVs.\n")
            }
        }
    }

    # number of columns
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
            ncols <- length(unique(data[, D])) - 1
        }
        if (treat.type == "continuous") {
            ncols <- length(D.sample)
        }
    }

    treat.info <- list()
    treat.info[["treat.type"]] <- treat.type
    treat.info[["other.treat"]] <- other.treat
    treat.info[["all.treat"]] <- all.treat
    treat.info[["base"]] <- base
    treat.info[["D.sample"]] <- D.sample
    treat.info[["ncols"]] <- ncols

    # diff values
    if (is.null(diff.values) == TRUE) {
        diff.values.plot <- NULL
        diff.values <- quantile(data[, X], probs = c(0.25, 0.5, 0.75))
        difference.name <- c("50% vs 25%", "75% vs 50%", "75% vs 25%")
    } else {
        if (estimator == "binning") {
            stop("Can't calculate the difference using binning estimator...")
        }

        if (percentile == TRUE) {
            diff.pc <- diff.values
            diff.values <- quantile(data[, X], probs = diff.values)
        }

        diff.values.plot <- diff.values

        if (length(diff.values) == 3) {
            if (percentile == FALSE) {
                difference.name <- c(
                    paste0(round(diff.values[2], 3), " vs ", round(diff.values[1], 3)),
                    paste0(round(diff.values[3], 3), " vs ", round(diff.values[2], 3)),
                    paste0(round(diff.values[3], 3), " vs ", round(diff.values[1], 3))
                )
            }
            if (percentile == TRUE) {
                difference.name <- c(
                    paste0(round(100 * diff.pc[2], 3), "%", " vs ", round(100 * diff.pc[1], 3), "%"),
                    paste0(round(100 * diff.pc[3], 3), "%", " vs ", round(100 * diff.pc[2], 3), "%"),
                    paste0(round(100 * diff.pc[3], 3), "%", " vs ", round(100 * diff.pc[1], 3), "%")
                )
            }
        }
        if (length(diff.values) == 2) {
            if (percentile == FALSE) {
                difference.name <- c(paste0(round(diff.values[2], 3), " vs ", round(diff.values[1], 3)))
            }
            if (percentile == TRUE) {
                difference.name <- c(paste0(round(100 * diff.pc[2], 3), "%", " vs ", round(100 * diff.pc[1], 3), "%"))
            }
        }
    }

    diff.info <- list()
    diff.info[["diff.values.plot"]] <- diff.values.plot
    diff.info[["diff.values"]] <- diff.values
    diff.info[["difference.name"]] <- difference.name

    # Fixed effects
    ## change fixed effect variable to factor
    if (is.null(FE) == FALSE) {
        if (length(FE) == 1) {
            data[, FE] <- as.numeric(as.factor(data[, FE]))
        } else {
            data[, FE] <- sapply(data[, FE], function(vec) {
                as.numeric(as.factor(vec))
            })
        }
    }
    if (is.null(cl) == FALSE) {
        data[, cl] <- as.numeric(as.factor(data[, cl]))
    }
    if (is.null(time) == FALSE) {
        data[, time] <- as.numeric(as.factor(data[, time]))
    }

    if (estimator == "linear") {
        output <- interflex.linear(
            data = data,
            Y = Y, # outcome
            D = D, # treatment indicator
            X = X, # moderator
            Z = Z, # covariates
            full.moderate = full.moderate, # fully moderated model
            Z.X = Z.X, # fully moderated terms
            FE = FE,
            IV = IV,
            weights = weights, # weighting variable
            # CI = CI,
            neval = neval,
            X.eval = X.eval,
            method = method,
            vartype = vartype,
            vcov.type = vcov.type,
            time = time,
            pairwise = pairwise,
            nboots = nboots,
            nsimu = nsimu,
            parallel = parallel,
            cores = cores,
            cl = cl, # variable to be clustered on
            # predict = predict,
            Z.ref = Z.ref, # same length as Z, set the value of Z when estimating marginal effects/predicted value
            treat.info = treat.info,
            diff.info = diff.info,
            figure = figure,
            CI = CI,
            subtitles = subtitles,
            show.subtitles = show.subtitles,
            Xdistr = Xdistr, # c("density","histogram","none")
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
            interval = interval,
            file = file,
            ncols = ncols,
            pool = pool,
            color = color,
            legend.title = legend.title,
            show.all = show.all,
            scale = scale,
            height = height,
            width = width
        )
    }
    if (estimator == "binning") {
        output <- interflex.binning(
            data = data,
            Y = Y, # outcome
            D = D, # treatment indicator
            X = X, # moderator
            treat.info = treat.info,
            Z = Z, # covariates
            FE = FE,
            IV = IV,
            full.moderate = full.moderate, # fully moderated model
            Z.X = Z.X, # fully moderated terms
            weights = weights, # weighting variable
            neval = neval,
            method = method, ## "probit"; "logit"; "poisson"; "nbinom"
            nbins = nbins,
            cutoffs = cutoffs,
            vartype = vartype, ## variance type "simu"; "bootstrap"; "delta"
            vcov.type = vcov.type, ## "homoscedastic"; "robust"; "cluster"; "pcse"
            time = time,
            pairwise = pairwise,
            nboots = nboots,
            nsimu = nsimu,
            parallel = parallel,
            cores = cores,
            cl = cl, # variable to be clustered on
            # predict = predict,
            Z.ref = Z.ref, # same length as Z, set the value of Z when estimating marginal effects/predicted value
            wald = wald,
            figure = figure,
            CI = CI,
            bin.labs = bin.labs,
            order = order,
            subtitles = subtitles,
            show.subtitles = show.subtitles,
            Xdistr = Xdistr, # c("density","histogram","none")
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
            interval = interval,
            file = file,
            ncols = ncols,
            pool = pool,
            color = color,
            legend.title = legend.title,
            show.all = show.all,
            scale = scale,
            height = height,
            width = width
        )
    }
    if (estimator == "kernel") {
        output <- interflex.kernel(
            data = data,
            Y = Y, # outcome
            D = D, # treatment indicator
            X = X, # moderator
            treat.info = treat.info,
            diff.info = diff.info,
            bw = bw,
            grid = grid,
            kfold = kfold,
            metric = metric,
            Z = Z, # covariates
            FE = FE, # fixed effects
            IV = IV, # instrumental variables
            full.moderate = full.moderate, # fully moderated model
            # Z.X = Z.X, # fully moderated terms
            weights = weights, # weighting variable
            neval = neval,
            X.eval = X.eval,
            method = method, ## "probit"; "logit"; "poisson"; "nbinom"
            CI = CI,
            vartype = vartype,
            nboots = nboots,
            parallel = parallel,
            cores = cores,
            cl = cl, # variable to be clustered on
            # predict = predict,
            Z.ref = Z.ref, # same length as Z, set the value of Z when estimating marginal effects/predicted value
            figure = figure,
            order = order,
            subtitles = subtitles,
            show.subtitles = show.subtitles,
            Xdistr = Xdistr, # c("density","histogram","none")
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
            interval = interval,
            file = file,
            ncols = ncols,
            pool = pool,
            color = color,
            legend.title = legend.title,
            show.all = show.all,
            scale = scale,
            height = height,
            width = width
        )
    }

    if (estimator == "gam") {
        if (treat.type == "discrete") {
            stop("gam can only be used when the treatment is continuous.\n")
        }
        output <- interflex.gam(
            data = data,
            Y = Y,
            D = D,
            X = X,
            Z = Z,
            Z.ref = Z.ref,
            weights = weights,
            full.moderate = full.moderate,
            FE = FE,
            CI = CI,
            method = method,
            k = gam.k,
            angle = angle,
            Ylabel = Ylabel,
            Dlabel = Dlabel,
            Xlabel = Xlabel,
            n.grid = neval,
            file = file,
            scale = scale,
            height = height,
            width = width
        )
    }

    if (estimator == "raw") {
        output <- interflex.raw(
            data = data,
            Y = Y,
            D = D,
            X = X,
            treat.info = treat.info,
            order = order,
            subtitles = subtitles,
            weights = weights,
            nbins = nbins,
            cutoffs = cutoffs,
            span = span,
            main = main,
            Ylabel = Ylabel,
            Dlabel = Dlabel,
            Xlabel = Xlabel,
            theme.bw = theme.bw,
            show.grid = show.grid,
            cex.main = cex.main,
            cex.lab = cex.lab,
            cex.axis = cex.axis,
            ncols = ncols,
            file = file,
            scale = scale,
            height = height,
            width = width,
            box.pos = box.pos,
            xlim = xlim,
            ylim = ylim
        )
    }
    if (estimator == "grf") {
        output <- interflex.grf(
            data = data,
            Y = Y, # outcome
            D = D, # treatment indicator
            X = X, # moderator
            Z = Z, # covariates
            weights = weights, # weighting variable
            num.trees = grf.num.trees,
            treat.info = treat.info,
            diff.info = diff.info,
            figure = figure,
            CI = CI,
            subtitles = subtitles,
            show.subtitles = show.subtitles,
            Xdistr = Xdistr, # c("density","histogram","none")
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
            interval = interval,
            file = file,
            ncols = ncols,
            pool = pool,
            color = color,
            legend.title = legend.title,
            show.all = show.all,
            scale = scale,
            height = height,
            width = width
        )
    }
    if ((estimator == "DML") | (estimator == "dml")) {
        output <- interflex.DML(
            data = data,
            Y = Y, # outcome
            D = D, # treatment indicator
            X = X, # moderator
            Z = Z, # covariates
            weights = weights, # weighting variable
            model.y = model.y,
            param.y = param.y,
            param.grid.y = param.grid.y,
            scoring.y = scoring.y,
            model.t = model.t,
            param.t = param.t,
            param.grid.t = param.grid.t,
            scoring.t = scoring.t,
            CV = CV,
            n.folds = n.folds,
            n.jobs = n.jobs,
            cf.n.folds = cf.n.folds,
            gate = gate,
            treat.info = treat.info,
            diff.info = diff.info,
            figure = figure,
            CI = CI,
            subtitles = subtitles,
            show.subtitles = show.subtitles,
            Xdistr = Xdistr, # c("density","histogram","none")
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
            interval = interval,
            file = file,
            ncols = ncols,
            pool = pool,
            color = color,
            legend.title = legend.title,
            show.all = show.all,
            scale = scale,
            height = height,
            width = width
        )
    }

    return(output)
}
