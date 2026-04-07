## Input Check
interflex <- function(estimator, # "linear", "kernel", "binning" , "gam", "raw", "dml", "lasso"
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
                      parallel = TRUE,
                      cores = NULL,
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
                      signal    = c("aipw"),
                      estimand  = c("ATE"),
                      basis.type          = c("bspline"),
                      include.interactions = TRUE,
                      poly.degree         = 2,
                      spline.df           = 4,
                      spline.degree       = 2,
                      lambda.seq         = NULL,
                      reduce.dimension   = c("kernel"),
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
                      theme.bw = TRUE,
                      show.grid = TRUE,
                      show.uniform.CI = TRUE,
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
                      box.pos = "down",
                      verbose = TRUE) {
    ##################################################### CHECK ##############################################################
    ## in case data is in tibble format
    if (!is.data.frame(data)) {
        data <- as.data.frame(data)
    }
    n <- dim(data)[1]

    # Reset per-call warning flags
    options(interflex.uniform_ci_warned = FALSE)

    estimator <- tolower(estimator)

    ## estimator
    if (!estimator %in% c("linear", "binning", "kernel", "gam", "raw", "grf", "dml", "lasso")) {
        stop("estimator must be one of the following: raw, linear, binning, kernel, gam, grf, dml or lasso.\n")
    }

    ## gate validation
    if (isTRUE(gate)) {
        if (estimator == "binning") {
            stop("gate = TRUE is not supported with estimator = 'binning'.\n")
        }
        n_unique_X <- length(unique(data[, X]))
        if (n_unique_X > 10) {
            stop("gate = TRUE requires a discrete moderator X. The variable '", X, "' has ", n_unique_X, " unique values, which suggests it is continuous.\n")
        }
        ## For kernel + gate: evaluate only at observed X levels
        if (estimator == "kernel") {
            X.eval <- sort(unique(data[[X]]))
            neval <- length(X.eval)
        }
    }

    ## Y
    if (!is.character(Y)) {
        stop("\"Y\" is not a string.")
    } else {
        Y <- Y[1]
    }
    if (!Y %in% colnames(data)) {
        f <- sprintf("%s is not in the data. \n", Y)
        stop(f)
    }


    ## D
    if (!is.character(D)) {
        stop("\"D\" is not a string.")
    } else {
        D <- D[1]
    }
    if (!D %in% colnames(data)) {
        f <- sprintf("%s is not in the data. \n", D)
        stop(f)
    }

    ## X
    if (!is.character(X)) {
        stop("\"X\" is not a string.")
    } else {
        X <- X[1]
    }
    if (!X %in% colnames(data)) {
        f <- sprintf("%s is not in the data. \n", X)
        stop(f)
    }

    ## Z
    if (!is.null(Z)) {
        for (i in 1:length(Z)) {
            if (!is.character(Z[i])) {
                stop("Some element in \"Z\" is not a string.")
            }
            if (!Z[i] %in% colnames(data)) {
                f <- sprintf("%s is not in the data. \n", Z[i])
                stop(f)
            }
        }
    }

    ## FE
    if (!is.null(FE)) {
        if (method != "linear") {
            stop("glm models with fixed effects are not supported in this version.")
        }
        requireNamespace("lfe")
        for (i in 1:length(FE)) {
            if (!is.character(FE[i])) {
                stop("Some element in \"FE\" is not a string.")
            }
            if (!FE[i] %in% colnames(data)) {
                f <- sprintf("%s is not in the data. \n", FE[i])
                stop(f)
            }
        }
    }

    ## IV
    if (!is.null(IV)) {
        if (method != "linear") {
            stop("glm models using instrumental variables are not supported in this version.")
        }
        if (vcov.type == "pcse") {
            stop("instrumental regression doesn't support pcse variances.\n")
        }
        requireNamespace("AER")
        for (i in 1:length(IV)) {
            if (!is.character(IV[i])) {
                stop("Some element in \"IV\" is not a string.")
            }
            if (!IV[i] %in% colnames(data)) {
                f <- sprintf("%s is not in the data. \n", IV[i])
                stop(f)
            }
        }
    }

    ## fully moderated model
    if (!is.logical(full.moderate) & !is.numeric(full.moderate)) {
        stop("\"full.moderate\" is not a logical flag.")
    }

    ## Weights
    if (!is.null(weights)) {
        if (!is.character(weights)) {
            stop("\"weights\" is not a string.")
        }
        if (!weights %in% colnames(data)) {
            f <- sprintf("%s is not in the data. \n", weights)
            stop(f)
        }
    }

    # cl
    if (!is.null(cl)) {
        if (!is.character(cl)) {
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
    if (!is.logical(na.rm) && !is.numeric(na.rm)) {
        stop("\"na.rm\" is not a logical flag.")
    }

    ## Xunif
    if (!is.logical(Xunif) && !is.numeric(Xunif)) {
        stop("\"Xunif\" is not a logical flag.")
    }

    ## CI
    if (!is.logical(CI) && !is.numeric(CI)) {
        stop("\"CI\" is not a logical flag.")
    }

    ## method
    if (is.null(method)) {
        method <- "linear"
    }

    if (!method %in% c("linear", "logit", "probit", "poisson", "nbinom")) {
        stop("\"method\" must be one of the following: \"linear\",\"logit\",\"probit\",\"poisson\",\"nbinom\".")
    }

    ## vartype
    if (is.null(vartype)) {
        vartype <- "simu"
    }
    if (!vartype %in% c("simu", "delta", "bootstrap")) {
        stop("\"vartype\" must be one of the following: \"simu\",\"delta\".")
    }

    ## vcov.type
    if (is.null(vcov.type)) {
        vcov.type <- "robust"
    }
    if (!vcov.type %in% c("robust", "homoscedastic", "cluster", "pcse")) {
        stop("\"vcov.type\" must be one of the following: \"robust\",\"homoscedastic\",\"cluster\",\"pcse\".")
    }

    if (vcov.type == "pcse") {
        if (!is.null(FE)) {
            vcov.type <- "cluster"
            warning("Fixed-effect models do not allow panel corrected standard errors; changed to clustered standard errors.")
        }
        if (method != "linear") {
            stop("panel corrected standard error can only be used in linear model.")
        }
        if (is.null(cl) | is.null(time)) {
            stop("\"cl\" and \"time\" should be specified when using panel corrected standard error.")
        }
        if (!is.logical(pairwise) & !is.numeric(pairwise)) {
            stop("\"pairwise\" should be a logical flag.")
        }
    }

    if (vcov.type == "cluster" & is.null(cl)) {
        stop("\"cl\" should be specified when using clustered standard error.")
    }

    # nboots
    if (!is.null(nboots)) {
        if (!is.numeric(nboots)) {
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
    if (!is.null(nsimu)) {
        if (!is.numeric(nsimu)) {
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
    if (!is.logical(parallel) & !is.numeric(parallel)) {
        stop("\"paralell\" is not a logical flag.")
    }

    # Cores: if not supplied, use min(available - 2, 8) to avoid hogging system resources
    cores_auto <- is.null(cores)
    if (cores_auto) {
        cores <- max(1L, min(parallelly::availableCores(omit = 2L), 8L))
    } else {
        if (!is.numeric(cores) || cores[1] %% 1 != 0 || cores[1] <= 0) {
            stop("\"cores\" is not a positive integer.")
        }
        cores <- cores[1]
    }
    # Only print parallel message when parallel computing will actually be used
    uses_parallel <- estimator %in% c("lasso", "dml", "grf") ||
        (estimator %in% c("kernel", "linear", "binning") && vartype == "bootstrap" && parallel)
    if (verbose && uses_parallel) {
        avail <- parallelly::availableCores()
        cat(sprintf("Parallel computing: using %d of %d available cores.\n", cores, avail))
        if (cores_auto) {
            cat("To change: set cores = <n> in interflex().\n")
        }
    }


    ## check missing values
    vars <- c(Y, D, X, Z, cl, time, weights, FE, IV)
    if (na.rm) {
        data <- na.omit(data[, vars])
    } else {
        if (sum(is.na(data[, vars])) > 0) {
            stop("Have Missing values. Try option na.rm = TRUE\n")
        }
    }

    # Z.ref
    if (!is.null(Z.ref)) {
        if (length(Z.ref) != length(Z)) {
            stop("\"Z.ref\" should have the same length as \"Z\".")
        }
    }

    # D.ref
    if (!is.null(D.ref)) {
        if (!is.numeric(D.ref)) {
            stop("\"D.ref\" is not a numeric variable.")
        }
        if (length(D.ref) > 9) {
            stop("Too many values in \"D.ref\".")
        }
    }

    # figure
    if (!is.logical(figure) & !is.numeric(figure)) {
        stop("\"figure\" is not a logical flag.")
    }

    # show.subtitles
    if (!is.null(show.subtitles)) {
        if (!is.logical(show.subtitles) & !is.numeric(show.subtitles)) {
            stop("\"show.subtitles\" is not a logical flag.")
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
        Ylabel <- Y
    } else {
        if (!is.character(Ylabel)) {
            stop("\"Ylabel\" is not a string.")
        } else {
            Ylabel <- Ylabel[1]
        }
    }

    # Dlabel
    if (is.null(Dlabel)) {
        Dlabel <- D
    } else {
        if (!is.character(Dlabel)) {
            stop("\"Dlabel\" is not a string.")
        } else {
            Dlabel <- Dlabel[1]
        }
    }

    # Xlabel
    if (is.null(Xlabel)) {
        Xlabel <- X
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

    ## xlim ylim
    ## Auto-trim: when xlim is not specified and X is continuous,
    ## clip tails that lack treatment variation (sparse/uninformative regions).
    ## Uses 1st and 99th percentiles as candidates, but only clips a tail
    ## if the observations beyond the cutoff lack treatment variation
    ## (e.g., all treated or all control) or have fewer than 10 observations.
    if (is.null(xlim) && !isTRUE(gate)) {
        x_vals <- data[[X]]
        d_vals <- data[[D]]
        n_unique_x <- length(unique(x_vals))
        if (n_unique_x > 10) {
            q01 <- as.numeric(quantile(x_vals, 0.01, na.rm = TRUE))
            q99 <- as.numeric(quantile(x_vals, 0.99, na.rm = TRUE))
            .tail_has_variation <- function(d_sub) {
                if (length(d_sub) < 10) return(FALSE)
                if (is.numeric(d_sub) && length(unique(d_sub)) <= 5) {
                    # binary/discrete D: need at least 2 distinct values
                    return(length(unique(d_sub)) >= 2)
                }
                # continuous D: need nonzero variance
                return(sd(d_sub, na.rm = TRUE) > 1e-10)
            }
            left_idx <- which(x_vals < q01)
            right_idx <- which(x_vals > q99)
            xlim_lo <- if (length(left_idx) > 0 && !.tail_has_variation(d_vals[left_idx])) q01 else min(x_vals, na.rm = TRUE)
            xlim_hi <- if (length(right_idx) > 0 && !.tail_has_variation(d_vals[right_idx])) q99 else max(x_vals, na.rm = TRUE)
            if (xlim_lo > min(x_vals, na.rm = TRUE) || xlim_hi < max(x_vals, na.rm = TRUE)) {
                xlim <- c(xlim_lo, xlim_hi)
            }
        }
    }
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

    # interval
    if (!is.null(interval)) {
        if (!is.numeric(interval)) {
            stop("Some element in \"interval\" is not numeric.")
        }
    }

    # file
    if (!is.null(file)) {
        if (!is.character(file)) {
            stop("Wrong file name.")
        }
    }

    # ncols
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

    ## pool
    if (!is.logical(pool) & !is.numeric(pool)) {
        stop("\"pool\" is not a logical flag.")
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

    ## legend.title
    if (!is.null(legend.title)) {
        legend.title <- as.character(legend.title)[1]
    }

    ## diff.values
    if (!is.null(diff.values)) {
        if (!is.numeric(diff.values)) {
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
    if (!is.logical(percentile) & !is.numeric(percentile)) {
        stop("\"percentile\" is not a logical flag.")
    }

    if (percentile) {
        for (a in diff.values) {
            if (a < 0 | a > 1) {
                stop("Elements in \"diff.values\" should be between 0 and 1 when percentile.")
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
    if (is.null(metric)) {
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
    if (Xunif) {
        x <- data[, X]
        data[, X] <- rank(x, ties.method = "average") / length(x) * 100
        Xlabel <- paste(Xlabel, "(Percentile)")
    }

    # Covariates & Z.ref
    to_dummy_var <- c()
    k <- 1
    Z.origin <- Z

    for (a in Z) {
        if (is.factor(data[, a])) {
            to_dummy_var <- c(to_dummy_var, a)
            unique.a <- unique(as.character(data[, a]))
            if (!is.null(Z.ref)) {
                if (!(Z.ref[k] %in% unique.a)) {
                    stop("A value in \"Z.ref\" doesn't exist in the data...")
                }
            }
        } else if (is.character(data[, a])) {
            stop("\"Z\" should be numeric or factorial.")
        } else {
            if (!is.null(Z.ref)) {
                suppressWarnings(
                    ww <- as.numeric(Z.ref[k])
                )

                if (is.na(ww)) {
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

        if (!is.null(Z.ref)) {
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
        if (!is.null(Z.ref)) {
            Z.ref <- data[1, Z]
            data <- data[2:dim(data)[1], ]
        }
    }

    if (is.null(Z.ref)) {
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
    if (full.moderate & estimator != "kernel") {
        for (sub.Z in Z) {
            Z.X <- c(Z.X, paste0(sub.Z, ".X"))
            data[, paste0(sub.Z, ".X")] <- data[, sub.Z] * data[, X]
        }
    }
    if (full.moderate) {
        cat("Use fully moderated model.\n")
    }

    ## treat
    if (is.null(treat.type)) {
        if (is.numeric(data[, D])) {
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
        if (is.null(base)) {
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

        if (!is.null(order)) {
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

        if (!is.null(subtitles)) {
            if (length(subtitles) != length(other.treat) & estimator != "raw") {
                stop("The number of elements in \"subtitles\" should be m-1(m is the number of different treatment arms).\n")
            }
            if (length(subtitles) != length(all.treat) & estimator == "raw") {
                stop("The number of elements in \"subtitles\" should be the number of different treatment arms.\n")
            }
        }

        if (!is.logical(show.subtitles) & !is.numeric(show.subtitles) & !is.null(show.subtitles)) {
            stop("\"show.subtitles\" is not a logical flag.")
        }

        if (!is.null(IV)) {
            if (length(IV) < length(other.treat)) {
                stop("Insufficient number of IVs.\n")
            }
        }
    }

    if (treat.type == "continuous") {
        base <- NULL
        other.treat <- NULL
        if (!is.numeric(data[, D])) {
            stop("\"D\" is not a numeric variable")
        }
        if (is.null(D.ref)) {
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
            if (!is.numeric(D.ref)) {
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
        if (!is.null(IV)) {
            if (length(IV) < 1) {
                stop("Insufficient number of IVs.\n")
            }
        }
    }

    # number of columns
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
    if (is.null(diff.values)) {
        diff.values.plot <- NULL
        if(length(unique(data[,X]))<5){
            diff.values <- NULL
            difference.name <- NULL          
        }else{
            diff.values <- quantile(data[, X], probs = c(0.25, 0.5, 0.75))
            difference.name <- c("50% vs 25%", "75% vs 50%", "75% vs 25%")                
        }
    } else {
        if (estimator == "binning") {
            stop("Can't calculate the difference using binning estimator...")
        }

        if (percentile) {
            diff.pc <- diff.values
            diff.values <- quantile(data[, X], probs = diff.values)
        }

        diff.values.plot <- diff.values

        if (length(diff.values) == 3) {
            if (!percentile) {
                difference.name <- c(
                    paste0(round(diff.values[2], 3), " vs ", round(diff.values[1], 3)),
                    paste0(round(diff.values[3], 3), " vs ", round(diff.values[2], 3)),
                    paste0(round(diff.values[3], 3), " vs ", round(diff.values[1], 3))
                )
            }
            if (percentile) {
                difference.name <- c(
                    paste0(round(100 * diff.pc[2], 3), "%", " vs ", round(100 * diff.pc[1], 3), "%"),
                    paste0(round(100 * diff.pc[3], 3), "%", " vs ", round(100 * diff.pc[2], 3), "%"),
                    paste0(round(100 * diff.pc[3], 3), "%", " vs ", round(100 * diff.pc[1], 3), "%")
                )
            }
        }
        if (length(diff.values) == 2) {
            if (!percentile) {
                difference.name <- c(paste0(round(diff.values[2], 3), " vs ", round(diff.values[1], 3)))
            }
            if (percentile) {
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
    if (!is.null(FE)) {
        if (length(FE) == 1) {
            data[, FE] <- as.numeric(as.factor(data[, FE]))
        } else {
            data[, FE] <- sapply(data[, FE], function(vec) {
                as.numeric(as.factor(vec))
            })
        }
    }
    if (!is.null(cl)) {
        data[, cl] <- as.numeric(as.factor(data[, cl]))
    }
    if (!is.null(time)) {
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
            show.uniform.CI = show.uniform.CI,
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
            width = width,
            gate = gate
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
            show.uniform.CI = show.uniform.CI,
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
            show.uniform.CI = show.uniform.CI,
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
            width = width,
            verbose = verbose
        )
        ## When gate = TRUE, copy kernel estimates to g.est and re-plot with by.group
        if (isTRUE(gate) && !is.null(output$est.kernel)) {
            output$g.est <- output$est.kernel
            if (figure) {
                class(output) <- "interflex"
                output$figure <- plot.interflex(output,
                    by.group = TRUE, CI = CI,
                    Xdistr = Xdistr, main = main,
                    Ylabel = Ylabel, Dlabel = Dlabel, Xlabel = Xlabel,
                    xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
                    theme.bw = theme.bw, show.grid = show.grid,
                    cex.main = cex.main, cex.sub = cex.sub,
                    cex.lab = cex.lab, cex.axis = cex.axis,
                    show.uniform.CI = show.uniform.CI)
            }
        }
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
            Y = Y,
            D = D,
            X = X,
            Z = Z,
            weights = weights,
            num.trees = grf.num.trees,
            treat.info = treat.info,
            diff.info = diff.info,
            figure = figure,
            show.uniform.CI = show.uniform.CI,
            CI = CI,
            subtitles = subtitles,
            show.subtitles = show.subtitles,
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
            interval = interval,
            file = file,
            ncols = ncols,
            pool = pool,
            color = color,
            legend.title = legend.title,
            show.all = show.all,
            scale = scale,
            height = height,
            width = width,
            gate = gate
        )
    }
    if (estimator == "dml") {
        output <- interflex.dml(
            data = data,
            Y = Y,
            D = D,
            X = X,
            Z = Z,
            weights = weights,
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
            show.uniform.CI = show.uniform.CI,
            CI = CI,
            subtitles = subtitles,
            show.subtitles = show.subtitles,
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

    if (estimator == "lasso") {
        # When gate = TRUE, basis expansion on discrete X is inappropriate
        if (isTRUE(gate) && basis.type != "none") {
            basis.type <- "none"
        }
        if(isTRUE(gate) || length(unique(data[,X]))<5){
            output <- interflex.lasso_discrete(
                data = data,
                Y = Y,
                D = D,
                X = X,
                Z = Z,
                FE = FE,
                weights = NULL,
                B = nboots,
                alpha = 0.05,
                model.y = 'lasso',
                model.t = 'lasso',
                signal = signal,
                estimand = estimand,
                neval = neval,
                basis.type = basis.type,
                include.interactions = include.interactions,
                poly.degree = poly.degree,
                spline.df = spline.df,
                spline.degree = spline.degree,
                lambda.seq = lambda.seq,
                cores = cores,
                verbose = TRUE,
                treat.info = treat.info,
                diff.info = diff.info,
                figure = figure,
                CI = CI,
                subtitles = subtitles,
                show.subtitles = show.subtitles,
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
        else{
            output <- interflex.lasso(
                data = data,
                Y = Y,
                D = D,
                X = X,
                Z = Z,
                FE = FE,
                weights = NULL,
                B = nboots,
                alpha = 0.05,
                model.y = 'lasso',
                model.t = 'lasso',
                signal = signal,
                estimand = estimand,
                neval = neval,
                basis.type = basis.type,
                include.interactions = include.interactions,
                poly.degree = poly.degree,
                spline.df = spline.df,
                spline.degree = spline.degree,
                lambda.seq = lambda.seq,
                reduce.dimension = reduce.dimension,
                cores = cores,
                verbose = TRUE,
                treat.info = treat.info,
                diff.info = diff.info,
                figure = figure,
                show.uniform.CI = show.uniform.CI,
                CI = CI,
                subtitles = subtitles,
                show.subtitles = show.subtitles,
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

    }

    return(output)
}
