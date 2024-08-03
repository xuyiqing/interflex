interflex.linear <- function(data,
                             Y, # outcome
                             D, # treatment indicator
                             X, # moderator
                             treat.info,
                             diff.info,
                             Z = NULL, # covariates
                             weights = NULL, # weighting variable
                             full.moderate = FALSE, # fully moderated model
                             Z.X = NULL, # fully moderated terms
                             FE = NULL, # fixed effects
                             IV = NULL,
                             neval = 50,
                             X.eval = NULL,
                             method = "linear", ## "probit"; "logit"; "poisson"; "nbinom"
                             vartype = "simu", ## variance type "simu"; "bootstrap"; "delta"
                             vcov.type = "robust", ## "homoscedastic"; "robust"; "cluster"; "pcse"
                             time = NULL,
                             pairwise = TRUE,
                             nboots = 200,
                             nsimu = 1000,
                             parallel = TRUE,
                             cores = 4,
                             cl = NULL, # variable to be clustered on
                             # predict = FALSE,
                             Z.ref = NULL, # same length as Z, set the value of Z when estimating marginal effects/predicted value
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
    WEIGHTS <- NULL
    n <- dim(data)[1]
    diff.values.plot <- diff.info[["diff.values.plot"]]
    diff.values <- diff.info[["diff.values"]]
    difference.name <- diff.info[["difference.name"]]

    treat.type <- treat.info[["treat.type"]]
    if (treat.type == "discrete") {
        other.treat <- treat.info[["other.treat"]]
        other.treat.origin <- names(other.treat)
        names(other.treat.origin) <- other.treat
        all.treat <- treat.info[["all.treat"]]
        all.treat.origin <- names(all.treat)
        names(all.treat.origin) <- all.treat
        base <- treat.info[["base"]]
    }
    if (treat.type == "continuous") {
        D.sample <- treat.info[["D.sample"]]
        label.name <- names(D.sample)
        # names(label.name) <- D.sample
    }
    ncols <- treat.info[["ncols"]]

    if (is.null(cl) == FALSE) { ## find clusters
        clusters <- unique(data[, cl])
        id.list <- split(1:n, data[, cl])
    }

    if (is.null(FE) == TRUE) {
        use_fe <- 0
    } else {
        use_fe <- 1
    }

    # pcse
    if (vcov.type == "pcse") {
        data.cl <- data[, cl]
        data.time <- data[, time]
    }

    # evaluation points
    X.eval0 <- X.eval
    X.eval <- seq(min(data[, X]), max(data[, X]), length.out = neval)
    X.eval <- sort(c(X.eval, X.eval0))
    neval <- length(X.eval)

    # construct the formula
    formula <- paste0(Y, "~", X)
    if (treat.type == "discrete") {
        for (char in other.treat) {
            data[, paste0("D", ".", char)] <- as.numeric(data[, D] == char)
            data[, paste0("DX", ".", char)] <- as.numeric(data[, D] == char) * data[, X]
            formula <- paste0(formula, "+", paste0("D", ".", char), "+", paste0("DX", ".", char))
        }
    } else {
        data[, "DX"] <- data[, D] * data[, X]
        formula <- paste0(formula, "+", D, "+DX")
    }

    if (is.null(Z) == FALSE) {
        formula <- paste0(formula, "+", paste0(Z, collapse = "+"))
        if (full.moderate == TRUE) {
            formula <- paste0(formula, "+", paste0(Z.X, collapse = "+"))
        }
    }

    if (use_fe == 1) {
        formula <- paste0(formula, "|", paste0(FE, collapse = "+"))
        if (vcov.type == "cluster") {
            formula <- paste0(formula, "| 0 |", paste0(cl, collapse = "+"))
        }
    }

    # iv formula
    if (!is.null(IV)) {
        if (use_fe == 0) {
            # mod2 <- ivreg(Y~D+X+D:X+Z1|W+X+W:X+Z1,data=s1)
            formula.iv <- X
            for (sub.iv in IV) {
                data[, paste0(X, ".", sub.iv)] <- data[, sub.iv] * data[, X]
                formula.iv <- paste0(formula.iv, "+", sub.iv)
                formula.iv <- paste0(formula.iv, "+", paste0(X, ".", sub.iv))
            }
            if (is.null(Z) == FALSE) {
                formula.iv <- paste0(formula.iv, "+", paste0(Z, collapse = "+"))
                if (full.moderate == TRUE) {
                    formula.iv <- paste0(formula.iv, "+", paste0(Z.X, collapse = "+"))
                }
            }
            formula <- paste0(formula, "|", formula.iv)
        }
        if (use_fe == 1) {
            # mod2 <- felm(Y~X+Z1+Z2|unit+year|(D|DX ~ W+WX)|0,data = s1)
            formula <- paste0(Y, "~", X)
            formula.iv <- ""
            name.update <- c(X)
            if (is.null(Z) == FALSE) {
                formula <- paste0(formula, "+", paste0(Z, collapse = "+"))
                name.update <- c(name.update, Z)
                if (full.moderate == TRUE) {
                    formula <- paste0(formula, "+", paste0(Z.X, collapse = "+"))
                    name.update <- c(name.update, Z.X)
                }
            }
            if (treat.type == "discrete") {
                for (char in other.treat) {
                    formula.iv <- paste0(formula.iv, "|", paste0("D", ".", char), "|", paste0("DX", ".", char))
                    name.update <- c(name.update, paste0("D", ".", char), paste0("DX", ".", char))
                }
            } else {
                formula.iv <- paste0(formula.iv, "|", D, "|DX")
                name.update <- c(name.update, D, "DX")
            }
            formula.iv <- substring(formula.iv, 2)
            formula.iv <- paste0(formula.iv, " ~ ")

            for (sub.iv in IV) {
                data[, paste0(X, ".", sub.iv)] <- data[, sub.iv] * data[, X]
                formula.iv <- paste0(formula.iv, "+", sub.iv)
                formula.iv <- paste0(formula.iv, "+", paste0(X, ".", sub.iv))
            }

            formula <- paste0(formula, "|", paste0(FE, collapse = "+"), "|")
            formula.iv <- paste0("(", formula.iv, ")")
            formula <- paste0(formula, formula.iv)
            if (vcov.type == "cluster") {
                formula <- paste0(formula, "|", paste0(cl, collapse = "+"))
            }
        }
    }

    formula <- as.formula(formula)
    if (is.null(weights) == TRUE) {
        w <- rep(1, n)
    } else {
        w <- data[, weights]
    }
    data[["WEIGHTS"]] <- w

    # model fit
    if (method == "linear") {
        if (is.null(IV)) {
            if (use_fe == 0) {
                suppressWarnings(
                    model <- glm(formula, data = data, weights = WEIGHTS)
                )
            }
            if (use_fe == 1) {
                suppressWarnings(
                    model <- felm(formula, data = data, weights = w)
                )
                model$converged <- TRUE
            }
        }
        if (!is.null(IV)) {
            if (use_fe == 0) {
                suppressWarnings(
                    model <- ivreg(formula, data = data, weights = WEIGHTS)
                )
                model$converged <- TRUE
            }
            if (use_fe == 1) {
                suppressWarnings(
                    model <- felm(formula, data = data, weights = w)
                )
                model$converged <- TRUE
            }
        }
    }
    if (method == "logit") {
        suppressWarnings(
            model <- glm(formula, data = data, family = binomial(link = "logit"), weights = WEIGHTS)
        )
    }
    if (method == "probit") {
        suppressWarnings(
            model <- glm(formula, data = data, family = binomial(link = "probit"), weights = WEIGHTS)
        )
    }
    if (method == "poisson") {
        suppressWarnings(
            model <- glm(formula, data = data, family = poisson, weights = WEIGHTS)
        )
    }
    if (method == "nbinom") {
        suppressWarnings(
            model <- glm.nb(formula, data = data, weights = WEIGHTS)
        )
    }

    if (use_fe == 0) {
        if (model$converged == FALSE) {
            stop("Linear estimator can't converge.")
        }
    }

    model.coef <- coef(model)
    if (!is.null(IV)) {
        if (use_fe == 1) {
            names(model.coef) <- name.update
        }
    }

    if (use_fe == 0) {
        if (vcov.type == "homoscedastic") {
            model.vcov <- vcov(model)
        }
        if (vcov.type == "robust") {
            model.vcov <- vcovHC(model, type = "HC1")
        }
        if (vcov.type == "cluster") {
            model.vcov <- vcovCluster(model, cluster = data[, cl])
        }
        if (vcov.type == "pcse") {
            model.vcov <- pcse(model, pairwise = pairwise, groupN = data.cl, groupT = data.time)$vcov
            rownames(model.vcov)[1] <- "(Intercept)"
            colnames(model.vcov)[1] <- "(Intercept)"
        }
    } else { # fe vcov
        if (vcov.type == "homoscedastic") {
            model.vcov <- vcov(model, type = "iid")
        }
        if (vcov.type == "robust") {
            model.vcov <- vcov(model, type = "robust")
        }
        if (vcov.type == "cluster") {
            model.vcov <- vcov(model, type = "cluster")
        }
    }

    if (!is.null(IV)) {
        if (use_fe == 1) {
            rownames(model.vcov) <- colnames(model.vcov) <- name.update
        }
    }

    model.df <- model$df.residual
    model.coef[which(is.na(model.coef))] <- 0
    model.vcov[which(is.na(model.vcov))] <- 0
    model.vcov.all <- matrix(0, nrow = length(model.coef), ncol = length(model.coef))
    colnames(model.vcov.all) <- names(model.coef)
    rownames(model.vcov.all) <- names(model.coef)
    for (a1 in rownames(model.vcov)) {
        for (a2 in colnames(model.vcov)) {
            model.vcov.all[a1, a2] <- model.vcov[a1, a2]
        }
    }

    if (isSymmetric.matrix(model.vcov.all, tol = 1e-6) == FALSE) {
        warning(paste0("Option vcov.type==", vcov.type, "leads to unstable standard error in linear estimator, by default use homoscedastic standard error as an alternative.\n"))
        model.vcov.all <- vcov(model)
        if (!is.null(IV)) {
            if (use_fe == 1) {
                rownames(model.vcov.all) <- colnames(model.vcov.all) <- name.update
            }
        }
        model.vcov.all[which(is.na(model.vcov.all))] <- 0
    }
    model.vcov <- model.vcov.all


    ## Function A
    # 1, estimate treatment effects/marginal effects given model.coef
    # 2, estimate difference of treatment effects/marginal effects at different values of the moderator
    # 3, input: model.coef;X.eval;char(discrete)/D.ref(continuous);diff.values(default to NULL);
    # 4, output: marginal effects/treatment effects/E.pred/E.base/diff.estimate

    gen.general.TE <- function(model.coef, char = NULL, D.ref = NULL, diff.values = NULL) {
        if (is.null(char) == TRUE) {
            treat.type <- "continuous"
        }
        if (is.null(D.ref) == TRUE) {
            treat.type <- "discrete"
        }

        gen.TE <- function(model.coef, X.eval) {
            neval <- length(X.eval)
            if (treat.type == "discrete") {
                link.1 <- model.coef["(Intercept)"] + X.eval * model.coef[X] + 1 * model.coef[paste0("D.", char)] + X.eval * model.coef[paste0("DX.", char)]
                link.0 <- model.coef["(Intercept)"] + X.eval * model.coef[X]
                if (is.null(Z) == FALSE) {
                    for (a in Z) {
                        target.Z <- Z.ref[a]
                        link.1 <- link.1 + target.Z * model.coef[a]
                        link.0 <- link.0 + target.Z * model.coef[a]
                        if (full.moderate == TRUE) {
                            link.1 <- link.1 + target.Z * model.coef[paste0(a, ".X")] * X.eval
                            link.0 <- link.0 + target.Z * model.coef[paste0(a, ".X")] * X.eval
                        }
                    }
                }
                if (method == "linear") {
                    TE <- link.1 - link.0
                    E.pred <- link.1
                    E.base <- link.0
                }

                if (method == "logit") {
                    E.pred <- E.prob.1 <- exp(link.1) / (1 + exp(link.1))
                    E.base <- E.prob.0 <- exp(link.0) / (1 + exp(link.0))
                    TE <- E.prob.1 - E.prob.0
                }

                if (method == "probit") {
                    E.pred <- E.prob.1 <- pnorm(link.1, 0, 1)
                    E.base <- E.prob.0 <- pnorm(link.0, 0, 1)
                    TE <- E.prob.1 - E.prob.0
                }

                if (method == "poisson" | method == "nbinom") {
                    E.pred <- E.y.1 <- exp(link.1)
                    E.base <- E.y.0 <- exp(link.0)
                    TE <- E.y.1 - E.y.0
                }
                names(link.1) <- rep(paste0("link.1.", char), neval)
                names(link.0) <- rep(paste0("link.0.", base), neval)
                names(TE) <- rep(paste0("TE.", char), neval)
                names(E.pred) <- rep(paste0("Predict.", char), neval)
                names(E.base) <- rep(paste0("Predict.", base), neval)
                gen.TE.output <- list(
                    TE = TE, E.pred = E.pred, E.base = E.base,
                    link.1 = link.1, link.0 = link.0
                )
            }

            if (treat.type == "continuous") {
                link <- model.coef["(Intercept)"] + X.eval * model.coef[X] + model.coef[D] * D.ref + model.coef["DX"] * X.eval * D.ref
                if (is.null(Z) == FALSE) {
                    for (a in Z) {
                        target.Z <- Z.ref[a]
                        link <- link + target.Z * model.coef[a]
                        if (full.moderate == TRUE) {
                            link <- link + target.Z * model.coef[paste0(a, ".X")] * X.eval
                        }
                    }
                }
                if (method == "logit") {
                    ME <- exp(link) / (1 + exp(link))^2 * (model.coef[D] + model.coef["DX"] * X.eval)
                    E.pred <- exp(link) / (1 + exp(link))
                }
                if (method == "probit") {
                    ME <- (model.coef[D] + model.coef["DX"] * X.eval) * dnorm(link)
                    E.pred <- pnorm(link, 0, 1)
                }
                if (method == "linear") {
                    ME <- model.coef[D] + model.coef["DX"] * X.eval
                    E.pred <- link
                }
                if (method == "poisson" | method == "nbinom") {
                    ME <- exp(link) * (model.coef[D] + model.coef["DX"] * X.eval)
                    E.pred <- exp(link)
                }
                names(link) <- rep("link", neval)
                names(ME) <- rep(paste0("ME.", names(D.sample)[D.sample == D.ref]), neval)
                names(E.pred) <- rep(paste0("Predict.", names(D.sample)[D.sample == D.ref]), neval)
                gen.TE.output <- list(ME = ME, E.pred = E.pred, link = link)
            }
            return(gen.TE.output)
        }

        gen.TE.fe <- function(model.coef, X.eval) {
            neval <- length(X.eval)
            if (treat.type == "discrete") {
                TE <- model.coef[paste0("D.", char)] + X.eval * model.coef[paste0("DX.", char)]
                names(TE) <- rep(paste0("TE.", char), neval)
                E.pred <- rep(0, neval) # doesn't return prediction value for fixed effects
                E.base <- rep(0, neval)
                gen.TE.output <- list(TE = TE, E.pred = E.pred, E.base = E.base, link.1 = E.pred, link.0 = E.base)
            }
            if (treat.type == "continuous") {
                ME <- model.coef[D] + model.coef["DX"] * X.eval
                names(ME) <- rep(paste0("ME.", names(D.sample)[D.sample == D.ref]), neval)
                E.pred <- rep(0, neval)
                gen.TE.output <- list(ME = ME, E.pred = E.pred, link = E.pred)
            }
            return(gen.TE.output)
        }

        if (use_fe == 0) {
            gen.TE.output <- gen.TE(model.coef = model.coef, X.eval = X.eval)
        }
        if (use_fe == 1) {
            gen.TE.output <- gen.TE.fe(model.coef = model.coef, X.eval = X.eval)
        }

        if (is.null(diff.values) == FALSE) {
            if (use_fe == 0) {
                gen.TE.diff <- gen.TE(model.coef = model.coef, X.eval = diff.values)
            }
            if (use_fe == 1) {
                gen.TE.diff <- gen.TE.fe(model.coef = model.coef, X.eval = diff.values)
            }

            if (treat.type == "discrete") {
                target <- gen.TE.diff$TE
            }
            if (treat.type == "continuous") {
                target <- gen.TE.diff$ME
            }
            if (length(diff.values) == 2) {
                difference.temp <- c(target[2] - target[1])
            }
            if (length(diff.values) == 3) {
                difference.temp <- c(
                    target[2] - target[1],
                    target[3] - target[2],
                    target[3] - target[1]
                )
            }

            if (treat.type == "discrete") {
                names(difference.temp) <- paste0(char, ".", difference.name)
            }

            if (treat.type == "continuous") {
                names(difference.temp) <- paste0(names(D.sample)[D.sample == D.ref], ".", difference.name)
            }
        } else {
            difference.temp <- NULL
        }
        if (treat.type == "discrete") {
            return(list(
                TE = gen.TE.output$TE,
                E.pred = gen.TE.output$E.pred,
                E.base = gen.TE.output$E.base,
                link.1 = gen.TE.output$link.1,
                link.0 = gen.TE.output$link.0,
                diff.estimate = difference.temp
            ))
        }

        if (treat.type == "continuous") {
            return(list(
                ME = gen.TE.output$ME,
                E.pred = gen.TE.output$E.pred,
                link = gen.TE.output$link,
                diff.estimate = difference.temp
            ))
        }
    }

    ## Function B delta method
    # 1, estimate sd of treatment effects/marginal effects using delta method
    # 2, estimate sd of predicted values using delta method
    # 3, estimate sd of diff.estimate using delta method
    # 4, estimate vcov of ME/TE using delta method
    # 5, estimate sd of linear predictor using delta method
    # 6, input: model.coef; model.vcov; char(discrete)/D.ref(continuous);diff.values
    # 7, output: sd of TE/ME; sd of diff.values; vcov of ME/TE
    gen.delta.TE <- function(model.coef, model.vcov, char = NULL, D.ref = NULL, diff.values = NULL, vcov = FALSE) {
        if (is.null(char) == TRUE) {
            treat.type <- "continuous"
            flag <- 1
        }
        if (is.null(D.ref) == TRUE) {
            treat.type <- "discrete"
            if (char == base) {
                flag <- 0
            } else {
                flag <- 1
            }
        }

        # sd for TE/ME
        gen.sd.fe <- function(x, to.diff = FALSE) {
            if (treat.type == "discrete") {
                target.slice <- c(paste0("D.", char), paste0("DX.", char))
                vec.1 <- c(1, x)
                vec.0 <- c(0, 0)
                vec <- vec.1 - vec.0
                temp.vcov.matrix <- model.vcov[target.slice, target.slice]
                if (to.diff == TRUE) {
                    return(list(vec = vec, temp.vcov.matrix = temp.vcov.matrix))
                }
                delta.sd <- sqrt((t(vec) %*% temp.vcov.matrix %*% vec)[1, 1])
                return(delta.sd)
            }

            if (treat.type == "continuous") {
                target.slice <- c(D, "DX")
                vec <- c(1, x)
                temp.vcov.matrix <- model.vcov[target.slice, target.slice]
                if (to.diff == TRUE) {
                    return(list(vec = vec, temp.vcov.matrix = temp.vcov.matrix))
                }

                delta.sd <- sqrt((t(vec) %*% temp.vcov.matrix %*% vec)[1, 1])
                return(delta.sd)
            }
        }

        gen.sd <- function(x, to.diff = FALSE) {
            if (use_fe == TRUE) {
                return(gen.sd.fe(x = x, to.diff = to.diff))
            }

            if (treat.type == "discrete") {
                link.1 <- model.coef["(Intercept)"] + x * model.coef[X] + 1 * model.coef[paste0("D.", char)] + x * model.coef[paste0("DX.", char)]
                link.0 <- model.coef["(Intercept)"] + x * model.coef[X]
                if (is.null(Z) == FALSE) {
                    for (a in Z) {
                        target.Z <- Z.ref[a]
                        link.1 <- link.1 + target.Z * model.coef[a]
                        link.0 <- link.0 + target.Z * model.coef[a]
                        if (full.moderate == TRUE) {
                            link.1 <- link.1 + target.Z * model.coef[paste0(a, ".X")] * x
                            link.0 <- link.0 + target.Z * model.coef[paste0(a, ".X")] * x
                        }
                    }
                }

                if (is.null(Z) == FALSE) {
                    if (full.moderate == FALSE) {
                        vec.1 <- c(1, x, 1, x, Z.ref)
                        vec.0 <- c(1, x, 0, 0, Z.ref)
                        target.slice <- c("(Intercept)", X, paste0("D.", char), paste0("DX.", char), Z)
                    }
                    if (full.moderate == TRUE) {
                        vec.1 <- c(1, x, 1, x, Z.ref, x * Z.ref)
                        vec.0 <- c(1, x, 0, 0, Z.ref, x * Z.ref)
                        target.slice <- c("(Intercept)", X, paste0("D.", char), paste0("DX.", char), Z, Z.X)
                    }
                } else {
                    vec.1 <- c(1, x, 1, x)
                    vec.0 <- c(1, x, 0, 0)
                    target.slice <- c("(Intercept)", X, paste0("D.", char), paste0("DX.", char))
                }
                temp.vcov.matrix <- model.vcov[target.slice, target.slice]
                if (method == "logit") {
                    vec <- vec.1 * exp(link.1) / (1 + exp(link.1))^2 - vec.0 * exp(link.0) / (1 + exp(link.0))^2
                }
                if (method == "probit") {
                    vec <- vec.1 * dnorm(link.1) - vec.0 * dnorm(link.0)
                }
                if (method == "poisson" | method == "nbinom") {
                    vec <- vec.1 * exp(link.1) - vec.0 * exp(link.0)
                }
                if (method == "linear") {
                    vec <- vec.1 - vec.0
                }
                if (to.diff == TRUE) {
                    return(list(vec = vec, temp.vcov.matrix = temp.vcov.matrix))
                }

                delta.sd <- sqrt((t(vec) %*% temp.vcov.matrix %*% vec)[1, 1])
                return(delta.sd)
            }

            if (treat.type == "continuous") {
                link <- model.coef["(Intercept)"] + x * model.coef[X] + model.coef[D] * D.ref + model.coef["DX"] * x * D.ref
                if (is.null(Z) == FALSE) {
                    for (a in Z) {
                        target.Z <- Z.ref[a]
                        link <- link + target.Z * model.coef[a]
                        if (full.moderate == TRUE) {
                            link <- link + x * target.Z * model.coef[paste0(a, ".X")]
                        }
                    }
                }

                if (is.null(Z) == FALSE) {
                    if (full.moderate == FALSE) {
                        vec1 <- c(1, x, D.ref, D.ref * x, Z.ref)
                        vec0 <- c(0, 0, 1, x, rep(0, length(Z)))
                        target.slice <- c("(Intercept)", X, D, "DX", Z)
                    }
                    if (full.moderate == TRUE) {
                        vec1 <- c(1, x, D.ref, D.ref * x, Z.ref, Z.ref * x)
                        vec0 <- c(0, 0, 1, x, rep(0, 2 * length(Z)))
                        target.slice <- c("(Intercept)", X, D, "DX", Z, Z.X)
                    }
                } else {
                    vec1 <- c(1, x, D.ref, D.ref * x)
                    vec0 <- c(0, 0, 1, x)
                    target.slice <- c("(Intercept)", X, D, "DX")
                }

                temp.vcov.matrix <- model.vcov[target.slice, target.slice]

                if (method == "logit") {
                    vec <- -(model.coef[D] + x * model.coef["DX"]) * (exp(link) - exp(-link)) / (2 + exp(link) + exp(-link))^2 * vec1 + exp(link) / (1 + exp(link))^2 * vec0
                }
                if (method == "probit") {
                    vec <- dnorm(link) * vec0 - (model.coef[D] + x * model.coef["DX"]) * link * dnorm(link) * vec1
                }
                if (method == "poisson" | method == "nbinom") {
                    vec <- (model.coef[D] + x * model.coef["DX"]) * exp(link) * vec1 + exp(link) * vec0
                }
                if (method == "linear") {
                    vec <- vec0
                }

                if (to.diff == TRUE) {
                    return(list(vec = vec, temp.vcov.matrix = temp.vcov.matrix))
                }

                delta.sd <- sqrt((t(vec) %*% temp.vcov.matrix %*% vec)[1, 1])
                return(delta.sd)
            }
        }

        if (flag == 1) {
            TE.sd <- c(sapply(X.eval, function(x) gen.sd(x)))
        } else {
            TE.sd <- NULL
        }

        if (treat.type == "discrete" & is.null(TE.sd) == FALSE) {
            names(TE.sd) <- rep(paste0("sd.", char), neval)
        }

        if (treat.type == "continuous") {
            names(TE.sd) <- rep(paste0("sd.", names(D.sample)[D.sample == D.ref]), neval)
        }

        # sd for linear predictor
        gen.link.sd <- function(x, base = FALSE) {
            if (treat.type == "discrete") {
                if (is.null(Z) == FALSE) {
                    if (base == FALSE) {
                        vec <- c(1, x, 1, x, Z.ref)
                        target.slice <- c("(Intercept)", X, paste0("D.", char), paste0("DX.", char), Z)
                    } else {
                        vec <- c(1, x, Z.ref)
                        target.slice <- c("(Intercept)", X, Z)
                    }

                    if (full.moderate == TRUE) {
                        vec <- c(vec, Z.ref * x)
                        target.slice <- c(target.slice, Z.X)
                    }
                } else {
                    if (base == FALSE) {
                        vec <- c(1, x, 1, x)
                        target.slice <- c("(Intercept)", X, paste0("D.", char), paste0("DX.", char))
                    } else {
                        vec <- c(1, x)
                        target.slice <- c("(Intercept)", X)
                    }
                }
                temp.vcov.matrix <- model.vcov[target.slice, target.slice]
                link.sd <- sqrt((t(vec) %*% temp.vcov.matrix %*% vec)[1, 1])
                return(link.sd)
            }

            if (treat.type == "continuous") {
                if (is.null(Z) == FALSE) {
                    vec <- c(1, x, D.ref, D.ref * x, Z.ref)
                    target.slice <- c("(Intercept)", X, D, "DX", Z)
                    if (full.moderate == TRUE) {
                        vec <- c(vec, Z.ref * x)
                        target.slice <- c(target.slice, Z.X)
                    }
                } else {
                    vec <- c(1, x, D.ref, D.ref * x)
                    target.slice <- c("(Intercept)", X, D, "DX")
                }

                temp.vcov.matrix <- model.vcov[target.slice, target.slice]
                link.sd <- sqrt((t(vec) %*% temp.vcov.matrix %*% vec)[1, 1])
                return(link.sd)
            }
        }

        gen.link.sd.fe <- function(x, base = FALSE) {
            return(0)
        }

        if (use_fe == FALSE) {
            if (treat.type == "discrete") {
                if (char == base) {
                    link.sd <- c(sapply(X.eval, function(x) gen.link.sd(x, base = TRUE)))
                } else {
                    link.sd <- c(sapply(X.eval, function(x) gen.link.sd(x)))
                }
            } else {
                link.sd <- c(sapply(X.eval, function(x) gen.link.sd(x)))
            }
        } else {
            link.sd <- rep(0, length(X.eval))
        }


        # sd for predicted value
        gen.predict.sd <- function(x, base = FALSE) {
            if (treat.type == "discrete") {
                if (base == FALSE) {
                    link <- model.coef["(Intercept)"] + x * model.coef[X] + 1 * model.coef[paste0("D.", char)] + x * model.coef[paste0("DX.", char)]
                }
                if (base == TRUE) {
                    link <- model.coef["(Intercept)"] + x * model.coef[X]
                }
                if (is.null(Z) == FALSE) {
                    for (a in Z) {
                        target.Z <- Z.ref[a]
                        link <- link + target.Z * model.coef[a]
                        if (full.moderate == TRUE) {
                            link <- link + target.Z * model.coef[paste0(a, ".X")] * x
                        }
                    }
                }

                if (is.null(Z) == FALSE) {
                    if (base == FALSE) {
                        vec <- c(1, x, 1, x, Z.ref)
                        target.slice <- c("(Intercept)", X, paste0("D.", char), paste0("DX.", char), Z)
                    } else {
                        vec <- c(1, x, Z.ref)
                        target.slice <- c("(Intercept)", X, Z)
                    }

                    if (full.moderate == TRUE) {
                        vec <- c(vec, Z.ref * x)
                        target.slice <- c(target.slice, Z.X)
                    }
                } else {
                    if (base == FALSE) {
                        vec <- c(1, x, 1, x)
                        target.slice <- c("(Intercept)", X, paste0("D.", char), paste0("DX.", char))
                    } else {
                        vec <- c(1, x)
                        target.slice <- c("(Intercept)", X)
                    }
                }

                temp.vcov.matrix <- model.vcov[target.slice, target.slice]
                if (method == "logit") {
                    vec <- vec * exp(link) / (1 + exp(link))^2
                }
                if (method == "probit") {
                    vec <- vec * dnorm(link)
                }
                if (method == "poisson" | method == "nbinom") {
                    vec <- vec * exp(link)
                }
                if (method == "linear") {
                    vec <- vec
                }
                predict.sd <- sqrt((t(vec) %*% temp.vcov.matrix %*% vec)[1, 1])
                return(predict.sd)
            }

            if (treat.type == "continuous") {
                link <- model.coef["(Intercept)"] + x * model.coef[X] + model.coef[D] * D.ref + model.coef["DX"] * x * D.ref
                if (is.null(Z) == FALSE) {
                    for (a in Z) {
                        target.Z <- Z.ref[a]
                        link <- link + target.Z * model.coef[a]
                        if (full.moderate == TRUE) {
                            link <- link + target.Z * model.coef[paste0(a, ".X")] * x
                        }
                    }
                }

                if (is.null(Z) == FALSE) {
                    vec <- c(1, x, D.ref, D.ref * x, Z.ref)
                    target.slice <- c("(Intercept)", X, D, "DX", Z)
                    if (full.moderate == TRUE) {
                        vec <- c(vec, Z.ref * x)
                        target.slice <- c(target.slice, Z.X)
                    }
                } else {
                    vec <- c(1, x, D.ref, D.ref * x)
                    target.slice <- c("(Intercept)", X, D, "DX")
                }

                temp.vcov.matrix <- model.vcov[target.slice, target.slice]
                if (method == "logit") {
                    vec <- vec * exp(link) / (1 + exp(link))^2
                }
                if (method == "probit") {
                    vec <- vec * dnorm(link)
                }
                if (method == "poisson" | method == "nbinom") {
                    vec <- vec * exp(link)
                }
                if (method == "linear") {
                    vec <- vec
                }
                predict.sd <- sqrt((t(vec) %*% temp.vcov.matrix %*% vec)[1, 1])
                return(predict.sd)
            }
        }

        gen.predict.sd.fe <- function(x, base = FALSE) {
            return(0)
        }

        if (use_fe == FALSE) {
            if (treat.type == "discrete") {
                if (char == base) {
                    predict.sd <- c(sapply(X.eval, function(x) gen.predict.sd(x, base = TRUE)))
                } else {
                    predict.sd <- c(sapply(X.eval, function(x) gen.predict.sd(x)))
                }
            } else {
                predict.sd <- c(sapply(X.eval, function(x) gen.predict.sd(x)))
            }
        } else {
            predict.sd <- rep(0, length(X.eval))
        }

        if (treat.type == "discrete") {
            names(predict.sd) <- rep(paste0("predict.sd.", char), neval)
        }
        if (treat.type == "continuous") {
            names(predict.sd) <- rep(paste0("predict.sd.", names(D.sample)[D.sample == D.ref]), neval)
        }

        # sd for diff.estimates
        if (is.null(diff.values) == FALSE & flag == 1) {
            if (length(diff.values) == 2) {
                vec.list2 <- gen.sd(x = diff.values[2], to.diff = TRUE)
                vec.list1 <- gen.sd(x = diff.values[1], to.diff = TRUE)
                vec1 <- vec.list1$vec
                vec2 <- vec.list2$vec
                vec <- vec2 - vec1
                temp.vcov.matrix <- vec.list1$temp.vcov.matrix
                delta.diff.sd <- c(sqrt((t(vec) %*% temp.vcov.matrix %*% vec)[1, 1]))
            }

            if (length(diff.values) == 3) {
                vec.list3 <- gen.sd(x = diff.values[3], to.diff = TRUE)
                vec.list2 <- gen.sd(x = diff.values[2], to.diff = TRUE)
                vec.list1 <- gen.sd(x = diff.values[1], to.diff = TRUE)
                vec1 <- vec.list1$vec
                vec2 <- vec.list2$vec
                vec3 <- vec.list3$vec
                vec21 <- vec2 - vec1
                vec32 <- vec3 - vec2
                vec31 <- vec3 - vec1
                temp.vcov.matrix <- vec.list1$temp.vcov.matrix

                delta.diff.sd <- c(
                    sqrt((t(vec21) %*% temp.vcov.matrix %*% vec21)[1, 1]),
                    sqrt((t(vec32) %*% temp.vcov.matrix %*% vec32)[1, 1]),
                    sqrt((t(vec31) %*% temp.vcov.matrix %*% vec31)[1, 1])
                )
            }

            if (treat.type == "discrete") {
                names(delta.diff.sd) <- paste0("sd.", char, ".", difference.name)
            }

            if (treat.type == "continuous") {
                names(delta.diff.sd) <- paste0("sd.", names(D.sample)[D.sample == D.ref], ".", difference.name)
            }
        } else {
            delta.diff.sd <- NULL
        }

        # vcov matrix
        if (flag == 1 & vcov == TRUE) {
            gen.cov.element <- function(x1, x2) {
                vec.list1 <- gen.sd(x1, to.diff = TRUE)
                vec.list2 <- gen.sd(x2, to.diff = TRUE)
                vec1 <- vec.list1$vec
                vec2 <- vec.list2$vec
                temp.vcov.matrix <- vec.list1$temp.vcov.matrix
                cov.element <- (t(vec1) %*% temp.vcov.matrix %*% vec2)[1, 1]
                return(cov.element)
            }
            gen.cov.vec <- function(colvec, x0) {
                output <- sapply(colvec, function(xx) {
                    gen.cov.element(xx, x0)
                })
                return(output)
            }
            TE.sd.vcov <- as.matrix(sapply(X.eval, function(xx) {
                gen.cov.vec(X.eval, xx)
            }))
        } else {
            TE.sd.vcov <- NULL
        }

        # output
        if (treat.type == "discrete") {
            return(list(
                TE.sd = TE.sd,
                predict.sd = predict.sd,
                link.sd = link.sd,
                sd.diff.estimate = delta.diff.sd,
                TE.vcov = TE.sd.vcov
            ))
        }
        if (treat.type == "continuous") {
            return(list(
                ME.sd = TE.sd,
                predict.sd = predict.sd,
                link.sd = link.sd,
                sd.diff.estimate = delta.diff.sd,
                ME.vcov = TE.sd.vcov
            ))
        }
    }


    ## Function C ATE(weighted average)
    # 1, estimate average treatment effects(ATE)/average marginal effects(AME)
    # 2, estimate sd of ATE/AME using delta method
    # 3, input: data(weights); model.coef; model.vcov; char(discrete); delta(Default to FALSE)
    # 4, output: ATE for char/AME; sd for ATE/AME

    gen.ATE.fe <- function(data, model.coef, model.vcov, char = NULL, delta = FALSE) {
        if (treat.type == "discrete") {
            sub.data <- data[data[, D] == char, ]
            weights <- data[data[, D] == char, "WEIGHTS"]
            TE <- model.coef[paste0("D.", char)] + sub.data[, X] * model.coef[paste0("DX.", char)]
            ATE <- weighted.mean(TE, weights, na.rm = TRUE)
            if (delta == FALSE) {
                return(ATE)
            }
        }

        if (treat.type == "continuous") {
            weights <- data[, "WEIGHTS"]
            ME <- model.coef[D] + model.coef["DX"] * data[, X]
            AME <- weighted.mean(ME, weights, na.rm = TRUE)
            if (delta == FALSE) {
                return(AME)
            }
        }
        if (delta == TRUE) {
            gen.ATE.sd.vec <- function(index) {
                if (treat.type == "discrete") {
                    x <- sub.data[index, X]
                    vec.1 <- c(1, x)
                    vec.0 <- c(0, 0)
                    vec <- vec.1 - vec.0
                    target.slice <- c(paste0("D.", char), paste0("DX.", char))
                    temp.vcov.matrix <- model.vcov[target.slice, target.slice]
                    return(vec)
                }

                if (treat.type == "continuous") {
                    vec <- vec0 <- c(1, data[index, X])
                    target.slice <- c(D, "DX")
                    temp.vcov.matrix <- model.vcov[target.slice, target.slice]
                    return(vec)
                }
            }

            if (treat.type == "discrete") {
                index.all <- c(1:dim(sub.data)[1])
                vec.all <- sapply(index.all, function(x) gen.ATE.sd.vec(x))
                vec.mean <- apply(vec.all, 1, function(x) weighted.mean(x, weights))
                target.slice <- c(paste0("D.", char), paste0("DX.", char))
                temp.vcov.matrix <- model.vcov[target.slice, target.slice]
                ATE.sd <- sqrt((t(vec.mean) %*% temp.vcov.matrix %*% vec.mean)[1, 1])
                return(list(ATE = ATE, sd = ATE.sd))
            }

            if (treat.type == "continuous") {
                index.all <- c(1:dim(data)[1])
                vec.all <- sapply(index.all, function(x) gen.ATE.sd.vec(x))
                vec.mean <- apply(vec.all, 1, function(x) weighted.mean(x, weights))
                target.slice <- c(D, "DX")
                temp.vcov.matrix <- model.vcov[target.slice, target.slice]
                AME.sd <- sqrt((t(vec.mean) %*% temp.vcov.matrix %*% vec.mean)[1, 1])
                return(list(AME = AME, sd = AME.sd))
            }
        }
    }


    gen.ATE <- function(data, model.coef, model.vcov, char = NULL, delta = FALSE) {
        if (use_fe == TRUE) {
            return(gen.ATE.fe(data = data, model.coef = model.coef, model.vcov = model.vcov, char = char, delta = delta))
        }

        if (treat.type == "discrete") {
            sub.data <- data[data[, D] == char, ]
            weights <- data[data[, D] == char, "WEIGHTS"]
            link.1 <- model.coef["(Intercept)"] + sub.data[, X] * model.coef[X] +
                model.coef[paste0("D.", char)] + sub.data[, X] * model.coef[paste0("DX.", char)]
            link.0 <- model.coef["(Intercept)"] + sub.data[, X] * model.coef[X]
            if (is.null(Z) == FALSE) {
                for (a in Z) {
                    target.Z <- sub.data[, a]
                    link.1 <- link.1 + target.Z * model.coef[a]
                    link.0 <- link.0 + target.Z * model.coef[a]
                    if (full.moderate == TRUE) {
                        link.1 <- link.1 + target.Z * model.coef[paste0(a, ".X")] * sub.data[, X]
                        link.0 <- link.0 + target.Z * model.coef[paste0(a, ".X")] * sub.data[, X]
                    }
                }
            }
            if (method == "linear") {
                TE <- link.1 - link.0
            }

            if (method == "logit") {
                E.pred <- E.prob.1 <- exp(link.1) / (1 + exp(link.1))
                E.base <- E.prob.0 <- exp(link.0) / (1 + exp(link.0))
                TE <- E.prob.1 - E.prob.0
            }

            if (method == "probit") {
                E.pred <- E.prob.1 <- pnorm(link.1, 0, 1)
                E.base <- E.prob.0 <- pnorm(link.0, 0, 1)
                TE <- E.prob.1 - E.prob.0
            }

            if (method == "poisson" | method == "nbinom") {
                E.pred <- E.y.1 <- exp(link.1)
                E.base <- E.y.0 <- exp(link.0)
                TE <- E.y.1 - E.y.0
            }
            ATE <- weighted.mean(TE, weights, na.rm = TRUE)
            if (delta == FALSE) {
                return(ATE)
            }
        }

        if (treat.type == "continuous") {
            weights <- data[, "WEIGHTS"]
            link <- model.coef["(Intercept)"] + data[, X] * model.coef[X] + model.coef[D] * data[, D] + model.coef["DX"] * data[, X] * data[, D]
            if (is.null(Z) == FALSE) {
                for (a in Z) {
                    target.Z <- data[, a]
                    link <- link + target.Z * model.coef[a]
                    if (full.moderate == TRUE) {
                        link <- link + target.Z * model.coef[paste0(a, ".X")] * data[, X]
                    }
                }
            }
            # return(link)
            if (method == "logit") {
                ME <- exp(link) / (1 + exp(link))^2 * (model.coef[D] + model.coef["DX"] * data[, X])
            }
            if (method == "probit") {
                ME <- (model.coef[D] + model.coef["DX"] * data[, X]) * dnorm(link)
            }
            if (method == "linear") {
                ME <- model.coef[D] + model.coef["DX"] * data[, X]
            }
            if (method == "poisson" | method == "nbinom") {
                ME <- exp(link) * (model.coef[D] + model.coef["DX"] * data[, X])
            }
            AME <- weighted.mean(ME, weights, na.rm = TRUE)
            if (delta == FALSE) {
                return(AME)
            }
        }

        if (delta == TRUE) {
            gen.ATE.sd.vec <- function(index) {
                if (treat.type == "discrete") {
                    x <- sub.data[index, X]
                    link.1 <- model.coef["(Intercept)"] + x * model.coef[X] + 1 * model.coef[paste0("D.", char)] + x * model.coef[X] * model.coef[paste0("DX.", char)]
                    link.0 <- model.coef["(Intercept)"] + x * model.coef[X]
                    if (is.null(Z) == FALSE) {
                        for (a in Z) {
                            target.Z <- sub.data[index, a]
                            link.1 <- link.1 + target.Z * model.coef[a]
                            link.0 <- link.0 + target.Z * model.coef[a]
                            if (full.moderate == TRUE) {
                                link.1 <- link.1 + target.Z * model.coef[paste0(a, ".X")] * x
                                link.0 <- link.0 + target.Z * model.coef[paste0(a, ".X")] * x
                            }
                        }
                    }
                    if (is.null(Z) == FALSE) {
                        vec.1 <- c(1, x, 1, x, as.matrix(sub.data[index, Z]))
                        vec.0 <- c(1, x, 0, 0, as.matrix(sub.data[index, Z]))
                        target.slice <- c("(Intercept)", X, paste0("D.", char), paste0("DX.", char), Z)
                        if (full.moderate == TRUE) {
                            vec.1 <- c(vec.1, as.matrix(sub.data[index, Z] * x))
                            vec.0 <- c(vec.0, as.matrix(sub.data[index, Z] * x))
                            target.slice <- c("(Intercept)", X, paste0("D.", char), paste0("DX.", char), Z, Z.X)
                        }
                    } else {
                        vec.1 <- c(1, x, 1, x)
                        vec.0 <- c(1, x, 0, 0)
                        target.slice <- c("(Intercept)", X, paste0("D.", char), paste0("DX.", char))
                    }
                    temp.vcov.matrix <- model.vcov[target.slice, target.slice]
                    if (method == "logit") {
                        vec <- vec.1 * exp(link.1) / (1 + exp(link.1))^2 - vec.0 * exp(link.0) / (1 + exp(link.0))^2
                    }
                    if (method == "probit") {
                        vec <- vec.1 * dnorm(link.1) - vec.0 * dnorm(link.0)
                    }
                    if (method == "poisson" | method == "nbinom") {
                        vec <- vec.1 * exp(link.1) - vec.0 * exp(link.0)
                    }
                    if (method == "linear") {
                        vec <- vec.1 - vec.0
                    }
                    return(vec)
                }

                if (treat.type == "continuous") {
                    link <- model.coef["(Intercept)"] + data[index, X] * model.coef[X] + model.coef[D] * data[index, D] + model.coef["DX"] * data[index, X] * data[index, D]
                    if (is.null(Z) == FALSE) {
                        for (a in Z) {
                            target.Z <- data[index, a]
                            link <- link + target.Z * model.coef[a]
                            if (full.moderate == TRUE) {
                                link <- link + target.Z * model.coef[paste0(a, ".X")] * data[index, X]
                            }
                        }
                    }

                    if (is.null(Z) == FALSE) {
                        vec1 <- c(1, data[index, X], data[index, D], data[index, D] * data[index, X], as.matrix(data[index, Z]))
                        vec0 <- c(0, 0, 1, data[index, X], rep(0, length(Z)))
                        target.slice <- c("(Intercept)", X, D, "DX", Z)
                        if (full.moderate == TRUE) {
                            vec1 <- c(vec1, as.matrix(data[index, Z] * data[index, X]))
                            vec0 <- c(vec0, rep(0, length(Z)))
                            target.slice <- c(target.slice, Z.X)
                        }
                    } else {
                        vec1 <- c(1, data[index, X], data[index, D], data[index, D] * data[index, X])
                        vec0 <- c(0, 0, 1, data[index, X])
                        target.slice <- c("(Intercept)", X, D, "DX")
                    }


                    temp.vcov.matrix <- model.vcov[target.slice, target.slice]

                    if (method == "logit") {
                        vec <- -(model.coef[D] + data[index, X] * model.coef["DX"]) * (exp(link) - exp(-link)) / (2 + exp(link) + exp(-link))^2 * vec1 + exp(link) / (1 + exp(link))^2 * vec0
                    }
                    if (method == "probit") {
                        vec <- dnorm(link) * vec0 - (model.coef[D] + data[index, X] * model.coef["DX"]) * link * dnorm(link) * vec1
                    }
                    if (method == "poisson" | method == "nbinom") {
                        vec <- (model.coef[D] + data[index, X] * model.coef["DX"]) * exp(link) * vec1 + exp(link) * vec0
                    }
                    if (method == "linear") {
                        vec <- vec0
                    }
                    return(vec)
                }
            }

            if (treat.type == "discrete") {
                index.all <- c(1:dim(sub.data)[1])
                vec.all <- sapply(index.all, function(x) gen.ATE.sd.vec(x))
                vec.mean <- apply(vec.all, 1, function(x) weighted.mean(x, weights))
                if (is.null(Z) == FALSE) {
                    target.slice <- c("(Intercept)", X, paste0("D.", char), paste0("DX.", char), Z)
                    if (full.moderate == TRUE) {
                        target.slice <- c(target.slice, Z.X)
                    }
                } else {
                    target.slice <- c("(Intercept)", X, paste0("D.", char), paste0("DX.", char))
                }
                temp.vcov.matrix <- model.vcov[target.slice, target.slice]

                ATE.sd <- sqrt((t(vec.mean) %*% temp.vcov.matrix %*% vec.mean)[1, 1])
                return(list(ATE = ATE, sd = ATE.sd))
            }

            if (treat.type == "continuous") {
                index.all <- c(1:dim(data)[1])
                vec.all <- sapply(index.all, function(x) gen.ATE.sd.vec(x))
                vec.mean <- apply(vec.all, 1, function(x) weighted.mean(x, weights))
                if (is.null(Z) == FALSE) {
                    target.slice <- c("(Intercept)", X, D, "DX", Z)
                    if (full.moderate == TRUE) {
                        target.slice <- c(target.slice, Z.X)
                    }
                } else {
                    target.slice <- c("(Intercept)", X, D, "DX")
                }

                temp.vcov.matrix <- model.vcov[target.slice, target.slice]
                AME.sd <- sqrt((t(vec.mean) %*% temp.vcov.matrix %*% vec.mean)[1, 1])
                return(list(AME = AME, sd = AME.sd))
            }
        }
    }


    # Part 1. Simu Vartype
    # generate TE/ME(x,mean,sd,0.025,0.975); matrix form; by group/D.ref
    # generate predict(x,mean,sd,0.025,0.975); matrix form; by group/D.ref
    # generate vcov of TE/ME; matrix form; by group/D.ref
    # generate diff.values of TE/ME(mean,sd,0.025,0.975); matrix form; by group/D.ref
    if (vartype == "simu") {
        # simu
        M <- nsimu
        simu.coef <- rmvnorm(M, model.coef, model.vcov)

        if (treat.type == "discrete") {
            TE.output.all.list <- list()
            pred.output.all.list <- list()
            link.output.all.list <- list()
            diff.output.all.list <- list()
            TE.vcov.list <- list()
            ATE.output.list <- list()
            for (char in other.treat) {
                gen.general.TE.output <- gen.general.TE(model.coef = model.coef, char = char, diff.values = diff.values)
                TE.output <- gen.general.TE.output$TE
                E.pred.output <- gen.general.TE.output$E.pred
                E.base.output <- gen.general.TE.output$E.base
                link.1.output <- gen.general.TE.output$link.1
                link.0.output <- gen.general.TE.output$link.0
                diff.estimate.output <- gen.general.TE.output$diff.estimate
                ATE.estimate <- gen.ATE(data = data, model.coef = model.coef, model.vcov = model.vcov, char = char)
                # simu
                one.simu.output <- function(model.coef, char = char, diff.values = diff.values) {
                    one.simu.TE.output <- gen.general.TE(model.coef = model.coef, char = char, diff.values = diff.values)
                    output1 <- one.simu.TE.output$TE
                    names(output1) <- rep("TE", length(output1))
                    output2 <- one.simu.TE.output$E.pred
                    names(output2) <- rep("pred", length(output2))
                    output3 <- one.simu.TE.output$E.base
                    names(output3) <- rep("base", length(output3))

                    output3a <- one.simu.TE.output$link.1
                    names(output3a) <- rep("link.1", length(output3a))
                    output3b <- one.simu.TE.output$link.0
                    names(output3b) <- rep("link.0", length(output3b))

                    output4 <- one.simu.TE.output$diff.estimate
                    names(output4) <- rep("diff", length(output4))
                    output5 <- gen.ATE(data = data, model.coef = model.coef, model.vcov = model.vcov, char = char)
                    output5 <- c(output5)
                    names(output5) <- c("ATE")
                    output.all <- c(output1, output2, output3, output3a, output3b, output4, output5)
                    return(output.all)
                }
                all.simu.matrix <- apply(simu.coef, 1, function(x) one.simu.output(model.coef = x, char = char, diff.values = diff.values))
                TE.simu.matrix <- all.simu.matrix[rownames(all.simu.matrix) == "TE", ]
                pred.simu.matrix <- all.simu.matrix[rownames(all.simu.matrix) == "pred", ]
                base.simu.matrix <- all.simu.matrix[rownames(all.simu.matrix) == "base", ]
                link.1.simu.matrix <- all.simu.matrix[rownames(all.simu.matrix) == "link.1", ]
                link.0.simu.matrix <- all.simu.matrix[rownames(all.simu.matrix) == "link.0", ]
                diff.simu.matrix <- all.simu.matrix[rownames(all.simu.matrix) == "diff", ]
                ATE.simu.matrix <- matrix(all.simu.matrix[rownames(all.simu.matrix) == "ATE", ], nrow = 1)

                if (length(diff.values) == 2) {
                    diff.simu.matrix <- matrix(diff.simu.matrix, nrow = 1)
                }
                if (length(diff.values) == 3) {
                    diff.simu.matrix <- as.matrix(diff.simu.matrix)
                }


                TE.simu.sd <- apply(TE.simu.matrix, 1, sd, na.rm = TRUE)
                pred.simu.sd <- apply(pred.simu.matrix, 1, sd, na.rm = TRUE)
                base.simu.sd <- apply(base.simu.matrix, 1, sd, na.rm = TRUE)
                link.1.simu.sd <- apply(link.1.simu.matrix, 1, sd, na.rm = TRUE)
                link.0.simu.sd <- apply(link.0.simu.matrix, 1, sd, na.rm = TRUE)
                diff.simu.sd <- apply(diff.simu.matrix, 1, sd, na.rm = TRUE)
                ATE.simu.sd <- apply(ATE.simu.matrix, 1, sd, na.rm = TRUE)

                TE.simu.CI <- t(apply(TE.simu.matrix, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
                pred.simu.CI <- t(apply(pred.simu.matrix, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
                base.simu.CI <- t(apply(base.simu.matrix, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
                link.1.simu.CI <- t(apply(link.1.simu.matrix, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
                link.0.simu.CI <- t(apply(link.0.simu.matrix, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
                diff.simu.CI <- t(apply(diff.simu.matrix, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
                ATE.simu.CI <- matrix(t(apply(ATE.simu.matrix, 1, quantile, c(0.025, 0.975), na.rm = TRUE)), nrow = 1)

                if (length(diff.values) == 2) {
                    diff.simu.CI <- matrix(diff.simu.CI, nrow = 1)
                }
                if (length(diff.values) == 3) {
                    diff.simu.CI <- as.matrix(diff.simu.CI)
                }

                TE.simu.vcov <- cov(t(TE.simu.matrix), use = "na.or.complete")
                rownames(TE.simu.vcov) <- NULL
                colnames(TE.simu.vcov) <- NULL

                TE.output.all <- cbind(X.eval, TE.output, TE.simu.sd, TE.simu.CI[, 1], TE.simu.CI[, 2])
                colnames(TE.output.all) <- c("X", "TE", "sd", "lower CI(95%)", "upper CI(95%)")
                rownames(TE.output.all) <- NULL
                TE.output.all.list[[other.treat.origin[char]]] <- TE.output.all

                pred.output.all <- cbind(X.eval, E.pred.output, pred.simu.sd, pred.simu.CI[, 1], pred.simu.CI[, 2])
                colnames(pred.output.all) <- c("X", "E(Y)", "sd", "lower CI(95%)", "upper CI(95%)")
                rownames(pred.output.all) <- NULL
                pred.output.all.list[[other.treat.origin[char]]] <- pred.output.all

                link.1.output.all <- cbind(X.eval, link.1.output, link.1.simu.sd, link.1.simu.CI[, 1], link.1.simu.CI[, 2])
                colnames(link.1.output.all) <- c("X", "E(Y)", "sd", "lower CI(95%)", "upper CI(95%)")
                rownames(link.1.output.all) <- NULL
                link.output.all.list[[other.treat.origin[char]]] <- link.1.output.all

                TE.vcov.list[[other.treat.origin[char]]] <- TE.simu.vcov

                z.value <- diff.estimate.output / diff.simu.sd
                p.value <- 2 * pnorm(-abs(z.value))
                diff.output.all <- cbind(diff.estimate.output, diff.simu.sd, z.value, p.value, diff.simu.CI[, 1], diff.simu.CI[, 2])
                colnames(diff.output.all) <- c("diff.estimate", "sd", "z-value", "p-value", "lower CI(95%)", "upper CI(95%)")
                rownames(diff.output.all) <- difference.name
                diff.output.all.list[[other.treat.origin[char]]] <- diff.output.all

                ATE.z.value <- ATE.estimate / ATE.simu.sd
                ATE.p.value <- 2 * pnorm(-abs(ATE.z.value))
                ATE.output <- c(ATE.estimate, ATE.simu.sd, ATE.z.value, ATE.p.value, ATE.simu.CI[, 1], ATE.simu.CI[, 2])
                names(ATE.output) <- c("ATE", "sd", "z-value", "p-value", "lower CI(95%)", "upper CI(95%)")
                ATE.output.list[[other.treat.origin[char]]] <- ATE.output
            }

            # base
            base.output.all <- cbind(X.eval, E.base.output, base.simu.sd, base.simu.CI[, 1], base.simu.CI[, 2])
            colnames(base.output.all) <- c("X", "E(Y)", "sd", "lower CI(95%)", "upper CI(95%)")
            rownames(base.output.all) <- NULL
            pred.output.all.list[[all.treat.origin[base]]] <- base.output.all

            link.0.output.all <- cbind(X.eval, link.0.output, link.0.simu.sd, link.0.simu.CI[, 1], link.0.simu.CI[, 2])
            colnames(link.0.output.all) <- c("X", "E(Y)", "sd", "lower CI(95%)", "upper CI(95%)")
            rownames(link.0.output.all) <- NULL
            link.output.all.list[[all.treat.origin[base]]] <- link.0.output.all
        }

        if (treat.type == "continuous") {
            ME.output.all.list <- list()
            pred.output.all.list <- list()
            link.output.all.list <- list()
            diff.output.all.list <- list()
            ME.vcov.list <- list()
            AME.output <- list()
            k <- 1
            for (D.ref in D.sample) {
                gen.general.ME.output <- gen.general.TE(model.coef = model.coef, D.ref = D.ref, diff.values = diff.values)
                ME.output <- gen.general.ME.output$ME
                E.pred.output <- gen.general.ME.output$E.pred
                link.output <- gen.general.ME.output$link
                diff.estimate.output <- gen.general.ME.output$diff.estimate
                # AME.estimate <- gen.ATE(data=data,model.coef=model.coef,model.vcov=model.vcov)
                # simu
                one.simu.output2 <- function(model.coef, D.ref = D.ref, diff.values = diff.values) {
                    one.simu.ME.output <- gen.general.TE(model.coef = model.coef, D.ref = D.ref, diff.values = diff.values)
                    output1 <- one.simu.ME.output$ME
                    names(output1) <- rep("ME", length(output1))
                    output2 <- one.simu.ME.output$E.pred
                    names(output2) <- rep("pred", length(output2))

                    output3 <- one.simu.ME.output$link
                    names(output3) <- rep("link", length(output3))

                    output4 <- one.simu.ME.output$diff.estimate
                    names(output4) <- rep("diff", length(output4))
                    # output5 <- gen.ATE(data=data,model.coef=model.coef,model.vcov=model.vcov)
                    # output5 <- c(output5)
                    # names(output5) <- c("AME")
                    # output.all <- c(output1,output2,output4,output5)
                    output.all <- c(output1, output2, output3, output4)
                    return(output.all)
                }
                all.simu.matrix <- apply(simu.coef, 1, function(x) one.simu.output2(model.coef = x, D.ref = D.ref, diff.values = diff.values))
                ME.simu.matrix <- all.simu.matrix[rownames(all.simu.matrix) == "ME", ]
                pred.simu.matrix <- all.simu.matrix[rownames(all.simu.matrix) == "pred", ]
                link.simu.matrix <- all.simu.matrix[rownames(all.simu.matrix) == "link", ]
                diff.simu.matrix <- all.simu.matrix[rownames(all.simu.matrix) == "diff", ]
                # AME.simu.matrix <- matrix(all.simu.matrix[rownames(all.simu.matrix)=="AME",],nrow=1)

                if (length(diff.values) == 2) {
                    diff.simu.matrix <- matrix(diff.simu.matrix, nrow = 1)
                }
                if (length(diff.values) == 3) {
                    diff.simu.matrix <- as.matrix(diff.simu.matrix)
                }


                ME.simu.sd <- apply(ME.simu.matrix, 1, sd)
                pred.simu.sd <- apply(pred.simu.matrix, 1, sd)
                link.simu.sd <- apply(link.simu.matrix, 1, sd)
                diff.simu.sd <- apply(diff.simu.matrix, 1, sd)
                # AME.simu.sd <- apply(AME.simu.matrix, 1, sd)

                ME.simu.CI <- t(apply(ME.simu.matrix, 1, quantile, c(0.025, 0.975)))
                pred.simu.CI <- t(apply(pred.simu.matrix, 1, quantile, c(0.025, 0.975)))
                link.simu.CI <- t(apply(link.simu.matrix, 1, quantile, c(0.025, 0.975)))
                diff.simu.CI <- t(apply(diff.simu.matrix, 1, quantile, c(0.025, 0.975)))
                # AME.simu.CI <- matrix(t(apply(AME.simu.matrix, 1, quantile, c(0.025,0.975))),nrow=1)

                if (length(diff.values) == 2) {
                    diff.simu.CI <- matrix(diff.simu.CI, nrow = 1)
                }
                if (length(diff.values) == 3) {
                    diff.simu.CI <- as.matrix(diff.simu.CI)
                }

                ME.simu.vcov <- cov(t(ME.simu.matrix), use = "na.or.complete")
                rownames(ME.simu.vcov) <- NULL
                colnames(ME.simu.vcov) <- NULL

                ME.output.all <- cbind(X.eval, ME.output, ME.simu.sd, ME.simu.CI[, 1], ME.simu.CI[, 2])
                colnames(ME.output.all) <- c("X", "ME", "sd", "lower CI(95%)", "upper CI(95%)")
                rownames(ME.output.all) <- NULL
                ME.output.all.list[[label.name[k]]] <- ME.output.all

                pred.output.all <- cbind(X.eval, E.pred.output, pred.simu.sd, pred.simu.CI[, 1], pred.simu.CI[, 2])
                colnames(pred.output.all) <- c("X", "E(Y)", "sd", "lower CI(95%)", "upper CI(95%)")
                rownames(pred.output.all) <- NULL
                pred.output.all.list[[label.name[k]]] <- pred.output.all

                link.output.all <- cbind(X.eval, link.output, link.simu.sd, link.simu.CI[, 1], link.simu.CI[, 2])
                colnames(link.output.all) <- c("X", "E(Y)", "sd", "lower CI(95%)", "upper CI(95%)")
                rownames(link.output.all) <- NULL
                link.output.all.list[[label.name[k]]] <- link.output.all

                ME.vcov.list[[label.name[k]]] <- ME.simu.vcov

                z.value <- diff.estimate.output / diff.simu.sd
                p.value <- 2 * pnorm(-abs(z.value))
                diff.output.all <- cbind(diff.estimate.output, diff.simu.sd, z.value, p.value, diff.simu.CI[, 1], diff.simu.CI[, 2])
                colnames(diff.output.all) <- c("diff.estimate", "sd", "z-value", "p-value", "lower CI(95%)", "upper CI(95%)")
                rownames(diff.output.all) <- difference.name
                diff.output.all.list[[label.name[k]]] <- diff.output.all

                # AME.z.value <- AME.estimate/AME.simu.sd
                # AME.p.value <- 2*pnorm(-abs(AME.z.value))
                # AME.output <- c(AME.estimate,AME.simu.sd,AME.z.value,AME.p.value,AME.simu.CI[,1],AME.simu.CI[,2])
                # names(AME.output) <- c("AME","sd","z-value","p-value","lower CI(95%)","upper CI(95%)")

                k <- k + 1
            }
            AME.estimate <- gen.ATE(data = data, model.coef = model.coef, model.vcov = model.vcov)
            all.simu.AME <- apply(simu.coef, 1, function(x) gen.ATE(model.coef = x, data = data, model.vcov = model.vcov))
            AME.simu.CI <- quantile(all.simu.AME, c(0.025, 0.975))
            AME.simu.sd <- sd(all.simu.AME)
            AME.z.value <- AME.estimate / AME.simu.sd
            AME.p.value <- 2 * pnorm(-abs(AME.z.value))
            AME.output <- c(AME.estimate, AME.simu.sd, AME.z.value, AME.p.value, AME.simu.CI[1], AME.simu.CI[2])
            names(AME.output) <- c("AME", "sd", "z-value", "p-value", "lower CI(95%)", "upper CI(95%)")
        }
    }


    # Part 2. Delta Vartype
    # generate TE/ME(x,mean,sd,0.025,0.975); matrix form; by group/D.ref
    # generate predict(x,mean,sd,0.025,0.975); matrix form; by group/D.ref
    # generate vcov of TE/ME; matrix form; by group/D.ref
    # generate diff.values of TE/ME(mean,sd,0.025,0.975); matrix form; by group/D.ref
    if (vartype == "delta") {
        crit <- abs(qt(.025, df = model.df))
        if (treat.type == "discrete") {
            TE.output.all.list <- list()
            pred.output.all.list <- list()
            link.output.all.list <- list()
            diff.output.all.list <- list()
            TE.vcov.list <- list()
            ATE.output.list <- list()
            for (char in other.treat) {
                gen.general.TE.output <- gen.general.TE(model.coef = model.coef, char = char, diff.values = diff.values)
                TE.output <- gen.general.TE.output$TE
                E.pred.output <- gen.general.TE.output$E.pred
                E.base.output <- gen.general.TE.output$E.base
                link.1.output <- gen.general.TE.output$link.1
                link.0.output <- gen.general.TE.output$link.0

                diff.estimate.output <- gen.general.TE.output$diff.estimate
                ATE.estimate.list <- gen.ATE(data = data, model.coef = model.coef, model.vcov = model.vcov, char = char, delta = TRUE)
                ATE.estimate <- ATE.estimate.list$ATE
                ATE.estimate.sd <- ATE.estimate.list$sd

                # delta
                gen.delta.TE.output <- gen.delta.TE(
                    model.coef = model.coef, model.vcov = model.vcov,
                    char = char, diff.values = diff.values, vcov = TRUE
                )
                TE.output.sd <- gen.delta.TE.output$TE.sd
                E.pred.sd <- gen.delta.TE.output$predict.sd
                link.1.sd <- gen.delta.TE.output$link.sd
                diff.estimate.sd <- gen.delta.TE.output$sd.diff.estimate
                TE.delta.vcov <- gen.delta.TE.output$TE.vcov
                colnames(TE.delta.vcov) <- NULL
                rownames(TE.delta.vcov) <- NULL

                # save
                TE.vcov.list[[other.treat.origin[char]]] <- TE.delta.vcov
                TE.uniform.q <- calculate_delta_uniformCI(TE.delta.vcov,alpha=0.05,N=2000)

                TE.output.all <- cbind(X.eval, TE.output, TE.output.sd, TE.output - crit * TE.output.sd, TE.output + crit * TE.output.sd, TE.output - TE.uniform.q * TE.output.sd, TE.output + TE.uniform.q * TE.output.sd)
                colnames(TE.output.all) <- c("X", "ME", "sd", "lower CI(95%)", "upper CI(95%)","lower uniform CI(95%)", "upper uniform CI(95%)")
                rownames(TE.output.all) <- NULL
                TE.output.all.list[[other.treat.origin[char]]] <- TE.output.all

                pred.output.all <- cbind(X.eval, E.pred.output, E.pred.sd, E.pred.output - crit * E.pred.sd, E.pred.output + crit * E.pred.sd)
                colnames(pred.output.all) <- c("X", "E(Y)", "sd", "lower CI(95%)", "upper CI(95%)")
                rownames(pred.output.all) <- NULL
                pred.output.all.list[[other.treat.origin[char]]] <- pred.output.all

                link.output.all <- cbind(X.eval, link.1.output, link.1.sd, link.1.output - crit * link.1.sd, link.1.output + crit * link.1.sd)
                colnames(link.output.all) <- c("X", "E(Y)", "sd", "lower CI(95%)", "upper CI(95%)")
                rownames(link.output.all) <- NULL
                link.output.all.list[[other.treat.origin[char]]] <- link.output.all

                z.value <- diff.estimate.output / diff.estimate.sd
                p.value <- 2 * pnorm(-abs(z.value))
                diff.output.all <- cbind(diff.estimate.output, diff.estimate.sd, z.value, p.value, diff.estimate.output - crit * diff.estimate.sd, diff.estimate.output + crit * diff.estimate.sd)
                colnames(diff.output.all) <- c("diff.estimate", "sd", "z-value", "p-value", "lower CI(95%)", "upper CI(95%)")
                rownames(diff.output.all) <- difference.name
                diff.output.all.list[[other.treat.origin[char]]] <- diff.output.all

                ATE.z.value <- ATE.estimate / ATE.estimate.sd
                ATE.p.value <- 2 * pnorm(-abs(ATE.z.value))
                ATE.output <- c(ATE.estimate, ATE.estimate.sd, ATE.z.value, ATE.p.value, ATE.estimate - crit * ATE.estimate.sd, ATE.estimate + crit * ATE.estimate.sd)
                names(ATE.output) <- c("ATE", "sd", "z-value", "p-value", "lower CI(95%)", "upper CI(95%)")
                ATE.output.list[[other.treat.origin[char]]] <- ATE.output
            }

            # base
            gen.delta.TE.base <- gen.delta.TE(
                model.coef = model.coef, model.vcov = model.vcov,
                char = base, diff.values = diff.values
            )


            E.pred.sd <- gen.delta.TE.base$predict.sd
            base.output.all <- cbind(X.eval, E.base.output, E.pred.sd, E.base.output - crit * E.pred.sd, E.base.output + crit * E.pred.sd)
            colnames(base.output.all) <- c("X", "E(Y)", "sd", "lower CI(95%)", "upper CI(95%)")
            rownames(base.output.all) <- NULL
            pred.output.all.list[[all.treat.origin[base]]] <- base.output.all

            link.0.sd <- gen.delta.TE.base$link.sd
            link.output.all <- cbind(X.eval, link.0.output, link.0.sd, link.0.output - crit * link.0.sd, link.0.output + crit * link.0.sd)
            colnames(link.output.all) <- c("X", "E(Y)", "sd", "lower CI(95%)", "upper CI(95%)")
            rownames(link.output.all) <- NULL
            link.output.all.list[[all.treat.origin[base]]] <- link.output.all
        }


        if (treat.type == "continuous") {
            ME.output.all.list <- list()
            pred.output.all.list <- list()
            link.output.all.list <- list()
            diff.output.all.list <- list()
            ME.vcov.list <- list()
            AME.output <- list()
            k <- 1
            for (D.ref in D.sample) {
                gen.general.ME.output <- gen.general.TE(model.coef = model.coef, D.ref = D.ref, diff.values = diff.values)
                ME.output <- gen.general.ME.output$ME
                E.pred.output <- gen.general.ME.output$E.pred
                link.output <- gen.general.ME.output$link
                diff.estimate.output <- gen.general.ME.output$diff.estimate

                # delta
                gen.delta.ME.output <- gen.delta.TE(
                    model.coef = model.coef, model.vcov = model.vcov,
                    D.ref = D.ref, diff.values = diff.values, vcov = TRUE
                )

                ME.output.sd <- gen.delta.ME.output$ME.sd
                E.pred.sd <- gen.delta.ME.output$predict.sd
                link.sd <- gen.delta.ME.output$link.sd
                diff.estimate.sd <- gen.delta.ME.output$sd.diff.estimate
                ME.delta.vcov <- gen.delta.ME.output$ME.vcov
                colnames(ME.delta.vcov) <- NULL
                rownames(ME.delta.vcov) <- NULL

                # save
                ME.vcov.list[[label.name[k]]] <- ME.delta.vcov
                ME.uniform.q <- calculate_delta_uniformCI(ME.delta.vcov,alpha=0.05,N=2000)
                

                ME.output.all <- cbind(X.eval, ME.output, ME.output.sd, ME.output - crit * ME.output.sd, ME.output + crit * ME.output.sd, ME.output - ME.uniform.q * ME.output.sd, ME.output + ME.uniform.q * ME.output.sd)
                colnames(ME.output.all) <- c("X", "ME", "sd", "lower CI(95%)", "upper CI(95%)","lower uniform CI(95%)", "upper uniform CI(95%)")
                rownames(ME.output.all) <- NULL
                ME.output.all.list[[label.name[k]]] <- ME.output.all

                pred.output.all <- cbind(X.eval, E.pred.output, E.pred.sd, E.pred.output - crit * E.pred.sd, E.pred.output + crit * E.pred.sd)
                colnames(pred.output.all) <- c("X", "E(Y)", "sd", "lower CI(95%)", "upper CI(95%)")
                rownames(pred.output.all) <- NULL
                pred.output.all.list[[label.name[k]]] <- pred.output.all

                link.output.all <- cbind(X.eval, link.output, link.sd, link.output - crit * link.sd, link.output + crit * link.sd)
                colnames(link.output.all) <- c("X", "E(Y)", "sd", "lower CI(95%)", "upper CI(95%)")
                rownames(link.output.all) <- NULL
                link.output.all.list[[label.name[k]]] <- link.output.all

                z.value <- diff.estimate.output / diff.estimate.sd
                p.value <- 2 * pnorm(-abs(z.value))
                diff.output.all <- cbind(diff.estimate.output, diff.estimate.sd, z.value, p.value, diff.estimate.output - crit * diff.estimate.sd, diff.estimate.output + crit * diff.estimate.sd)
                colnames(diff.output.all) <- c("diff.estimate", "sd", "z-value", "p-value", "lower CI(95%)", "upper CI(95%)")
                rownames(diff.output.all) <- difference.name
                diff.output.all.list[[label.name[k]]] <- diff.output.all
                k <- k + 1
            }

            AME.estimate.list <- gen.ATE(data = data, model.coef = model.coef, model.vcov = model.vcov, delta = TRUE)
            AME.estimate <- AME.estimate.list$AME
            AME.estimate.sd <- AME.estimate.list$sd
            AME.z.value <- AME.estimate / AME.estimate.sd
            AME.p.value <- 2 * pnorm(-abs(AME.z.value))
            AME.output <- c(AME.estimate, AME.estimate.sd, AME.z.value, AME.p.value, AME.estimate - crit * AME.estimate.sd, AME.estimate + crit * AME.estimate.sd)
            names(AME.output) <- c("AME", "sd", "z-value", "p-value", "lower CI(95%)", "upper CI(95%)")
        }
    }


    # Part 3. Bootstrap
    # generate TE/ME(x,mean,sd,0.025,0.975); matrix form; by group/D.ref
    # generate predict(x,mean,sd,0.025,0.975); matrix form; by group/D.ref
    # generate vcov of TE/ME; matrix form; by group/D.ref
    # generate diff.values of TE/ME(mean,sd,0.025,0.975); matrix form; by group/D.ref
    if (vartype == "bootstrap") {
        if (treat.type == "discrete") {
            all.length <- neval * length(other.treat) + # TE
                neval * length(all.treat) + # pred
                neval * length(all.treat) + # link
                length(other.treat) * length(difference.name) + # diff
                length(other.treat) # ATE
        }

        if (treat.type == "continuous") {
            all.length <- neval * length(label.name) + # ME
                neval * length(label.name) + # pred
                neval * length(label.name) + # link
                length(label.name) * length(difference.name) + 1 # diff&AME
        }

        one.boot <- function() {
            if (is.null(cl) == TRUE) {
                smp <- sample(1:n, n, replace = TRUE)
            } else { ## block bootstrap
                cluster.boot <- sample(clusters, length(clusters), replace = TRUE)
                smp <- unlist(id.list[match(cluster.boot, clusters)])
            }
            data.boot <- data[smp, ]

            boot.out <- matrix(NA, nrow = all.length, ncol = 0)

            ## check input...
            if (treat.type == "discrete") {
                if (length(unique(data.boot[, D])) != length(unique(data[, D]))) {
                    return(boot.out)
                }
            }
            if (is.null(IV)) {
                if (use_fe == 0) {
                    if (method == "linear") {
                        suppressWarnings(
                            model.boot <- glm(formula, data = data.boot, weights = WEIGHTS)
                        )
                    }
                    if (method == "logit") {
                        suppressWarnings(
                            model.boot <- glm(formula, data = data.boot, family = binomial(link = "logit"), weights = WEIGHTS)
                        )
                    }
                    if (method == "probit") {
                        suppressWarnings(
                            model.boot <- glm(formula, data = data.boot, family = binomial(link = "probit"), weights = WEIGHTS)
                        )
                    }
                    if (method == "poisson") {
                        suppressWarnings(
                            model.boot <- glm(formula, data = data.boot, family = poisson, weights = WEIGHTS)
                        )
                    }
                    if (method == "nbinom") {
                        suppressWarnings(
                            model.boot <- glm.nb(formula, data = data.boot, weights = WEIGHTS)
                        )
                    }
                    ## check converge...
                    if (model.boot$converged == FALSE) {
                        return(boot.out)
                    }
                }
                if (use_fe == TRUE) {
                    w.boot <- data.boot[, "WEIGHTS"]
                    model.boot <- felm(formula, data = data.boot, weights = w.boot)
                }
            }

            if (!is.null(IV)) {
                if (use_fe == 0) {
                    suppressWarnings(
                        model.boot <- ivreg(formula, data = data.boot, weights = WEIGHTS)
                    )
                    model.boot$converged <- TRUE
                }
                if (use_fe == 1) {
                    w.boot <- data.boot[, "WEIGHTS"]
                    suppressWarnings(
                        model.boot <- felm(formula, data = data.boot, weights = w.boot)
                    )
                    model.boot$converged <- TRUE
                }
            }

            coef.boot <- coef(model.boot)
            if (!is.null(IV)) {
                if (use_fe == 1) {
                    names(coef.boot) <- name.update
                }
            }
            boot.one.round <- c()

            if (treat.type == "discrete") {
                for (char in other.treat) {
                    gen.general.TE.output <- gen.general.TE(model.coef = coef.boot, char = char, diff.values = diff.values)
                    TE.output <- gen.general.TE.output$TE
                    names(TE.output) <- rep(paste0("TE.", char), neval)
                    E.pred.output <- gen.general.TE.output$E.pred
                    names(E.pred.output) <- rep(paste0("pred.", char), neval)
                    E.base.output <- gen.general.TE.output$E.base

                    link.1.output <- gen.general.TE.output$link.1
                    names(link.1.output) <- rep(paste0("link.", char), neval)
                    link.0.output <- gen.general.TE.output$link.0

                    diff.estimate.output <- c(gen.general.TE.output$diff.estimate)
                    names(diff.estimate.output) <- rep(paste0("diff.", char), length(difference.name))
                    ATE.estimate <- c(gen.ATE(data = data.boot, model.coef = coef.boot, model.vcov = NULL, char = char))
                    names(ATE.estimate) <- paste0("ATE.", char)
                    boot.one.round <- c(boot.one.round, TE.output, E.pred.output, link.1.output, diff.estimate.output, ATE.estimate)
                }
                names(E.base.output) <- rep(paste0("pred.", base), neval)
                names(link.0.output) <- rep(paste0("link.", base), neval)
                boot.one.round <- c(boot.one.round, E.base.output, link.0.output)
            }


            if (treat.type == "continuous") {
                k <- 1
                for (D.ref in D.sample) {
                    gen.general.ME.output <- gen.general.TE(model.coef = coef.boot, D.ref = D.ref, diff.values = diff.values)
                    ME.output <- gen.general.ME.output$ME
                    names(ME.output) <- rep(paste0("ME.", label.name[k]), neval)
                    E.pred.output <- gen.general.ME.output$E.pred
                    names(E.pred.output) <- rep(paste0("pred.", label.name[k]), neval)
                    link.output <- gen.general.ME.output$link
                    names(link.output) <- rep(paste0("link.", label.name[k]), neval)
                    diff.estimate.output <- c(gen.general.ME.output$diff.estimate)
                    names(diff.estimate.output) <- rep(paste0("diff.", label.name[k]), length(difference.name))
                    boot.one.round <- c(boot.one.round, ME.output, E.pred.output, link.output, diff.estimate.output)
                    k <- k + 1
                }
                AME.estimate <- c(gen.ATE(data = data.boot, model.coef = coef.boot, model.vcov = NULL))
                names(AME.estimate) <- c("AME")
                boot.one.round <- c(boot.one.round, AME.estimate)
            }

            boot.out <- cbind(boot.out, boot.one.round)
            rownames(boot.out) <- names(boot.one.round)
            return(boot.out)
        }

        if (parallel == TRUE) {
            requireNamespace("doParallel")
            ## require(iterators)
            maxcores <- detectCores()
            cores <- min(maxcores, cores)
            pcl <- future::makeClusterPSOCK(cores)
            doParallel::registerDoParallel(pcl)
            cat("Parallel computing with", cores, "cores...\n")

            suppressWarnings(
                bootout <- foreach(
                    i = 1:nboots, .combine = cbind,
                    .export = c("one.boot"), .packages = c("lfe", "AER"),
                    .inorder = FALSE
                ) %dopar% {
                    one.boot()
                }
            )
            suppressWarnings(stopCluster(pcl))
            cat("\r")
        } else {
            bootout <- matrix(NA, all.length, 0)
            for (i in 1:nboots) {
                tempdata <- one.boot()
                if (is.null(tempdata) == FALSE) {
                    bootout <- cbind(bootout, tempdata)
                }
                if (i %% 50 == 0) cat(i) else cat(".")
            }
            cat("\r")
        }

        if (treat.type == "discrete") {
            TE.output.all.list <- list()
            pred.output.all.list <- list()
            link.output.all.list <- list()
            diff.output.all.list <- list()
            TE.vcov.list <- list()
            ATE.output.list <- list()
            for (char in other.treat) {
                gen.general.TE.output <- gen.general.TE(model.coef = model.coef, char = char, diff.values = diff.values)
                TE.output <- gen.general.TE.output$TE
                E.pred.output <- gen.general.TE.output$E.pred
                E.base.output <- gen.general.TE.output$E.base
                link.1.output <- gen.general.TE.output$link.1
                link.0.output <- gen.general.TE.output$link.0

                diff.estimate.output <- gen.general.TE.output$diff.estimate
                ATE.estimate <- gen.ATE(data = data, model.coef = model.coef, model.vcov = model.vcov, char = char)

                TE.boot.matrix <- bootout[rownames(bootout) == paste0("TE.", char), ]
                pred.boot.matrix <- bootout[rownames(bootout) == paste0("pred.", char), ]
                base.boot.matrix <- bootout[rownames(bootout) == paste0("pred.", base), ]
                link.1.boot.matrix <- bootout[rownames(bootout) == paste0("link.", char), ]
                link.0.boot.matrix <- bootout[rownames(bootout) == paste0("link.", base), ]

                diff.boot.matrix <- bootout[rownames(bootout) == paste0("diff.", char), ]
                ATE.boot.matrix <- matrix(bootout[rownames(bootout) == paste0("ATE.", char), ], nrow = 1)

                if (length(diff.values) == 2) {
                    diff.boot.matrix <- matrix(diff.boot.matrix, nrow = 1)
                }
                if (length(diff.values) == 3) {
                    diff.boot.matrix <- as.matrix(diff.boot.matrix)
                }

                TE.boot.sd <- apply(TE.boot.matrix, 1, sd, na.rm = TRUE)
                pred.boot.sd <- apply(pred.boot.matrix, 1, sd, na.rm = TRUE)
                base.boot.sd <- apply(base.boot.matrix, 1, sd, na.rm = TRUE)
                link.1.boot.sd <- apply(link.1.boot.matrix, 1, sd, na.rm = TRUE)
                link.0.boot.sd <- apply(link.0.boot.matrix, 1, sd, na.rm = TRUE)
                diff.boot.sd <- apply(diff.boot.matrix, 1, sd, na.rm = TRUE)
                ATE.boot.sd <- apply(ATE.boot.matrix, 1, sd, na.rm = TRUE)

                TE.boot.CI <- t(apply(TE.boot.matrix, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
                pred.boot.CI <- t(apply(pred.boot.matrix, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
                base.boot.CI <- t(apply(base.boot.matrix, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
                link.1.boot.CI <- t(apply(link.1.boot.matrix, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
                link.0.boot.CI <- t(apply(link.0.boot.matrix, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
                diff.boot.CI <- t(apply(diff.boot.matrix, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
                ATE.boot.CI <- matrix(t(apply(ATE.boot.matrix, 1, quantile, c(0.025, 0.975), na.rm = TRUE)), nrow = 1)


                TE.boot.uniform.CI <- calculate_uniform_quantiles(TE.boot.matrix,0.05)
                uniform.coverage <- TE.boot.uniform.CI$coverage
                TE.boot.uniform.CI <- TE.boot.uniform.CI$Q_j

                pred.boot.uniform.CI <- calculate_uniform_quantiles(pred.boot.matrix,0.05)
                pred.boot.uniform.CI <- pred.boot.uniform.CI$Q_j

                link.1.boot.uniform.CI <- calculate_uniform_quantiles(link.1.boot.matrix,0.05)
                link.1.boot.uniform.CI <- link.1.boot.uniform.CI$Q_j

                if (length(diff.values) == 2) {
                    diff.boot.CI <- matrix(diff.boot.CI, nrow = 1)
                }
                if (length(diff.values) == 3) {
                    diff.boot.CI <- as.matrix(diff.boot.CI)
                }

                TE.boot.vcov <- cov(t(TE.boot.matrix), use = "na.or.complete")
                rownames(TE.boot.vcov) <- NULL
                colnames(TE.boot.vcov) <- NULL

                TE.output.all <- cbind(X.eval, TE.output, TE.boot.sd, TE.boot.CI[, 1], TE.boot.CI[, 2],TE.boot.uniform.CI[,1],TE.boot.uniform.CI[,2])
                colnames(TE.output.all) <- c("X", "TE", "sd", "lower CI(95%)", "upper CI(95%)","lower uniform CI(95%)", "upper uniform CI(95%)")
                rownames(TE.output.all) <- NULL
                TE.output.all.list[[other.treat.origin[char]]] <- TE.output.all

                pred.output.all <- cbind(X.eval, E.pred.output, pred.boot.sd, pred.boot.CI[, 1], pred.boot.CI[, 2], pred.boot.uniform.CI[,1], pred.boot.uniform.CI[,2])
                colnames(pred.output.all) <- c("X", "E(Y)", "sd", "lower CI(95%)", "upper CI(95%)","lower uniform CI(95%)", "upper uniform CI(95%)")
                rownames(pred.output.all) <- NULL
                pred.output.all.list[[other.treat.origin[char]]] <- pred.output.all

                link.output.all <- cbind(X.eval, link.1.output, link.1.boot.sd, link.1.boot.CI[, 1], link.1.boot.CI[, 2], link.1.boot.uniform.CI[,1], link.1.boot.uniform.CI[,2])
                colnames(link.output.all) <- c("X", "E(Y)", "sd", "lower CI(95%)", "upper CI(95%)","lower uniform CI(95%)", "upper uniform CI(95%)")
                rownames(link.output.all) <- NULL
                link.output.all.list[[other.treat.origin[char]]] <- link.output.all


                TE.vcov.list[[other.treat.origin[char]]] <- TE.boot.vcov

                z.value <- diff.estimate.output / diff.boot.sd
                p.value <- 2 * pnorm(-abs(z.value))
                diff.output.all <- cbind(
                    diff.estimate.output, diff.boot.sd,
                    z.value, p.value, diff.boot.CI[, 1], diff.boot.CI[, 2]
                )
                colnames(diff.output.all) <- c("diff.estimate", "sd", "z-value", "p-value", "lower CI(95%)", "upper CI(95%)")
                rownames(diff.output.all) <- difference.name
                diff.output.all.list[[other.treat.origin[char]]] <- diff.output.all

                ATE.z.value <- ATE.estimate / ATE.boot.sd
                ATE.p.value <- 2 * pnorm(-abs(ATE.z.value))
                ATE.output <- c(ATE.estimate, ATE.boot.sd, ATE.z.value, ATE.p.value, ATE.boot.CI[, 1], ATE.boot.CI[, 2])
                names(ATE.output) <- c("ATE", "sd", "z-value", "p-value", "lower CI(95%)", "upper CI(95%)")
                ATE.output.list[[other.treat.origin[char]]] <- ATE.output
            }
            # base
            base.boot.uniform.CI <- calculate_uniform_quantiles(base.boot.matrix,0.05)
            base.boot.uniform.CI <- base.boot.uniform.CI$Q_j

            link.0.boot.uniform.CI <- calculate_uniform_quantiles(link.0.boot.matrix,0.05)
            link.0.boot.uniform.CI <- link.0.boot.uniform.CI$Q_j

            # base
            base.output.all <- cbind(X.eval, E.base.output, base.boot.sd, base.boot.CI[, 1], base.boot.CI[, 2],base.boot.uniform.CI[,1],base.boot.uniform.CI[,2])
            colnames(base.output.all) <- c("X", "E(Y)", "sd", "lower CI(95%)", "upper CI(95%)","lower uniform CI(95%)", "upper uniform CI(95%)")
            rownames(base.output.all) <- NULL
            pred.output.all.list[[all.treat.origin[base]]] <- base.output.all

            link.0.output.all <- cbind(X.eval, link.0.output, link.0.boot.sd, link.0.boot.CI[, 1], link.0.boot.CI[, 2],link.0.boot.uniform.CI[,1], link.0.boot.uniform.CI[,2])
            colnames(link.0.output.all) <- c("X", "E(Y)", "sd", "lower CI(95%)", "upper CI(95%)","lower uniform CI(95%)", "upper uniform CI(95%)")
            rownames(link.0.output.all) <- NULL
            link.output.all.list[[all.treat.origin[base]]] <- link.0.output.all
        }

        if (treat.type == "continuous") {
            ME.output.all.list <- list()
            pred.output.all.list <- list()
            link.output.all.list <- list()
            diff.output.all.list <- list()
            ME.vcov.list <- list()
            AME.output <- list()
            k <- 1
            for (D.ref in D.sample) {
                label <- label.name[k]
                gen.general.ME.output <- gen.general.TE(model.coef = model.coef, D.ref = D.ref, diff.values = diff.values)
                ME.output <- gen.general.ME.output$ME
                E.pred.output <- gen.general.ME.output$E.pred
                link.output <- gen.general.ME.output$link
                diff.estimate.output <- gen.general.ME.output$diff.estimate

                ME.boot.matrix <- bootout[rownames(bootout) == paste0("ME.", label), ]
                pred.boot.matrix <- bootout[rownames(bootout) == paste0("pred.", label), ]
                link.boot.matrix <- bootout[rownames(bootout) == paste0("link.", label), ]
                diff.boot.matrix <- bootout[rownames(bootout) == paste0("diff.", label), ]
                if (length(diff.values) == 2) {
                    diff.boot.matrix <- matrix(diff.boot.matrix, nrow = 1)
                }
                if (length(diff.values) == 3) {
                    diff.boot.matrix <- as.matrix(diff.boot.matrix)
                }

                ME.boot.vcov <- cov(t(ME.boot.matrix), use = "na.or.complete")
                rownames(ME.boot.vcov) <- NULL
                colnames(ME.boot.vcov) <- NULL

                ME.boot.sd <- apply(ME.boot.matrix, 1, sd, na.rm = TRUE)
                pred.boot.sd <- apply(pred.boot.matrix, 1, sd, na.rm = TRUE)
                link.boot.sd <- apply(link.boot.matrix, 1, sd, na.rm = TRUE)
                diff.boot.sd <- apply(diff.boot.matrix, 1, sd, na.rm = TRUE)

                ME.boot.CI <- t(apply(ME.boot.matrix, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
                pred.boot.CI <- t(apply(pred.boot.matrix, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
                link.boot.CI <- t(apply(link.boot.matrix, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
                diff.boot.CI <- t(apply(diff.boot.matrix, 1, quantile, c(0.025, 0.975), na.rm = TRUE))


                ME.boot.uniform.CI <- calculate_uniform_quantiles(ME.boot.matrix,0.05)
                uniform.coverage <- ME.boot.uniform.CI$coverage
                ME.boot.uniform.CI <- ME.boot.uniform.CI$Q_j

                pred.boot.uniform.CI <- calculate_uniform_quantiles(pred.boot.matrix,0.05)
                pred.boot.uniform.CI <- pred.boot.uniform.CI$Q_j

                link.boot.uniform.CI <- calculate_uniform_quantiles(link.boot.matrix,0.05)
                link.boot.uniform.CI <- link.boot.uniform.CI$Q_j


                ME.output.all <- cbind(X.eval, ME.output, ME.boot.sd, ME.boot.CI[, 1], ME.boot.CI[, 2], ME.boot.uniform.CI[,1], ME.boot.uniform.CI[,2])
                colnames(ME.output.all) <- c("X", "ME", "sd", "lower CI(95%)", "upper CI(95%)","lower uniform CI(95%)", "upper uniform CI(95%)")
                rownames(ME.output.all) <- NULL
                ME.output.all.list[[label]] <- ME.output.all

                pred.output.all <- cbind(X.eval, E.pred.output, pred.boot.sd, pred.boot.CI[, 1], pred.boot.CI[, 2], pred.boot.uniform.CI[, 1], pred.boot.uniform.CI[, 2])
                colnames(pred.output.all) <- c("X", "E(Y)", "sd", "lower CI(95%)", "upper CI(95%)","lower uniform CI(95%)", "upper uniform CI(95%)")
                rownames(pred.output.all) <- NULL
                pred.output.all.list[[label]] <- pred.output.all

                link.output.all <- cbind(X.eval, link.output, link.boot.sd, link.boot.CI[, 1], link.boot.CI[, 2],link.boot.uniform.CI[,1],link.boot.uniform.CI[,2])
                colnames(link.output.all) <- c("X", "E(Y)", "sd", "lower CI(95%)", "upper CI(95%)","lower uniform CI(95%)", "upper uniform CI(95%)")
                rownames(link.output.all) <- NULL
                link.output.all.list[[label]] <- link.output.all

                ME.vcov.list[[label]] <- ME.boot.vcov

                z.value <- diff.estimate.output / diff.boot.sd
                p.value <- 2 * pnorm(-abs(z.value))
                diff.output.all <- cbind(
                    diff.estimate.output, diff.boot.sd,
                    z.value, p.value, diff.boot.CI[, 1], diff.boot.CI[, 2]
                )
                colnames(diff.output.all) <- c("diff.estimate", "sd", "z-value", "p-value", "lower CI(95%)", "upper CI(95%)")
                rownames(diff.output.all) <- difference.name
                diff.output.all.list[[label]] <- diff.output.all

                k <- k + 1
            }
            AME.estimate <- gen.ATE(data = data, model.coef = model.coef, model.vcov = model.vcov)
            AME.boot.matrix <- matrix(bootout[rownames(bootout) == "AME", ], nrow = 1)
            AME.boot.sd <- apply(AME.boot.matrix, 1, sd)
            AME.boot.CI <- matrix(t(apply(AME.boot.matrix, 1, quantile, c(0.025, 0.975), na.rm = TRUE)), nrow = 1)
            AME.z.value <- AME.estimate / AME.boot.sd
            AME.p.value <- 2 * pnorm(-abs(AME.z.value))
            AME.output <- c(AME.estimate, AME.boot.sd, AME.z.value, AME.p.value, AME.boot.CI[, 1], AME.boot.CI[, 2])
            names(AME.output) <- c("AME", "sd", "z-value", "p-value", "lower CI(95%)", "upper CI(95%)")
        }
    }


    if (TRUE) { # density or histogram
        if (treat.type == "discrete") { ## discrete D
            # density
            if (is.null(weights) == TRUE) {
                de <- density(data[, X])
            } else {
                suppressWarnings(de <- density(data[, X], weights = data[, "WEIGHTS"]))
            }

            treat_den <- list()
            for (char in all.treat) {
                de.name <- paste0("den.", char)
                if (is.null(weights) == TRUE) {
                    de.tr <- density(data[data[, D] == char, X])
                } else {
                    suppressWarnings(de.tr <- density(data[data[, D] == char, X],
                        weights = data[data[, D] == char, "WEIGHTS"]
                    ))
                }
                treat_den[[all.treat.origin[char]]] <- de.tr
            }

            # histogram
            if (is.null(weights) == TRUE) {
                hist.out <- hist(data[, X], breaks = 80, plot = FALSE)
            } else {
                suppressWarnings(hist.out <- hist(data[, X], data[, "WEIGHTS"],
                    breaks = 80, plot = FALSE
                ))
            }
            n.hist <- length(hist.out$mids)

            # count the number of treated
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

        if (treat.type == "continuous") { ## continuous D
            if (is.null(weights) == TRUE) {
                de <- density(data[, X])
            } else {
                suppressWarnings(de <- density(data[, X], weights = data[, "WEIGHTS"]))
            }
            if (is.null(weights) == TRUE) {
                hist.out <- hist(data[, X], breaks = 80, plot = FALSE)
            } else {
                suppressWarnings(hist.out <- hist(data[, X], data[, "WEIGHTS"],
                    breaks = 80, plot = FALSE
                ))
            }
            de.co <- de.tr <- NULL
        }
    }


    ## L kurtosis
    requireNamespace("Lmoments")
    Lkurtosis <- Lmoments(data[, X], returnobject = TRUE)$ratios[4]
    tests <- list(
        treat.type = treat.type,
        X.Lkurtosis = sprintf(Lkurtosis, fmt = "%#.3f")
    )

    ## Output

    if (treat.type == "discrete") {
        for (char in other.treat.origin) {
            new.diff.est <- as.data.frame(diff.output.all.list[[char]])
            for (i in 1:dim(new.diff.est)[1]) {
                new.diff.est[i, ] <- sprintf(new.diff.est[i, ], fmt = "%#.3f")
            }
            diff.output.all.list[[char]] <- new.diff.est

            outss <- data.frame(lapply(ATE.output.list[[char]], function(x) t(data.frame(x))))
            colnames(outss) <- c("ATE", "sd", "z-value", "p-value", "lower CI(95%)", "upper CI(95%)")
            rownames(outss) <- c("Average Treatment Effect")
            outss[1, ] <- sprintf(outss[1, ], fmt = "%#.3f")
            ATE.output.list[[char]] <- outss
        }


        final.output <- list(
            diff.info = diff.info,
            treat.info = treat.info,
            est.lin = TE.output.all.list,
            pred.lin = pred.output.all.list,
            link.lin = link.output.all.list,
            diff.estimate = diff.output.all.list,
            vcov.matrix = TE.vcov.list,
            Avg.estimate = ATE.output.list,
            Xlabel = Xlabel,
            Dlabel = Dlabel,
            Ylabel = Ylabel,
            de = de,
            de.tr = treat_den, # density
            hist.out = hist.out,
            count.tr = treat.hist,
            tests = tests,
            estimator = "linear",
            model.linear = model,
            use.fe = use_fe
        )
    }

    if (treat.type == "continuous") {
        for (label in label.name) {
            new.diff.est <- as.data.frame(diff.output.all.list[[label]])
            for (i in 1:dim(new.diff.est)[1]) {
                new.diff.est[i, ] <- sprintf(new.diff.est[i, ], fmt = "%#.3f")
            }
            diff.output.all.list[[label]] <- new.diff.est
        }

        outss <- data.frame(lapply(AME.output, function(x) t(data.frame(x))))
        colnames(outss) <- c("AME", "sd", "z-value", "p-value", "lower CI(95%)", "upper CI(95%)")
        rownames(outss) <- c("Average Marginal Effect")
        outss[1, ] <- sprintf(outss[1, ], fmt = "%#.3f")
        AME.output <- outss


        final.output <- list(
            diff.info = diff.info,
            treat.info = treat.info,
            est.lin = ME.output.all.list,
            pred.lin = pred.output.all.list,
            link.lin = link.output.all.list,
            diff.estimate = diff.output.all.list,
            vcov.matrix = ME.vcov.list,
            Avg.estimate = AME.output,
            Xlabel = Xlabel,
            Dlabel = Dlabel,
            Ylabel = Ylabel,
            de = de, # density
            de.tr = de.tr,
            hist.out = hist.out,
            count.tr = NULL,
            tests = tests,
            estimator = "linear",
            model.linear = model,
            use.fe = use_fe
        )
    }

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
