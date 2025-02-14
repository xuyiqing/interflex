interflex.kernel <- function(data,
                             Y, # outcome
                             D, # treatment indicator
                             X, # moderator
                             treat.info,
                             diff.info,
                             bw = NULL,
                             kfold = 10,
                             grid = 30,
                             metric = "MSE",
                             Z = NULL, # covariates
                             FE = NULL, # fixed effects
                             IV = NULL, # instrumental variables
                             weights = NULL, # weighting variable
                             full.moderate = FALSE, # fully moderated model
                             # Z.X = NULL, # fully moderated terms
                             neval = 50,
                             X.eval = NULL,
                             method = "linear", ## "probit"; "logit"; "poisson"; "nbinom"
                             CI = TRUE,
                             vartype = "bootstrap", ## "delta"; "bootstrap"
                             nboots = 200,
                             parallel = FALSE,
                             cores = 4,
                             cl = NULL, # variable to be clustered on
                             # predict = FALSE,
                             Z.ref = NULL, # same length as Z, set the value of Z when estimating marginal effects/predicted value
                             figure = TRUE,
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
    uniform.coverage <- NULL
    n <- dim(data)[1]

    binary <- FALSE
    if (length(unique(data[, Y])) == 2) {
        if ((0 %in% unique(data[, Y])) & (1 %in% unique(data[, Y]))) {
            binary <- TRUE
        }
    }

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

    # X.eval
    X.eval0 <- X.eval
    X.eval <- seq(min(data[, X]), max(data[, X]), length.out = neval)
    X.eval <- sort(c(X.eval, X.eval0))
    neval <- length(X.eval)
    # X.eval <- seq(min(data[,X]), max(data[,X]), length.out = neval)

    if (treat.type == "discrete") {
        for (char in other.treat) {
            data[, paste0("D", ".", char)] <- as.numeric(data[, D] == char)
        }
    }

    if (is.null(weights) == TRUE) {
        w <- rep(1, n)
    } else {
        w <- data[, weights]
    }

    data[, "WEIGHTS"] <- w

    # Xdensity
    suppressWarnings(Xdensity <- density(data[, X], weights = w))
    # fixed effects
    if (is.null(FE) == TRUE) {
        use_fe <- 0
    } else {
        use_fe <- 1
    }

    ## Function A: weighted glm
    # generate estimate of coef at a specific value(x) of the moderator
    # input:x, data, bw, weights, Xdensity
    # weights: vector

    wls.iv.fe <- function(x, data, bw, weights, Xdensity) {
        data.touse <- data
        xx <- data.touse[, "delta.x"] <- data.touse[, X] - x
        use.variable <- c(Y)
        n.coef <- 1

        # construct endogenous variables
        endogenous.var <- c()
        if (treat.type == "discrete") {
            for (char in other.treat) {
                data.touse[, paste0("D.delta.x", ".", char)] <- data.touse[, paste0("D", ".", char)] * data.touse[, "delta.x"]
                use.variable <- c(use.variable, c(paste0("D", ".", char), paste0("D.delta.x", ".", char)))
                endogenous.var <- c(endogenous.var, c(paste0("D", ".", char), paste0("D.delta.x", ".", char)))
                n.coef <- n.coef + 2
            }
        }
        if (treat.type == "continuous") {
            data.touse[, "D.delta.x"] <- data.touse[, D] * data.touse[, "delta.x"]
            use.variable <- c(use.variable, c(D, "D.delta.x"))
            endogenous.var <- c(endogenous.var, c(D, "D.delta.x"))
            n.coef <- n.coef + 2
        }

        # construct weight
        temp_density <- Xdensity$y[which.min(abs(Xdensity$x - x))]
        density.mean <- exp(mean(log(Xdensity$y)))
        bw.adapt <- bw * (1 + log(max(Xdensity$y) / temp_density))
        # bw.adapt <- bw*sqrt(density.mean/temp_density)
        w <- dnorm(data.touse[, "delta.x"] / bw.adapt) * weights
        data.touse[, "WEIGHTS"] <- w

        if (max(data.touse[, "WEIGHTS"]) == 0) {
            result <- rep(NA, 1 + n.coef)
            return(list(
                result = result, model.vcov = NULL,
                model.df = NULL, data.touse = NULL
            ))
        }

        formula <- paste0(use.variable[1], "~", paste0(Z, collapse = "+"), "|", paste0(FE, collapse = "+"), "|", paste0(endogenous.var, collapse = "+"), "~", paste0(excluded.iv, collapse = "+"))
        fe_res <- feols(fml = as.formula(formula), data = data.touse, weights = w, vcov = "hetero")

        if (typeof(fe_res) != "list") {
            result <- rep(NA, 1 + n.coef)
            names(result) <- c("x0", "(Intercept)", endogenous.var, Z)
            return(list(
                result = result, model.vcov = NULL,
                model.df = NULL, data.touse = NULL
            ))
        }
        result <- c(x, mean(fe_res$sumFE), coef(fe_res))
        result[which(is.nan(result))] <- 0

        model.vcov.original <- vcov(fe_res, vcov = "hetero")
        model.vcov <- cbind(0, rbind(0, model.vcov.original))
        rownames(model.vcov) <- c("(Intercept)", rownames(model.vcov.original))
        colnames(model.vcov) <- c("(Intercept)", colnames(model.vcov.original))

        return(list(
            result = result, model.vcov = model.vcov,
            model.df = degrees_freedom(fe_res, type = "k"), data.touse = data.touse
        ))
    }

    wls.fe <- function(x, data, bw, weights, Xdensity) {
        data.touse <- data
        data.touse[, "delta.x"] <- data.touse[, X] - x
        use.variable <- c(Y, "delta.x")
        n.coef <- 1

        if (treat.type == "discrete") {
            for (char in other.treat) {
                data.touse[, paste0("D.delta.x", ".", char)] <- data.touse[, paste0("D", ".", char)] * data.touse[, "delta.x"]
                use.variable <- c(use.variable, c(paste0("D", ".", char), paste0("D.delta.x", ".", char)))
                n.coef <- n.coef + 2
            }
        }

        if (treat.type == "continuous") {
            data.touse[, "D.delta.x"] <- data.touse[, D] * data.touse[, "delta.x"]
            use.variable <- c(use.variable, c(D, "D.delta.x"))
            n.coef <- n.coef + 2
        }

        if (is.null(Z) == FALSE) {
            use.variable <- c(use.variable, Z)
            n.coef <- n.coef + length(Z)
            if (full.moderate == TRUE) {
                for (a in Z) {
                    data.touse[, paste0(a, ".delta.x")] <- data.touse[, a] * data.touse[, "delta.x"]
                    use.variable <- c(use.variable, paste0(a, ".delta.x"))
                    n.coef <- n.coef + 1
                }
            }
        }

        temp_density <- Xdensity$y[which.min(abs(Xdensity$x - x))]
        bw.adapt <- bw * (1 + log(max(Xdensity$y) / temp_density))
        w <- dnorm(data.touse[, "delta.x"] / bw.adapt) * weights
        data.touse[, "WEIGHTS"] <- w

        if (max(data.touse[, "WEIGHTS"]) == 0) {
            result <- rep(NA, 1 + n.coef)
            return(list(result = result, model.vcov = NULL, model.df = NULL, data.touse = NULL))
        }

        formula <- paste0(use.variable[1], "~", paste0(use.variable[2:length(use.variable)], collapse = "+"), "|", paste0(FE, collapse = "+"))
        fe_res <- feols(fml = as.formula(formula), data = data.touse, weights = w, vcov = "hetero")

        if (typeof(fe_res) != "list") {
            result <- rep(NA, 1 + n.coef)
            names(result) <- c("x0", "(Intercept)", use.variable[2:length(use.variable)])
            return(list(result = result, model.vcov = NULL, model.df = NULL, data.touse = NULL))
        }

        result <- c(x, mean(fe_res$sumFE), coef(fe_res))
        result[which(is.nan(result))] <- 0
        names(result) <- c("x0", "(Intercept)", use.variable[2:length(use.variable)])

        model.vcov.original <- vcov(fe_res, vcov = "hetero")
        model.vcov <- cbind(0, rbind(0, model.vcov.original))
        rownames(model.vcov) <- c("(Intercept)", rownames(model.vcov.original))
        colnames(model.vcov) <- c("(Intercept)", colnames(model.vcov.original))

        return(list(
            result = result, model.vcov = model.vcov,
            model.df = degrees_freedom(fe_res, type = "k"), data.touse = data.touse
        ))
    }

    wls.iv <- function(x, data, bw, weights, Xdensity) {
        data.touse <- data
        data.touse[, "delta.x"] <- data.touse[, X] - x
        formula <- paste0(Y, "~delta.x")
        all.var.name <- c("(Intercept)", "delta.x")
        n.coef <- 2
        if (treat.type == "discrete") {
            for (char in other.treat) {
                data.touse[, paste0("D.delta.x", ".", char)] <- data.touse[, paste0("D", ".", char)] * data.touse[, "delta.x"]
                formula <- paste0(formula, "+", paste0("D", ".", char), "+", paste0("D.delta.x", ".", char))
                all.var.name <- c(all.var.name, paste0("D", ".", char), paste0("D.delta.x", ".", char))
                n.coef <- n.coef + 2
            }
        }
        if (treat.type == "continuous") {
            data.touse[, "D.delta.x"] <- data.touse[, D] * data.touse[, "delta.x"]
            formula <- paste0(formula, "+", D, "+", "D.delta.x")
            all.var.name <- c(all.var.name, D, "D.delta.x")
            n.coef <- n.coef + 2
        }
        if (is.null(Z) == FALSE) {
            formula <- paste0(formula, "+", paste0(Z, collapse = "+"))
            all.var.name <- c(all.var.name, Z)
            n.coef <- n.coef + length(Z)
            if (full.moderate == TRUE) {
                for (a in Z) {
                    data.touse[, paste0(a, ".delta.x")] <- data.touse[, a] * data.touse[, "delta.x"]
                    formula <- paste0(formula, "+", paste0(a, ".delta.x"))
                    all.var.name <- c(all.var.name, paste0(a, ".delta.x"))
                    n.coef <- n.coef + 1
                }
            }
        }
        formula.iv <- "delta.x"
        for (sub.iv in IV) {
            data.touse[, paste0("delta.x.", sub.iv)] <- data.touse[, sub.iv] * data.touse[, "delta.x"]
            formula.iv <- paste0(formula.iv, "+", sub.iv)
            formula.iv <- paste0(formula.iv, "+", paste0("delta.x.", sub.iv))
        }
        # if(is.null(Z)==FALSE){
        # 	formula.iv <- paste0(formula.iv,"+",paste0(Z,collapse="+"))
        # 	if(full.moderate==TRUE){
        # 		formula.iv <- paste0(formula.iv,"+",paste0(Z.X,collapse="+"))
        # 	}
        # }
        if (is.null(Z) == FALSE) {
            formula.iv <- paste0(formula.iv, "+", paste0(Z, collapse = "+"))
            if (full.moderate == TRUE) {
                for (a in Z) {
                    formula.iv <- paste0(formula.iv, "+", paste0(a, ".delta.x"))
                }
            }
        }
        formula <- paste0(formula, "|", formula.iv)
        formula <- as.formula(formula)
        temp_density <- Xdensity$y[which.min(abs(Xdensity$x - x))]
        density.mean <- exp(mean(log(Xdensity$y)))
        bw.adapt <- bw * (1 + log(max(Xdensity$y) / temp_density))
        # bw.adapt <- bw*sqrt(density.mean/temp_density)
        w <- dnorm(data.touse[, "delta.x"] / bw.adapt) * weights
        data.touse[, "WEIGHTS"] <- w
        if (max(data.touse[, "WEIGHTS"]) == 0) {
            result <- rep(NA, 1 + n.coef)
            names(result) <- c("x0", all.var.name)
            return(list(
                result = result, model.vcov = NULL,
                model.df = NULL, data.touse = NULL
            ))
        }
        suppressWarnings( # correct
            iv.reg <- tryCatch(
                ivreg(formula, data = data.touse, weights = WEIGHTS),
                error = function(e) {
                    return("error")
                }
            )
        )
        if (typeof(iv.reg) != "list") {
            result <- rep(NA, 1 + n.coef)
            names(result) <- c("x0", all.var.name)
            return(list(
                result = result, model.vcov = NULL,
                model.df = NULL, data.touse = NULL
            ))
        }

        result <- c(x, iv.reg$coef)
        names(result) <- c("x0", names(iv.reg$coef))
        result[which(is.na(result))] <- 0
        return(result)
        return(list(
            result = result, model.vcov = vcov(iv.reg, type = "H2"),
            model.df = iv.reg$df.residual, data.touse = data.touse
        ))
    }

    wls.nofe <- function(x, data, bw, weights, Xdensity) {
        data.touse <- data
        data.touse[, "delta.x"] <- data.touse[, X] - x
        formula <- paste0(Y, "~delta.x")
        all.var.name <- c("(Intercept)", "delta.x")
        n.coef <- 2
        if (treat.type == "discrete") {
            for (char in other.treat) {
                data.touse[, paste0("D.delta.x", ".", char)] <- data.touse[, paste0("D", ".", char)] * data.touse[, "delta.x"]
                formula <- paste0(formula, "+", paste0("D", ".", char), "+", paste0("D.delta.x", ".", char))
                all.var.name <- c(all.var.name, paste0("D", ".", char), paste0("D.delta.x", ".", char))
                n.coef <- n.coef + 2
            }
        }
        if (treat.type == "continuous") {
            data.touse[, "D.delta.x"] <- data.touse[, D] * data.touse[, "delta.x"]
            formula <- paste0(formula, "+", D, "+", "D.delta.x")
            all.var.name <- c(all.var.name, D, "D.delta.x")
            n.coef <- n.coef + 2
        }
        if (is.null(Z) == FALSE) {
            formula <- paste0(formula, "+", paste0(Z, collapse = "+"))
            all.var.name <- c(all.var.name, Z)
            n.coef <- n.coef + length(Z)
            if (full.moderate == TRUE) {
                for (a in Z) {
                    data.touse[, paste0(a, ".delta.x")] <- data.touse[, a] * data.touse[, "delta.x"]
                    formula <- paste0(formula, "+", paste0(a, ".delta.x"))
                    all.var.name <- c(all.var.name, paste0(a, ".delta.x"))
                    n.coef <- n.coef + 1
                }
            }
        }

        formula <- as.formula(formula)
        
        temp_density <- Xdensity$y[which.min(abs(Xdensity$x - x))]
        bw.adapt <- bw * (1 + log(max(Xdensity$y) / temp_density))
        w <- dnorm(data.touse[, "delta.x"] / bw.adapt) * weights
        if (0 %in% w) {
            w <- w + min(w[w != 0])
        }
        # density.mean <- exp(mean(log(Xdensity$y)))
        # bw.adapt <- bw * sqrt(density.mean/temp_density)

        data.touse[, "WEIGHTS"] <- w
        if (max(data.touse[, "WEIGHTS"]) == 0) {
            result <- rep(NA, 1 + n.coef)
            names(result) <- c("x0", all.var.name)
            return(result)
        }

        if (method == "linear") {
            suppressWarnings( # correct
                glm.reg <- tryCatch(
                    glm(formula, data = data.touse, weights = WEIGHTS),
                    error = function(e) {
                        return("error")
                    }
                )
            )
        }

        if (method == "logit") {
            suppressWarnings( # correct
                glm.reg <- tryCatch(
                    glm(formula, data = data.touse, weights = WEIGHTS, family = binomial(link = "logit")),
                    error = function(e) {
                        return("error")
                    }
                )
            )
        }

        if (method == "probit") {
            suppressWarnings( # correct
                glm.reg <- tryCatch(
                    glm(formula, data = data.touse, weights = WEIGHTS, family = binomial(link = "probit")),
                    error = function(e) {
                        return("error")
                    }
                )
            )
        }

        if (method == "poisson") {
            suppressWarnings( # correct
                glm.reg <- tryCatch(
                    glm(formula, data = data.touse, weights = WEIGHTS, family = poisson),
                    error = function(e) {
                        return("error")
                    }
                )
            )
        }

        if (method == "nbinom") {
            suppressWarnings( # correct
                glm.reg <- tryCatch(
                    glm.nb(formula, data = data.touse, weights = WEIGHTS, control = glm.control(epsilon = 1e-5, maxit = 200)),
                    error = function(e) {
                        return("error")
                    }
                )
            )
        }

        if (typeof(glm.reg) != "list") {
            result <- rep(NA, 1 + n.coef)
            names(result) <- c(
                "x0",
                names(glm.reg$coef)
            )
            return(list(result = result, model.vcov = NULL, model.df = NULL, data.touse = NULL))
        }

        glm.reg.summary <- summary(glm.reg, robust = "HC2")
        glm.reg.vcov <- vcovHC(glm.reg, type = "HC2")
        glm.reg.df <- glm.reg$df.residual
        if (glm.reg$converged == FALSE) {
            result <- rep(NA, 1 + length(glm.reg$coef))
            names(result) <- c(
                "x0",
                names(glm.reg$coef)
            )
            return(list(result = result, model.vcov = NULL, model.df = NULL, data.touse = NULL))
        }

        if (glm.reg$converged == TRUE) {
            result <- c(
                x,
                glm.reg.summary$coef[, 1]
            )
            names(result) <- c(
                "x0",
                names(glm.reg$coef)
            )
            result[which(is.na(result))] <- 0
            return(list(result = result, model.vcov = glm.reg.vcov, model.df = glm.reg.df, data.touse = data.touse))
        }
    }

    wls <- function(x, data, bw, weights, Xdensity) {
        if (use_fe == TRUE & is.null(IV)) {
            return(wls.fe(x = x, data = data, bw = bw, weights = weights, Xdensity = Xdensity))
        }
        if (use_fe == FALSE & is.null(IV)) {
            return(wls.nofe(x = x, data = data, bw = bw, weights = weights, Xdensity = Xdensity))
        }
        if (use_fe == FALSE & !is.null(IV)) {
            return(wls.iv(x = x, data = data, bw = bw, weights = weights, Xdensity = Xdensity))
        }
        if (use_fe == TRUE & !is.null(IV)) {
            return(wls.iv.fe(x = x, data = data, bw = bw, weights = weights, Xdensity = Xdensity))
        }
    }

    # Function # Cross-Validation
    if (is.null(bw) == TRUE) {
        CV <- TRUE
        cat("Cross-validating bandwidth ... \n")
        if (length(grid) == 1) {
            rangeX <- max(data[, X]) - min(data[, X])
            bw.grid <- exp(seq(log(rangeX / 50), log(rangeX), length.out = grid))
        } else {
            bw.grid <- grid
        }
    } else {
        CV <- FALSE
    }
    if (CV == TRUE) {
        # weights: name of a variable
        getError.CV <- function(train, test, bw, neval, weights_name, Xdensity) {
            if (is.null(weights_name) == TRUE) {
                w.touse.cv <- rep(1, dim(train)[1])
            } else {
                w.touse.cv <- train[, weights_name]
            }

            if (use_fe == TRUE) {
                fe_res <- feols(fml = as.formula(paste0(Y, " ~ 1 | ", paste(FE, collapse = "+"))), data = train, weights = w.touse.cv)
                FEvalues <- fixef(fe_res)
                FEnumbers <- fe_res$fixef_sizes
                FE_coef <- matrix(0, nrow = sum(FEnumbers), ncol = 1)
                rowname <- c()
                fe_index_name <- c()
                for (fe in FE) {
                    for (i in 1:FEnumbers[[fe]]) {
                        FE_coef[i, 1] <- FEvalues[[fe]][i]
                        rowname <- c(rowname, paste0(fe, ".", names(FEvalues[[fe]])[i]))
                        fe_index_name <- c(fe_index_name, fe)
                    }
                }
                rownames(FE_coef) <- rowname
                train[, Y] <- fe_res$residuals
            }

            X.eval.cv <- seq(min(train[, X]), max(train[, X]), length.out = neval)

            # demean -> use wls without fixed effects#
            if (is.null(IV)) {
                coef.grid.cv <- c()
                for (x in X.eval) {
                    coef.grid.cv <- rbind(coef.grid.cv, wls.nofe(x = x, data = train, bw = bw, weights = w.touse.cv, Xdensity = Xdensity)$result)
                }
            } else {
                coef.grid.cv <- t(sapply(X.eval.cv, function(x) wls.iv(x = x, data = train, bw = bw, weights = w.touse.cv, Xdensity = Xdensity)))
            }
            coef.grid.cv <- na.omit(coef.grid.cv)
            eff.eval.point <- dim(coef.grid.cv)[1]
            X.eval.cv <- coef.grid.cv[, "x0"]
            if (dim(coef.grid.cv)[1] <= neval / 2) {
                output <- rep(NA, 5)
                names(output) <- c("Num.Eff.Points", "Cross Entropy", "AUC", "MSE", "MAE")
                return(output)
            }

            esCoef.cv <- function(x) { ## obtain the coefficients for x[i]
                Xnew.cv <- abs(X.eval.cv - x)
                d1 <- min(Xnew.cv)
                label1 <- which.min(Xnew.cv)
                Xnew.cv[label1] <- Inf
                d2 <- min(Xnew.cv) ## distance between x[i] and the second nearest x in training set
                label2 <- which.min(Xnew.cv)
                if (d1 == 0) {
                    func <- coef.grid.cv[label1, ]
                } else if (d2 == 0) {
                    func <- coef.grid.cv[label2, ]
                } else { ## weighted average
                    func1 <- coef.grid.cv[label1, ]
                    func2 <- coef.grid.cv[label2, ]
                    func <- (func1 * d2 + func2 * d1) / (d1 + d2)
                }
                return(func)
            }

            # Test
            ## limit test set in the support of the train dataset
            test <- test[which(test[, X] >= min(train[, X]) & test[, X] <= max(train[, X])), ]

            add_FE <- rep(0, dim(test)[1])
            if (use_fe == TRUE) {
                # addictive FE
                add_FE <- matrix(0, nrow = dim(test)[1], ncol = length(FE))
                colnames(add_FE) <- FE
                for (fe in FE) {
                    add_FE[, fe] <- 0
                    fe_name <- paste0(fe, ".", test[, fe])
                    find_FE_index <- which(fe_name %in% rownames(FE_coef))
                    not_find_FE_index <- which(!fe_name %in% rownames(FE_coef))
                    add_FE[find_FE_index, fe] <- FE_coef[fe_name[find_FE_index], ]
                    add_FE[not_find_FE_index, fe] <- mean(FE_coef[which(fe_index_name == fe), ])
                }
                add_FE <- rowSums(add_FE)
                add_FE <- add_FE + mean(fe_res$sumFE)
            }

            if (dim(test)[1] < 3) {
                output <- rep(NA, 5)
                names(output) <- c("Num.Eff.Points", "Cross Entropy", "AUC", "MSE", "MAE")
                return(output)
            }

            Knn <- t(sapply(test[, X], esCoef.cv))
            link <- Knn[, "(Intercept)"]

            if (treat.type == "discrete") {
                for (char in other.treat) {
                    link <- link + Knn[, paste0("D.", char)] * as.numeric(test[, D] == char) + Knn[, paste0("D.delta.x", ".", char)] * (test[, X] - Knn[, "x0"]) * as.numeric(test[, D] == char)
                }
                link <- link + Knn[, "delta.x"] * (test[, X] - Knn[, "x0"])
            }

            if (treat.type == "continuous") {
                link <- link + Knn[, D] * test[, D] + Knn[, "delta.x"] * (test[, X] - Knn[, "x0"]) + Knn[, "D.delta.x"] * (test[, X] - Knn[, "x0"]) * test[, D]
            }

            if (is.null(Z) == FALSE) {
                for (a in Z) {
                    link <- link + test[, a] * Knn[, a]
                    if (full.moderate == TRUE) {
                        for (a in Z) {
                            link <- link + Knn[, paste0(a, ".delta.x")] * (test[, X] - Knn[, "x0"]) * test[, a]
                        }
                    }
                }
            }

            if (method == "logit") {
                E.pred <- exp(link) / (1 + exp(link))
            }

            if (method == "probit") {
                E.pred <- pnorm(link, 0, 1)
            }

            if (method == "linear") {
                E.pred <- link

                if (use_fe == TRUE) {
                    E.pred <- E.pred + add_FE
                }

                if (binary == TRUE) {
                    E.pred[which(E.pred > 1)] <- 1
                    E.pred[which(E.pred < 0)] <- 0
                }
            }

            if (method == "poisson" | method == "nbinom") {
                E.pred <- exp(link)
            }

            true <- test[which(is.na(E.pred) == FALSE), Y]
            E.pred <- na.omit(E.pred)


            if (length(unique(true)) <= 1 | length(E.pred) < 3) { # all 0 or all 1 or few observations(auc)
                output <- rep(NA, 5)
                names(output) <- c("Num.Eff.Points", "Cross Entropy", "AUC", "MSE", "MAE")
                return(output)
            }

            # cross entropy
            if (binary == TRUE) {
                cross.entropy <- logLoss(true, E.pred)
                suppressMessages(
                    roc.y <- roc(true, E.pred, levels = c(0, 1))
                )
                # auc
                auc <- roc.y$auc
            } else {
                cross.entropy <- NA
                auc <- NA
            }

            # mse L2
            mse <- mean((true - E.pred)^2)
            # mse L1
            mae <- mean(abs(true - E.pred))
            output <- c(cross.entropy, auc, mse, mae)

            if (method == "poisson" | method == "nbinom") {
                mse.link <- mean((log(true + 1) - log(E.pred + 1))^2)
                mae.link <- mean(abs(log(true + 1) - log(E.pred + 1)))
                output <- c(mse.link, mae.link, mse, mae)
            }

            output <- c(eff.eval.point, output)
            names(output) <- c("Num.Eff.Points", "Cross Entropy", "AUC", "MSE", "MAE")
            return(output)
        }
        fold <- createFolds(factor(data[, D]), k = kfold, list = FALSE)
        # kfold <- min(n,kfold)
        # cat("#folds =",kfold)
        # cat("\n")
        # fold <- c(0:(n-1))%%kfold + 1
        # fold <- sample(fold, n, replace = FALSE)

        cv.new <- function(bw, neval) {
            error <- matrix(NA, kfold, 5)
            for (j in 1:kfold) { # K-fold CV
                testid <- which(fold == j)
                train <- data[-testid, ]
                test <- data[testid, ]
                error[j, ] <- getError.CV(train = train, test = test, bw = bw, neval = neval, weights_name = "WEIGHTS", Xdensity = Xdensity)
            }

            colnames(error) <- c("Num.Eff.Points", "cross entropy", "auc", "L2", "L1")
            for (pp in colnames(error)) {
                if (any(is.na(error[, pp])) == TRUE & all(is.na(error[, pp])) == FALSE) {
                    if (pp != "auc") {
                        error[is.na(error[, pp]), pp] <- max(error[, pp], na.rm = T)
                    } else {
                        error[is.na(error[, pp]), pp] <- min(error[, pp], na.rm = T)
                    }
                }
            }

            output <- c(bw, apply(error, 2, mean, na.rm = TRUE))
            names(output) <- c("Num.Eff.Points", "bw", "cross entropy", "auc", "L2", "L1")
            return(output)
        }

        ## test
        # try <- cv.new(0.02789474,neval=neval)
        # return(try)
        ##

        ## -----------------------------------------------------------------------------------
        if (parallel == TRUE) {
            requireNamespace("doParallel")
            ## require(iterators)
            maxcores <- detectCores()
            cores <- min(maxcores, cores)
            pcl <- future::makeClusterPSOCK(cores)
            doParallel::registerDoParallel(pcl)
            cat("Parallel computing with", cores, "cores...\n")
            Error <- suppressWarnings(foreach(
                bw = bw.grid, .combine = rbind,
                .packages = c("ModelMetrics", "pROC", "MASS", "AER"),
                .inorder = FALSE
            ) %dopar% {
                cv.output.sub <- try(cv.new(bw, neval = neval),silent = TRUE)
                if('try-error' %in% class(cv.output.sub)){
                    return(NA)
                }else{
                    return(cv.output.sub)
                }
            })

            suppressWarnings(stopCluster(pcl))
            # return(Error)
        } else {
            Error <- matrix(NA, length(bw.grid), 6)
            for (i in 1:length(bw.grid)) {
                suppressWarnings(cv.output.sub <- try(cv.new(bw = bw.grid[i], neval = neval), silent = FALSE))
                if('try-error' %in% class(cv.output.sub)){
                    Error[i, ] <- NA
                }else{
                    Error[i, ] <- cv.output.sub
                }
                cat(".")
            }
        }

        colnames(Error) <- c("bw", "Num.Eff.Points", "Cross Entropy", "AUC", "MSE", "MAE")
        rownames(Error) <- NULL

        if (binary == FALSE) {
            if (method == "poisson" | method == "nbinom") {
                colnames(Error) <- c("bw", "Num.Eff.Points", "MSE.link", "MAE.link", "MSE", "MAE")
            } else {
                Error <- Error[, c("bw", "Num.Eff.Points", "MSE", "MAE")]
            }
            if (metric == "MSE") {
                bw <- bw.grid[which.min(Error[, "MSE"] / Error[, "Num.Eff.Points"])]
            }
            if (metric == "MAE") {
                bw <- bw.grid[which.min(Error[, "MAE"] / Error[, "Num.Eff.Points"])]
            }
        }
        if (binary == TRUE) {
            if (metric == "MSE") {
                bw <- bw.grid[which.min(Error[, "MSE"] / Error[, "Num.Eff.Points"])]
            }
            if (metric == "MAE") {
                bw <- bw.grid[which.min(Error[, "MAE"] / Error[, "Num.Eff.Points"])]
            }
            if (metric == "Cross Entropy") {
                bw <- bw.grid[which.min(Error[, "Cross Entropy"] / Error[, "Num.Eff.Points"])]
            }
            if (metric == "AUC") {
                bw <- bw.grid[which.max(Error[, "AUC"] * Error[, "Num.Eff.Points"])]
            }
        }
        cat(paste0("Optimal bw=", round(bw, 4), ".\n"))
    } else {
        Error <- NULL
    }
    # Core Estimation, gen grid points

    count <- 1
    results <- list()
    coef.grid <- c()
    model.vcovs <- list()
    model.dfs <- c()
    for (x in X.eval) {
        results[[count]] <- wls(x = x, data = data, bw = bw, weights = w, Xdensity = Xdensity)
        coef.grid <- rbind(coef.grid, results[[count]]$result)
        model.vcovs[[count]] <- results[[count]]$model.vcov
        model.dfs <- c(model.dfs, results[[count]]$model.df)
        count <- count + 1
    }

    coef.grid <- na.omit(coef.grid)
    if (dim(coef.grid)[1] <= neval / 2) {
        warning("Inappropriate bandwidth.\n")
    }
    if (dim(coef.grid)[1] <= 3) {
        stop("Inappropriate bandwidth.")
    }
    X.eval <- coef.grid[, "x0"]
    neval <- length(X.eval)

    cat(paste0("Number of evaluation points:", neval, "\n"))

    gen.sd <- function(result, char = NULL, D.ref=NULL, to.diff = FALSE) {
        coef.grid <- result$result
        x_prev <- coef.grid["x0"]
        model.vcov <- result$model.vcov
        data.touse <- result$data.touse
        x <- data.touse[which.min(abs(data.touse[[X]] - x_prev)), "delta.x"]
        if (treat.type == "discrete") {
            link.1 <- coef.grid["(Intercept)"] + x * coef.grid["delta.x"] + 1 * coef.grid[paste0("D.", char)] + x * coef.grid[paste0("D.delta.x.", char)]
            link.0 <- coef.grid["(Intercept)"] + x * coef.grid["delta.x"]
            if (is.null(Z) == FALSE) {
                for (a in Z) {
                    target.Z <- Z.ref[a]
                    link.1 <- link.1 + target.Z * coef.grid[a]
                    link.0 <- link.0 + target.Z * coef.grid[a]
                }
            }

            if (is.null(Z) == FALSE) {
                if (full.moderate == FALSE) {
                    vec.1 <- c(1, x, 1, x, Z.ref)
                    vec.0 <- c(1, x, 0, 0, Z.ref)
                    target.slice <- c("(Intercept)", "delta.x", paste0("D.", char), paste0("D.delta.x.", char), Z)
                } else {
                    vec.1 <- c(1, x, 1, x, Z.ref, x * Z.ref)
                    vec.0 <- c(1, x, 0, 0, Z.ref, x * Z.ref)
                    target.slice <- c("(Intercept)", "delta.x", paste0("D.", char), paste0("D.delta.x.", char), Z, paste0(Z, ".delta.x"))
                }
            } else {
                vec.1 <- c(1, x, 1, x)
                vec.0 <- c(1, x, 0, 0)
                target.slice <- c("(Intercept)", "delta.x", paste0("D.", char), paste0("D.delta.x.", char))
            }
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
            temp.vcov.matrix <- model.vcov[target.slice, target.slice]
            if (to.diff == TRUE) {
                return(list(vec = vec, temp.vcov.matrix = temp.vcov.matrix))
            }
            delta.sd <- sqrt((t(vec) %*% temp.vcov.matrix %*% vec)[1, 1])
            return(delta.sd)
        }

        if (treat.type == "continuous") {
            link <- coef.grid["(Intercept)"] + x * coef.grid[X] + coef.grid[D] * D.ref + coef.grid["D.delta.x"] * x * D.ref
            if (is.null(Z) == FALSE) {
                for (a in Z) {
                    target.Z <- Z.ref[a]
                    link <- link + target.Z * coef.grid[a]
                }
            }

            if (is.null(Z) == FALSE) {
                if (full.moderate == FALSE) {
                    vec1 <- c(1, x, D.ref, D.ref * x, Z.ref)
                    vec0 <- c(0, 0, 1, x, rep(0, length(Z)))
                    target.slice <- c("(Intercept)", "delta.x", D, "D.delta.x", Z)
                } else {
                    vec1 <- c(1, x, D.ref, D.ref * x, Z.ref, x * Z.ref)
                    vec0 <- c(0, 0, 1, x, rep(0, length(Z)), rep(0, length(Z)))
                    target.slice <- c("(Intercept)", "delta.x", D, "D.delta.x", Z, paste0(Z, ".delta.x"))
                }
            } else {
                vec1 <- c(1, x, D.ref, D.ref * x)
                vec0 <- c(0, 0, 1, x)
                target.slice <- c("(Intercept)", "delta.x", D, "D.delta.x")
            }

            temp.vcov.matrix <- model.vcov[target.slice, target.slice]

            if (method == "logit") {
                vec <- -(coef.grid[D] + x * coef.grid["D.delta.x"]) * (exp(link) - exp(-link)) / (2 + exp(link) + exp(-link))^2 * vec1 + exp(link) / (1 + exp(link))^2 * vec0
            }
            if (method == "probit") {
                vec <- dnorm(link) * vec0 - (coef.grid[D] + x * coef.grid["D.delta.x"]) * link * dnorm(link) * vec1
            }
            if (method == "poisson" | method == "nbinom") {
                vec <- (coef.grid[D] + x * coef.grid["D.delta.x"]) * exp(link) * vec1 + exp(link) * vec0
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

    ## Function B: estimate TE/ME; E.predict(E.base);
    # 1, estimate treatment effects/marginal effects given model.coef;
    # 2, estimate E.pred given treat/D;
    # 3, input: coef.grid; char(discrete)/D.ref(continuous);
    # 4, output: marginal effects/treatment effects/E.pred/E.base
    gen.kernel.TE <- function(coef.grid, char = NULL, D.ref = NULL, base.flag=FALSE) {
        if (is.null(char) == TRUE) {
            treat.type <- "continuous"
        }
        if (is.null(D.ref) == TRUE) {
            treat.type == "discrete"
        }
        neval <- dim(coef.grid)[1]

        gen.link.sd <- function(result, base = FALSE) {
            coef.grid <- result$result
            x_prev <- coef.grid["x0"]
            model.vcov <- result$model.vcov
            data.touse <- result$data.touse
            x <- data.touse[which.min(abs(data.touse[[X]] - x_prev)), "delta.x"]
            if (treat.type == "discrete") {
                if (is.null(Z) == FALSE) {
                    if (base == FALSE) {
                        vec <- c(1, x, 1, x, Z.ref)
                        target.slice <- c("(Intercept)", "delta.x", paste0("D.", char), paste0("D.delta.x.", char), Z)
                    } else {
                        vec <- c(1, x, Z.ref)
                        target.slice <- c("(Intercept)", "delta.x", Z)
                    }
                } else {
                    if (base == FALSE) {
                        vec <- c(1, x, 1, x)
                        target.slice <- c("(Intercept)", "delta.x", paste0("D.", char), paste0("D.delta.x.", char))
                    } else {
                        vec <- c(1, x)
                        target.slice <- c("(Intercept)", "delta.x")
                    }
                }

                temp.vcov.matrix <- model.vcov[target.slice, target.slice]
                link.sd <- sqrt((t(vec) %*% temp.vcov.matrix %*% vec)[1, 1])
                return(link.sd)
            }

            if (treat.type == "continuous") {
                if (is.null(Z) == FALSE) {
                    vec <- c(1, x, D.ref, D.ref * x, Z.ref)
                    target.slice <- c("(Intercept)", "delta.x", D, "D.delta.x", Z)
                    if (full.moderate == TRUE) {
                        vec <- c(vec, Z.ref * x)
                        target.slice <- c(target.slice, paste0(Z, ".delta.x"))
                    }
                } else {
                    vec <- c(1, x, D.ref, D.ref * x)
                    target.slice <- c("(Intercept)", "delta.x", D, "D.delta.x")
                }

                temp.vcov.matrix <- model.vcov[target.slice, target.slice]
                link.sd <- sqrt((t(vec) %*% temp.vcov.matrix %*% vec)[1, 1])
                return(link.sd)
            }
        }

        gen.predict.sd <- function(result, base = FALSE) {
            coef.grid <- result$result
            x_prev <- coef.grid["x0"]
            model.vcov <- result$model.vcov
            data.touse <- result$data.touse
            x <- data.touse[which.min(abs(data.touse[[X]] - x_prev)), "delta.x"]
            if (treat.type == "discrete") {
                if (base == FALSE) {
                    link <- coef.grid["(Intercept)"] + x * coef.grid[X] + 1 * coef.grid[paste0("D.", char)] + x * coef.grid[paste0("D.delta.x.", char)]
                }
                if (base == TRUE) {
                    link <- coef.grid["(Intercept)"] + x * coef.grid[X]
                }
                if (is.null(Z) == FALSE) {
                    for (a in Z) {
                        target.Z <- Z.ref[a]
                        link <- link + target.Z * coef.grid[a]
                        if (full.moderate == TRUE) {
                            link <- link + target.Z * coef.grid[paste0(a, ".X")] * x
                        }
                    }
                }

                if (is.null(Z) == FALSE) {
                    if (base == FALSE) {
                        vec <- c(1, x, 1, x, Z.ref)
                        target.slice <- c("(Intercept)", "delta.x", paste0("D.", char), paste0("D.delta.x.", char), Z)
                    } else {
                        vec <- c(1, x, Z.ref)
                        target.slice <- c("(Intercept)", "delta.x", Z)
                    }
                } else {
                    if (base == FALSE) {
                        vec <- c(1, x, 1, x)
                        target.slice <- c("(Intercept)", "delta.x", paste0("D.", char), paste0("D.delta.x.", char))
                    } else {
                        vec <- c(1, x)
                        target.slice <- c("(Intercept)", "delta.x")
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
                link <- coef.grid["(Intercept)"] + x * coef.grid[X] + coef.grid[D] * D.ref + coef.grid["D.delta.x"] * x * D.ref
                if (is.null(Z) == FALSE) {
                    for (a in Z) {
                        target.Z <- Z.ref[a]
                        link <- link + target.Z * coef.grid[a]
                        if (full.moderate == TRUE) {
                            link <- link + target.Z * coef.grid[paste0(a, ".X")] * x
                        }
                    }
                }

                if (is.null(Z) == FALSE) {
                    vec <- c(1, x, D.ref, D.ref * x, Z.ref)
                    target.slice <- c("(Intercept)", "delta.x", D, "D.delta.x", Z)
                } else {
                    vec <- c(1, x, D.ref, D.ref * x)
                    target.slice <- c("(Intercept)", "delta.x", D, "D.delta.x")
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

        gen.TE <- function(coef.grid) {
            if (treat.type == "discrete") {
                link.1 <- coef.grid[, "(Intercept)"] + coef.grid[, paste0("D.", char)]
                link.0 <- coef.grid[, "(Intercept)"] + 0
                if (is.null(Z) == FALSE) {
                    for (a in Z) {
                        target.Z <- Z.ref[a]
                        link.1 <- link.1 + target.Z * coef.grid[, a]
                        link.0 <- link.0 + target.Z * coef.grid[, a]
                        # if(full.moderate==TRUE){
                        # 	link.1 <- link.1 + target.Z*coef.grid[,paste0(a,".X")]*coef.grid[,'x0']
                        # 	link.0 <- link.0 + target.Z*coef.grid[,paste0(a,".X")]*coef.grid[,'x0']
                        # }
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
                names(TE) <- rep(paste0("TE.", char), neval)
                names(E.pred) <- rep(paste0("Predict.", char), neval)
                names(E.base) <- rep(paste0("Predict.", base), neval)
                gen.TE.output <- list(TE = TE, E.pred = E.pred, E.base = E.base, link.1 = link.1, link.0 = link.0)
            }

            if (treat.type == "continuous") {
                link <- coef.grid[, "(Intercept)"] + D.ref * coef.grid[, D]
                if (is.null(Z) == FALSE) {
                    for (a in Z) {
                        target.Z <- Z.ref[a]
                        link <- link + target.Z * coef.grid[, a]
                        # if(full.moderate==TRUE){
                        # 	link <- link + target.Z*coef.grid[,paste0(a,".X")]*coef.grid[,'x0']
                        # }
                    }
                }
                if (method == "logit") {
                    ME <- exp(link) / (1 + exp(link))^2 * coef.grid[, D]
                    E.pred <- exp(link) / (1 + exp(link))
                }
                if (method == "probit") {
                    ME <- coef.grid[, D] * dnorm(link)
                    E.pred <- pnorm(link, 0, 1)
                }
                if (method == "linear") {
                    ME <- coef.grid[, D]
                    E.pred <- link
                }
                if (method == "poisson" | method == "nbinom") {
                    ME <- exp(link) * coef.grid[, D]
                    E.pred <- exp(link)
                }
                names(ME) <- rep(paste0("ME.", names(D.sample)[D.sample == D.ref]), neval)
                names(E.pred) <- rep(paste0("Predict.", names(D.sample)[D.sample == D.ref]), neval)
                gen.TE.output <- list(ME = ME, E.pred = E.pred, link = link)
            }
            return(gen.TE.output)
        }

        if (base.flag == FALSE) {
            TE.sd <- c(mapply(function(x) gen.sd(x, char = char, D.ref = D.ref), x = results))
        } else {
            TE.sd <- NULL
        }
        
        if (treat.type == "discrete") {
            if (char == base) {
                link.sd <- c(sapply(results, function(x) gen.link.sd(x, base = TRUE)))
                predict.sd <- c(sapply(results, function(x) gen.predict.sd(x, base = TRUE)))
            } else {
                link.sd <- c(sapply(results, function(x) gen.link.sd(x)))
                predict.sd <- c(sapply(results, function(x) gen.predict.sd(x)))
            }
            names(predict.sd) <- rep(paste0("predict.sd.", char), neval)
        } else {
            link.sd <- c(sapply(results, function(x) gen.link.sd(x)))
            predict.sd <- c(sapply(results, function(x) gen.predict.sd(x)))
            names(predict.sd) <- rep(paste0("predict.sd.", names(D.sample)[D.sample == D.ref]), neval)
        }
        if (treat.type == "discrete") {
            if (char == base) {
                return(list(
                    TE.sd = TE.sd,
                    predict.sd = predict.sd,
                    link.sd = link.sd,
                    model.df = model.dfs
                ))
            }
        }

        gen.TE.output <- gen.TE(coef.grid = coef.grid)
        if (treat.type == "discrete") {
            return(list(
                TE = gen.TE.output$TE,
                E.pred = gen.TE.output$E.pred,
                E.base = gen.TE.output$E.base,
                link.1 = gen.TE.output$link.1,
                link.0 = gen.TE.output$link.0,
                TE.sd = TE.sd,
                predict.sd = predict.sd,
                link.sd = link.sd,
                model.df = model.dfs
            ))
        }

        if (treat.type == "continuous") {
            return(list(
                ME = gen.TE.output$ME,
                E.pred = gen.TE.output$E.pred,
                link = gen.TE.output$link,
                ME.sd = TE.sd,
                predict.sd = predict.sd,
                link.sd = link.sd,
                model.df = model.dfs
            ))
        }
    }

    ## Function C: estimate difference of TE/ME at different values of the moderator
    # 1,	input: coef.grid; char/D.ref; diff.values
    # 2,	output: difference of TE/ME at different values of the moderator
    gen.kernel.difference <- function(coef.grid, diff.values, char = NULL, D.ref = NULL) {
        if (is.null(char) == TRUE) {
            treat.type <- "continuous"

            est.ME <- function(x) {
                Xnew <- abs(X.eval - x)
                d1 <- min(Xnew)
                label1 <- which.min(Xnew)
                Xnew[label1] <- Inf
                d2 <- min(Xnew)
                label2 <- which.min(Xnew)

                if (d1 == 0) {
                    link <- coef.grid[label1, "(Intercept)"] + coef.grid[label1, D] * D.ref
                    if (is.null(Z) == FALSE) {
                        for (a in Z) {
                            target.Z <- Z.ref[a]
                            link <- link + target.Z * coef.grid[label1, a]
                            # if(full.moderate==TRUE){
                            # 	link <- link + target.Z*coef.grid[label1,paste0(a,".X")]*x
                            # }
                        }
                    }
                    coef.grid.D <- coef.grid[label1, D]
                } else if (d2 == 0) {
                    link <- coef.grid[label2, "(Intercept)"] + coef.grid[label2, D] * D.ref
                    if (is.null(Z) == FALSE) {
                        for (a in Z) {
                            target.Z <- Z.ref[a]
                            link <- link + target.Z * coef.grid[label2, a]
                            # if(full.moderate==TRUE){
                            # 	link <- link + target.Z*coef.grid[label2,paste0(a,".X")]*x
                            # }
                        }
                    }
                    coef.grid.D <- coef.grid[label2, D]
                } else { ## weighted average
                    link.1 <- coef.grid[label1, "(Intercept)"] + coef.grid[label1, D] * D.ref +
                        coef.grid[label1, "D.delta.x"] * D.ref * (x - X.eval[label1]) +
                        coef.grid[label1, "delta.x"] * (x - X.eval[label1])

                    link.2 <- coef.grid[label2, "(Intercept)"] + coef.grid[label2, D] * D.ref +
                        coef.grid[label2, "D.delta.x"] * D.ref * (x - X.eval[label2]) +
                        coef.grid[label2, "delta.x"] * (x - X.eval[label2])

                    if (is.null(Z) == FALSE) {
                        for (a in Z) {
                            target.Z <- Z.ref[a]
                            link.1 <- link.1 + target.Z * coef.grid[label1, a]
                            link.2 <- link.2 + target.Z * coef.grid[label2, a]
                            if (full.moderate == TRUE) {
                                link.1 <- link.1 + target.Z * coef.grid[label1, paste0(a, ".delta.x")] * (x - X.eval[label1])
                                link.2 <- link.2 + target.Z * coef.grid[label2, paste0(a, ".delta.x")] * (x - X.eval[label2])
                            }
                        }
                    }

                    coef.grid.D.1 <- coef.grid[label1, D] + coef.grid[label1, "D.delta.x"] * (x - X.eval[label1])
                    coef.grid.D.2 <- coef.grid[label2, D] + coef.grid[label2, "D.delta.x"] * (x - X.eval[label2])
                    coef.grid.D <- (coef.grid.D.1 * d2 + coef.grid.D.2 * d1) / (d1 + d2)
                    link <- (link.1 * d2 + link.2 * d1) / (d1 + d2)
                }

                if (method == "logit") {
                    ME <- exp(link) / (1 + exp(link))^2 * coef.grid.D
                }
                if (method == "probit") {
                    ME <- coef.grid.D * dnorm(link)
                }
                if (method == "linear") {
                    ME <- coef.grid.D
                }
                if (method == "poisson" | method == "nbinom") {
                    ME <- exp(link) * coef.grid.D
                }

                return(ME)
            }
        }
        if (is.null(D.ref) == TRUE) {
            treat.type == "discrete"

            est.TE <- function(x) { ## estimate the treatment effects for an observation
                Xnew <- abs(X.eval - x)
                d1 <- min(Xnew)
                label1 <- which.min(Xnew)
                Xnew[label1] <- Inf
                d2 <- min(Xnew)
                label2 <- which.min(Xnew)
                if (d1 == 0) {
                    link.1 <- coef.grid[label1, "(Intercept)"] + coef.grid[label1, paste0("D.", char)]
                    link.0 <- coef.grid[label1, "(Intercept)"] + 0
                    if (is.null(Z) == FALSE) {
                        for (a in Z) {
                            target.Z <- Z.ref[a]
                            link.1 <- link.1 + target.Z * coef.grid[label1, a]
                            link.0 <- link.0 + target.Z * coef.grid[label1, a]
                            # if(full.moderate==TRUE){
                            # 	link.1 <- link.1 + target.Z*coef.grid[label1,paste0(a,".X")]*x
                            # 	link.0 <- link.0 + target.Z*coef.grid[label1,paste0(a,".X")]*x
                            # }
                        }
                    }
                } else if (d2 == 0) {
                    link.1 <- coef.grid[label2, "(Intercept)"] + coef.grid[label2, paste0("D.", char)]
                    link.0 <- coef.grid[label2, "(Intercept)"] + 0
                    if (is.null(Z) == FALSE) {
                        for (a in Z) {
                            target.Z <- Z.ref[a]
                            link.1 <- link.1 + target.Z * coef.grid[label2, a]
                            link.0 <- link.0 + target.Z * coef.grid[label2, a]
                            # if(full.moderate==TRUE){
                            # 	link.1 <- link.1 + target.Z*coef.grid[label2,paste0(a,".X")]*x
                            # 	link.0 <- link.0 + target.Z*coef.grid[label2,paste0(a,".X")]*x
                            # }
                        }
                    }
                } else { ## weighted average
                    link.1.1 <- coef.grid[label1, "(Intercept)"] + coef.grid[label1, paste0("D.", char)] +
                        (coef.grid[label1, paste0("D.delta.x.", char)] + coef.grid[label1, "delta.x"]) * (x - X.eval[label1])
                    link.1.2 <- coef.grid[label2, "(Intercept)"] + coef.grid[label2, paste0("D.", char)] +
                        (coef.grid[label2, paste0("D.delta.x.", char)] + coef.grid[label2, "delta.x"]) * (x - X.eval[label2])
                    link.0.1 <- coef.grid[label1, "(Intercept)"] + coef.grid[label1, "delta.x"] * (x - X.eval[label1])
                    link.0.2 <- coef.grid[label2, "(Intercept)"] + coef.grid[label2, "delta.x"] * (x - X.eval[label2])
                    if (is.null(Z) == FALSE) {
                        for (a in Z) {
                            target.Z <- Z.ref[a]
                            link.1.1 <- link.1.1 + target.Z * coef.grid[label1, a]
                            link.1.2 <- link.1.2 + target.Z * coef.grid[label2, a]
                            link.0.1 <- link.0.1 + target.Z * coef.grid[label1, a]
                            link.0.2 <- link.0.2 + target.Z * coef.grid[label2, a]
                            if (full.moderate == TRUE) {
                                link.1.1 <- link.1.1 + target.Z * coef.grid[label1, paste0(a, ".delta.x")] * (x - X.eval[label1])
                                link.1.2 <- link.1.2 + target.Z * coef.grid[label2, paste0(a, ".delta.x")] * (x - X.eval[label2])
                                link.0.1 <- link.0.1 + target.Z * coef.grid[label1, paste0(a, ".delta.x")] * (x - X.eval[label1])
                                link.0.2 <- link.0.2 + target.Z * coef.grid[label2, paste0(a, ".delta.x")] * (x - X.eval[label2])
                            }
                        }
                    }

                    link.1 <- c((link.1.1 * d2 + link.1.2 * d1) / (d1 + d2))
                    link.0 <- c((link.0.1 * d2 + link.0.2 * d1) / (d1 + d2))
                }

                if (method == "linear") {
                    TE <- link.1 - link.0
                }
                if (method == "logit") {
                    E.prob.1 <- exp(link.1) / (1 + exp(link.1))
                    E.prob.0 <- exp(link.0) / (1 + exp(link.0))
                    TE <- E.prob.1 - E.prob.0
                }
                if (method == "probit") {
                    E.prob.1 <- pnorm(link.1, 0, 1)
                    E.prob.0 <- pnorm(link.0, 0, 1)
                    TE <- E.prob.1 - E.prob.0
                }
                if (method == "poisson" | method == "nbinom") {
                    E.y.1 <- exp(link.1)
                    E.y.0 <- exp(link.0)
                    TE <- E.y.1 - E.y.0
                }
                return(TE)
            }
        }


        if (length(diff.values) == 2) {
            if (treat.type == "discrete") {
                difference <- c(est.TE(diff.values[2]) - est.TE(diff.values[1]))
            }
            if (treat.type == "continuous") {
                difference <- c(est.ME(diff.values[2]) - est.ME(diff.values[1]))
            }

            vec.list2 <- gen.sd(wls(x = diff.values[2], data = data, bw = bw, weights = w, Xdensity = Xdensity), char = char, D.ref = D.ref, to.diff = TRUE)
            vec.list1 <- gen.sd(wls(x = diff.values[1], data = data, bw = bw, weights = w, Xdensity = Xdensity), char = char, D.ref = D.ref, to.diff = TRUE)
            vec1 <- vec.list1$vec
            vec2 <- vec.list2$vec
            vec <- vec2 - vec1
            temp.vcov.matrix <- vec.list1$temp.vcov.matrix
            difference.sd <- c(sqrt((t(vec) %*% temp.vcov.matrix %*% vec)[1, 1]))
        }

        if (length(diff.values) == 3) {
            if (treat.type == "discrete") {
                difference1 <- est.TE(diff.values[2]) - est.TE(diff.values[1])
                difference2 <- est.TE(diff.values[3]) - est.TE(diff.values[2])
                difference3 <- est.TE(diff.values[3]) - est.TE(diff.values[1])
                difference <- c(difference1, difference2, difference3)
            }

            if (treat.type == "continuous") {
                difference1 <- est.ME(diff.values[2]) - est.ME(diff.values[1])
                difference2 <- est.ME(diff.values[3]) - est.ME(diff.values[2])
                difference3 <- est.ME(diff.values[3]) - est.ME(diff.values[1])
                difference <- c(difference1, difference2, difference3)
            }

            vec.list3 <- gen.sd(wls(x = diff.values[3], data = data, bw = bw, weights = w, Xdensity = Xdensity), char = char, D.ref = D.ref, to.diff = TRUE)
            vec.list2 <- gen.sd(wls(x = diff.values[2], data = data, bw = bw, weights = w, Xdensity = Xdensity), char = char, D.ref = D.ref, to.diff = TRUE)
            vec.list1 <- gen.sd(wls(x = diff.values[1], data = data, bw = bw, weights = w, Xdensity = Xdensity), char = char, D.ref = D.ref, to.diff = TRUE)
            vec1 <- vec.list1$vec
            vec2 <- vec.list2$vec
            vec3 <- vec.list3$vec
            vec21 <- vec2 - vec1
            vec32 <- vec3 - vec2
            vec31 <- vec3 - vec1
            temp.vcov.matrix <- vec.list1$temp.vcov.matrix

            difference.sd <- c(
                sqrt((t(vec21) %*% temp.vcov.matrix %*% vec21)[1, 1]),
                sqrt((t(vec32) %*% temp.vcov.matrix %*% vec32)[1, 1]),
                sqrt((t(vec31) %*% temp.vcov.matrix %*% vec31)[1, 1])
            )
        }

        if (treat.type == "discrete") {
            names(difference) <- paste0(char, ".", difference.name)
            names(difference.sd) <- paste0("sd.", char, ".", difference.name)
        }
        if (treat.type == "continuous") {
            names(difference) <- paste0(names(D.sample)[D.sample == D.ref], ".", difference.name)
            names(difference.sd) <- paste0("sd.", names(D.sample)[D.sample == D.ref], ".", difference.name)
        }
        return(list(difference = difference, difference.sd = difference.sd))
    }

    ## Function D: estimate ATE/AME
    gen.ATE <- function(data, coef.grid, model.vcovs, char = NULL) {
        if (is.null(char) == TRUE) {
            treat.type <- "continuous"
            weights <- data[, "WEIGHTS"]
        } else {
            treat.type <- "discrete"
            which.index <- which(data[, D] == char)
            sub.data <- data[which.index, ]
            weights <- data[which.index, "WEIGHTS"]
        }

        gen.ATE.sub <- function(index) {
            if (treat.type == "discrete") {
                x <- sub.data[index, X]
                Xnew <- abs(X.eval - x)
                d1 <- min(Xnew)
                label1 <- which.min(Xnew)
                Xnew[label1] <- Inf
                d2 <- min(Xnew)
                label2 <- which.min(Xnew)
                if (d1 == 0) {
                    link.1 <- coef.grid[label1, "(Intercept)"] + coef.grid[label1, paste0("D.", char)]
                    link.0 <- coef.grid[label1, "(Intercept)"] + 0
                    if (is.null(Z) == FALSE) {
                        for (a in Z) {
                            target.Z <- sub.data[index, a]
                            link.1 <- link.1 + target.Z * coef.grid[label1, a]
                            link.0 <- link.0 + target.Z * coef.grid[label1, a]
                            # if(full.moderate==TRUE){
                            # 	link.1 <- link.1 + target.Z*coef.grid[label1,paste0(a,".X")]*x
                            # 	link.0 <- link.0 + target.Z*coef.grid[label1,paste0(a,".X")]*x
                            # }
                        }
                    }
                } else if (d2 == 0) {
                    link.1 <- coef.grid[label2, "(Intercept)"] + coef.grid[label2, paste0("D.", char)]
                    link.0 <- coef.grid[label2, "(Intercept)"] + 0
                    if (is.null(Z) == FALSE) {
                        for (a in Z) {
                            target.Z <- sub.data[index, a]
                            link.1 <- link.1 + target.Z * coef.grid[label2, a]
                            link.0 <- link.0 + target.Z * coef.grid[label2, a]
                            # if(full.moderate==TRUE){
                            # 	link.1 <- link.1 + target.Z*coef.grid[label2,paste0(a,".X")]*x
                            # 	link.0 <- link.0 + target.Z*coef.grid[label2,paste0(a,".X")]*x
                            # }
                        }
                    }
                } else { ## weighted average
                    link.1.1 <- coef.grid[label1, "(Intercept)"] + coef.grid[label1, paste0("D.", char)] +
                        (coef.grid[label1, paste0("D.delta.x.", char)] + coef.grid[label1, "delta.x"]) * (x - X.eval[label1])
                    link.1.2 <- coef.grid[label2, "(Intercept)"] + coef.grid[label2, paste0("D.", char)] +
                        (coef.grid[label2, paste0("D.delta.x.", char)] + coef.grid[label2, "delta.x"]) * (x - X.eval[label2])
                    link.0.1 <- coef.grid[label1, "(Intercept)"] + coef.grid[label1, "delta.x"] * (x - X.eval[label1])
                    link.0.2 <- coef.grid[label2, "(Intercept)"] + coef.grid[label2, "delta.x"] * (x - X.eval[label2])
                    if (is.null(Z) == FALSE) {
                        for (a in Z) {
                            target.Z <- sub.data[index, a]
                            link.1.1 <- link.1.1 + target.Z * coef.grid[label1, a]
                            link.1.2 <- link.1.2 + target.Z * coef.grid[label2, a]
                            link.0.1 <- link.0.1 + target.Z * coef.grid[label1, a]
                            link.0.2 <- link.0.2 + target.Z * coef.grid[label2, a]
                            if (full.moderate == TRUE) {
                                link.1.1 <- link.1.1 + coef.grid[label1, paste0(a, ".delta.x")] * (x - X.eval[label1])
                                link.1.2 <- link.1.2 + coef.grid[label2, paste0(a, ".delta.x")] * (x - X.eval[label2])
                                link.0.1 <- link.0.1 + coef.grid[label1, paste0(a, ".delta.x")] * (x - X.eval[label1])
                                link.0.2 <- link.0.2 + coef.grid[label2, paste0(a, ".delta.x")] * (x - X.eval[label2])
                            }
                        }
                    }
                    link.1 <- c((link.1.1 * d2 + link.1.2 * d1) / (d1 + d2))
                    link.0 <- c((link.0.1 * d2 + link.0.2 * d1) / (d1 + d2))
                }

                if (method == "linear") {
                    TE <- link.1 - link.0
                }
                if (method == "logit") {
                    E.prob.1 <- exp(link.1) / (1 + exp(link.1))
                    E.prob.0 <- exp(link.0) / (1 + exp(link.0))
                    TE <- E.prob.1 - E.prob.0
                }
                if (method == "probit") {
                    E.prob.1 <- pnorm(link.1, 0, 1)
                    E.prob.0 <- pnorm(link.0, 0, 1)
                    TE <- E.prob.1 - E.prob.0
                }
                if (method == "poisson" | method == "nbinom") {
                    E.y.1 <- exp(link.1)
                    E.y.0 <- exp(link.0)
                    TE <- E.y.1 - E.y.0
                }

                return(TE)
            }

            if (treat.type == "continuous") {
                x <- data[index, X]
                Xnew <- abs(X.eval - x)
                d1 <- min(Xnew)
                label1 <- which.min(Xnew)
                Xnew[label1] <- Inf
                d2 <- min(Xnew)
                label2 <- which.min(Xnew)

                if (d1 == 0) {
                    link <- coef.grid[label1, "(Intercept)"] + coef.grid[label1, D] * data[index, D]
                    if (is.null(Z) == FALSE) {
                        for (a in Z) {
                            target.Z <- data[index, a]
                            link <- link + target.Z * coef.grid[label1, a]
                            # if(full.moderate==TRUE){
                            # 	link <- link + target.Z*coef.grid[label1,paste0(a,".X")]*x
                            # }
                        }
                    }
                    coef.grid.D <- coef.grid[label1, D]
                } else if (d2 == 0) {
                    link <- coef.grid[label2, "(Intercept)"] + coef.grid[label2, D] * data[index, D]
                    if (is.null(Z) == FALSE) {
                        for (a in Z) {
                            target.Z <- data[index, a]
                            link <- link + target.Z * coef.grid[label2, a]
                            # if(full.moderate==TRUE){
                            # 	link <- link + target.Z*coef.grid[label2,paste0(a,".X")]*x
                            # }
                        }
                    }
                    coef.grid.D <- coef.grid[label2, D]
                } else { ## weighted average
                    link.1 <- coef.grid[label1, "(Intercept)"] + coef.grid[label1, D] * data[index, D] +
                        coef.grid[label1, "D.delta.x"] * data[index, D] * (x - X.eval[label1]) +
                        coef.grid[label1, "delta.x"] * (x - X.eval[label1])

                    link.2 <- coef.grid[label2, "(Intercept)"] + coef.grid[label2, D] * data[index, D] +
                        coef.grid[label2, "D.delta.x"] * data[index, D] * (x - X.eval[label2]) +
                        coef.grid[label2, "delta.x"] * (x - X.eval[label2])

                    if (is.null(Z) == FALSE) {
                        for (a in Z) {
                            target.Z <- data[index, a]
                            link.1 <- link.1 + target.Z * coef.grid[label1, a]
                            link.2 <- link.2 + target.Z * coef.grid[label2, a]
                            if (full.moderate == TRUE) {
                                link.1 <- link.1 + coef.grid[label1, paste0(a, ".delta.x")] * (x - X.eval[label1])
                                link.2 <- link.2 + coef.grid[label2, paste0(a, ".delta.x")] * (x - X.eval[label2])
                            }
                        }
                    }

                    coef.grid.D.1 <- coef.grid[label1, D] + coef.grid[label1, "D.delta.x"] * (x - X.eval[label1])
                    coef.grid.D.2 <- coef.grid[label2, D] + coef.grid[label2, "D.delta.x"] * (x - X.eval[label2])

                    link <- (link.1 * d2 + link.2 * d1) / (d1 + d2)
                    coef.grid.D <- (coef.grid.D.1 * d2 + coef.grid.D.2 * d1) / (d1 + d2)
                }
                if (method == "logit") {
                    ME <- exp(link) / (1 + exp(link))^2 * coef.grid.D
                }
                if (method == "probit") {
                    ME <- coef.grid.D * dnorm(link)
                }
                if (method == "linear") {
                    ME <- coef.grid.D
                }
                if (method == "poisson" | method == "nbinom") {
                    ME <- exp(link) * coef.grid.D
                }
                names(ME) <- NULL
                return(ME)
            }
        }

        gen.ATE.sd.vec <- function(index) {
            if (treat.type == "discrete") {
                x <- sub.data[index, X]
                link.1 <- coef.grid["(Intercept)"] + x * coef.grid[X] + 1 * coef.grid[paste0("D.", char)] + x * coef.grid[X] * coef.grid[paste0("D.delta.x.", char)]
                link.0 <- coef.grid["(Intercept)"] + x * coef.grid[X]
                if (is.null(Z) == FALSE) {
                    for (a in Z) {
                        target.Z <- sub.data[index, a]
                        link.1 <- link.1 + target.Z * coef.grid[a]
                        link.0 <- link.0 + target.Z * coef.grid[a]
                    }
                }
                if (is.null(Z) == FALSE) {
                    vec.1 <- c(1, x, 1, x, as.matrix(sub.data[index, Z]))
                    vec.0 <- c(1, x, 0, 0, as.matrix(sub.data[index, Z]))
                    target.slice <- c("(Intercept)", "delta.x", paste0("D.", char), paste0("D.delta.x.", char), Z)
                } else {
                    vec.1 <- c(1, x, 1, x)
                    vec.0 <- c(1, x, 0, 0)
                    target.slice <- c("(Intercept)", "delta.x", paste0("D.", char), paste0("D.delta.x.", char))
                }
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
                link <- coef.grid["(Intercept)"] + data[index, X] * coef.grid[X] + coef.grid[D] * data[index, D] + coef.grid["D.delta.x"] * data[index, X] * data[index, D]
                if (is.null(Z) == FALSE) {
                    for (a in Z) {
                        target.Z <- data[index, a]
                        link <- link + target.Z * coef.grid[a]
                        if (full.moderate == TRUE) {
                            link <- link + target.Z * coef.grid[paste0(a, ".X")] * data[index, X]
                        }
                    }
                }

                if (is.null(Z) == FALSE) {
                    vec1 <- c(1, data[index, X], data[index, D], data[index, D] * data[index, X], as.matrix(data[index, Z]))
                    vec0 <- c(0, 0, 1, data[index, X], rep(0, length(Z)))
                    target.slice <- c("(Intercept)", X, D, "D.delta.x", Z)
                } else {
                    vec1 <- c(1, data[index, X], data[index, D], data[index, D] * data[index, X])
                    vec0 <- c(0, 0, 1, data[index, X])
                    target.slice <- c("(Intercept)", X, D, "D.delta.x")
                }

                if (method == "logit") {
                    vec <- -(coef.grid[D] + data[index, X] * coef.grid["D.delta.x"]) * (exp(link) - exp(-link)) / (2 + exp(link) + exp(-link))^2 * vec1 + exp(link) / (1 + exp(link))^2 * vec0
                }
                if (method == "probit") {
                    vec <- dnorm(link) * vec0 - (coef.grid[D] + data[index, X] * coef.grid["D.delta.x"]) * link * dnorm(link) * vec1
                }
                if (method == "poisson" | method == "nbinom") {
                    vec <- (coef.grid[D] + data[index, X] * coef.grid["D.delta.x"]) * exp(link) * vec1 + exp(link) * vec0
                }
                if (method == "linear") {
                    vec <- vec0
                }
                return(vec)
            }
        }

        if (treat.type == "discrete") {
            index.all <- c(1:dim(sub.data)[1])
            TE.all.real <- sapply(index.all, function(x) gen.ATE.sub(x))
            ATE <- weighted.mean(TE.all.real, weights, na.rm = TRUE)
            vec.all <- sapply(index.all, function(x) gen.ATE.sd.vec(x))
            vec.mean <- apply(vec.all, 1, function(x) weighted.mean(x, weights))
            if (is.null(Z) == FALSE) {
                target.slice <- c("(Intercept)", "delta.x", paste0("D.", char), paste0("D.delta.x.", char), Z)
            } else {
                target.slice <- c("(Intercept)", "delta.x", paste0("D.", char), paste0("D.delta.x.", char))
            }
            temp.vcov.matrix.mean <- Reduce("+", model.vcovs) / length(model.vcovs)
            temp.vcov.matrix <- temp.vcov.matrix.mean[target.slice, target.slice]
            ATE.sd <- sqrt((t(vec.mean) %*% temp.vcov.matrix %*% vec.mean)[1, 1])
            return(list(ATE = ATE, sd = ATE.sd))
        }

        if (treat.type == "continuous") {
            index.all <- c(1:dim(data)[1])
            ME.all.real <- sapply(index.all, function(x) gen.ATE.sub(x))
            names(ME.all.real) <- NULL
            AME <- weighted.mean(ME.all.real, weights, na.rm = TRUE)
            vec.all <- sapply(index.all, function(x) gen.ATE.sd.vec(x))
            vec.mean <- apply(vec.all, 1, function(x) weighted.mean(x, weights))
            if (is.null(Z) == FALSE) {
                target.slice <- c("(Intercept)", "delta.x", D, "D.delta.x", Z)
            } else {
                target.slice <- c("(Intercept)", "delta.x", D, "D.delta.x")
            }

            temp.vcov.matrix.mean <- Reduce("+", model.vcovs) / length(model.vcovs)
            temp.vcov.matrix <- temp.vcov.matrix.mean[target.slice, target.slice]
            AME.sd <- sqrt((t(vec.mean) %*% temp.vcov.matrix %*% vec.mean)[1, 1])
            return(list(AME = AME, sd = AME.sd))
        }
    }

    all.output.noCI <- list()
    if (treat.type == "discrete") {
        for (char in other.treat) {
            gen.TE.output <- gen.kernel.TE(coef.grid = coef.grid, char = char)
            gen.diff.output <- gen.kernel.difference(
                coef.grid = coef.grid,
                diff.values = diff.values,
                char = char
            )
            gen.ATE.output <- gen.ATE(data = data, coef.grid = coef.grid, model.vcovs = model.vcovs, char = char)

            all.output.noCI[[char]] <- list(
                TE = gen.TE.output,
                diff = gen.diff.output,
                ATE = gen.ATE.output
            )
        }
        gen.TE.output.base <- gen.kernel.TE(coef.grid = coef.grid, char = base, base.flag = TRUE)
    }

    if (treat.type == "continuous") {
        k <- 1
        for (D.ref in D.sample) {
            gen.ME.output <- gen.kernel.TE(coef.grid = coef.grid, D.ref = D.ref)

            gen.diff.output <- gen.kernel.difference(
                coef.grid = coef.grid,
                diff.values = diff.values,
                D.ref = D.ref
            )
            all.output.noCI[[label.name[k]]] <- list(
                ME = gen.ME.output,
                diff = gen.diff.output
            )
            k <- k + 1
        }

        AME.estimate <- gen.ATE(coef.grid = coef.grid, model.vcovs = model.vcovs, data = data)
        all.output.noCI[["AME"]] <- AME.estimate
    }

    if (CI == TRUE) {
        if (vartype == "bootstrap") {
            # Part1: Bootstrap
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
                if (treat.type == "discrete") {
                    if (length(unique(data.boot[, D])) != length(unique(data[, D]))) {
                        return(boot.out)
                    }
                }

                if (is.null(weights) == TRUE) {
                    w.touse <- rep(1, dim(data.boot)[1])
                } else {
                    w.touse <- data.boot[, weights]
                }

                # Xdensity
                suppressWarnings(Xdensity.boot <- density(data.boot[, X], weights = w.touse))
                coef.grid.boot <- c()
                for (x in X.eval) {
                    coef.grid.boot <- rbind(coef.grid.boot, wls(x = x, data = data.boot, bw = bw, weights = w.touse, Xdensity = Xdensity.boot)$result)
                }

                boot.one.round <- c()
                if (treat.type == "discrete") {
                    for (char in other.treat) {
                        gen.TE.output <- gen.kernel.TE(coef.grid = coef.grid.boot, char = char)
                        gen.diff.output <- gen.kernel.difference(
                            coef.grid = coef.grid.boot,
                            diff.values = diff.values,
                            char = char
                        )

                        gen.ATE.output <- gen.ATE(coef.grid = coef.grid.boot, data = data.boot, model.vcovs = model.vcovs, char = char)

                        TE.output <- gen.TE.output$TE
                        names(TE.output) <- rep(paste0("TE.", char), neval)

                        E.pred.output <- gen.TE.output$E.pred
                        names(E.pred.output) <- rep(paste0("pred.", char), neval)
                        E.base.output <- gen.TE.output$E.base

                        link.output <- gen.TE.output$link.1
                        names(link.output) <- rep(paste0("link.", char), neval)
                        link0.output <- gen.TE.output$link.0

                        diff.estimate.output <- c(gen.diff.output$difference)
                        names(diff.estimate.output) <- rep(paste0("diff.", char), length(difference.name))

                        ATE.estimate <- c(gen.ATE.output$ATE)
                        names(ATE.estimate) <- paste0("ATE.", char)
                        boot.one.round <- c(boot.one.round, TE.output, E.pred.output, link.output, diff.estimate.output, ATE.estimate)
                    }
                    names(E.base.output) <- rep(paste0("pred.", base), neval)
                    boot.one.round <- c(boot.one.round, E.base.output)

                    names(link0.output) <- rep(paste0("link.", base), neval)
                    boot.one.round <- c(boot.one.round, link0.output)
                }

                if (treat.type == "continuous") {
                    k <- 1
                    for (D.ref in D.sample) {
                        gen.kernel.ME.output <- gen.kernel.TE(coef.grid = coef.grid.boot, D.ref = D.ref)
                        ME.output <- gen.kernel.ME.output$ME
                        names(ME.output) <- rep(paste0("ME.", label.name[k]), neval)
                        E.pred.output <- gen.kernel.ME.output$E.pred
                        names(E.pred.output) <- rep(paste0("pred.", label.name[k]), neval)

                        link.output <- gen.kernel.ME.output$link
                        names(link.output) <- rep(paste0("link.", label.name[k]), neval)

                        gen.diff.output <- gen.kernel.difference(
                            coef.grid = coef.grid.boot,
                            diff.values = diff.values,
                            D.ref = D.ref
                        )

                        diff.estimate.output <- c(gen.diff.output$difference)
                        names(diff.estimate.output) <- rep(paste0("diff.", label.name[k]), length(difference.name))
                        boot.one.round <- c(boot.one.round, ME.output, E.pred.output, link.output, diff.estimate.output)
                        k <- k + 1
                    }
                    AME.estimate <- c(gen.ATE(coef.grid = coef.grid.boot, model.vcovs = model.vcovs, data = data.boot)$AME)
                    names(AME.estimate) <- c("AME")
                    boot.one.round <- c(boot.one.round, AME.estimate)
                }
                boot.out <- cbind(boot.out, boot.one.round)
                rownames(boot.out) <- names(boot.one.round)
                colnames(boot.out) <- NULL

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
                        .export = c("one.boot"), .packages = c("MASS", "AER"),
                        .inorder = FALSE
                    ) %dopar% {
                        output.all <- try(one.boot(),silent = TRUE)
                        if('try-error' %in% class(output.all)){
                            return(NA)
                        }
                        else{
                            return(output.all)
                        }
                    }
                )
                suppressWarnings(stopCluster(pcl))
                cat("\r")
            } else {
                bootout <- matrix(NA, all.length, 0)
                for (i in 1:nboots) {
                    suppressWarnings(tempdata <- try(one.boot(),silent = TRUE))
                    if('try-error' %in% class(tempdata)){
                        bootout <- cbind(bootout, NA)
                    }
                    else{
                        bootout <- cbind(bootout, tempdata)
                    }                    
                    #if (is.null(tempdata) == FALSE) {
                    #    bootout <- cbind(bootout, tempdata)
                    #}
                    if (i %% 50 == 0) cat(i) else cat(".")
                }
                cat("\r")
            }

            if (treat.type == "discrete") {
                TE.output.all.list <- list()
                pred.output.all.list <- list()
                diff.output.all.list <- list()
                TE.vcov.list <- list()
                ATE.output.list <- list()
                link.output.all.list <- list()
                for (char in other.treat) {
                    gen.general.TE.output <- all.output.noCI[[char]]$TE
                    TE.output <- gen.general.TE.output$TE
                    E.pred.output <- gen.general.TE.output$E.pred
                    E.base.output <- gen.general.TE.output$E.base
                    link.output <- gen.general.TE.output$link.1
                    link0.output <- gen.general.TE.output$link.0
                    diff.estimate.output <- all.output.noCI[[char]]$diff$difference
                    ATE.estimate <- all.output.noCI[[char]]$ATE$ATE

                    TE.boot.matrix <- bootout[rownames(bootout) == paste0("TE.", char), ]
                    pred.boot.matrix <- bootout[rownames(bootout) == paste0("pred.", char), ]
                    base.boot.matrix <- bootout[rownames(bootout) == paste0("pred.", base), ]
                    link.boot.matrix <- bootout[rownames(bootout) == paste0("link.", char), ]
                    link0.boot.matrix <- bootout[rownames(bootout) == paste0("link.", base), ]
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
                    link.boot.sd <- apply(link.boot.matrix, 1, sd, na.rm = TRUE)
                    link0.boot.sd <- apply(link0.boot.matrix, 1, sd, na.rm = TRUE)
                    diff.boot.sd <- apply(diff.boot.matrix, 1, sd, na.rm = TRUE)
                    ATE.boot.sd <- apply(ATE.boot.matrix, 1, sd, na.rm = TRUE)

                    TE.boot.CI <- t(apply(TE.boot.matrix, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
                    pred.boot.CI <- t(apply(pred.boot.matrix, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
                    base.boot.CI <- t(apply(base.boot.matrix, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
                    link.boot.CI <- t(apply(link.boot.matrix, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
                    link0.boot.CI <- t(apply(link0.boot.matrix, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
                    diff.boot.CI <- t(apply(diff.boot.matrix, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
                    ATE.boot.CI <- matrix(t(apply(ATE.boot.matrix, 1, quantile, c(0.025, 0.975), na.rm = TRUE)), nrow = 1)

                    TE.boot.uniform.CI <- calculate_uniform_quantiles(TE.boot.matrix,0.05)
                    uniform.coverage <- TE.boot.uniform.CI$coverage
                    TE.boot.uniform.CI <- TE.boot.uniform.CI$Q_j

                    pred.boot.uniform.CI <- calculate_uniform_quantiles(pred.boot.matrix,0.05)
                    pred.boot.uniform.CI <- pred.boot.uniform.CI$Q_j

                    link.boot.uniform.CI <- calculate_uniform_quantiles(link.boot.matrix,0.05)
                    link.boot.uniform.CI <- link.boot.uniform.CI$Q_j

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

                    link.output.all <- cbind(X.eval, link.output, link.boot.sd, link.boot.CI[, 1], link.boot.CI[, 2], link.boot.uniform.CI[,1], link.boot.uniform.CI[,2])
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

                link0.boot.uniform.CI <- calculate_uniform_quantiles(link0.boot.matrix,0.05)
                link0.boot.uniform.CI <- link0.boot.uniform.CI$Q_j

                # base
                base.output.all <- cbind(X.eval, E.base.output, base.boot.sd, base.boot.CI[, 1], base.boot.CI[, 2],base.boot.uniform.CI[,1],base.boot.uniform.CI[,2])
                colnames(base.output.all) <- c("X", "E(Y)", "sd", "lower CI(95%)", "upper CI(95%)","lower uniform CI(95%)", "upper uniform CI(95%)")
                rownames(base.output.all) <- NULL
                pred.output.all.list[[all.treat.origin[base]]] <- base.output.all

                link0.output.all <- cbind(X.eval, link0.output, link0.boot.sd, link0.boot.CI[, 1], link0.boot.CI[, 2],link0.boot.uniform.CI[,1], link0.boot.uniform.CI[,2])
                colnames(link0.output.all) <- c("X", "E(Y)", "sd", "lower CI(95%)", "upper CI(95%)","lower uniform CI(95%)", "upper uniform CI(95%)")
                rownames(link0.output.all) <- NULL
                link.output.all.list[[all.treat.origin[base]]] <- link0.output.all
            }


            if (treat.type == "continuous") {
                ME.output.all.list <- list()
                pred.output.all.list <- list()
                diff.output.all.list <- list()
                ME.vcov.list <- list()
                link.output.all.list <- list()
                k <- 1
                for (D.ref in D.sample) {
                    label <- label.name[k]
                    gen.general.ME.output <- all.output.noCI[[label]]$ME
                    ME.output <- gen.general.ME.output$ME
                    E.pred.output <- gen.general.ME.output$E.pred
                    link.output <- gen.general.ME.output$link
                    diff.estimate.output <- all.output.noCI[[label]]$diff$difference

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
                AME.estimate <- all.output.noCI$AME$AME
                AME.boot.matrix <- matrix(bootout[rownames(bootout) == "AME", ], nrow = 1)
                AME.boot.sd <- apply(AME.boot.matrix, 1, sd, na.rm = TRUE)
                AME.boot.CI <- matrix(t(apply(AME.boot.matrix, 1, quantile, c(0.025, 0.975), na.rm = TRUE)), nrow = 1)
                AME.z.value <- AME.estimate / AME.boot.sd
                AME.p.value <- 2 * pnorm(-abs(AME.z.value))
                AME.output <- c(AME.estimate, AME.boot.sd, AME.z.value, AME.p.value, AME.boot.CI[, 1], AME.boot.CI[, 2])
                names(AME.output) <- c("AME", "sd", "z-value", "p-value", "lower CI(95%)", "upper CI(95%)")
            }
        }

        if (vartype == "delta") {
            crit <- abs(qt(.025, df = model.dfs))
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

            if (treat.type == "discrete") {
                TE.output.all.list <- list()
                pred.output.all.list <- list()
                diff.output.all.list <- list()
                TE.vcov.list <- list()
                ATE.output.list <- list()
                link.output.all.list <- list()
                for (char in other.treat) {
                    gen.general.TE.output <- all.output.noCI[[char]]$TE
                    TE.output <- gen.general.TE.output$TE
                    E.pred.output <- gen.general.TE.output$E.pred
                    E.base.output <- gen.general.TE.output$E.base
                    link.output <- gen.general.TE.output$link.1
                    link0.output <- gen.general.TE.output$link.0
                    diff.estimate.output <- all.output.noCI[[char]]$diff
                    ATE.estimate <- all.output.noCI[[char]]$ATE

                    TE.delta.sd <- gen.general.TE.output$TE.sd
                    TE.output.all <- cbind(X.eval, TE.output, TE.delta.sd, TE.output - crit * TE.delta.sd, TE.output + crit * TE.delta.sd)
                    colnames(TE.output.all) <- c("X", "TE", "sd", "lower CI(95%)", "upper CI(95%)")
                    rownames(TE.output.all) <- NULL
                    TE.output.all.list[[other.treat.origin[char]]] <- TE.output.all

                    pred.delta.sd <- gen.general.TE.output$predict.sd
                    pred.output.all <- cbind(X.eval, E.pred.output, pred.delta.sd, E.pred.output - crit * pred.delta.sd, E.pred.output + crit * pred.delta.sd)
                    colnames(pred.output.all) <- c("X", "E(Y)", "sd", "lower CI(95%)", "upper CI(95%)")
                    rownames(pred.output.all) <- NULL
                    pred.output.all.list[[other.treat.origin[char]]] <- pred.output.all

                    link.delta.sd <- gen.general.TE.output$link.sd
                    link.output.all <- cbind(X.eval, link.output, link.delta.sd, link.output - crit * link.delta.sd, link.output + crit * link.delta.sd)
                    colnames(link.output.all) <- c("X", "E(Y)", "sd", "lower CI(95%)", "upper CI(95%)")
                    rownames(link.output.all) <- NULL
                    link.output.all.list[[other.treat.origin[char]]] <- link.output.all

                    TE.vcov.list[[other.treat.origin[char]]] <- NULL

                    diff.estimate.value <- diff.estimate.output$difference
                    diff.delta.sd <- diff.estimate.output$difference.sd

                    z.value <- diff.estimate.value / diff.delta.sd
                    p.value <- 2 * pnorm(-abs(z.value))
                    diff.output.all <- cbind(
                        diff.estimate.value, diff.delta.sd,
                        z.value, p.value, diff.estimate.value - crit[1] * diff.delta.sd, diff.estimate.value + crit[1] * diff.delta.sd
                    )

                    colnames(diff.output.all) <- c("diff.estimate", "sd", "z-value", "p-value", "lower CI(95%)", "upper CI(95%)")
                    rownames(diff.output.all) <- difference.name

                    diff.output.all.list[[other.treat.origin[char]]] <- diff.output.all

                    ATE.estimate.value <- ATE.estimate$ATE
                    ATE.delta.sd <- ATE.estimate$sd
                    ATE.z.value <- ATE.estimate.value / ATE.delta.sd
                    ATE.p.value <- 2 * pnorm(-abs(ATE.z.value))
                    ATE.output <- c(ATE.estimate.value, ATE.delta.sd, ATE.z.value, ATE.p.value, ATE.estimate.value - crit[1] * ATE.delta.sd, ATE.estimate.value + crit[1] * ATE.delta.sd)
                    names(ATE.output) <- c("ATE", "sd", "z-value", "p-value", "lower CI(95%)", "upper CI(95%)")
                    ATE.output.list[[other.treat.origin[char]]] <- ATE.output
                }
                # base
                base.delta.sd <- gen.TE.output.base$predict.sd
                base.output.all <- cbind(X.eval, E.base.output, base.delta.sd, E.base.output - crit * base.delta.sd, E.base.output + crit * base.delta.sd)
                colnames(base.output.all) <- c("X", "E(Y)", "sd", "lower CI(95%)", "upper CI(95%)")
                rownames(base.output.all) <- NULL

                pred.output.all.list[[all.treat.origin[base]]] <- base.output.all

                link0.delta.sd <- gen.TE.output.base$link.sd
                link0.output.all <- cbind(X.eval, link0.output, link0.delta.sd, link0.output - crit * link0.delta.sd, link0.output + crit * link0.delta.sd)
                colnames(link0.output.all) <- c("X", "E(Y)", "sd", "lower CI(95%)", "upper CI(95%)")
                rownames(link0.output.all) <- NULL
                link.output.all.list[[all.treat.origin[base]]] <- link0.output.all
            }


            if (treat.type == "continuous") {
                ME.output.all.list <- list()
                pred.output.all.list <- list()
                diff.output.all.list <- list()
                ME.vcov.list <- list()
                link.output.all.list <- list()
                k <- 1
                for (D.ref in D.sample) {
                    label <- label.name[k]
                    gen.general.ME.output <- all.output.noCI[[label]]$ME
                    ME.output <- gen.general.ME.output$ME
                    E.pred.output <- gen.general.ME.output$E.pred
                    link.output <- gen.general.ME.output$link
                    diff.estimate.output <- all.output.noCI[[label]]$diff

                    ME.delta.sd <- gen.general.ME.output$ME.sd
                    ME.output.all <- cbind(X.eval, ME.output, ME.delta.sd, ME.output - crit * ME.delta.sd, ME.output + crit * ME.delta.sd)
                    colnames(ME.output.all) <- c("X", "ME", "sd", "lower CI(95%)", "upper CI(95%)")
                    rownames(ME.output.all) <- NULL
                    ME.output.all.list[[label]] <- ME.output.all

                    pred.delta.sd <- gen.general.ME.output$predict.sd
                    pred.output.all <- cbind(X.eval, E.pred.output, pred.delta.sd, E.pred.output - crit * pred.delta.sd, E.pred.output + crit * pred.delta.sd)
                    colnames(pred.output.all) <- c("X", "E(Y)", "sd", "lower CI(95%)", "upper CI(95%)")
                    rownames(pred.output.all) <- NULL
                    pred.output.all.list[[label]] <- pred.output.all

                    link.delta.sd <- gen.general.ME.output$link.sd
                    link.output.all <- cbind(X.eval, link.output, link.delta.sd, link.output - crit * link.delta.sd, link.output + crit * link.delta.sd)
                    colnames(link.output.all) <- c("X", "E(Y)", "sd", "lower CI(95%)", "upper CI(95%)")
                    rownames(link.output.all) <- NULL
                    link.output.all.list[[label]] <- link.output.all

                    ME.vcov.list[[label]] <- NULL

                    diff.estimate.value <- diff.estimate.output$difference
                    diff.delta.sd <- diff.estimate.output$difference.sd
                    z.value <- diff.estimate.value / diff.delta.sd
                    p.value <- 2 * pnorm(-abs(z.value))
                    diff.output.all <- cbind(
                        diff.estimate.value, diff.delta.sd,
                        z.value, p.value, diff.estimate.value - crit[1] * diff.delta.sd, diff.estimate.value + crit[1] * diff.delta.sd
                    )
                    colnames(diff.output.all) <- c("diff.estimate", "sd", "z-value", "p-value", "lower CI(95%)", "upper CI(95%)")
                    rownames(diff.output.all) <- difference.name
                    diff.output.all.list[[label]] <- diff.output.all

                    k <- k + 1
                }
                AME.estimate.value <- all.output.noCI$AME$AME
                AME.delta.sd <- all.output.noCI$AME$sd
                AME.z.value <- AME.estimate.value / AME.delta.sd
                AME.p.value <- 2 * pnorm(-abs(AME.z.value))
                AME.output <- c(AME.estimate.value, AME.delta.sd, AME.z.value, AME.p.value, AME.estimate.value - crit[1] * AME.delta.sd, AME.estimate.value + crit[1] * AME.delta.sd)
                names(AME.output) <- c("AME", "sd", "z-value", "p-value", "lower CI(95%)", "upper CI(95%)")
            }
        }
    }

    if (CI == FALSE) {
        if (treat.type == "discrete") {
            TE.output.all.list <- list()
            pred.output.all.list <- list()
            diff.output.all.list <- list()
            link.output.all.list <- list()
            TE.vcov.list <- NULL
            ATE.output.list <- list()
            for (char in other.treat) {
                gen.general.TE.output <- all.output.noCI[[char]]$TE
                TE.output <- gen.general.TE.output$TE
                E.pred.output <- gen.general.TE.output$E.pred
                E.base.output <- gen.general.TE.output$E.base
                link.output <- gen.general.TE.output$link.1
                link0.output <- gen.general.TE.output$link.0
                diff.estimate.output <- all.output.noCI[[char]]$diff$difference
                ATE.estimate <- all.output.noCI[[char]]$ATE$ATE

                TE.output.all <- cbind(X.eval, TE.output)
                colnames(TE.output.all) <- c("X", "TE")
                rownames(TE.output.all) <- NULL
                TE.output.all.list[[other.treat.origin[char]]] <- TE.output.all

                pred.output.all <- cbind(X.eval, E.pred.output)
                colnames(pred.output.all) <- c("X", "E(Y)")
                rownames(pred.output.all) <- NULL
                pred.output.all.list[[other.treat.origin[char]]] <- pred.output.all

                link.output.all <- cbind(X.eval, link.output)
                colnames(link.output.all) <- c("X", "E(Y)")
                rownames(link.output.all) <- NULL
                link.output.all.list[[other.treat.origin[char]]] <- link.output.all

                diff.output.all <- cbind(diff.estimate.output)
                colnames(diff.output.all) <- c("diff.estimate")
                rownames(diff.output.all) <- difference.name
                diff.output.all.list[[other.treat.origin[char]]] <- diff.output.all

                ATE.output <- c(ATE.estimate)
                names(ATE.output) <- c("ATE")
                ATE.output.list[[other.treat.origin[char]]] <- ATE.output
            }

            # base
            base.output.all <- cbind(X.eval, E.base.output)
            colnames(base.output.all) <- c("X", "E(Y)")
            rownames(base.output.all) <- NULL
            pred.output.all.list[[all.treat.origin[base]]] <- base.output.all

            link0.output.all <- cbind(X.eval, link0.output)
            colnames(link0.output.all) <- c("X", "E(Y)")
            rownames(link0.output.all) <- NULL
            link.output.all.list[[all.treat.origin[base]]] <- link0.output.all
        }

        if (treat.type == "continuous") {
            ME.output.all.list <- list()
            pred.output.all.list <- list()
            diff.output.all.list <- list()
            ME.vcov.list <- NULL
            link.output.all.list <- list()
            k <- 1
            for (D.ref in D.sample) {
                label <- label.name[k]
                gen.general.ME.output <- all.output.noCI[[label]]$ME
                ME.output <- gen.general.ME.output$ME
                E.pred.output <- gen.general.ME.output$E.pred
                link.output <- gen.general.ME.output$link
                diff.estimate.output <- all.output.noCI[[label]]$diff$difference

                ME.output.all <- cbind(X.eval, ME.output)
                colnames(ME.output.all) <- c("X", "ME")
                rownames(ME.output.all) <- NULL
                ME.output.all.list[[label]] <- ME.output.all

                pred.output.all <- cbind(X.eval, E.pred.output)
                colnames(pred.output.all) <- c("X", "E(Y)")
                rownames(pred.output.all) <- NULL
                pred.output.all.list[[label]] <- pred.output.all

                link.output.all <- cbind(X.eval, link.output)
                colnames(link.output.all) <- c("X", "E(Y)")
                rownames(link.output.all) <- NULL
                link.output.all.list[[label]] <- link.output.all

                diff.output.all <- cbind(diff.estimate.output)
                colnames(diff.output.all) <- c("diff.estimate")
                rownames(diff.output.all) <- difference.name
                diff.output.all.list[[label]] <- diff.output.all

                k <- k + 1
            }
            AME.estimate <- all.output.noCI$AME$AME
            AME.output <- c(AME.estimate)
            names(AME.output) <- c("AME")
        }
    }


    # density or histogram
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

    ## Output
    if (treat.type == "discrete") {
        for (char in other.treat.origin) {
            if (CI == TRUE) {
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
        }

        final.output <- list(
            diff.info = diff.info,
            treat.info = treat.info,
            bw = bw,
            CV.output = Error,
            CI = CI,
            est.kernel = TE.output.all.list,
            uniform.coverage = uniform.coverage,
            pred.kernel = pred.output.all.list,
            link.kernel = link.output.all.list,
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
            estimator = "kernel",
            use.fe = use_fe
        )
    }

    if (treat.type == "continuous") {
        for (label in label.name) {
            if (CI == TRUE) {
                new.diff.est <- as.data.frame(diff.output.all.list[[label]])
                for (i in 1:dim(new.diff.est)[1]) {
                    new.diff.est[i, ] <- sprintf(new.diff.est[i, ], fmt = "%#.3f")
                }
                diff.output.all.list[[label]] <- new.diff.est
            }
        }


        if (CI == TRUE) {
            outss <- data.frame(lapply(AME.output, function(x) t(data.frame(x))))
            colnames(outss) <- c("AME", "sd", "z-value", "p-value", "lower CI(95%)", "upper CI(95%)")
            rownames(outss) <- c("Average Marginal Effect")
            outss[1, ] <- sprintf(outss[1, ], fmt = "%#.3f")
            AME.output <- outss
        }

        final.output <- list(
            diff.info = diff.info,
            treat.info = treat.info,
            bw = bw,
            CV.output = Error,
            CI = CI,
            est.kernel = ME.output.all.list,
            uniform.coverage = uniform.coverage,
            pred.kernel = pred.output.all.list,
            link.kernel = link.output.all.list,
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
            estimator = "kernel",
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

## This is the codes of "createFolds" in the package caret
"createFolds" <-
    function(y, k = 10, list = TRUE, returnTrain = FALSE) {
        if (class(y)[1] == "Surv") y <- y[, "time"]
        if (is.numeric(y)) {
            ## Group the numeric data based on their magnitudes
            ## and sample within those groups.

            ## When the number of samples is low, we may have
            ## issues further slicing the numeric data into
            ## groups. The number of groups will depend on the
            ## ratio of the number of folds to the sample size.
            ## At most, we will use quantiles. If the sample
            ## is too small, we just do regular unstratified
            ## CV
            cuts <- floor(length(y) / k)
            if (cuts < 2) cuts <- 2
            if (cuts > 5) cuts <- 5
            breaks <- unique(quantile(y, probs = seq(0, 1, length = cuts)))
            y <- cut(y, breaks, include.lowest = TRUE)
        }

        if (k < length(y)) {
            ## reset levels so that the possible levels and
            ## the levels in the vector are the same
            y <- factor(as.character(y))
            numInClass <- table(y)
            foldVector <- vector(mode = "integer", length(y))

            ## For each class, balance the fold allocation as far
            ## as possible, then resample the remainder.
            ## The final assignment of folds is also randomized.
            for (i in 1:length(numInClass)) {
                ## create a vector of integers from 1:k as many times as possible without
                ## going over the number of samples in the class. Note that if the number
                ## of samples in a class is less than k, nothing is produced here.
                min_reps <- numInClass[i] %/% k
                if (min_reps > 0) {
                    spares <- numInClass[i] %% k
                    seqVector <- rep(1:k, min_reps)
                    ## add enough random integers to get  length(seqVector) == numInClass[i]
                    if (spares > 0) seqVector <- c(seqVector, sample(1:k, spares))
                    ## shuffle the integers for fold assignment and assign to this classes's data
                    foldVector[which(y == names(numInClass)[i])] <- sample(seqVector)
                } else {
                    ## Here there are less records in the class than unique folds so
                    ## randomly sprinkle them into folds.
                    foldVector[which(y == names(numInClass)[i])] <- sample(1:k, size = numInClass[i])
                }
            }
        } else {
            foldVector <- seq(along = y)
        }

        if (list) {
            out <- split(seq(along = y), foldVector)
            names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))), sep = "")
            if (returnTrain) out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
        } else {
            out <- foldVector
        }
        out
    }
