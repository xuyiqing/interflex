interflex.aipw <- function(data,
                           Y,
                           D,
                           X,
                           treat.info,
                           diff.info,
                           Z = NULL,
                           FE = NULL,
                           signal = c("outcome", "ipw", "aipw"),
                           weights = NULL,
                           B = 200,
                           alpha = 0.05,
                           model.y = "lasso",
                           model.t = "lasso",
                           model.ps = "lasso",
                           basis.type = c("polynomial", "bspline", "none"),
                           include.interactions = TRUE,
                           poly.degree = 2,
                           spline.df = 4,
                           spline.degree = 2,
                           lambda.seq = NULL,
                           reduce.dimension = c("bspline", "kernel"),
                           bw = NULL,
                           best.span = NULL,
                           x.eval = NULL,
                           verbose = TRUE,
                           figure = TRUE,
                           CI = CI,
                           order = NULL,
                           subtitles = NULL,
                           show.subtitles = NULL,
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
  diff.values.plot <- diff.info[["diff.values.plot"]]
  treat.type <- treat.info[["treat.type"]]
  if (treat.type == "discrete") {
    other.treat <- treat.info[["other.treat"]]
    other.treat.origin <- names(other.treat)
    names(other.treat.origin) <- other.treat
    all.treat <- treat.info[["all.treat"]]
    all.treat.origin <- names(all.treat)
    names(all.treat.origin) <- all.treat
  }
  if (treat.type == "continuous") {
    D.sample <- treat.info[["D.sample"]]
    label.name <- names(D.sample)
    # names(label.name) <- D.sample
  }

  n <- dim(data)[1]
  if (is.null(weights) == TRUE) {
    w <- rep(1, n)
  } else {
    w <- data[, weights]
  }
  data[["WEIGHTS"]] <- w

  if (TRUE) {
    if (treat.type == "discrete") {
      if (is.null(weights) == TRUE) {
        de <- density(data[, X])
      } else {
        suppressWarnings(de <- density(data[, X], weights = data[, "WEIGHTS"]))
      }

      treat_den <- list()
      for (char in all.treat) {
        if (is.null(weights) == TRUE) {
          de.tr <- density(data[data[, D] == char, X])
        } else {
          suppressWarnings(de.tr <- density(data[data[, D] == char, X], weights = data[data[, D] == char, "WEIGHTS"]))
        }
        treat_den[[all.treat.origin[char]]] <- de.tr
      }

      if (is.null(weights) == TRUE) {
        hist.out <- hist(data[, X], breaks = 80, plot = FALSE)
      } else {
        suppressWarnings(hist.out <- hist(data[, X], data[, "WEIGHTS"],
          breaks = 80, plot = FALSE
        ))
      }
      n.hist <- length(hist.out$mids)

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
      de.tr <- NULL
    }
  }

  treat.base <- treat.info[["base"]]

  # reticulate::use_virtualenv("r-reticulate", required = TRUE)
  # use_condaenv(condaenv = "r-reticulate", required = TRUE)

  TE.output.all.list <- list()
  if (treat.type == "discrete") {
    for (char in other.treat) {
      data_part <- data[data[[D]] %in% c(treat.base, char), ]
      data_part[data_part[[D]] == treat.base, D] <- 0L
      data_part[data_part[[D]] == char, D] <- 1L
      result <- bootstrapCME(
        data_part,
        Y = Y,
        D = D,
        X = X,
        Z = Z,
        FE = FE,
        signal = signal,
        B = B,
        alpha = alpha,
        outcome_model_type = model.y,
        ps_model_type = model.ps,
        basis_type = basis.type,
        include_interactions = include.interactions,
        poly_degree = poly.degree,
        spline_df = spline.df,
        spline_degree = spline.degree,
        lambda_seq = lambda.seq,
        reduce.dimension = reduce.dimension,
        best_span = best.span,
        x.eval = x.eval,
        verbose = verbose
      )
      TE.output.all <- data.frame(result$results, check.names = FALSE)
      TE.output.all.list[[other.treat.origin[char]]] <- TE.output.all
    }
  } else if (treat.type == "continuous") {
    result <- bootstrapCME_PLR(
      data,
      Y = Y,
      D = D,
      X = X,
      Z = Z,
      B = B,
      alpha = alpha,
      outcome_model_type = model.y,
      treatment_model_type = model.t,
      basis_type = basis.type,
      include_interactions = include.interactions,
      poly_degree = poly.degree,
      spline_df = spline.df,
      spline_degree = spline.degree,
      lambda_seq = lambda.seq,
      reduce.dimension = reduce.dimension,
      bw = bw,
      x.eval = x.eval,
      verbose = verbose
    )
    TE.output.all <- data.frame(result$results, check.names = FALSE)

    k <- 1
    for (d_ref in D.sample) {
      label <- label.name[k]
      TE.output.all.list[[label]] <- TE.output.all
      k <- k + 1
    }
  }
  if (treat.type == "discrete") {
    final.output <- list(
      diff.info = diff.info,
      treat.info = treat.info,
      est.aipw = TE.output.all.list,
      Xlabel = Xlabel,
      Dlabel = Dlabel,
      Ylabel = Ylabel,
      de = de,
      hist.out = hist.out,
      de.tr = treat_den,
      count.tr = treat.hist,
      estimator = "aipw"
    )
  } else if (treat.type == "continuous") {
    final.output <- list(
      diff.info = diff.info,
      treat.info = treat.info,
      est.aipw = TE.output.all.list,
      Xlabel = Xlabel,
      Dlabel = Dlabel,
      Ylabel = Ylabel,
      de = de,
      hist.out = hist.out,
      de.tr = de.tr,
      count.tr = NULL,
      estimator = "aipw"
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
