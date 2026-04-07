interflex.lasso <- function(data, 
                         Y,
                         D,
                         X,
                         treat.info,
                         diff.info,
                         Z = NULL,
                         FE = NULL,
                         
                         # For binary D only:
                         signal    = c("outcome", "ipw", "aipw"),
                         estimand  = c("ATE", "ATT"),
                         
                         weights = NULL,
                         B       = 200,
                         alpha   = 0.05,
                         
                         model.y  = "lasso",
                         model.t  = "lasso",
                         neval = 50,
                         basis.type          = c("polynomial", "bspline", "none"),
                         include.interactions = TRUE,
                         poly.degree         = 2,
                         spline.df           = 4,
                         spline.degree       = 2,
                         
                         lambda.seq         = NULL,
                         reduce.dimension   = c("bspline", "kernel"),
                         bw                 = NULL,
                         best.span          = NULL,
                         x.eval             = NULL,
                         cores              = 8,
                         verbose            = TRUE,
                         figure             = TRUE,
                         show.uniform.CI    = TRUE,
                         CI                 = TRUE,
                         order              = NULL,
                         subtitles          = NULL,
                         show.subtitles     = NULL,
                         Xdistr             = "histogram",
                         main               = NULL,
                         Ylabel             = NULL,
                         Dlabel             = NULL,
                         Xlabel             = NULL,
                         xlab               = NULL,
                         ylab               = NULL,
                         xlim               = NULL,
                         ylim               = NULL,
                         theme.bw           = TRUE,
                         show.grid          = TRUE,
                         cex.main           = NULL,
                         cex.sub            = NULL,
                         cex.lab            = NULL,
                         cex.axis           = NULL,
                         interval           = NULL,
                         file               = NULL,
                         ncols              = NULL,
                         pool               = FALSE,
                         color              = NULL,
                         legend.title       = NULL,
                         show.all           = FALSE,
                         scale              = 1.1,
                         height             = 7,
                         width              = 10) {
  basis.type        <- match.arg(basis.type)
  reduce.dimension  <- match.arg(reduce.dimension)
  signal            <- match.arg(signal)
  estimand          <- match.arg(estimand)
  
  diff.values.plot <- diff.info[["diff.values.plot"]]
  ti <- .extract_treat_info(treat.info)
  treat.type <- ti$treat.type
  all.treat <- all.treat.origin <- NULL
  if (treat.type == "discrete") {
    other.treat <- ti$other.treat
    other.treat.origin <- ti$other.treat.origin
    all.treat <- ti$all.treat
    all.treat.origin <- ti$all.treat.origin
  }
  if (treat.type == "continuous") {
    D.sample <- ti$D.sample
    label.name <- ti$label.name
  }
  
  n <- nrow(data)
  if (is.null(weights)) {
    w <- rep(1, n)
  } else {
    w <- data[[weights]]
  }
  data[["WEIGHTS"]] <- w
  
  # Compute density/histogram of X and (if discrete) by treatment levels
  dens <- .compute_density(data, X, D, weights, treat.type, all.treat, all.treat.origin)
  de <- dens$de
  treat_den <- dens$treat_den
  hists <- .compute_histograms(data, X, D, weights, treat.type, all.treat, all.treat.origin)
  hist.out <- hists$hist.out
  treat.hist <- hists$treat.hist
  if (treat.type == "continuous") {
    de.tr <- NULL
  }
  
  treat.base        <- treat.info[["base"]]
  TE.output.all.list <- list()
  fit_full.list <- list()

  # Set up parallel backend once for all treatment arms
  pcfg <- .parallel_config(B, cores)
  if (pcfg$use_parallel) {
    .setup_parallel(cores)
    on.exit(future::plan(future::sequential), add = TRUE)
  }

  # -------------------------------
  # Branch 1: Discrete D / Binary
  # -------------------------------
  if (treat.type == "discrete") {
    # Report message
    if (verbose) {
      cat(">> Treatment is binary/discrete.\n")
      cat("   Signal = ", signal, "\n", sep = "")
      cat("   Estimand = ", estimand, "\n", sep = "")
      cat("   Basis extension = ", basis.type, "\n", sep = "")
      cat("   Outcome model = ", model.y, "; PS model = ", model.t, "\n", sep = "")
    }
    for (char in other.treat) {
      data_part <- data[data[[D]] %in% c(treat.base, char), ]
      data_part[data_part[[D]] == treat.base, D] <- 0L
      data_part[data_part[[D]] == char,      D] <- 1L
      data_part[, D] <- as.numeric(data_part[[D]])
      
      result <- bootstrapCME(
        data_part,
        Y                    = Y,
        D                    = D,
        X                    = X,
        Z                    = Z,
        FE                   = FE,
        signal               = signal,
        estimand             = estimand,
        B                    = B,
        alpha                = alpha,
        outcome_model_type   = model.y,
        ps_model_type        = model.t,
        basis_type           = basis.type,
        include_interactions = include.interactions,
        poly_degree          = poly.degree,
        spline_df            = spline.df,
        spline_degree        = spline.degree,
        lambda_seq           = lambda.seq,
        reduce.dimension     = reduce.dimension,
        best_span            = best.span,
        x.eval               = x.eval,
        neval = neval,
        CI = CI,
        cores = cores,
        parallel_ready       = pcfg$use_parallel,
        verbose              = verbose
      )

      TE.output.all <- data.frame(result$results, check.names = FALSE)
      TE.output.all.list[[ other.treat.origin[char] ]] <- TE.output.all
      fit_full.list[[ other.treat.origin[char] ]] <- result$fit_full
    }
  }
  
  # ---------------------------------
  # Branch 2: Continuous D 
  # ---------------------------------
  else if (treat.type == "continuous") {
    # Report message
    if (verbose) {
      cat(">> Treatment is continuous.\n")
      cat("   Outcome model = ", model.y, "; Treatment model = ", model.t, "\n", sep = "")
      cat("   Basis extension = ", basis.type, "\n", sep = "")
    }
    result <- bootstrapCME_PLR(
      data                  = data,
      Y                     = Y,
      D                     = D,
      X                     = X,
      Z                     = Z,
      FE                    = FE,
      B                     = B,
      alpha                 = alpha,
      outcome_model_type    = model.y,
      treatment_model_type  = model.t,
      basis_type            = basis.type,
      include_interactions  = include.interactions,
      poly_degree           = poly.degree,
      spline_df             = spline.df,
      spline_degree         = spline.degree,
      lambda_seq            = lambda.seq,
      reduce.dimension      = reduce.dimension,
      bw                    = bw,
      x.eval                = x.eval,
      neval = neval,
      CI = CI,
      cores = cores,
      parallel_ready        = pcfg$use_parallel,
      verbose               = verbose
    )
    TE.output.all <- data.frame(result$results, check.names = FALSE)
    
    # replicate same plot for each D.sample label
    label <- label.name[1]
    TE.output.all.list[[ label ]] <- TE.output.all
    fit_full.list[[label]] <- result$fit_full
  }
  
  # -------------------------------
  # Assemble final output list
  # -------------------------------
  if (treat.type == "discrete") {
    final.output <- list(
      diff.info  = diff.info,
      treat.info = treat.info,
      est.lasso    = TE.output.all.list,
      Xlabel     = Xlabel,
      Dlabel     = Dlabel,
      Ylabel     = Ylabel,
      de         = de,
      hist.out   = hist.out,
      de.tr      = treat_den,
      count.tr   = treat.hist,
      fit_full = fit_full.list,
      estimator  = "lasso"
    )
  } else if (treat.type == "continuous") {
    final.output <- list(
      diff.info  = diff.info,
      treat.info = treat.info,
      est.lasso    = TE.output.all.list,
      Xlabel     = Xlabel,
      Dlabel     = Dlabel,
      Ylabel     = Ylabel,
      de         = de,
      hist.out   = hist.out,
      de.tr      = de.tr,
      count.tr   = NULL,
      fit_full = fit_full.list,
      estimator  = "lasso"
    )
  }
  
  # -------------------------------
  # Plotting (unchanged)
  # -------------------------------
  if (figure) {
    class(final.output) <- "interflex"
    figure.output <- plot.interflex(
      x             = final.output,
      order         = order,
      subtitles     = subtitles,
      show.subtitles= show.subtitles,
      CI            = CI,
      diff.values   = diff.values.plot,
      Xdistr        = Xdistr,
      main          = main,
      Ylabel        = Ylabel,
      Dlabel        = Dlabel,
      Xlabel        = Xlabel,
      xlab          = xlab,
      ylab          = ylab,
      xlim          = xlim,
      ylim          = ylim,
      theme.bw      = theme.bw,
      show.grid     = show.grid,
      cex.main      = cex.main,
      cex.sub       = cex.sub,
      cex.lab       = cex.lab,
      cex.axis      = cex.axis,
      interval      = interval,
      file          = file,
      ncols         = ncols,
      pool          = pool,
      legend.title  = legend.title,
      color         = color,
      show.all      = show.all,
      show.uniform.CI = show.uniform.CI,
      scale         = scale,
      height        = height,
      width         = width
    )
    final.output <- c(final.output, list(figure = figure.output))
  }
  
  class(final.output) <- "interflex"
  return(final.output)
}
