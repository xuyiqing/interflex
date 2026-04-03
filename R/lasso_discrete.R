interflex.lasso_discrete <- function(
  data,
  Y,
  D,
  X,
  treat.info,
  diff.info,
  Z = NULL,
  FE = NULL,

  signal    = c("outcome", "ipw", "aipw"),
  estimand  = c("ATE", "ATT"),

  weights = NULL,
  B       = 200,
  alpha   = 0.05,

  model.y  = "lasso",
  model.t  = "lasso",
  neval = 50,
  basis.type           = c("polynomial", "bspline", "none"),
  include.interactions = FALSE,
  poly.degree          = 2,
  spline.df            = 4,
  spline.degree        = 2,

  lambda.seq           = NULL,
  cores              = 8,

  verbose            = TRUE,
  figure             = TRUE,
  CI                 = NULL,
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
  theme.bw           = FALSE,
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
  width              = 10
) {
  basis.type       <- match.arg(basis.type)
  signal           <- match.arg(signal)
  estimand         <- match.arg(estimand)

  diff.values.plot <- diff.info[["diff.values.plot"]]
  ti <- .extract_treat_info(treat.info)
  treat.type <- ti$treat.type

  # weights
  n <- nrow(data)
  if (is.null(weights)) {
    w <- rep(1, n)
  } else {
    w <- data[[weights]]
  }
  data[["WEIGHTS"]] <- w

  # X-distribution (density and histogram)
  dens <- .compute_density(data, X, D, weights, treat.type)
  de <- dens$de
  hists <- .compute_histograms(data, X, D, weights, treat.type)
  hist.out <- hists$hist.out

  TE.output.all.list <- list()

  # Set up parallel backend once for all treatment arms
  pcfg <- .parallel_config(B, cores)
  if (pcfg$use_parallel) {
    .setup_parallel(cores)
    on.exit(future::plan(future::sequential), add = TRUE)
  }

  # 1) DISCRETE TREATMENT ---------------------------------------------------
  if (treat.type == "discrete") {
    other.treat <- ti$other.treat
    other.treat.origin <- ti$other.treat.origin
    all.treat <- ti$all.treat
    all.treat.origin <- ti$all.treat.origin

    if (verbose) {
      cat(">> Treatment is discrete/binary.\n")
      cat("   Signal     =", signal,   "\n")
      cat("   Estimand   =", estimand, "\n")
      cat("   Basis type =", basis.type,"\n")
      cat("   Outcome model   =", model.y, ";",
          "Propensity model =", model.t, "\n")
    }

    for (char in other.treat) {
      data_part <- subset(data, data[[D]] %in% c(treat.info[["base"]], char))
      data_part[data_part[[D]] == treat.info[["base"]], D] <- 0L
      data_part[data_part[[D]] == char,            D] <- 1L

      result <- bootstrapGTE(
        data                  = data_part,
        Y                     = Y,
        D                     = D,
        X                     = X,
        Z                     = Z,
        FE                    = FE,
        estimand              = estimand,
        signal                = signal,
        B                     = B,
        alpha                 = alpha,
        outcome_model_type    = model.y,
        ps_model_type         = model.t,
        basis_type            = basis.type,
        include_interactions  = include.interactions,
        poly_degree           = poly.degree,
        spline_df             = spline.df,
        spline_degree         = spline.degree,
        lambda_seq            = lambda.seq,
        CI = CI,
        cores = cores,
        parallel_ready        = pcfg$use_parallel,
        verbose               = verbose
      )

      TE.output.all <- data.frame(result$results, check.names = FALSE)
      TE.output.all.list[[ other.treat.origin[char] ]] <- TE.output.all
    }
  }

  # 2) CONTINUOUS TREATMENT -------------------------------------------------
  else if (treat.type == "continuous") {
    D.sample <- ti$D.sample
    label.name <- ti$label.name

    if (verbose) {
      cat(">> Treatment is continuous.\n")
      cat("   Outcome model   =", model.y, ";",
          "Treatment model =", model.t, "\n")
      cat("   Basis type =", basis.type, "\n")
    }

    result <- bootstrapGATE_PLR(
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
      CI = CI,
      cores = cores,
      parallel_ready        = pcfg$use_parallel,
      verbose               = verbose
    )

    TE.output.all <- data.frame(result$results, check.names = FALSE)
    k <- 1
    for (d_ref in D.sample) {
      lbl <- label.name[k]
      TE.output.all.list[[ lbl ]] <- TE.output.all
      k <- k + 1
    }
  }

  # 3) Assemble final output -----------------------------------------------
  final.output <- list(
    diff.info    = diff.info,
    treat.info   = treat.info,
    est.lasso    = TE.output.all.list,
    Xlabel       = Xlabel,
    Dlabel       = Dlabel,
    Ylabel       = Ylabel,
    de           = de,
    hist.out     = hist.out,
    de.tr        = if (treat.type=="discrete") NULL else NULL,
    count.tr     = if (treat.type=="discrete") NULL else NULL,
    estimator    = "lasso"
  )

  # 4) Plotting (same as interflex.lasso) -----------------------------------
  if (figure) {
    class(final.output) <- "interflex"
    fig <- plot.interflex(
      x               = final.output,
      order           = order,
      subtitles       = subtitles,
      show.subtitles  = show.subtitles,
      CI              = CI,
      diff.values     = diff.values.plot,
      Xdistr          = Xdistr,
      main            = main,
      Ylabel          = Ylabel,
      Dlabel          = Dlabel,
      Xlabel          = Xlabel,
      xlab            = xlab,
      ylab            = ylab,
      xlim            = xlim,
      ylim            = ylim,
      theme.bw        = theme.bw,
      show.grid       = show.grid,
      cex.main        = cex.main,
      cex.sub         = cex.sub,
      cex.lab         = cex.lab,
      cex.axis        = cex.axis,
      interval        = interval,
      file            = file,
      ncols           = ncols,
      pool            = pool,
      legend.title    = legend.title,
      color           = color,
      show.all        = show.all,
      scale           = scale,
      height          = height,
      width           = width
    )
    final.output$figure <- fig
  }

  class(final.output) <- "interflex"
  final.output
}
