# PAD-001 regression: prediction grid must be restricted to user-supplied
# xlim for continuous-treatment plots, so the visible curve actually ends at
# (lo, hi) and the .pad_xlim(mult=0.04) coord_cartesian pad becomes visible
# whitespace. See test-spec.md §1 (B1-B6), §2-§3.
#
# Pipeline isolation: this file does NOT depend on rendering, quarto, or
# snapshots. It builds interflex objects on bundled toy datasets, calls
# plot.interflex(..., show.all = TRUE), and inspects ggplot_build() output.

skip_if_not_installed("ggplot2")

library(interflex)
library(ggplot2)

# --- Helpers (mirror those in test-plot-limits.R) ----------------------------

load_toy <- function(name) {
    e <- new.env()
    ok <- tryCatch({ data(list = name, envir = e); TRUE },
                   error = function(e) FALSE,
                   warning = function(w) FALSE)
    if (!ok || !exists(name, envir = e)) {
        rdata <- system.file("data", "interflex.RData", package = "interflex")
        if (!nzchar(rdata)) rdata <- file.path("..", "..", "data", "interflex.RData")
        if (file.exists(rdata)) load(rdata, envir = e)
    }
    if (!exists(name, envir = e)) skip(paste("toy dataset", name, "not available"))
    get(name, envir = e)
}

get_inner_ggplots <- function(out, ...) {
    # plot.interflex is not exported from the namespace, so resolve it via :::
    # which works whether the package is installed or loaded via devtools.
    plot_fn <- getFromNamespace("plot.interflex", "interflex")
    pg <- plot_fn(out, show.all = TRUE, ...)
    if (inherits(pg, "ggplot")) return(list(pg))
    stopifnot(is.list(pg), length(pg) >= 1L)
    for (g in pg) stopifnot(inherits(g, "ggplot"))
    pg
}

coord_limits <- function(p) {
    co <- p$coordinates
    stopifnot(inherits(co, "CoordCartesian"))
    list(x = co$limits$x, y = co$limits$y)
}

# Extract per-layer x range for GeomLine / GeomPath / GeomRibbon layers
# from a built ggplot. Returns a list of named lists.
layer_x_range <- function(g) {
    b <- ggplot2::ggplot_build(g)
    out <- list()
    for (i in seq_along(g$layers)) {
        l <- g$layers[[i]]
        df <- b$data[[i]]
        if (is.null(df)) next
        geom_class <- class(l$geom)[1]
        if (!geom_class %in% c("GeomLine", "GeomPath", "GeomRibbon")) next
        xs <- c()
        if ("x" %in% names(df))    xs <- c(xs, df$x)
        if ("xmin" %in% names(df)) xs <- c(xs, df$xmin)
        if ("xmax" %in% names(df)) xs <- c(xs, df$xmax)
        xs <- xs[is.finite(xs)]
        if (length(xs) == 0) next
        out[[length(out) + 1L]] <- list(
            i = i, geom = geom_class,
            min = min(xs), max = max(xs)
        )
    }
    out
}

# --- Per-estimator scenarios -------------------------------------------------
#
# Toy dataset choice: use s2 (continuous D, 200 rows, X in [-0.4, 6.2]). Pick
# xlim = c(0, 5) so both endpoints are strictly interior to the X range and a
# pre-fix run (which spans the full data range) would FAIL B2 on the high side
# (data x extends to ~6.2 > 5) AND B1 (data x reaches ~-0.4 < 0).

CONT_DS    <- "s2"
LO         <- 0
HI         <- 5
SPAN       <- HI - LO                  # 5
EPS        <- 1e-9 * max(1, SPAN)      # ~5e-9
PAD        <- 0.04 * SPAN              # 0.20
EXP_COORD  <- c(LO - PAD, HI + PAD)    # c(-0.20, 5.20)

assert_xlim_grid <- function(g, label) {
    layers <- layer_x_range(g)
    expect_true(length(layers) >= 1L)
    for (lr in layers) {
        ctx <- paste(label, "layer", lr$i, lr$geom)
        # B1: layer x stays inside [lo, hi]
        expect_gte(lr$min, LO - EPS, label = paste(ctx, "min"))
        expect_lte(lr$max, HI + EPS, label = paste(ctx, "max"))
        # B2: layer x reaches the supplied endpoints
        expect_lte(lr$min, LO + EPS, label = paste(ctx, "min reaches lo"))
        expect_gte(lr$max, HI - EPS, label = paste(ctx, "max reaches hi"))
    }
    # B3: coord pad is exactly 4% of span on each side
    lims <- coord_limits(g)
    expect_false(is.null(lims$x))
    expect_equal(lims$x[1], EXP_COORD[1], tolerance = 1e-9)
    expect_equal(lims$x[2], EXP_COORD[2], tolerance = 1e-9)
}

assert_default_xlim_null <- function(g, label) {
    lims <- coord_limits(g)
    expect_null(lims$x, label = paste(label, "default coord$xlim"))
}

# --- Linear ------------------------------------------------------------------

test_that("PAD-001 linear: layer x grid is restricted to user xlim", {
  skip_on_cran()
    s <- load_toy(CONT_DS)
    out <- interflex(estimator = "linear", data = s,
                     Y = "Y", D = "D", X = "X",
                     treat.type = "continuous", na.rm = TRUE,
                     xlim = c(LO, HI))
    inner <- get_inner_ggplots(out)
    expect_true(length(inner) >= 1L)
    for (g in inner) assert_xlim_grid(g, "linear")
})

test_that("PAD-001 linear: default (no xlim) leaves coord$xlim NULL", {
  skip_on_cran()
    s <- load_toy(CONT_DS)
    out <- interflex(estimator = "linear", data = s,
                     Y = "Y", D = "D", X = "X",
                     treat.type = "continuous", na.rm = TRUE)
    inner <- get_inner_ggplots(out)
    for (g in inner) assert_default_xlim_null(g, "linear default")
})

# --- Kernel ------------------------------------------------------------------

test_that("PAD-001 kernel: layer x grid is restricted to user xlim", {
  skip_on_cran()
    s <- load_toy(CONT_DS)
    out <- tryCatch(
        interflex(estimator = "kernel", data = s,
                  Y = "Y", D = "D", X = "X",
                  treat.type = "continuous", na.rm = TRUE,
                  nboots = 2, vartype = "delta",
                  xlim = c(LO, HI)),
        error = function(e) { skip(paste("kernel build failed:", conditionMessage(e))) })
    inner <- get_inner_ggplots(out)
    expect_true(length(inner) >= 1L)
    for (g in inner) assert_xlim_grid(g, "kernel")
})

test_that("PAD-001 kernel: default (no xlim) leaves coord$xlim NULL", {
  skip_on_cran()
    s <- load_toy(CONT_DS)
    out <- tryCatch(
        interflex(estimator = "kernel", data = s,
                  Y = "Y", D = "D", X = "X",
                  treat.type = "continuous", na.rm = TRUE,
                  nboots = 2, vartype = "delta"),
        error = function(e) { skip(paste("kernel build failed:", conditionMessage(e))) })
    inner <- get_inner_ggplots(out)
    for (g in inner) assert_default_xlim_null(g, "kernel default")
})

# --- Binning -----------------------------------------------------------------

test_that("PAD-001 binning: layer x grid is restricted to user xlim", {
  skip_on_cran()
    s <- load_toy(CONT_DS)
    out <- tryCatch(
        interflex(estimator = "binning", data = s,
                  Y = "Y", D = "D", X = "X",
                  treat.type = "continuous", na.rm = TRUE,
                  nbins = 3, nboots = 2,
                  xlim = c(LO, HI)),
        error = function(e) { skip(paste("binning build failed:", conditionMessage(e))) })
    inner <- get_inner_ggplots(out)
    expect_true(length(inner) >= 1L)
    for (g in inner) assert_xlim_grid(g, "binning")
})

test_that("PAD-001 binning: default (no xlim) leaves coord$xlim NULL", {
  skip_on_cran()
    s <- load_toy(CONT_DS)
    out <- tryCatch(
        interflex(estimator = "binning", data = s,
                  Y = "Y", D = "D", X = "X",
                  treat.type = "continuous", na.rm = TRUE,
                  nbins = 3, nboots = 2),
        error = function(e) { skip(paste("binning build failed:", conditionMessage(e))) })
    inner <- get_inner_ggplots(out)
    for (g in inner) assert_default_xlim_null(g, "binning default")
})

# --- Lasso (PO-lasso, continuous PLR) ----------------------------------------

test_that("PAD-001 lasso: layer x grid is restricted to user xlim", {
  skip_on_cran()
    skip_if_not_installed("glmnet")
    s <- load_toy(CONT_DS)
    out <- suppressWarnings(suppressMessages(
        interflex(estimator = "lasso", data = s,
                  Y = "Y", D = "D", X = "X",
                  treat.type = "continuous", na.rm = TRUE,
                  nboots = 2,
                  xlim = c(LO, HI))))
    inner <- get_inner_ggplots(out)
    expect_true(length(inner) >= 1L)
    for (g in inner) assert_xlim_grid(g, "lasso")
})

test_that("PAD-001 lasso: default (no xlim) leaves coord$xlim NULL", {
  skip_on_cran()
    skip_if_not_installed("glmnet")
    s <- load_toy(CONT_DS)
    out <- suppressWarnings(suppressMessages(
        interflex(estimator = "lasso", data = s,
                  Y = "Y", D = "D", X = "X",
                  treat.type = "continuous", na.rm = TRUE,
                  nboots = 2)))
    inner <- get_inner_ggplots(out)
    for (g in inner) assert_default_xlim_null(g, "lasso default")
})

# --- DML ---------------------------------------------------------------------

test_that("PAD-001 DML: layer x grid is restricted to user xlim", {
  skip_on_cran()
    skip_if_not_installed("DoubleML")
    skip_if_not_installed("mlr3")
    s <- load_toy(CONT_DS)
    out <- suppressWarnings(suppressMessages(
        interflex(estimator = "dml", data = s,
                  Y = "Y", D = "D", X = "X",
                  treat.type = "continuous", na.rm = TRUE,
                  model.y = "linear", model.t = "linear",
                  xlim = c(LO, HI))))
    inner <- get_inner_ggplots(out)
    expect_true(length(inner) >= 1L)
    for (g in inner) assert_xlim_grid(g, "DML")
})

# --- §2.6 .pad_xlim default mult is still 0.04 -------------------------------

test_that(".pad_xlim default mult is 0.04 (PAD-001 invariant)", {
  skip_on_cran()
    expect_equal(interflex:::.pad_xlim(c(0, 10)), c(-0.4, 10.4))
})

# --- §2.5 Discrete-treatment plots are unaffected (B6) -----------------------

test_that("PAD-001 discrete-treatment plot with explicit xlim still pads", {
  skip_on_cran()
    # s1: D in {0,1}, X in [-0.4, 6.2]. Pick interior xlim so the test
    # exercises an explicit user xlim on a discrete plot.
    s1 <- load_toy("s1")
    user_xlim <- c(0, 5)
    span <- diff(user_xlim)
    out <- interflex(estimator = "linear", data = s1,
                     Y = "Y", D = "D", X = "X", na.rm = TRUE,
                     xlim = user_xlim)
    inner <- get_inner_ggplots(out)
    for (g in inner) {
        lims <- coord_limits(g)
        expect_false(is.null(lims$x))
        expect_equal(lims$x[1], user_xlim[1] - 0.04 * span, tolerance = 1e-9)
        expect_equal(lims$x[2], user_xlim[2] + 0.04 * span, tolerance = 1e-9)
    }
})

# --- §3 Edge cases -----------------------------------------------------------

# --- §2.8 PASS 2 — Path 2 (plot-time xlim) ----------------------------------
#
# Path 2 fits interflex WITHOUT xlim, then calls plot(out, xlim = c(lo, hi)).
# B7: B1 + B3 must hold on g2. B2 is relaxed to "within one grid step".
# B8: idempotence — Path1+Path2 with same xlim yields identical layer x range
#     to Path1-only; narrower plot-time xlim shrinks further.
# B9: default plot(out) (no xlim anywhere) yields NULL coord x and full range.
# B10: discrete unchanged under plot-time xlim.

plot_fn_pi <- function() getFromNamespace("plot.interflex", "interflex")

# Relaxed B2 tolerance: max gap is data_range / (neval - 1), default neval=50
relaxed_gap <- function(out_obj, span_xlim) {
    # Try to find the moderator data range on the fitted object.
    xvec <- NULL
    if (!is.null(out_obj$data) && !is.null(out_obj$X) &&
        is.character(out_obj$X) && out_obj$X %in% names(out_obj$data)) {
        xvec <- out_obj$data[[out_obj$X]]
    }
    if (is.null(xvec) && !is.null(out_obj$Xvec)) xvec <- out_obj$Xvec
    if (is.null(xvec)) return(span_xlim)   # permissive fallback
    xvec <- xvec[is.finite(xvec)]
    if (length(xvec) < 2L) return(span_xlim)
    diff(range(xvec)) / (50 - 1)
}

assert_xlim_grid_path2 <- function(g, label, out_obj) {
    layers <- layer_x_range(g)
    expect_true(length(layers) >= 1L)
    gap <- relaxed_gap(out_obj, SPAN)
    for (lr in layers) {
        ctx <- paste("Path2", label, "layer", lr$i, lr$geom)
        # B1: stays inside [lo, hi]
        expect_gte(lr$min, LO - EPS, label = paste(ctx, "min<lo"))
        expect_lte(lr$max, HI + EPS, label = paste(ctx, "max>hi"))
        # B2 relaxed: reaches within one grid step
        expect_lte(lr$min, LO + gap + EPS,
                   label = paste(ctx, "min did not reach near lo"))
        expect_gte(lr$max, HI - gap - EPS,
                   label = paste(ctx, "max did not reach near hi"))
    }
    # B3: coord pad still exactly 4% of span
    lims <- coord_limits(g)
    expect_false(is.null(lims$x))
    expect_equal(lims$x[1], EXP_COORD[1], tolerance = 1e-9)
    expect_equal(lims$x[2], EXP_COORD[2], tolerance = 1e-9)
}

# B7 — Path 2 per-estimator -------------------------------------------------

test_that("PAD-001 B7 linear Path 2: plot-time xlim restricts layer grid", {
  skip_on_cran()
    s <- load_toy(CONT_DS)
    out <- interflex(estimator = "linear", data = s,
                     Y = "Y", D = "D", X = "X",
                     treat.type = "continuous", na.rm = TRUE)
    pf <- plot_fn_pi()
    g2 <- pf(out, show.all = TRUE, xlim = c(LO, HI))
    inner <- if (inherits(g2, "ggplot")) list(g2) else g2
    expect_true(length(inner) >= 1L)
    for (g in inner) assert_xlim_grid_path2(g, "linear", out)
})

test_that("PAD-001 B7 kernel Path 2: plot-time xlim restricts layer grid", {
  skip_on_cran()
    s <- load_toy(CONT_DS)
    out <- tryCatch(
        interflex(estimator = "kernel", data = s,
                  Y = "Y", D = "D", X = "X",
                  treat.type = "continuous", na.rm = TRUE,
                  nboots = 2, vartype = "delta"),
        error = function(e) { skip(paste("kernel build failed:", conditionMessage(e))) })
    pf <- plot_fn_pi()
    g2 <- pf(out, show.all = TRUE, xlim = c(LO, HI))
    inner <- if (inherits(g2, "ggplot")) list(g2) else g2
    expect_true(length(inner) >= 1L)
    for (g in inner) assert_xlim_grid_path2(g, "kernel", out)
})

test_that("PAD-001 B7 binning Path 2: plot-time xlim restricts layer grid", {
  skip_on_cran()
    s <- load_toy(CONT_DS)
    out <- tryCatch(
        interflex(estimator = "binning", data = s,
                  Y = "Y", D = "D", X = "X",
                  treat.type = "continuous", na.rm = TRUE,
                  nbins = 3, nboots = 2),
        error = function(e) { skip(paste("binning build failed:", conditionMessage(e))) })
    pf <- plot_fn_pi()
    g2 <- pf(out, show.all = TRUE, xlim = c(LO, HI))
    inner <- if (inherits(g2, "ggplot")) list(g2) else g2
    expect_true(length(inner) >= 1L)
    for (g in inner) assert_xlim_grid_path2(g, "binning", out)
})

test_that("PAD-001 B7 lasso Path 2: plot-time xlim restricts layer grid", {
  skip_on_cran()
    skip_if_not_installed("glmnet")
    s <- load_toy(CONT_DS)
    out <- tryCatch(
        interflex(estimator = "lasso", data = s,
                  Y = "Y", D = "D", X = "X",
                  treat.type = "continuous", na.rm = TRUE),
        error = function(e) { skip(paste("lasso build failed:", conditionMessage(e))) })
    pf <- plot_fn_pi()
    g2 <- tryCatch(pf(out, show.all = TRUE, xlim = c(LO, HI)),
                   error = function(e) { skip(paste("lasso plot failed:", conditionMessage(e))) })
    inner <- if (inherits(g2, "ggplot")) list(g2) else g2
    expect_true(length(inner) >= 1L)
    for (g in inner) assert_xlim_grid_path2(g, "lasso", out)
})

test_that("PAD-001 B7 DML Path 2: plot-time xlim restricts layer grid", {
  skip_on_cran()
    skip_if_not_installed("DoubleML")
    skip_if_not_installed("mlr3")
    s <- load_toy(CONT_DS)
    out <- tryCatch(
        interflex(estimator = "DML", data = s,
                  Y = "Y", D = "D", X = "X",
                  treat.type = "continuous", na.rm = TRUE),
        error = function(e) { skip(paste("DML build failed:", conditionMessage(e))) })
    pf <- plot_fn_pi()
    g2 <- tryCatch(pf(out, show.all = TRUE, xlim = c(LO, HI)),
                   error = function(e) { skip(paste("DML plot failed:", conditionMessage(e))) })
    inner <- if (inherits(g2, "ggplot")) list(g2) else g2
    expect_true(length(inner) >= 1L)
    for (g in inner) assert_xlim_grid_path2(g, "DML", out)
})

# B8 — Idempotence: Path1 and Path1+Path2(same) yield identical layer x ------

test_that("PAD-001 B8 linear: Path 1 + Path 2 idempotence and narrowing", {
  skip_on_cran()
    s <- load_toy(CONT_DS)
    out1 <- interflex(estimator = "linear", data = s,
                      Y = "Y", D = "D", X = "X",
                      treat.type = "continuous", na.rm = TRUE,
                      xlim = c(LO, HI))
    pf <- plot_fn_pi()
    g_a <- pf(out1, show.all = TRUE)
    g_b <- pf(out1, show.all = TRUE, xlim = c(LO, HI))
    g_c <- pf(out1, show.all = TRUE, xlim = c(-0.5, 0.5))
    inner_a <- if (inherits(g_a, "ggplot")) list(g_a) else g_a
    inner_b <- if (inherits(g_b, "ggplot")) list(g_b) else g_b
    inner_c <- if (inherits(g_c, "ggplot")) list(g_c) else g_c
    expect_equal(length(inner_a), length(inner_b))
    # g_a and g_b must be bit-identical on layer x ranges
    for (i in seq_along(inner_a)) {
        la <- layer_x_range(inner_a[[i]])
        lb <- layer_x_range(inner_b[[i]])
        expect_equal(length(la), length(lb))
        for (j in seq_along(la)) {
            expect_equal(la[[j]]$min, lb[[j]]$min, tolerance = 1e-12)
            expect_equal(la[[j]]$max, lb[[j]]$max, tolerance = 1e-12)
        }
    }
    # g_c narrows further to [-0.5, 0.5]
    for (g in inner_c) {
        layers <- layer_x_range(g)
        expect_true(length(layers) >= 1L)
        for (lr in layers) {
            expect_gte(lr$min, -0.5 - 1e-9)
            expect_lte(lr$max,  0.5 + 1e-9)
        }
        lims <- coord_limits(g)
        expect_equal(lims$x[1], -0.5 - 0.04, tolerance = 1e-9)
        expect_equal(lims$x[2],  0.5 + 0.04, tolerance = 1e-9)
    }
})

# B9 — default plot(out) with no xlim anywhere -------------------------------

test_that("PAD-001 B9 linear: default plot(out) leaves coord$xlim NULL and full range", {
  skip_on_cran()
    s <- load_toy(CONT_DS)
    out <- interflex(estimator = "linear", data = s,
                     Y = "Y", D = "D", X = "X",
                     treat.type = "continuous", na.rm = TRUE)
    pf <- plot_fn_pi()
    g <- pf(out, show.all = TRUE)
    inner <- if (inherits(g, "ggplot")) list(g) else g
    for (gi in inner) {
        lims <- coord_limits(gi)
        expect_null(lims$x)
        # Layer x should span ~full data range, not a restricted window.
        layers <- layer_x_range(gi)
        if (length(layers) >= 1L) {
            data_rng <- range(s$X, na.rm = TRUE)
            for (lr in layers) {
                # Allow a tiny slack — full-data grid
                expect_lte(lr$min, data_rng[1] + 1e-6)
                expect_gte(lr$max, data_rng[2] - 1e-6)
            }
        }
    }
})

# B10 — Discrete unchanged under plot-time xlim ------------------------------

test_that("PAD-001 B10 discrete: plot-time xlim still pads, no layer filter", {
  skip_on_cran()
    s1 <- load_toy("s1")
    out <- interflex(estimator = "linear", data = s1,
                     Y = "Y", D = "D", X = "X", na.rm = TRUE)
    pf <- plot_fn_pi()
    g <- pf(out, show.all = TRUE, xlim = c(0, 5))
    inner <- if (inherits(g, "ggplot")) list(g) else g
    for (gi in inner) {
        lims <- coord_limits(gi)
        expect_false(is.null(lims$x))
        expect_equal(lims$x[1], 0 - 0.04 * 5, tolerance = 1e-9)
        expect_equal(lims$x[2], 5 + 0.04 * 5, tolerance = 1e-9)
    }
})

test_that("PAD-001 narrow xlim entirely inside data: linear", {
  skip_on_cran()
    s <- load_toy(CONT_DS)
    narrow_lo <- 1
    narrow_hi <- 2
    out <- interflex(estimator = "linear", data = s,
                     Y = "Y", D = "D", X = "X",
                     treat.type = "continuous", na.rm = TRUE,
                     xlim = c(narrow_lo, narrow_hi))
    inner <- get_inner_ggplots(out)
    span <- narrow_hi - narrow_lo
    eps <- 1e-9 * max(1, span)
    for (g in inner) {
        layers <- layer_x_range(g)
        expect_true(length(layers) >= 1L)
        for (lr in layers) {
            expect_gte(lr$min, narrow_lo - eps)
            expect_lte(lr$max, narrow_hi + eps)
            expect_lte(lr$min, narrow_lo + eps)
            expect_gte(lr$max, narrow_hi - eps)
        }
    }
})
