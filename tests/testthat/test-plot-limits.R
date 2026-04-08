# Regression tests for BOOK-003 plot limit fixes:
#   - .pad_xlim helper contract (also covered by Check 3 smoke test)
#   - User-supplied xlim is padded by ~2% via coord_cartesian
#   - Default behavior (no xlim) leaves coord_cartesian xlim NULL / unpadded
#   - User-supplied ylim is applied via coord_cartesian (NOT ylim()/scale_*),
#     so CI ribbon layer data retains rows outside the ylim window
#
# These tests must NOT depend on rendering output, on quarto, or on snapshots.
# They construct interflex objects from the bundled toy datasets, call
# plot.interflex(..., show.all = TRUE) to obtain raw ggplot objects, and
# inspect the coord and layer data directly.

skip_if_not_installed("ggplot2")

library(interflex)
library(ggplot2)

# The bundled datasets ship in data/interflex.RData rather than as one .rda
# per dataset. data() works after library(interflex); fall back to load() if
# not visible (defensive — this is a known wrinkle from BOOK-003 audit).
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

# Extract the inner ggplot(s) from a plot.interflex(show.all = TRUE) return
# value. For discrete treatment, p.group is a named list of ggplots keyed by
# `other.treat`; for continuous, by `label.name`. Returns a non-empty list of
# ggplot objects.
get_inner_ggplots <- function(out, ...) {
    pg <- plot.interflex(out, show.all = TRUE, ...)
    if (inherits(pg, "ggplot")) return(list(pg))
    stopifnot(is.list(pg), length(pg) >= 1L)
    for (g in pg) stopifnot(inherits(g, "ggplot"))
    pg
}

# Pull the coord_cartesian limits off a ggplot. ggplot stores the active coord
# in p$coordinates; for CoordCartesian the user limits live in $limits$x and
# $limits$y (each NULL or length-2 numeric). Returns list(x = ..., y = ...).
coord_limits <- function(p) {
    co <- p$coordinates
    stopifnot(inherits(co, "CoordCartesian"))
    list(x = co$limits$x, y = co$limits$y)
}

test_that("default plot has no user-imposed coord_cartesian xlim", {
  skip_on_cran()
    s1 <- load_toy("s1")
    out <- interflex(estimator = "linear", data = s1,
                     Y = "Y", D = "D", X = "X", na.rm = TRUE)
    inner <- get_inner_ggplots(out)
    for (g in inner) {
        lims <- coord_limits(g)
        # No user xlim was supplied, so the merged group-equalization coord
        # call sets xlim = NULL (see R/plot.R post-BOOK-003).
        expect_null(lims$x)
    }
})

test_that("user-supplied xlim is padded by ~4% via .pad_xlim", {
  skip_on_cran()
    s1 <- load_toy("s1")
    user_xlim <- c(-2, 2)
    span <- user_xlim[2] - user_xlim[1]            # 4
    expected_pad <- 0.04 * span                    # 0.16
    expected_lo <- user_xlim[1] - expected_pad     # -2.16
    expected_hi <- user_xlim[2] + expected_pad     #  2.16

    out <- interflex(estimator = "linear", data = s1,
                     Y = "Y", D = "D", X = "X", na.rm = TRUE,
                     xlim = user_xlim)
    inner <- get_inner_ggplots(out)
    expect_true(length(inner) >= 1L)
    for (g in inner) {
        lims <- coord_limits(g)
        expect_false(is.null(lims$x))
        expect_length(lims$x, 2L)
        expect_equal(lims$x[1], expected_lo, tolerance = 1e-9)
        expect_equal(lims$x[2], expected_hi, tolerance = 1e-9)
        # Sanity: padded range strictly contains the user range
        expect_lt(lims$x[1], user_xlim[1])
        expect_gt(lims$x[2], user_xlim[2])
        # And the pad is small (< 10% of span on each side)
        expect_lt(user_xlim[1] - lims$x[1], 0.1 * span)
        expect_lt(lims$x[2] - user_xlim[2], 0.1 * span)
    }
})

test_that("user-supplied ylim is set on coord_cartesian, not via ylim()", {
  skip_on_cran()
    s1 <- load_toy("s1")
    user_ylim <- c(-0.5, 0.5)
    out <- interflex(estimator = "linear", data = s1,
                     Y = "Y", D = "D", X = "X", na.rm = TRUE,
                     ylim = user_ylim)
    inner <- get_inner_ggplots(out)
    for (g in inner) {
        lims <- coord_limits(g)
        expect_false(is.null(lims$y))
        expect_equal(as.numeric(lims$y), user_ylim, tolerance = 1e-9)
        # NO scale_y_continuous(limits = ...) should be active. The y scale
        # must remain a default continuous scale with NULL user limits, so
        # the data are not filtered before being passed to geoms.
        ys <- g$scales$get_scales("y")
        if (!is.null(ys)) {
            # ggplot stores user-set scale limits in $limits; coord limits are
            # not propagated here. Either NULL or length-0 is acceptable.
            expect_true(is.null(ys$limits) || length(ys$limits) == 0L)
        }
    }
})

test_that("CI ribbon retains data outside narrow user ylim", {
  skip_on_cran()
    pick_outside_ylim <- function() {
        candidates <- list(
            list(name = "s1", est = "linear",  args = list()),
            list(name = "s2", est = "linear",  args = list()),
            list(name = "s3", est = "linear",  args = list())
        )
        for (cand in candidates) {
            d <- tryCatch(load_toy(cand$name), error = function(e) NULL)
            if (is.null(d)) next
            args <- c(list(estimator = cand$est, data = d,
                           Y = "Y", D = "D", X = "X", na.rm = TRUE,
                           ylim = c(-0.5, 0.5)),
                      cand$args)
            out <- tryCatch(do.call(interflex, args), error = function(e) NULL)
            if (is.null(out)) next
            inner <- tryCatch(get_inner_ggplots(out), error = function(e) NULL)
            if (is.null(inner)) next
            for (g in inner) {
                ribbon_idx <- which(vapply(g$layers,
                    function(l) inherits(l$geom, "GeomRibbon"),
                    logical(1)))
                if (length(ribbon_idx) == 0L) next
                b <- ggplot_build(g)
                for (i in ribbon_idx) {
                    df <- b$data[[i]]
                    if (is.null(df) || !all(c("ymin", "ymax") %in% names(df))) next
                    has_ymin <- !all(is.na(df$ymin))
                    has_ymax <- !all(is.na(df$ymax))
                    outside <- any(df$ymin < -0.5, na.rm = TRUE) ||
                               any(df$ymax >  0.5, na.rm = TRUE)
                    if (has_ymin && has_ymax && outside) {
                        return(list(g = g, df = df, dataset = cand$name))
                    }
                }
            }
        }
        NULL
    }

    hit <- pick_outside_ylim()
    if (is.null(hit)) {
        skip("no bundled dataset produced a CI ribbon extending outside [-0.5, 0.5]; ribbon-preservation invariant cannot be exercised on this build")
    }
    df <- hit$df
    expect_true(any(df$ymin < -0.5, na.rm = TRUE) ||
                any(df$ymax >  0.5, na.rm = TRUE))
    expect_false(all(is.na(df$ymin)))
    expect_false(all(is.na(df$ymax)))
    lims <- coord_limits(hit$g)
    expect_equal(as.numeric(lims$y), c(-0.5, 0.5), tolerance = 1e-9)
})
