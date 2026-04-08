## Minimal fast smoke tests that run on CRAN.
##
## Everything heavy (bootstrap inference, kernel bandwidth CV, DML, lasso,
## GATE, discrete-outcome GLMs, xlim grid suites) is gated with
## skip_on_cran() in its own file so it runs only on the developer machine
## or on CI. This file exists so CRAN still sees interflex fit, plot, and
## the raw wrapper exercise every estimator-independent code path at least
## once, without adding more than a few seconds to R CMD check.

test_that("linear estimator: minimal fit + figure", {
  data(interflex, envir = environment())
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "linear", data = s1,
              Y = "Y", D = "D", X = "X", Z = "Z1",
              CI = FALSE, vartype = "delta")
  ))
  expect_s3_class(out, "interflex")
  expect_false(is.null(out$figure))
  expect_false(is.null(out$est.lin))
})

test_that("raw estimator: minimal fit + figure", {
  data(interflex, envir = environment())
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "raw", data = s1,
              Y = "Y", D = "D", X = "X", Z = "Z1",
              theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
  expect_identical(out$type, "raw")
  expect_false(is.null(out$figure))
})

test_that("plot.interflex returns a ggplot on a linear fit", {
  data(interflex, envir = environment())
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "linear", data = s1,
              Y = "Y", D = "D", X = "X", Z = "Z1",
              CI = FALSE, vartype = "delta")
  ))
  p <- plot(out, xlim = c(-1, 1), ylim = c(-5, 5))
  expect_true(inherits(p, "ggplot"))
})
