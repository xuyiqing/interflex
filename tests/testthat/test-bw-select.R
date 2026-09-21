# bw.select / bw.se.mult / bw.ess.min: the bandwidth selection rules of the
# kernel estimator, exposed through interflex() and passed through to the
# kernel smoothing step of the lasso estimator with a continuous treatment.

.bw_select_data <- function(n = 300) {
  set.seed(20260920)
  x <- rnorm(n, 3, 1)
  d <- rnorm(n, 3, 1)
  z <- rnorm(n, 3, 1)
  y <- 5 - 4 * x - 9 * d + 3 * x * d + z + 2 * rnorm(n)
  data.frame(Y = y, D = d, X = x, Z1 = z)
}

test_that("kernel estimator: cv.1se never picks a smaller bandwidth than cv.min on the same folds", {
  skip_on_cran()
  s <- .bw_select_data()
  grid <- c(0.5, 1, 2)
  set.seed(1)
  out.min <- suppressMessages(interflex(
    estimator = "kernel", data = s, Y = "Y", D = "D", X = "X", Z = "Z1",
    grid = grid, kfold = 5, CI = FALSE, parallel = FALSE, figure = FALSE, verbose = FALSE
  ))
  set.seed(1)
  out.1se <- suppressMessages(interflex(
    estimator = "kernel", data = s, Y = "Y", D = "D", X = "X", Z = "Z1",
    grid = grid, kfold = 5, bw.select = "cv.1se", CI = FALSE, parallel = FALSE,
    figure = FALSE, verbose = FALSE
  ))
  expect_equal(out.min$bw.select, "cv.min")
  expect_equal(out.1se$bw.select, "cv.1se")
  expect_true(out.min$bw %in% grid)
  expect_true(out.1se$bw %in% grid)
  expect_gte(out.1se$bw, out.min$bw)
  cv <- as.data.frame(out.1se$CV.output)
  expect_equal(cv$bw[cv$CV.selected], out.1se$bw)
})

test_that("kernel estimator: ess-guarded rules run and record their settings", {
  skip_on_cran()
  s <- .bw_select_data()
  set.seed(2)
  out <- suppressMessages(suppressWarnings(interflex(
    estimator = "kernel", data = s, Y = "Y", D = "D", X = "X", Z = "Z1",
    grid = c(0.5, 1, 2), kfold = 5, bw.select = "cv.1se.ess", bw.se.mult = 2, bw.ess.min = 40,
    CI = FALSE, parallel = FALSE, figure = FALSE, verbose = FALSE
  )))
  expect_equal(out$bw.select, "cv.1se.ess")
  expect_equal(out$bw.selection.settings$bw.se.mult, 2)
  expect_equal(out$bw.selection.settings$bw.ess.min, 40)
  expect_true(is.finite(out$bw) && out$bw > 0)
})

test_that("invalid selector settings are rejected up front", {
  skip_on_cran()
  s <- .bw_select_data()
  expect_error(
    interflex(estimator = "kernel", data = s, Y = "Y", D = "D", X = "X", bw.select = "banana",
              CI = FALSE, parallel = FALSE, figure = FALSE, verbose = FALSE),
    "bw.select"
  )
  expect_error(
    interflex(estimator = "kernel", data = s, Y = "Y", D = "D", X = "X", bw.se.mult = -1,
              CI = FALSE, parallel = FALSE, figure = FALSE, verbose = FALSE),
    "bw.se.mult"
  )
  expect_error(
    interflex(estimator = "kernel", data = s, Y = "Y", D = "D", X = "X", bw.ess.min = 0,
              CI = FALSE, parallel = FALSE, figure = FALSE, verbose = FALSE),
    "bw.ess.min"
  )
})

test_that("lasso estimator with continuous treatment: bw and bw.select reach the kernel smoother", {
  skip_on_cran()
  s <- .bw_select_data()
  set.seed(3)
  out.1se <- suppressMessages(suppressWarnings(interflex(
    estimator = "lasso", data = s, Y = "Y", D = "D", X = "X", Z = "Z1",
    bw.select = "cv.1se", CI = FALSE, figure = FALSE
  )))
  expect_s3_class(out.1se, "interflex")
  expect_equal(out.1se$bw.select, "cv.1se")
  expect_equal(out.1se$fit_full[[1]]$bw.select, "cv.1se")
  expect_true(is.finite(out.1se$bw) && out.1se$bw > 0)

  set.seed(3)
  out.fixed <- suppressMessages(suppressWarnings(interflex(
    estimator = "lasso", data = s, Y = "Y", D = "D", X = "X", Z = "Z1",
    bw = 1, CI = FALSE, figure = FALSE
  )))
  expect_equal(out.fixed$bw, 1)
  expect_equal(out.fixed$fit_full[[1]]$bw, 1)
})
