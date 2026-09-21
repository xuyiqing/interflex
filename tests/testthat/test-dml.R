test_that("DML estimator works with binary treatment (linear + linear)", {
  skip_on_cran()
  data(interflex, envir = environment())
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "DML", data = app_hma2015,
              Y = "totangry", D = "threat", X = "pidentity",
              Z = c("issuestr2", "knowledge", "educ", "male", "age10"),
              model.y = "linear", model.t = "linear",
              theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
  expect_true(!is.null(out$figure))
  expect_true(!is.null(out$est.dml))
})

test_that("DML estimator works with rf models", {
  skip_on_cran()
  skip_if_not_installed("ranger")
  data(interflex, envir = environment())
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "DML", data = app_hma2015,
              Y = "totangry", D = "threat", X = "pidentity",
              Z = c("issuestr2", "knowledge", "educ", "male", "age10"),
              model.y = "rf", model.t = "rf",
              theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
  expect_true(!is.null(out$est.dml))
})

test_that("DML estimator works with continuous treatment", {
  skip_on_cran()
  set.seed(1234)
  n <- 200
  d2 <- rnorm(n, 3, 1)
  x <- rnorm(n, 3, 1)
  z <- rnorm(n, 3, 1)
  e <- rnorm(n, 0, 1)
  y2 <- 5 - 4 * x - 9 * d2 + 3 * x * d2 + 1 * z + 2 * e
  s2 <- data.frame(Y = y2, D = d2, X = x, Z1 = z)

  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "DML", data = s2,
              Y = "Y", D = "D", X = "X", Z = "Z1",
              treat.type = "continuous",
              model.y = "linear", model.t = "linear",
              theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
  expect_true(!is.null(out$est.dml))
})

test_that("DML estimator result structure is correct", {
  skip_on_cran()
  data(interflex, envir = environment())
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "DML", data = app_hma2015,
              Y = "totangry", D = "threat", X = "pidentity",
              Z = c("issuestr2", "knowledge", "educ", "male", "age10"),
              model.y = "linear", model.t = "linear",
              theme.bw = TRUE)
  ))
  expect_true(is.list(out$est.dml))
  expect_true(length(out$est.dml) >= 1)
  est_df <- out$est.dml[[1]]
  expect_true(is.data.frame(est_df) || is.matrix(est_df))
})

test_that("boosting learners default to a silent LightGBM (verbose = -1)", {
  skip_on_cran()
  skip_if_not_installed("lightgbm")
  skip_if_not_installed("mlr3extralearners")
  expect_equal(interflex:::.map_boosting_params(NULL)$verbose, -1L)
  expect_equal(interflex:::.map_boosting_params(list())$verbose, -1L)
  expect_equal(interflex:::.map_boosting_params(list(max_iter = 10))$verbose, -1L)
  expect_equal(interflex:::.map_boosting_params(list(verbose = 1L))$verbose, 1L)
})

test_that("DML with hgb learners prints no LightGBM warnings", {
  skip_on_cran()
  skip_if_not_installed("lightgbm")
  skip_if_not_installed("mlr3extralearners")
  set.seed(20260921)
  n <- 300
  x <- rnorm(n); z <- rnorm(n)
  d <- rbinom(n, 1, plogis(0.5 * x))
  y <- 1 + x + 2 * d + d * x + z + rnorm(n)
  dat <- data.frame(Y = y, D = d, X = x, Z1 = z)
  msgs <- character(0)
  out <- capture.output(
    withCallingHandlers(
      fit <- suppressWarnings(interflex(estimator = "dml", data = dat, Y = "Y", D = "D", X = "X", Z = "Z1",
                                        model.y = "hgb", model.t = "hgb", CI = FALSE, figure = FALSE)),
      message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
    ),
    type = "output"
  )
  expect_s3_class(fit, "interflex")
  expect_false(any(grepl("[LightGBM]", c(out, msgs), fixed = TRUE)))
})
