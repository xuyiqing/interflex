test_that("lasso (AIPW) estimator works with binary treatment", {
  skip_on_cran()
  data(interflex, envir = environment())
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "lasso", data = app_hma2015,
              Y = "totangry", D = "threat", X = "pidentity",
              Z = c("issuestr2", "knowledge", "educ", "male", "age10"),
              nboots = 20, theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
  expect_true(!is.null(out$figure))
  expect_true(!is.null(out$est.lasso))
})

test_that("lasso estimator works with signal = outcome", {
  skip_on_cran()
  data(interflex, envir = environment())
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "lasso", data = app_hma2015,
              Y = "totangry", D = "threat", X = "pidentity",
              Z = c("issuestr2", "knowledge", "educ", "male", "age10"),
              signal = "outcome", nboots = 20, theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
})

test_that("lasso estimator works with signal = ipw", {
  skip_on_cran()
  data(interflex, envir = environment())
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "lasso", data = app_hma2015,
              Y = "totangry", D = "threat", X = "pidentity",
              Z = c("issuestr2", "knowledge", "educ", "male", "age10"),
              signal = "ipw", nboots = 20, theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
})

test_that("lasso estimator works with estimand = ATT", {
  skip_on_cran()
  data(interflex, envir = environment())
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "lasso", data = app_hma2015,
              Y = "totangry", D = "threat", X = "pidentity",
              Z = c("issuestr2", "knowledge", "educ", "male", "age10"),
              estimand = "ATT", nboots = 20, theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
})

test_that("lasso estimator works with reduce.dimension = bspline", {
  skip_on_cran()
  data(interflex, envir = environment())
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "lasso", data = app_hma2015,
              Y = "totangry", D = "threat", X = "pidentity",
              Z = c("issuestr2", "knowledge", "educ", "male", "age10"),
              reduce.dimension = "bspline", nboots = 20, theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
})

test_that("lasso estimator works with basis.type = polynomial", {
  skip_on_cran()
  data(interflex, envir = environment())
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "lasso", data = app_hma2015,
              Y = "totangry", D = "threat", X = "pidentity",
              Z = c("issuestr2", "knowledge", "educ", "male", "age10"),
              basis.type = "polynomial", nboots = 20, theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
})

test_that("lasso estimator works with basis.type = none", {
  skip_on_cran()
  data(interflex, envir = environment())
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "lasso", data = app_hma2015,
              Y = "totangry", D = "threat", X = "pidentity",
              Z = c("issuestr2", "knowledge", "educ", "male", "age10"),
              basis.type = "none", nboots = 20, theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
})

test_that("lasso (PO) estimator works with continuous treatment", {
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
    interflex(estimator = "lasso", data = s2,
              Y = "Y", D = "D", X = "X", Z = "Z1",
              nboots = 20, theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
  expect_true(!is.null(out$est.lasso))
})
