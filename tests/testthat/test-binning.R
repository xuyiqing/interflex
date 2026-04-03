test_that("binning estimator works with binary treatment", {
  data(interflex, envir = environment())
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "binning", data = app_hma2015,
              Y = "totangry", D = "threat", X = "pidentity",
              Z = c("issuestr2", "knowledge", "educ", "male", "age10"),
              theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
  expect_true(!is.null(out$figure))
  expect_true(!is.null(out$est.bin))
  expect_true(!is.null(out$tests))
  expect_true(is.numeric(out$tests$p.wald))
})

test_that("binning estimator works with full.moderate", {
  data(interflex, envir = environment())
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "binning", data = app_hma2015,
              Y = "totangry", D = "threat", X = "pidentity",
              Z = c("issuestr2", "knowledge", "educ", "male", "age10"),
              theme.bw = TRUE, full.moderate = TRUE)
  ))
  expect_s3_class(out, "interflex")
  expect_true(!is.null(out$est.bin))
})

test_that("binning estimator works with nbins option", {
  data(interflex, envir = environment())
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "binning", data = app_hma2015,
              Y = "totangry", D = "threat", X = "pidentity",
              Z = c("issuestr2", "knowledge", "educ", "male", "age10"),
              nbins = 5, theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
})

test_that("binning estimator works with Xunif", {
  data(interflex, envir = environment())
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "binning", data = app_hma2015,
              Y = "totangry", D = "threat", X = "pidentity",
              Z = c("issuestr2", "knowledge", "educ", "male", "age10"),
              Xunif = TRUE, theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
})

test_that("binning estimator works with continuous treatment", {
  set.seed(1234)
  n <- 200
  d2 <- rnorm(n, 3, 1)
  x <- rnorm(n, 3, 1)
  z <- rnorm(n, 3, 1)
  e <- rnorm(n, 0, 1)
  y2 <- 5 - 4 * x - 9 * d2 + 3 * x * d2 + 1 * z + 2 * e
  s2 <- data.frame(Y = y2, D = d2, X = x, Z1 = z)

  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "binning", data = s2,
              Y = "Y", D = "D", X = "X", Z = "Z1",
              theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
  expect_true(!is.null(out$est.bin))
})
