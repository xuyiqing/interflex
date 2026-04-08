test_that("binning estimator works with binary treatment", {
  skip_on_cran()
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
  ## p.wald is currently stored as a formatted character string for display;
  ## assert presence and parseability rather than numeric type.
  expect_false(is.null(out$tests$p.wald))
  expect_false(is.na(suppressWarnings(as.numeric(out$tests$p.wald))))
})

test_that("binning estimator works with full.moderate", {
  skip_on_cran()
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
  skip_on_cran()
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
  skip_on_cran()
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
    interflex(estimator = "binning", data = s2,
              Y = "Y", D = "D", X = "X", Z = "Z1",
              theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
  expect_true(!is.null(out$est.bin))
})
