test_that("kernel estimator works with binary treatment", {
  skip_on_cran()
  data(interflex, envir = environment())
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "kernel", data = app_hma2015,
              Y = "totangry", D = "threat", X = "pidentity",
              Z = c("issuestr2", "knowledge", "educ", "male", "age10"),
              bw = 0.916, theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
  expect_true(!is.null(out$figure))
  expect_true(!is.null(out$est.kernel))
})

test_that("kernel estimator works with bootstrap vartype", {
  skip_on_cran()
  data(interflex, envir = environment())
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "kernel", data = app_hma2015,
              Y = "totangry", D = "threat", X = "pidentity",
              Z = c("issuestr2", "knowledge", "educ", "male", "age10"),
              vartype = "bootstrap", nboots = 20, bw = 0.916,
              theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
  expect_true(!is.null(out$est.kernel))
})

test_that("kernel estimator works with continuous treatment", {
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
    interflex(estimator = "kernel", data = s2,
              Y = "Y", D = "D", X = "X", Z = "Z1",
              bw = 1.0, theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
  expect_true(!is.null(out$est.kernel))
})

test_that("kernel estimator works with full.moderate", {
  skip_on_cran()
  set.seed(1234)
  n <- 200
  d1 <- sample(c(0, 1), n, replace = TRUE)
  x <- rnorm(n, 3, 1)
  z <- rnorm(n, 3, 1)
  e <- rnorm(n, 0, 1)
  y1 <- 5 - 4 * x - 9 * d1 + 3 * x * d1 + 1 * z + 2 * e
  s1 <- data.frame(Y = y1, D = d1, X = x, Z1 = z)

  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "kernel", data = s1,
              Y = "Y", D = "D", X = "X", Z = "Z1",
              bw = 1.0, full.moderate = TRUE, theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
})
