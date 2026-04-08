test_that("fixed effects work with binning estimator", {
  skip_on_cran()
  set.seed(1234)
  n <- 500
  d4 <- sample(c(0, 1), n, replace = TRUE)
  x4 <- runif(n, -3, 3)
  z4 <- rnorm(n, 3, 1)
  alpha <- 20 * rep(rnorm(n / 10), each = 10)
  xi <- rep(rnorm(10), n / 10)
  y4 <- d4 * (x4^2 - 2.5) + (1 - d4) * (-1 * x4^2 + 2.5) + 1 * z4 +
    alpha + xi + 2 * rnorm(n, 0, 1)
  s4 <- data.frame(D = d4, X = x4, Y = y4, Z1 = z4,
                   unit = rep(1:(n / 10), each = 10),
                   year = rep(1:10, n / 10))
  s4$wgt <- 1

  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "binning", data = s4,
              Y = "Y", D = "D", X = "X", Z = "Z1",
              FE = c("unit", "year"), cl = "unit",
              weights = "wgt", theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
  expect_true(!is.null(out$figure))
})

test_that("fixed effects work with kernel estimator", {
  skip_on_cran()
  set.seed(1234)
  n <- 200
  d4 <- sample(c(0, 1), n, replace = TRUE)
  x4 <- runif(n, -3, 3)
  z4 <- rnorm(n, 3, 1)
  alpha <- 20 * rep(rnorm(n / 10), each = 10)
  xi <- rep(rnorm(10), n / 10)
  y4 <- d4 * (x4^2 - 2.5) + (1 - d4) * (-1 * x4^2 + 2.5) + 1 * z4 +
    alpha + xi + 2 * rnorm(n, 0, 1)
  s4 <- data.frame(D = d4, X = x4, Y = y4, Z1 = z4,
                   unit = rep(1:(n / 10), each = 10),
                   year = rep(1:10, n / 10))
  s4$wgt <- 1

  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "kernel", data = s4,
              Y = "Y", D = "D", X = "X", Z = "Z1",
              FE = c("unit", "year"), cl = "unit",
              bw = 0.62, weights = "wgt", theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
})

test_that("multiple treatment arms work with binning estimator", {
  skip_on_cran()
  data(interflex)
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "binning", data = s5,
              Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"),
              theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
  expect_true(!is.null(out$figure))
})

test_that("multiple treatment arms work with base group option", {
  skip_on_cran()
  data(interflex)
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "binning", data = s5,
              Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"),
              base = "B", theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
})

test_that("multiple treatment arms work with linear estimator", {
  skip_on_cran()
  data(interflex)
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "linear", data = s5,
              Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"),
              vartype = "delta", theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
  expect_true(!is.null(out$est.lin))
})

test_that("predict method works on binning output", {
  skip_on_cran()
  data(interflex)
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "binning", data = s5,
              Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"),
              theme.bw = TRUE)
  ))
  pred <- suppressWarnings(suppressMessages(
    predict(out, theme.bw = TRUE)
  ))
  expect_true(!is.null(pred))
})

test_that("predict method works on linear output", {
  skip_on_cran()
  data(interflex)
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "linear", data = s5,
              Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"),
              vartype = "delta", theme.bw = TRUE)
  ))
  pred <- suppressWarnings(suppressMessages(
    predict(out, theme.bw = TRUE)
  ))
  expect_true(!is.null(pred))
})

test_that("inter.test works on linear output", {
  skip_on_cran()
  data(interflex)
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "linear", data = s5,
              Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"),
              vartype = "delta", theme.bw = TRUE)
  ))
  result <- suppressWarnings(suppressMessages(
    inter.test(out, diff.values = c(-2, 0, 2))
  ))
  expect_true(!is.null(result))
})

test_that("inter.test works with percentile = TRUE", {
  skip_on_cran()
  data(interflex)
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "linear", data = s5,
              Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"),
              vartype = "delta", theme.bw = TRUE)
  ))
  result <- suppressWarnings(suppressMessages(
    inter.test(out, diff.values = c(0.25, 0.5, 0.75), percentile = TRUE)
  ))
  expect_true(!is.null(result))
})

test_that("diff.values option works in linear estimator", {
  skip_on_cran()
  data(interflex)
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "linear", data = s5,
              Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"),
              diff.values = c(-2, 0, 2),
              vartype = "delta", theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
  expect_true(!is.null(out$diff.estimate))
})
