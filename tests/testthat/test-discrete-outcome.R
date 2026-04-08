test_that("binning estimator works with logit method", {
  skip_on_cran()
  set.seed(110)
  n <- 500
  x <- runif(n, -3, 3)
  d <- sample(c(0, 1), n, replace = TRUE)
  Z1 <- runif(n, 0, 1)
  link <- -1 + d + x + d * x
  prob <- exp(link) / (1 + exp(link))
  y <- as.numeric(runif(n) < prob)
  df <- data.frame(Y = y, D = d, X = x, Z1 = Z1)

  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "binning", method = "logit", data = df,
              Y = "Y", D = "D", X = "X", Z = "Z1",
              na.rm = TRUE, theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
  expect_true(!is.null(out$figure))
  expect_true(!is.null(out$est.bin))
})

test_that("linear estimator works with logit method", {
  skip_on_cran()
  set.seed(110)
  n <- 500
  x <- runif(n, -3, 3)
  d <- sample(c(0, 1), n, replace = TRUE)
  Z1 <- runif(n, 0, 1)
  link <- -1 + d + x + d * x
  prob <- exp(link) / (1 + exp(link))
  y <- as.numeric(runif(n) < prob)
  df <- data.frame(Y = y, D = d, X = x, Z1 = Z1)

  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "linear", method = "logit", data = df,
              Y = "Y", D = "D", X = "X", Z = "Z1",
              na.rm = TRUE, theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
  expect_true(!is.null(out$est.lin))
})

test_that("linear estimator works with probit method", {
  skip_on_cran()
  set.seed(110)
  n <- 500
  x <- runif(n, -3, 3)
  d <- sample(c(0, 1), n, replace = TRUE)
  Z1 <- runif(n, 0, 1)
  link <- -1 + d + x + d * x
  prob <- exp(link) / (1 + exp(link))
  y <- as.numeric(runif(n) < prob)
  df <- data.frame(Y = y, D = d, X = x, Z1 = Z1)

  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "linear", method = "probit", data = df,
              Y = "Y", D = "D", X = "X", Z = "Z1",
              na.rm = TRUE, theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
  expect_true(!is.null(out$est.lin))
})

test_that("kernel estimator works with logit method", {
  skip_on_cran()
  set.seed(110)
  n <- 500
  x <- runif(n, -3, 3)
  d <- sample(c(0, 1), n, replace = TRUE)
  Z1 <- runif(n, 0, 1)
  link <- -1 + d + x + d * x
  prob <- exp(link) / (1 + exp(link))
  y <- as.numeric(runif(n) < prob)
  df <- data.frame(Y = y, D = d, X = x, Z1 = Z1)

  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "kernel", method = "logit", data = df,
              Y = "Y", D = "D", X = "X", Z = "Z1",
              na.rm = TRUE, bw = 1.0, theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
  expect_true(!is.null(out$est.kernel))
})

test_that("linear estimator works with poisson method (count outcome)", {
  skip_on_cran()
  set.seed(42)
  n <- 500
  x <- runif(n, -2, 2)
  d <- sample(c(0, 1), n, replace = TRUE)
  Z1 <- runif(n, 0, 1)
  lambda <- exp(1 + 0.5 * d + 0.3 * x + 0.2 * d * x + 0.1 * Z1)
  y <- rpois(n, lambda)
  df <- data.frame(Y = y, D = d, X = x, Z1 = Z1)

  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "linear", method = "poisson", data = df,
              Y = "Y", D = "D", X = "X", Z = "Z1",
              na.rm = TRUE, theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
  expect_true(!is.null(out$est.lin))
})

test_that("linear estimator works with nbinom method (count outcome)", {
  skip_on_cran()
  set.seed(42)
  n <- 500
  x <- runif(n, -2, 2)
  d <- sample(c(0, 1), n, replace = TRUE)
  Z1 <- runif(n, 0, 1)
  lambda <- exp(1 + 0.5 * d + 0.3 * x + 0.2 * d * x + 0.1 * Z1)
  y <- rpois(n, lambda)
  df <- data.frame(Y = y, D = d, X = x, Z1 = Z1)

  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "linear", method = "nbinom", data = df,
              Y = "Y", D = "D", X = "X", Z = "Z1",
              na.rm = TRUE, theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
  expect_true(!is.null(out$est.lin))
})
