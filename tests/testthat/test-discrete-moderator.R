test_that("DML with gate = TRUE works for binary treatment with discrete moderator", {
  skip_on_cran()
  set.seed(9999)
  n <- 500
  D <- sample(c(0, 1), n, replace = TRUE)
  X <- sample(1:5, size = n, replace = TRUE, prob = c(0.1, 0.3, 0.3, 0.2, 0.1))
  Z1 <- runif(n, min = -2, max = 2)
  Z2 <- rgamma(n, shape = 2, rate = 1)
  slopes_by_x <- c(`1` = 2, `2` = 4, `3` = 7, `4` = 3, `5` = 10)
  tau_i <- slopes_by_x[as.character(X)] + rnorm(n, mean = 0, sd = 0.5)
  e <- rnorm(n, mean = 0, sd = 2)
  Y <- 5 + tau_i * D + 1.5 * (Z1^2) - 2.0 * sqrt(Z2) + e
  s10 <- data.frame(Y = Y, D = D, X = X, Z1 = Z1, Z2 = Z2)

  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "dml", data = s10,
              Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"),
              gate = TRUE, theme.bw = TRUE,
              model.y = "linear", model.t = "linear")
  ))
  expect_s3_class(out, "interflex")
  expect_true(!is.null(out$figure))
  expect_true(!is.null(out$g.est.dml))
})

test_that("DML with gate = TRUE works for continuous treatment with discrete moderator", {
  skip_on_cran()
  set.seed(9999)
  n <- 500
  D <- rnorm(n, mean = 3, sd = 1)
  X <- sample(1:5, size = n, replace = TRUE, prob = c(0.1, 0.3, 0.3, 0.2, 0.1))
  Z1 <- runif(n, min = -2, max = 2)
  Z2 <- rgamma(n, shape = 2, rate = 1)
  slopes_by_x <- c(`1` = 2, `2` = 4, `3` = 7, `4` = 3, `5` = 10)
  tau_i <- slopes_by_x[as.character(X)] + rnorm(n, mean = 0, sd = 0.5)
  e <- rnorm(n, mean = 0, sd = 2)
  Y <- 5 + tau_i * D + 1.5 * (Z1^2) - 2.0 * sqrt(Z2) + e
  s10 <- data.frame(Y = Y, D = D, X = X, Z1 = Z1, Z2 = Z2)

  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "dml", data = s10,
              Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"),
              gate = TRUE, theme.bw = TRUE,
              model.y = "linear", model.t = "linear")
  ))
  expect_s3_class(out, "interflex")
  expect_true(!is.null(out$figure))
})

test_that("plot with by.group = TRUE works for discrete moderator", {
  skip_on_cran()
  set.seed(9999)
  n <- 500
  D <- sample(c(0, 1), n, replace = TRUE)
  X <- sample(1:5, size = n, replace = TRUE, prob = c(0.1, 0.3, 0.3, 0.2, 0.1))
  Z1 <- runif(n, min = -2, max = 2)
  Z2 <- rgamma(n, shape = 2, rate = 1)
  slopes_by_x <- c(`1` = 2, `2` = 4, `3` = 7, `4` = 3, `5` = 10)
  tau_i <- slopes_by_x[as.character(X)] + rnorm(n, mean = 0, sd = 0.5)
  e <- rnorm(n, mean = 0, sd = 2)
  Y <- 5 + tau_i * D + 1.5 * (Z1^2) - 2.0 * sqrt(Z2) + e
  s10 <- data.frame(Y = Y, D = D, X = X, Z1 = Z1, Z2 = Z2)

  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "dml", data = s10,
              Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"),
              gate = TRUE, theme.bw = TRUE,
              model.y = "linear", model.t = "linear")
  ))

  p <- suppressWarnings(suppressMessages(
    plot(out, by.group = TRUE, theme.bw = TRUE)
  ))
  expect_true(!is.null(p))
})
