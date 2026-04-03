test_that("raw estimator works with binary treatment", {
  data(interflex, envir = environment())
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "raw", data = app_hma2015,
              Y = "totangry", D = "threat", X = "pidentity",
              Z = c("issuestr2", "knowledge", "educ", "male", "age10"),
              theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
  expect_true(!is.null(out$figure))
})

test_that("raw estimator works with continuous treatment (simulated)", {
  set.seed(1234)
  n <- 200
  d2 <- rnorm(n, 3, 1)
  x <- rnorm(n, 3, 1)
  z <- rnorm(n, 3, 1)
  e <- rnorm(n, 0, 1)
  y2 <- 5 - 4 * x - 9 * d2 + 3 * x * d2 + 1 * z + 2 * e
  s2 <- data.frame(Y = y2, D = d2, X = x, Z1 = z)

  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "raw", data = s2,
              Y = "Y", D = "D", X = "X", Z = "Z1",
              theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
  expect_true(!is.null(out$figure))
})

test_that("raw estimator works with built-in s1 data", {
  data(interflex, envir = environment())
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "raw", data = s1,
              Y = "Y", D = "D", X = "X", Z = "Z1",
              theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
})
