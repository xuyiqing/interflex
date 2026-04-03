test_that("plot method works with Xdistr = density", {
  data(interflex)
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "binning", data = app_hma2015,
              Y = "totangry", D = "threat", X = "pidentity",
              Z = c("issuestr2", "knowledge", "educ", "male", "age10"),
              theme.bw = TRUE)
  ))
  p <- suppressWarnings(suppressMessages(
    plot(out, Xdistr = "density", theme.bw = TRUE)
  ))
  expect_true(!is.null(p))
})

test_that("plot method works with Xdistr = none", {
  data(interflex)
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "binning", data = app_hma2015,
              Y = "totangry", D = "threat", X = "pidentity",
              Z = c("issuestr2", "knowledge", "educ", "male", "age10"),
              theme.bw = TRUE)
  ))
  p <- suppressWarnings(suppressMessages(
    plot(out, Xdistr = "none", theme.bw = TRUE)
  ))
  expect_true(!is.null(p))
})

test_that("plot method works with bin.labs = FALSE", {
  data(interflex)
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "binning", data = app_hma2015,
              Y = "totangry", D = "threat", X = "pidentity",
              Z = c("issuestr2", "knowledge", "educ", "male", "age10"),
              theme.bw = TRUE)
  ))
  p <- suppressWarnings(suppressMessages(
    plot(out, bin.labs = FALSE, theme.bw = TRUE)
  ))
  expect_true(!is.null(p))
})

test_that("plot method works with show.grid = FALSE", {
  data(interflex)
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "linear", data = app_hma2015,
              Y = "totangry", D = "threat", X = "pidentity",
              Z = c("issuestr2", "knowledge", "educ", "male", "age10"),
              vartype = "delta", theme.bw = TRUE)
  ))
  p <- suppressWarnings(suppressMessages(
    plot(out, show.grid = FALSE, theme.bw = TRUE)
  ))
  expect_true(!is.null(p))
})

test_that("plot method works with order and subtitles for multi-arm", {
  data(interflex)
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "binning", data = s5,
              Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"),
              theme.bw = TRUE)
  ))
  p <- suppressWarnings(suppressMessages(
    plot(out, order = c("C", "B"),
         subtitles = c("Group C", "Group B"),
         theme.bw = TRUE)
  ))
  expect_true(!is.null(p))
})

test_that("plot method works with pool = TRUE for multi-arm", {
  data(interflex)
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "binning", data = s5,
              Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"),
              theme.bw = TRUE)
  ))
  p <- suppressWarnings(suppressMessages(
    plot(out, order = c("C", "B"),
         subtitles = c("Control Group", "Group C", "Group B"),
         pool = TRUE, theme.bw = TRUE)
  ))
  expect_true(!is.null(p))
})

test_that("plot method works with CI = FALSE", {
  data(interflex)
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "linear", data = s5,
              Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"),
              vartype = "delta", theme.bw = TRUE)
  ))
  p <- suppressWarnings(suppressMessages(
    plot(out, CI = FALSE, theme.bw = TRUE)
  ))
  expect_true(!is.null(p))
})

test_that("plot method works with ncols option", {
  data(interflex)
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "raw", data = s5,
              Y = "Y", D = "D", X = "X",
              ncols = 3, theme.bw = TRUE)
  ))
  expect_s3_class(out, "interflex")
})

test_that("plot method works with show.all = TRUE", {
  data(interflex)
  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "linear", data = s5,
              Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"),
              vartype = "delta", theme.bw = TRUE)
  ))
  p <- suppressWarnings(suppressMessages(
    plot(out, show.all = TRUE, theme.bw = TRUE)
  ))
  expect_true(is.list(p))
  expect_true(length(p) >= 1)
})
