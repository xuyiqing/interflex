# test-gate.R — Unified GATE tests across estimators
# Generated from test-spec.md (GATE-001)

# Shared DGP
make_gate_data <- function() {
  set.seed(1234)
  n <- 500
  X <- sample(1:4, n, replace = TRUE, prob = c(0.25, 0.25, 0.25, 0.25))
  Z1 <- rnorm(n)
  Z2 <- runif(n, -1, 1)
  D_binary <- rbinom(n, 1, plogis(0.3 * Z1 - 0.2 * Z2))
  D_continuous <- 2 + 0.5 * Z1 + rnorm(n)

  true_gate <- c(`1` = 2, `2` = 4, `3` = 6, `4` = 3)
  tau_i <- true_gate[as.character(X)]

  e <- rnorm(n, sd = 1.5)
  Y_binary <- 3 + tau_i * D_binary + 0.8 * Z1 - 0.5 * Z2 + e
  Y_continuous <- 3 + tau_i * D_continuous + 0.8 * Z1 - 0.5 * Z2 + e

  data.frame(Y_bin = Y_binary, Y_cont = Y_continuous,
             D_bin = D_binary, D_cont = D_continuous,
             X = X, Z1 = Z1, Z2 = Z2)
}

# --- T1. Linear + gate=TRUE + binary D ---
test_that("T1: Linear + gate=TRUE + binary D produces correct GATE", {
  skip_on_cran()
  dat <- make_gate_data()
  true_gate <- c(`1` = 2, `2` = 4, `3` = 6, `4` = 3)

  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "linear", data = dat,
              Y = "Y_bin", D = "D_bin", X = "X", Z = c("Z1", "Z2"),
              gate = TRUE, figure = FALSE, nboots = 50, CI = TRUE)
  ))

  expect_s3_class(out, "interflex")
  expect_true(!is.null(out$g.est))
  expect_true(is.list(out$g.est))
  expect_true(length(out$g.est) >= 1)

  # Check structure of each element
  for (i in seq_along(out$g.est)) {
    df <- out$g.est[[i]]
    expect_true(is.data.frame(df))
    expect_true("X" %in% names(df))
    expect_true("ME" %in% names(df))
    expect_equal(nrow(df), 4)
  }

  # Check GATE estimates within tolerance (abs < 3)
  df <- out$g.est[[1]]
  for (j in seq_len(nrow(df))) {
    x_val <- as.character(df$X[j])
    expect_true(abs(df$ME[j] - true_gate[x_val]) < 3,
                info = paste("GATE for X =", x_val, ": ME =", df$ME[j],
                             ", true =", true_gate[x_val]))
  }

  # Smooth CME curve still produced
  expect_true(!is.null(out$est.lin))
})

# --- T2. Linear + gate=TRUE + continuous D ---
test_that("T2: Linear + gate=TRUE + continuous D", {
  skip_on_cran()
  dat <- make_gate_data()

  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "linear", data = dat,
              Y = "Y_cont", D = "D_cont", X = "X", Z = c("Z1", "Z2"),
              gate = TRUE, figure = FALSE, nboots = 50, CI = TRUE)
  ))

  expect_true(!is.null(out$g.est))
  expect_true(is.list(out$g.est))
  for (i in seq_along(out$g.est)) {
    df <- out$g.est[[i]]
    expect_true("X" %in% names(df))
    expect_true("ME" %in% names(df))
  }
  expect_true(!is.null(out$est.lin))
})

# --- T3. GRF + gate=TRUE + binary D ---
test_that("T3: GRF + gate=TRUE + binary D", {
  skip_on_cran()
  dat <- make_gate_data()

  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "grf", data = dat,
              Y = "Y_bin", D = "D_bin", X = "X", Z = c("Z1", "Z2"),
              gate = TRUE, figure = FALSE)
  ))

  expect_true(!is.null(out$g.est))
  expect_true(is.list(out$g.est))
  for (i in seq_along(out$g.est)) {
    df <- out$g.est[[i]]
    expect_true(is.data.frame(df))
    expect_true("X" %in% names(df))
    expect_true("ME" %in% names(df))
    expect_equal(nrow(df), 4)
  }
  expect_true(!is.null(out$est.grf))
})

# --- T4. DML + gate=TRUE — regression test ---
test_that("T4: DML + gate=TRUE regression test", {
  skip_on_cran()
  dat <- make_gate_data()

  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "dml", data = dat,
              Y = "Y_bin", D = "D_bin", X = "X", Z = c("Z1", "Z2"),
              gate = TRUE, figure = FALSE,
              model.y = "linear", model.t = "linear")
  ))

  expect_true(!is.null(out$g.est))
  expect_true(is.list(out$g.est))
  expect_true(!is.null(out$g.est.dml))
  expect_true(is.list(out$g.est.dml))

  # g.est and g.est.dml should contain the same data
  for (i in seq_along(out$g.est)) {
    expect_equal(out$g.est[[i]], out$g.est.dml[[i]])
  }

  for (i in seq_along(out$g.est)) {
    df <- out$g.est[[i]]
    expect_true("X" %in% names(df))
    expect_true("ME" %in% names(df))
  }
})

# --- T5. Lasso + gate=TRUE ---
test_that("T5: Lasso + gate=TRUE explicit", {
  skip_on_cran()
  dat <- make_gate_data()

  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "lasso", data = dat,
              Y = "Y_bin", D = "D_bin", X = "X", Z = c("Z1", "Z2"),
              gate = TRUE, figure = FALSE, CI = TRUE, nboots = 50)
  ))

  expect_true(!is.null(out$g.est))
  expect_true(!is.null(out$est.lasso))
  expect_s3_class(out, "interflex")
})

# --- T6. Plot by.group=TRUE — linear ---
test_that("T6: Plot by.group=TRUE works for linear", {
  skip_on_cran()
  dat <- make_gate_data()

  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "linear", data = dat,
              Y = "Y_bin", D = "D_bin", X = "X", Z = c("Z1", "Z2"),
              gate = TRUE, figure = FALSE, nboots = 50, CI = TRUE)
  ))

  p <- suppressWarnings(suppressMessages(
    plot(out, by.group = TRUE)
  ))
  expect_true(!is.null(p))
})

# --- T7. Plot by.group=TRUE — grf ---
test_that("T7: Plot by.group=TRUE works for grf", {
  skip_on_cran()
  dat <- make_gate_data()

  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "grf", data = dat,
              Y = "Y_bin", D = "D_bin", X = "X", Z = c("Z1", "Z2"),
              gate = TRUE, figure = FALSE)
  ))

  p <- suppressWarnings(suppressMessages(
    plot(out, by.group = TRUE)
  ))
  expect_true(!is.null(p))
})

# --- T8. Plot by.group=TRUE — dml ---
test_that("T8: Plot by.group=TRUE works for dml", {
  skip_on_cran()
  dat <- make_gate_data()

  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "dml", data = dat,
              Y = "Y_bin", D = "D_bin", X = "X", Z = c("Z1", "Z2"),
              gate = TRUE, figure = FALSE,
              model.y = "linear", model.t = "linear")
  ))

  p <- suppressWarnings(suppressMessages(
    plot(out, by.group = TRUE)
  ))
  expect_true(!is.null(p))
})

# --- T9. Plot by.group=TRUE — lasso ---
test_that("T9: Plot by.group=TRUE works for lasso", {
  skip_on_cran()
  dat <- make_gate_data()

  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "lasso", data = dat,
              Y = "Y_bin", D = "D_bin", X = "X", Z = c("Z1", "Z2"),
              gate = TRUE, figure = FALSE, CI = TRUE, nboots = 50)
  ))

  p <- suppressWarnings(suppressMessages(
    plot(out, by.group = TRUE)
  ))
  expect_true(!is.null(p))
})

# --- E1. gate=TRUE + kernel -> error ---
test_that("E1: gate=TRUE + kernel errors", {
  skip_on_cran()
  dat <- make_gate_data()

  expect_error(
    interflex(estimator = "kernel", data = dat,
              Y = "Y_bin", D = "D_bin", X = "X", gate = TRUE),
    "gate = TRUE is not supported"
  )
})

# --- E2. gate=TRUE + binning -> error ---
test_that("E2: gate=TRUE + binning errors", {
  skip_on_cran()
  dat <- make_gate_data()

  expect_error(
    interflex(estimator = "binning", data = dat,
              Y = "Y_bin", D = "D_bin", X = "X", gate = TRUE),
    "gate = TRUE is not supported"
  )
})

# --- E3. gate=TRUE + continuous X -> error ---
test_that("E3: gate=TRUE + continuous X errors", {
  skip_on_cran()
  dat <- make_gate_data()
  dat_cont_x <- dat
  dat_cont_x$X <- rnorm(nrow(dat))

  expect_error(
    interflex(estimator = "linear", data = dat_cont_x,
              Y = "Y_bin", D = "D_bin", X = "X", gate = TRUE),
    "discrete moderator|continuous"
  )
})

# --- E4. gate=FALSE (default) — no g.est ---
test_that("E4: gate=FALSE produces no g.est", {
  skip_on_cran()
  dat <- make_gate_data()

  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "linear", data = dat,
              Y = "Y_bin", D = "D_bin", X = "X", figure = FALSE)
  ))

  expect_null(out$g.est)
})

# --- Property-based invariants ---
test_that("P1-P4: Property-based invariants for GATE output", {
  skip_on_cran()
  dat <- make_gate_data()

  out <- suppressWarnings(suppressMessages(
    interflex(estimator = "linear", data = dat,
              Y = "Y_bin", D = "D_bin", X = "X", Z = c("Z1", "Z2"),
              gate = TRUE, figure = FALSE, nboots = 50, CI = TRUE)
  ))

  for (i in seq_along(out$g.est)) {
    df <- out$g.est[[i]]

    # P1: Structure consistency
    expect_true(is.data.frame(df))
    expect_true(all(c("X", "ME") %in% names(df)))

    # P2: Level coverage
    expect_true(all(sort(unique(dat$X)) %in% df$X))

    # P3: Sign consistency (all true effects positive)
    expect_true(all(df$ME > 0), info = paste("Negative ME found:", paste(df$ME, collapse = ", ")))

    # P4: CI monotonicity
    if ("lower CI(95%)" %in% names(df) && "upper CI(95%)" %in% names(df)) {
      expect_true(all(df[["lower CI(95%)"]] <= df$ME))
      expect_true(all(df$ME <= df[["upper CI(95%)"]]))
    }
  }
})
