# Regression tests for the kernel estimator's adaptive-bandwidth normalizer
# fix (run KBW-20260917). These are builder's own tests, written from
# spec.md: they check that the implementation matches the spec, not the
# tester's independent scenarios.
#
# Conventions: skip_on_cran() first in every block; parallel = FALSE and
# figure = FALSE (no test here is about the figure); a seed is set before
# every fit; internal helpers are reached with interflex:::.

## Collects every warning whose message matches `pattern` raised while
## evaluating `expr`, without letting those warnings propagate (other
## warnings still propagate normally). Used only by the optional Malesky
## block (test 9) to check that no "dropped ..." warning is raised.
.collect_matching_warnings <- function(expr, pattern) {
  hits <- character(0)
  result <- withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl(pattern, conditionMessage(w))) {
        hits <<- c(hits, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    }
  )
  list(value = result, hits = hits)
}

test_that("kernel_nearest_index matches the which.min rule on random points, grid points, midpoints, and NA", {
  skip_on_cran()
  set.seed(20260917)
  x <- rlnorm(300)
  dens <- density(x)
  grid <- dens$x

  v <- runif(2000, min(grid) - 1, max(grid) + 1)
  got <- interflex:::.kernel_nearest_index(grid, v)
  want <- vapply(v, function(u) which.min(abs(grid - u)), 1L)
  expect_identical(got, want)

  got.grid <- interflex:::.kernel_nearest_index(grid, grid)
  want.grid <- vapply(grid, function(u) which.min(abs(grid - u)), 1L)
  expect_identical(got.grid, want.grid)

  mid <- (grid[-length(grid)] + grid[-1]) / 2
  got.mid <- interflex:::.kernel_nearest_index(grid, mid)
  want.mid <- vapply(mid, function(u) which.min(abs(grid - u)), 1L)
  expect_identical(got.mid, want.mid)

  expect_true(is.na(interflex:::.kernel_nearest_index(grid, NA_real_)))
})

test_that("kernel_prepare_density computes the weighted-geometric-mean g and the floor correctly", {
  skip_on_cran()
  set.seed(20260918)
  x <- rlnorm(60)

  # unit weights: g equals the unweighted geometric mean of f(X_i)
  w1 <- rep(1, 60)
  d1 <- suppressWarnings(density(x, weights = w1))
  prep1 <- interflex:::.kernel_prepare_density(d1, x, w1)
  idx1 <- interflex:::.kernel_nearest_index(prep1$x, x)
  fx1 <- prep1$y[idx1]
  manual.g1 <- exp(mean(log(fx1[fx1 > 0])))
  expect_equal(prep1$adapt.g, manual.g1, tolerance = 1e-12)

  # unequal positive weights: g equals the weighted formula
  set.seed(20260919)
  w2 <- runif(60, 0.5, 5)
  d2 <- suppressWarnings(density(x, weights = w2))
  prep2 <- interflex:::.kernel_prepare_density(d2, x, w2)
  idx2 <- interflex:::.kernel_nearest_index(prep2$x, x)
  fx2 <- prep2$y[idx2]
  ok2 <- fx2 > 0
  manual.g2 <- exp(sum(w2[ok2] * log(fx2[ok2])) / sum(w2[ok2]))
  expect_equal(prep2$adapt.g, manual.g2, tolerance = 1e-12)

  # scaling all weights by 7 scales adapt.g and dens$y by the same factor,
  # so the local bandwidth (which only depends on the ratio g / f) is unchanged
  d2b <- suppressWarnings(density(x, weights = w2 * 7))
  prep2b <- interflex:::.kernel_prepare_density(d2b, x, w2 * 7)
  for (x0 in x[1:5]) {
    h.a <- interflex:::.kernel_local_bw(x0, 1, prep2)
    h.b <- interflex:::.kernel_local_bw(x0, 1, prep2b)
    expect_equal(h.a, h.b, tolerance = 1e-12)
  }

  # an observation with weight 0 does not enter g
  x3 <- c(x, 500)
  w3 <- c(w1, 0)
  prep3 <- interflex:::.kernel_prepare_density(d1, x3, w3)
  expect_equal(prep3$adapt.g, prep1$adapt.g, tolerance = 1e-12)

  # adapt.floor equals the smallest positive dens$y
  expect_equal(prep1$adapt.floor, min(d1$y[d1$y > 0]))
})

test_that("kernel_local_bw matches the closed form, honours nearest-index ties, and rejects bw<=0/NA g", {
  skip_on_cran()
  set.seed(20260920)
  x <- rlnorm(80)
  w <- rep(1, 80)
  dens <- suppressWarnings(density(x, weights = w))
  prep <- interflex:::.kernel_prepare_density(dens, x, w)
  x0 <- x[10]
  f0 <- prep$y[interflex:::.kernel_nearest_index(prep$x, x0)]
  expect_equal(interflex:::.kernel_local_bw(x0, 1.7, prep), 1.7 * sqrt(prep$adapt.g / f0), tolerance = 1e-12)

  # hand-built density-like list with an exact zero, per spec.md section 10.3
  dl <- list(x = 0:9, y = c(5, 4, 0, 0, 0, 0, 0, 3, 2, 1))
  xs <- c(0, 1, 7, 8, 9)
  prep2 <- interflex:::.kernel_prepare_density(dl, xs, rep(1, 5))
  manual.g <- exp(mean(log(c(5, 4, 3, 2, 1))))
  expect_equal(prep2$adapt.g, manual.g, tolerance = 1e-12)
  expect_equal(prep2$adapt.floor, 1)

  bw0 <- 3
  # x0 = 4: nearest grid y is 0 -> floor used
  expect_equal(interflex:::.kernel_local_bw(4, bw0, prep2), bw0 * sqrt(manual.g / 1), tolerance = 1e-12)
  # x0 = 2.5: tie between two zero grid points -> floor used either way
  expect_equal(interflex:::.kernel_local_bw(2.5, bw0, prep2), bw0 * sqrt(manual.g / 1), tolerance = 1e-12)
  # x0 = 1.5: tie between index 2 (y=4) and index 3 (y=0) -> lower index wins, y=4
  expect_equal(interflex:::.kernel_local_bw(1.5, bw0, prep2), bw0 * sqrt(manual.g / 4), tolerance = 1e-12)

  na.dens <- list(x = 0:1, y = c(0, 0), adapt.g = NA_real_, adapt.floor = NA_real_)
  expect_true(is.na(interflex:::.kernel_local_bw(0, 1, na.dens)))
  expect_true(is.na(interflex:::.kernel_local_bw(4, 0, prep2)))
  expect_true(is.na(interflex:::.kernel_local_bw(4, -1, prep2)))
})

test_that("interflex(kernel) succeeds on a site-level moderator with heavy local ties (base code crashed here)", {
  skip_on_cran()
  set.seed(7)
  xs <- c(round(rlnorm(45, log(0.5), 0.4), 2), 4, 6.5, 6.5, 8.6, 8.6)
  m <- sample(5:15, 50, TRUE)
  X <- rep(xs, m)
  D <- rbinom(length(X), 1, 0.5)
  Y <- 1 + 0.2 * X + D * (0.5 - 0.1 * X) + rnorm(length(X))
  dat <- data.frame(Y = Y, D = D, X = X)

  for (bw in c(0.5, 2, 8)) {
    set.seed(7)
    out <- interflex(
      estimator = "kernel", Y = "Y", D = "D", X = "X", data = dat, bw = bw,
      CI = FALSE, figure = FALSE, parallel = FALSE, verbose = FALSE
    )
    expect_true(is.finite(out$bw))
    est <- out$est.kernel[[1]]
    expect_true(all(is.finite(est)))
    expect_gte(nrow(est), 45)
  }
})

test_that("degenerate far-apart clusters: unusable points are dropped with a clear warning, or the call stops", {
  skip_on_cran()
  set.seed(11)
  A <- data.frame(X = runif(60, 0, 1), D = 0)
  B <- data.frame(X = runif(60, 100, 101), D = 1)
  C <- data.frame(X = runif(120, 50, 51), D = rbinom(120, 1, 0.5))
  ABC <- rbind(A, B, C)
  ABC$Y <- 1 + 0.1 * ABC$X + ABC$D + rnorm(nrow(ABC))

  X.eval <- c(0.2, 0.5, 0.8, seq(50.1, 50.9, length.out = 5))

  set.seed(12)
  expect_warning(
    out <- interflex(
      estimator = "kernel", Y = "Y", D = "D", X = "X", data = ABC, bw = 0.05,
      X.eval = X.eval, diff.values = c(50.3, 50.7),
      CI = FALSE, figure = FALSE, parallel = FALSE, verbose = FALSE
    ),
    "dropped 3 of 8"
  )
  est <- out$est.kernel[[1]]
  expect_equal(nrow(est), 5)
  expect_true(all(is.finite(est)))
  expect_equal(sort(est[, "X"]), sort(seq(50.1, 50.9, length.out = 5)))
  expect_true(all(is.finite(out$diff.estimate[[1]][, "diff.estimate"])))

  set.seed(13)
  expect_warning(
    out2 <- interflex(
      estimator = "kernel", Y = "Y", D = "D", X = "X", data = ABC, bw = 0.05,
      X.eval = X.eval, diff.values = c(0.5, 50.5),
      CI = FALSE, figure = FALSE, parallel = FALSE, verbose = FALSE
    ),
    "not usable"
  )
  expect_true(all(is.na(out2$diff.estimate[[1]][, "diff.estimate"])))

  AB <- rbind(A, B)
  AB$Y <- 1 + 0.1 * AB$X + AB$D + rnorm(nrow(AB))
  set.seed(14)
  expect_error(
    interflex(
      estimator = "kernel", Y = "Y", D = "D", X = "X", data = AB, bw = 0.05,
      X.eval = c(0.2, 0.5, 0.8, 100.2, 100.5),
      CI = FALSE, figure = FALSE, parallel = FALSE, verbose = FALSE
    ),
    "Inappropriate bandwidth"
  )
})

test_that("bw is validated up front with a clear error", {
  skip_on_cran()
  set.seed(15)
  dat <- data.frame(X = runif(40), D = rbinom(40, 1, 0.5))
  dat$Y <- 1 + dat$X + dat$D + rnorm(40)

  bad.bw.values <- list(-1, 0, NA_real_, Inf, c(1, 2))
  for (bad.bw in bad.bw.values) {
    expect_error(
      interflex(
        estimator = "kernel", Y = "Y", D = "D", X = "X", data = dat, bw = bad.bw,
        CI = FALSE, figure = FALSE, parallel = FALSE, verbose = FALSE
      ),
      "bw"
    )
  }
})

test_that("cross-validation with no usable candidate bandwidth errors and returns nothing", {
  skip_on_cran()
  set.seed(16)
  A <- data.frame(X = runif(60, 0, 1), D = 0)
  B <- data.frame(X = runif(60, 100, 101), D = 1)
  C <- data.frame(X = runif(120, 50, 51), D = rbinom(120, 1, 0.5))
  ABC <- rbind(A, B, C)
  ABC$Y <- 1 + 0.1 * ABC$X + ABC$D + rnorm(nrow(ABC))

  set.seed(17)
  expect_error(
    interflex(
      estimator = "kernel", Y = "Y", D = "D", X = "X", data = ABC,
      grid = c(1e-8, 2e-8),
      CI = FALSE, figure = FALSE, parallel = FALSE, verbose = FALSE
    ),
    "Inappropriate bandwidth|Bandwidth selection failed"
  )
})

test_that("interflex() never leaves the global uniform_ci_warned option set", {
  skip_on_cran()
  set.seed(7)
  xs <- c(round(rlnorm(45, log(0.5), 0.4), 2), 4, 6.5, 6.5, 8.6, 8.6)
  m <- sample(5:15, 50, TRUE)
  X <- rep(xs, m)
  D <- rbinom(length(X), 1, 0.5)
  Y <- 1 + 0.2 * X + D * (0.5 - 0.1 * X) + rnorm(length(X))
  dat <- data.frame(Y = Y, D = D, X = X)

  set.seed(18)
  invisible(interflex(
    estimator = "kernel", Y = "Y", D = "D", X = "X", data = dat, bw = 2,
    vartype = "bootstrap", nboots = 10, CI = TRUE,
    figure = FALSE, parallel = FALSE, verbose = FALSE
  ))
  expect_null(getOption("interflex.uniform_ci_warned"))
})

test_that("Malesky raw-X kernel fits succeed at several bandwidths and under default CV (optional, real data)", {
  skip_on_cran()
  p <- Sys.getenv("INTERFLEX_MALESKY_CSV")
  skip_if(!nzchar(p) || !file.exists(p))

  d <- read.csv(p)
  expect_true(all(c("pci_id", "t2", "internet_users100", "d_question_count") %in% colnames(d)))

  for (bw in c(0.5, 2, 8, 50)) {
    set.seed(1234)
    res <- .collect_matching_warnings(
      interflex(
        estimator = "kernel", Y = "d_question_count", D = "t2", X = "internet_users100",
        data = d, bw = bw, CI = FALSE, figure = FALSE, parallel = FALSE, verbose = FALSE
      ),
      "dropped"
    )
    expect_length(res$hits, 0)
    out <- res$value
    expect_true(is.finite(out$bw))
    est <- out$est.kernel[[1]]
    expect_true(all(is.finite(est)))
    expect_equal(nrow(est), 50)
  }

  set.seed(1234)
  res.cv <- .collect_matching_warnings(
    interflex(
      estimator = "kernel", Y = "d_question_count", D = "t2", X = "internet_users100",
      data = d, CI = FALSE, figure = FALSE, parallel = FALSE, verbose = FALSE
    ),
    "dropped"
  )
  expect_length(res.cv$hits, 0)
  out.cv <- res.cv$value
  expect_true(is.finite(out.cv$bw))
  est.cv <- out.cv$est.kernel[[1]]
  expect_true(all(is.finite(est.cv)))
  expect_equal(nrow(est.cv), 50)
})
