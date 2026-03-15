# Internal utility functions for interflex
# Dot-prefixed to prevent auto-export via exportPattern("^[[:alpha:]]+")

.extract_treat_info <- function(treat.info) {
  treat.type <- treat.info[["treat.type"]]
  result <- list(treat.type = treat.type)

  if (treat.type == "discrete") {
    other.treat <- treat.info[["other.treat"]]
    other.treat.origin <- names(other.treat)
    names(other.treat.origin) <- other.treat
    all.treat <- treat.info[["all.treat"]]
    all.treat.origin <- names(all.treat)
    names(all.treat.origin) <- all.treat
    base <- treat.info[["base"]]

    result$other.treat <- other.treat
    result$other.treat.origin <- other.treat.origin
    result$all.treat <- all.treat
    result$all.treat.origin <- all.treat.origin
    result$base <- base
  }
  if (treat.type == "continuous") {
    D.sample <- treat.info[["D.sample"]]
    label.name <- names(D.sample)
    result$D.sample <- D.sample
    result$label.name <- label.name
  }
  if (!is.null(treat.info[["ncols"]])) {
    result$ncols <- treat.info[["ncols"]]
  }
  result
}

.compute_density <- function(data, X, D, weights, treat.type,
                             all.treat = NULL, all.treat.origin = NULL) {
  if (is.null(weights)) {
    de <- density(data[, X])
  } else {
    suppressWarnings(de <- density(data[, X], weights = data[, "WEIGHTS"]))
  }

  treat_den <- NULL
  if (treat.type == "discrete") {
    treat_den <- list()
    for (char in all.treat) {
      if (is.null(weights)) {
        de.tr <- density(data[data[, D] == char, X])
      } else {
        suppressWarnings(
          de.tr <- density(data[data[, D] == char, X],
                           weights = data[data[, D] == char, "WEIGHTS"])
        )
      }
      treat_den[[all.treat.origin[char]]] <- de.tr
    }
  }

  list(de = de, treat_den = treat_den)
}

.compute_histograms <- function(data, X, D, weights, treat.type,
                                all.treat = NULL, all.treat.origin = NULL) {
  if (is.null(weights)) {
    hist.out <- hist(data[, X], breaks = 80, plot = FALSE)
  } else {
    suppressWarnings(
      hist.out <- hist(data[, X], data[, "WEIGHTS"],
                       breaks = 80, plot = FALSE)
    )
  }

  treat.hist <- NULL
  if (treat.type == "discrete") {
    n.hist <- length(hist.out$mids)
    treat.hist <- list()
    for (char in all.treat) {
      count1 <- rep(0, n.hist)
      treat_index <- which(data[, D] == char)
      for (i in 1:n.hist) {
        count1[i] <- sum(data[treat_index, X] >= hist.out$breaks[i] &
                         data[treat_index, X] < hist.out$breaks[i + 1])
      }
      count1[n.hist] <- sum(data[treat_index, X] >= hist.out$breaks[n.hist] &
                            data[treat_index, X] <= hist.out$breaks[n.hist + 1])
      treat.hist[[all.treat.origin[char]]] <- count1
    }
  }

  list(hist.out = hist.out, treat.hist = treat.hist)
}
