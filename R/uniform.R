calculate_uniform_quantiles <- function(theta_matrix, alpha) {
  # drop columns with NA
  cols_with_na <- apply(theta_matrix, 2, function(col) any(is.na(col)))
  theta_matrix <- theta_matrix[, !cols_with_na]
  
  k <- nrow(theta_matrix)  # number of parameters
  N <- ncol(theta_matrix)  # number of bootstrap draws
  
  check_condition <- function(zeta) {
    # For each parameter (row), find the zeta and 1-zeta quantiles across the columns
    Q_j <- t(apply(theta_matrix, 1, quantile, c(zeta, 1-zeta)))
    # Check if each bootstrap sample for each parameter is within [Q_j[,1], Q_j[,2]]
    condition_mat <- theta_matrix >= Q_j[, 1] & theta_matrix <= Q_j[, 2]
    # We want to see for each column (bootstrap draw), do *all* parameters fall within?
    condition_mat <- apply(condition_mat, 2, all)
    # Return the proportion of bootstrap draws for which all parameters are in the band
    return(mean(condition_mat))
  }
  
  # Start bisection search bounds
  zeta_lower <- alpha / (2 * k)
  zeta_upper <- alpha / 2
  
  # Bisection search to find zeta_hat
  while (zeta_upper - zeta_lower > 1e-6) {
    zeta_mid <- (zeta_lower + zeta_upper) / 2
    # If coverage < 1-alpha, we need a bigger band => increase zeta
    if (check_condition(zeta_mid) < 1 - alpha) {
      zeta_upper <- zeta_mid
    } else {
      zeta_lower <- zeta_mid
    }
  }
  
  zeta_hat <- zeta_lower
  coverage <- check_condition(zeta_hat)
  
  if (zeta_hat == alpha / (2 * k)) {
    warning("Insufficient bootstrap samples for Bootstrapped Uniform CI. Using Bonferroni CI by default.\n")
  }
  
  Q_j <- t(apply(theta_matrix, 1, quantile, c(zeta_hat, 1 - zeta_hat)))
  return(list(
    Q_j       = Q_j,
    zeta_hat  = zeta_hat,
    coverage  = coverage
  ))
}


calculate_delta_uniformCI <- function(Sigma_hat,alpha=0.05,N=2000){
    k <- nrow(Sigma_hat)
    V <- mvrnorm(N, mu = rep(0, k), Sigma = Sigma_hat)
    Sigma_inv_sqrt <- 1/(diag(Sigma_hat)^{0.5})
    max_vals <- apply(V, 1, function(v) {
        max(abs(Sigma_inv_sqrt * v))
    })
    q <- quantile(max_vals,1-alpha)
    return(q)
}