#' Bernoulli test on the violations of the VaR
#'
#' @param et vector or matrix of violations. On the columns different VaR on the rows observations.
#' @param VaR vector of values at risk
#' @param ci confidence level for the test
#' @keywords diagnostic
#' @export
VaR_test <- function(et, VaR, ci = 0.01){
  e_mat <- 0
  if (is.vector(et)){
    e_mat <- matrix(et, ncol = 1)
  } else{
    stopifnot(is.matrix(et))
    e_mat <- et
  }
  # Number of observations
  n <- nrow(e_mat)
  # Vectorized
  VaR_hat <- apply(e_mat, 2, mean)
  # Test Statistics
  stat <- sqrt(n) * (apply(e_mat, 2, mean) - VaR) / sqrt(VaR * (1 - VaR))
  # p.values
  p.value <- 2 * pnorm(-abs(stat))
  # Critical values
  rejection_lev <- abs(qnorm(ci/2))
  tibble(
    n = n,
    ci = ci,
    stat = stat,
    VaR_hat = VaR_hat,
    VaR = VaR,
    rejection_lev = rejection_lev,
    p.value = p.value,
    H0 = ifelse(abs(stat) > rejection_lev, "Rejected", "Non-Rejected"),
  )
}

#' Violations of the VaR
#'
#' @param x vector, realized data
#' @param q vector or matrix, on the columns different VaR on the rows observations.
#' @keywords diagnostic
#' @export
VaR_viol <- function(x, q){
  stopifnot(is.vector(x))
  if (is.vector(q)){
    q_mat <- matrix(q, ncol = 1)
  } else{
    stopifnot(is.matrix(q))
    q_mat <- q
  }
  # Number of VaR
  n.VaR <- ncol(q_mat)
  # Number of observations
  n <- length(x)
  # Vectorized
  x_mat <- matrix(rep(x, n.VaR), nrow = n)
  # Violations
  e_mat <- ifelse(x_mat <= q_mat, 1, 0)
  e_mat
}
