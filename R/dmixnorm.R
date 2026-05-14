#' Gaussian mixture distribution
#' 
#' @description
#' `r lifecycle::badge("stable")` 
#' Density, distribution function, quantile function, and random generation for
#' a univariate Gaussian mixture.
#' 
#' @param x Numeric vector of quantiles.
#' @param q Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param n Number of observations.
#' @param mean Numeric vector of component means.
#' @param sd Numeric vector of component standard deviations.
#' @param alpha Numeric vector of component probabilities.
#' @param log.p Logical. If `TRUE`, probabilities are supplied or returned on
#'   the log scale.
#' @param log Logical. If `TRUE`, `dmixnorm()` returns log-densities.
#' @param lower.tail Logical. If `TRUE`, probabilities are \eqn{P[X \le x]};
#'   otherwise, \eqn{P[X > x]}.
#'
#' @return
#' - `dmixnorm()` returns a numeric vector of density values.
#' - `pmixnorm()` returns a numeric vector of probabilities.
#' - `qmixnorm()` returns a numeric vector of quantiles.
#' - `rmixnorm()` returns a tibble with simulated component values, component
#'   indicators, and the combined random draw.
#'
# @references Mixture Models [\href{https://en.wikipedia.org/wiki/Mixture_model}{W}].
#'
#' @examples
#' mean <- c(-1, 1)
#' sd <- c(0.5, 1)
#' alpha <- c(0.4, 0.6)
#' dmixnorm(c(-1, 0, 1), mean, sd, alpha)
#' pmixnorm(c(-1, 0, 1), mean, sd, alpha)
#' qmixnorm(c(0.25, 0.75), mean, sd, alpha)
#'
#' set.seed(1)
#' rmixnorm(3, mean, sd, alpha)
#'
#' @rdname dmixnorm
#' @name dmixnorm
#' @aliases dmixnorm
#' @aliases pmixnorm
#' @aliases qmixnorm
#' @aliases rmixnorm
#' @keywords distributions
#' @family distributions
#' @note Version 1.0.1.
#' @export
dmixnorm <- function(x, mean = rep(0, 2), sd = rep(1, 2), alpha = rep(1/2, 2), log = FALSE){
  # Number of observations
  N <- length(x)
  # Number of components
  K <- length(mean)
  # Initialize a matrix for the components
  f <- matrix(0, nrow = N, ncol = K)
  for(k in 1:K){
    f[,k] <- alpha[k]*dnorm(x, mean = mean[k], sd = sd[k])
  }
  # Density
  p <- rowSums(f)
  # Log-probability
  if (log) {
    p <- base::log(p)
  }
  return(p)
}

#' @rdname dmixnorm
#' @export
pmixnorm <- function(q, mean = rep(0, 2), sd = rep(1, 2), alpha = rep(1/2, 2), lower.tail = TRUE, log.p = FALSE){
  # Number of observations
  N <- length(q)
  # Number of components
  K <- length(mean)
  # Initialize a matrix for the components
  f <- matrix(0, nrow = N, ncol = K)
  for(k in 1:K){
    f[,k] <- alpha[k]*pnorm(q, mean = mean[k], sd = sd[k], lower.tail = lower.tail)
  }
  # Density
  p <- rowSums(f)
  # Log-probability
  if (log.p) {
    p <- base::log(p)
  }
  return(p)
}

#' @rdname dmixnorm
#' @export
qmixnorm <- function(p, mean = rep(0, 2), sd = rep(1, 2), alpha = rep(1/2, 2), lower.tail = TRUE, log.p = FALSE) {
  # Log-probabilities
  if (log.p) {
    p <- base::exp(p)
  }
  # Distribution function
  cdf <- function(x) pmixnorm(x, mean, sd, alpha, lower.tail = lower.tail)
  # Empirical quantile
  quantile_numeric <- Quantile(cdf, interval = c(min(mean) - max(sd)*10, max(mean) + max(sd)*10))
  # Quantiles
  q <- quantile_numeric(p)
  #q <- p
  #q[p <= 0] <- -Inf
  #q[p >= 1] <- Inf
  #q[p > 0 & p < 1] <- quantile_numeric(p[p > 0 & p < 1])
  return(q)
}

#' @rdname dmixnorm
#' @export
rmixnorm <- function(n, mean = rep(0, 3), sd = rep(1, 3), alpha = rep(1/3, 3)){
  # Number of components
  k <- length(mean)
  X <- matrix(NA, nrow = n, ncol = k)
  B <- matrix(0, nrow = n, ncol = k)
  index <- 1:n
  for(s in 1:k){
    # Simulated bernoulli
    if (s == k) {
      B[index,][,s] <- 1
    } else {
      B[index,][,s] <- rbinom(n, 1, alpha[s]/sum(alpha[s:k]))
    }
    # Simulated component
    X[B[,s] == 1, s] <- rnorm(sum(B[, s]), mean = mean[s], sd = sd[s])
    # Update number of remaining elements
    n <- n - sum(B[,s])
    # Update the remaining indexes
    index <- index[!(index %in% which(!is.na(X[,s])) )]
    # Substitue NA values with 0
    X[,s] <- ifelse(is.na(X[,s]), 0, X[,s])
  }
  colnames(X) <- paste0("X", 1:k)
  colnames(B) <- paste0("B", 1:k)
  sim <- dplyr::bind_cols(t = 1:nrow(X), X, B, X = rowSums(X))
  return(sim)
}
