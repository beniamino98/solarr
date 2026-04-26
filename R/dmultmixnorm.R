#' Density of a multinomial Gaussian mixture with 2 components each
#'
#' @param x matrix, on the rows the observations, on the columns the locations.
#' @param means matrix of means for each state.
#' @param Sigma List of variance-covariance matrices for each state.
#' @param probs vector of joint probabilities for each state.
#' @param log Logical, when `TRUE`, probabilities are returned as `log(p)`.
#' @param log.p Logical, when `TRUE`, probabilities p are given as `log(p)`.
#' @param lower.tail Logical, when `TRUE`, the default, the computed probabilities are \eqn{\mathbb{P}(X < x)}, otherwise \eqn{\mathbb{P}(X \ge x)}.
#' @examples
#' library(tidyverse)
#' x <- matrix(runif(200, -3, 3), ncol = 2)
#'
#' system.time(dmix2norm(x), gcFirst = F)
#' system.time(map(1:500, ~dmultmixnorm(x)), gcFirst = F)
#'
#' system.time(map(1:5, ~pmix2norm(x)), gcFirst = F)
#' system.time(map(1:5, ~pmultmixnorm(x)), gcFirst = F)
#'
#' qmix2norm(c(0.1, 0.2), x0 = c(-1, 1))
#' qmultmixnorm(c(0.1, 0.2))
#'
#' @rdname dmultmixnorm
#' @aliases dmultmixnorm
#' @aliases pmultmixnorm
#' @aliases qmultmixnorm
#' @aliases d2mixnorm
#' @aliases p2mixnorm
#' @aliases q2mixnorm
#' @keywords distributions
#' @details Version 1.0.0.
#' @export
dmultmixnorm <- function(x, means = matrix(0, 4, 2), Sigma = rep(list(diag(1, 2)), 4), probs = rep(1/4, 4), log = FALSE){
  # Number of observations
  N <- nrow(x)
  # Number of locations
  L <- ncol(x)
  # Number of components
  K <- 2^L
  # Compute density for each state
  f <- matrix(0, nrow = N, ncol = K)
  for(k in 1:K){
    f[, k] <- probs[k] * mvtnorm::dmvnorm(x, mean = means[k,], sigma = Sigma[[k]])
  }
  # Density
  p <- rowSums(f)

  # Log-probability
  if (log) {
    p <- base::log(p)
  }

  return(p)
}

#' Distribution of a multinomial Gaussian mixture with 2 components each
#'
#' @rdname dmultmixnorm
#' @export
pmultmixnorm <- function(x, means = matrix(0, 4, 2), Sigma = rep(list(diag(1, 2)), 4), probs = rep(1/4, 4), log = FALSE){
  # Number of observations
  N <- nrow(x)
  # Number of locations
  L <- ncol(x)
  # Number of components
  K <- 2^L
  # Lower limit
  lower <- rep(-Inf, L)
  # Helper
  Phi <- function(x, mean, sigma) purrr::map_dbl(1:nrow(x), ~mvtnorm::pmvnorm(lower = lower, upper = x[.x,], mean = mean, sigma = sigma))
  # Compute density for each state
  f <- matrix(0, nrow = N, ncol = K)
  for(k in 1:K){
    f[, k] <- probs[k] * Phi(x, means[k,], Sigma[[k]])
  }
  # Total density
  p <- rowSums(f)
  # Log-probability
  if (log) {
    p <- base::log(p)
  }
  return(p)
}

#' Quantile of a multinomial Gaussian mixture with 2 components each
#'
#' @rdname dmultmixnorm
#' @export
qmultmixnorm <- function(p, means = matrix(0, 4, 2), Sigma = rep(list(diag(1, 2)), 4), probs = rep(1/4, 4), x_init, log.p = FALSE, lower.tail = TRUE){
  # Log-probability
  if (log.p) {
    p <- exp(p)
  }
  # Lower tail
  if (!lower.tail) {
    p <- 1 - p
  }
  # Number of joint pairs
  N <- length(p)
  # Number of locations
  L <- ncol(means)
  # Loss function
  loss <- function(x, p_target = 0.05){
    x <- matrix(x, nrow = 1, ncol = L)
    cdf_x <- pmultmixnorm(x, means, Sigma, probs)
    sum((cdf_x - p_target)^2)
  }
  # Initial values
  if (missing(x_init)) {
    x_init <- matrix(0, nrow = N, ncol = L)
  }
  # Numerical quantiles
  col.names <- paste0("r", 1:L)
  x <- purrr:::map_df(1:N, ~setNames(optim(x_init[.x,], loss, p_target = p[.x])$par, col.names))
  x <- as.matrix(x)
  dimnames(x) <- NULL
  x
}
