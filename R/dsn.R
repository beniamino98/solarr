#' Skew-normal distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the skew-normal distribution.
#'
#' @param x Numeric vector of quantiles.
#' @param q Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param n Number of observations.
#' @param location Numeric location parameter.
#' @param scale Numeric scale parameter.
#' @param shape Numeric shape parameter controlling skewness.
#' @param log.p Logical. If `TRUE`, probabilities are supplied or returned on
#'   the log scale.
#' @param log Logical. If `TRUE`, `dsnorm()` returns log-densities.
#' @param lower.tail Logical. If `TRUE`, probabilities are \eqn{P[X \le x]};
#'   otherwise, \eqn{P[X > x]}.
#'
#' @return
#' - `dsnorm()` returns a numeric vector of density values.
#' - `psnorm()` returns a numeric vector of probabilities.
#' - `qsnorm()` returns a numeric vector of quantiles.
#' - `rsnorm()` returns a numeric vector of random draws.
#'
# @references Skewed Normal Distribution [\href{https://en.wikipedia.org/wiki/Skew_normal_distribution}{W}].
#'
#' @examples
#' dsnorm(c(-1, 0, 1), shape = 2)
#' psnorm(c(-1, 0, 1), shape = 2)
#' qsnorm(c(0.25, 0.75), shape = 2)
#'
#' set.seed(1)
#' rsnorm(3, shape = 2)
#' @name dsnorm
#' @rdname dsnorm
#' @aliases dsnorm
#' @aliases psnorm
#' @aliases qsnorm
#' @aliases rsnorm
#' @keywords distributions
#' @details Version 1.0.0.
#' @export
dsnorm <- function(x, location = 0, scale = 1, shape = 0, log = FALSE){
  # Standardized values
  z <- (x - location)/scale
  # Probabilities
  p <- 2*dnorm(z)*pnorm(shape*z)
  # Log-probabilities
  if (log) {
    p <- base::log(p)
  }
  return(p)
}


#' @rdname dsnorm
#' @export
psnorm <- function(q, location = 0, scale = 1, shape = 0, log.p = FALSE, lower.tail = TRUE){
  # Standardized values
  z <- (q - location)/scale
  # Distribution function
  cdf <- CDF(dsnorm, location = location, scale = scale, shape = shape, lower=-Inf, log = FALSE)
  # Cumulated probabilities
  p <- cdf(z)
  # Lower tail
  if (!lower.tail) {
    p <- 1 - p
  }
  # Log-probabilities
  if (log.p) {
    p <- base::log(p)
  }
  return(p)
}


#' @rdname dsnorm
#' @export
qsnorm <- function(p, location = 0, scale = 1, shape = 0, log.p = FALSE, lower.tail = TRUE) {
  # Distribution function
  cdf <- function(x) psnorm(x, location = location, scale = scale, shape = shape)
  # Quantile function
  quantile_numeric <- Quantile(cdf, interval = c(-location - scale*10, location + scale*10))
  # Quantiles
  q <- quantile_numeric(p, log.p = log.p, lower.tail = lower.tail)
  return(q)
}


#' @rdname dsnorm
#' @export
rsnorm <- function(n, location = 0, scale = 1, shape = 0){
  # Simulated grades
  u <- runif(n, min = 0, max = 1)
  # Quantiles
  q <- qsnorm(u, location, scale, shape, log.p = FALSE, lower.tail = TRUE)
  return(q)
}
