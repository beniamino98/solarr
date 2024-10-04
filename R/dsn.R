#' Skewed Normal random variable
#'
#' Skewed Normal density, distribution, quantile and random generator.
#'
#' @param x vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If `length(n) > 1`, the length is taken to be the number required.
#' @param location location parameter.
#' @param scale scale parameter.
#' @param shape skewness parameter.
#' @param log.p logical; if `TRUE`, probabilities p are given as `log(p)`.
#' @param log logical; if `TRUE`, probabilities are returned as `log(p)`.
#' @param lower.tail logical; if TRUE (default), probabilities are `P[X < x]` otherwise, `P[X > x]`.
#'
#' @references Skewed Normal Distribution [\href{https://en.wikipedia.org/wiki/Skew_normal_distribution}{W}].
#'
#' @examples
#' # Grid of points
#' x <- seq(-5, 5, 0.01)
#'
#' # Density function (right tailed)
#' p <- dsnorm(x, shape = 4.9)
#' plot(x, p, type = "l")
#'
#' # Density function (left tailed)
#' p <- dsnorm(x, shape = -4.9)
#' plot(x, p, type = "l")
#'
#' # Distribution function
#' p <- psnorm(x)
#' plot(x, p, type = "l")
#'
#' # Quantile function
#' dsnorm(0.1)
#' psnorm(qsnorm(0.9))
#'
#' # Random numbers
#' plot(rsnorm(100), type = "l")
#'
#' @name dsnorm
#' @rdname dsnorm
#' @aliases dsnorm
#' @aliases psnorm
#' @aliases qsnorm
#' @aliases rsnorm
#' @export
dsnorm <- function(x, location = 0, scale = 1, shape = 0, log = FALSE){
  # Standardized values
  z <- (x - location)/scale
  # Probabilities
  probs <- 2*dnorm(z)*pnorm(shape*z)
  # Log-probabilities
  if (log) {
    probs <- base::log(probs)
  }
  return(probs)
}


#' @rdname dsnorm
#' @export
psnorm <- function(x, location = 0, scale = 1, shape = 0, log.p = FALSE, lower.tail = TRUE){
  # Standardized values
  z <- (x - location)/scale
  # Distribution function
  cdf <- CDF(dsnorm, location = location, scale = scale, shape = shape, lower=-Inf, log = FALSE)
  # Cumulated probabilities
  probs <- cdf(z)

  # Lower tail
  if (!lower.tail) {
    probs <- 1 - probs
  }
  # Log-probabilities
  if (log.p) {
    probs <- base::log(probs)
  }
  return(probs)
}


#' @rdname dsnorm
#' @export
qsnorm <- function(p, location = 0, scale = 1, shape = 0, log.p = FALSE, lower.tail = TRUE) {
  # Distribution function
  cdf <- function(x) psnorm(x, location = location, scale = scale, shape = shape)
  # Quantile function
  qfunct <- Quantile(cdf, lower = -Inf)
  # Quantiles
  x <- qfunct(p, log.p = log.p, lower.tail = lower.tail)
  return(x)
}


#' @rdname dsnorm
#' @export
rsnorm <- function(n, location = 0, scale = 1, shape = 0){
  # Control
  if (length(n) > 1){
    n <- length(n)
  }
  u <- runif(n, min = 0, max = 1)
  qsnorm(u, location, scale, shape, log.p = FALSE, lower.tail = TRUE)
}

