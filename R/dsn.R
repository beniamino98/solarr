#' Skewed Normal
#'
#' Probability for a skewed normal random variable.
#'
#' @param x vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If `length(n) > 1`, the length is taken to be the number required.
#' @param mean vector of means.
#' @param sd vector of standard deviations.
#' @param skew vector of skewness.
#' @param log.p logical; if `TRUE`, probabilities p are given as `log(p)`.
#' @param log logical; if `TRUE`, probabilities are returned as `log(p)`.
#' @param lower.tail logical; if TRUE (default), probabilities are `P[X < x]` otherwise, `P[X > x]`.
#'
#' @examples
#' x <- seq(-5, 5, 0.01)
#' # Density function
#' p <- dsnorm(x, mean = 0, sd = 1)
#' plot(x, p, type = "l")
#' # Distribution function
#' p <- psnorm(x, mean = 0, sd = 1)
#' plot(x, p, type = "l")
#' # Quantile function
#' dsnorm(0.1)
#' psnorm(qsnorm(0.1))
#' # Random numbers
#' rsnorm(1000)
#' plot(rsnorm(1000), type = "l")
#'
#' @name snorm
#' @rdname snorm
#' @aliases dsnorm
#' @aliases psnorm
#' @aliases qsnorm
#' @aliases rsnorm
#' @export
dsnorm <- function(x, mean = 0, sd = 1, skew = 0, log = FALSE){
  z <- (x - mean)/sd
  p <- 2*dnorm(z)*pnorm(skew*z)

  if (log) {
    return(base::log(p))
  }
  return(p)
}


#' @rdname snorm
#' @export
psnorm <- function(x, mean = 0, sd = 1, skew = 0, log.p = FALSE, lower.tail = TRUE){
  z <- (x - mean)/sd
  p <- c()
  for(i in 1:length(z)){
    p[i] <- integrate(dsnorm, lower=-Inf, upper = z[i], skew = skew)$value
  }
  if (!lower.tail) {
    p <- 1 - p
  }

  if (log.p) {
    return(base::log(p))
  }
  return(p)
}


#' @rdname snorm
#' @export
qsnorm <- function(p, mean = 0, sd = 1, skew = 0, log.p = FALSE, lower.tail = TRUE) {

  if (log.p) {
    p <- exp(p)
  }

  loss <- function(x, p) {
    p_hat <- psnorm(x, mean, sd, skew, lower.tail = lower.tail)
    (p_hat - p)^2
  }
  x <- c()
  for(i in 1:length(p)){
    x[i] <- suppressWarnings(optim(par = mean, loss, p = p[i])$par)
  }
  return(x)
}


#' @rdname snorm
#' @export
rsnorm <- function(n, mean = 0, sd = 1, skew = 0){
  u <- runif(n, min = 0, max = 1)
  qsnorm(u, mean = mean, sd = sd, skew = skew)
}


