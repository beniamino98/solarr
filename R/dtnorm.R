#' Truncated normal distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' a normal distribution truncated to an interval.
#'
#' @param x Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param n Number of observations.
#' @param mean Numeric mean parameter.
#' @param sd Numeric standard deviation parameter.
#' @param a Numeric lower truncation bound.
#' @param b Numeric upper truncation bound.
#' @param log.p Logical. If `TRUE`, probabilities are supplied or returned on
#'   the log scale.
#' @param log Logical. If `TRUE`, `dtnorm()` returns log-densities.
#' @param lower.tail Logical. If `TRUE`, probabilities are \eqn{P[X \le x]};
#'   otherwise, \eqn{P[X > x]}.
#'
#' @return
#' - `dtnorm()` returns a numeric vector of density values.
#' - `ptnorm()` returns a numeric vector of probabilities.
#' - `qtnorm()` returns a numeric vector of quantiles.
#' - `rtnorm()` returns a numeric vector of random draws.
#'
#' @examples
#' dtnorm(c(-1, 0, 1), mean = 0, sd = 1, a = -2, b = 2)
#' ptnorm(c(-1, 0, 1), mean = 0, sd = 1, a = -2, b = 2)
#' qtnorm(c(0.25, 0.75), mean = 0, sd = 1, a = -2, b = 2)
#'
#' set.seed(1)
#' rtnorm(3, mean = 0, sd = 1, a = -2, b = 2)
#' @name dtnorm
#' @rdname dtnorm
#' @aliases dtnorm
#' @aliases ptnorm
#' @aliases qtnorm
#' @aliases rtnorm
#' 
#' @keywords distributions
#' @family distributions
#' @note Version 1.0.0.
#' @export
dtnorm <- function(x, mean = 0, sd = 1, a = -3, b = 3, log = FALSE){
  x[x < a | x > b] <- NA
  z <- (x - mean)/sd
  p <- (1/sd)*(dnorm(z)/(pnorm((b - mean)/sd) - pnorm((a - mean)/sd)))
  p[which(is.na(p))] <- 0

  if (log){
    return(base::log(p))
  }
  return(p)
}


#' @rdname dtnorm
#' @export
ptnorm <- function(x, mean = 0, sd = 1, a = -3, b = 3, log.p = FALSE, lower.tail = TRUE){
  z <- (x - mean)/sd
  p <- c()
  for(i in 1:length(z)){
    p[i] <- integrate(dtnorm, lower=-Inf, upper = z[i], a = a, b = b)$value
  }

  if (!lower.tail) {
    p <- 1 - p
  }

  if (log.p) {
    return(base::log(p))
  }
  return(p)
}


#' @export
#' @rdname dtnorm
qtnorm <- function(p, mean = 0, sd = 1, a = -3, b = 3, log.p = FALSE, lower.tail = TRUE) {

  if (log.p) {
    p <- exp(p)
  }

  loss <- function(x, p) {
    p_hat <- ptnorm(x, mean, sd, a = a, b = b, lower.tail = lower.tail)
    (p_hat - p)^2
  }
  x <- c()
  for(i in 1:length(p)){
    x[i] <- suppressWarnings(optim(par = (a+b)/2, loss, p = p[i])$par)
  }
  return(x)
}


#' @export
#' @rdname dtnorm
rtnorm <- function(n, mean = 0, sd = 1, a = -100, b = 100){
  u <- runif(n, min = 0, max = 1)
  x <- qtnorm(u, mean = mean, sd = sd, a = a, b = b)
  x[x < a] = a
  x[x > b] = b
  x
}

