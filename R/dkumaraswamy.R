#' Kumaraswamy distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the Kumaraswamy distribution.
#'
#' @param x Numeric vector of quantiles.
#' @param q Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param n Number of observations.
#' @param a Numeric shape parameter. Must be positive.
#' @param b Numeric shape parameter. Must be positive.
#' @param log Logical. If `TRUE`, `dkumaraswamy()` returns log-densities.
#' @param log.p Logical. If `TRUE`, probabilities are supplied or returned on
#'   the log scale.
#' @param lower.tail Logical. If `TRUE`, probabilities are \eqn{P[X \le x]};
#'   otherwise, \eqn{P[X > x]}.
#'
#' @return
#' - `dkumaraswamy()` returns a numeric vector of density values.
#' - `pkumaraswamy()` returns a numeric vector of probabilities.
#' - `qkumaraswamy()` returns a numeric vector of quantiles.
#' - `rkumaraswamy()` returns a numeric vector of random draws.
#'
# @references Kumaraswamy Distribution \href{https://en.wikipedia.org/wiki/Kumaraswamy_distribution}{W}.
#'
#' @examples
#' dkumaraswamy(c(0.25, 0.5), a = 2, b = 1.5)
#' pkumaraswamy(c(0.25, 0.5), a = 2, b = 1.5)
#' qkumaraswamy(c(0.25, 0.75), a = 2, b = 1.5)
#'
#' set.seed(1)
#' rkumaraswamy(3, a = 2, b = 1.5)
#' @name dkumaraswamy
#' @rdname dkumaraswamy
#' @aliases dkumaraswamy
#' @aliases pkumaraswamy
#' @aliases qkumaraswamy
#' @aliases rkumaraswamy
#' 
#' @keywords distributions
#' @family distributions
#' @note Version 1.0.0.
#' 
#' @export
dkumaraswamy <- function(x, a = 1, b = 1, log = FALSE){
  # Density function
  p <- a*b*x^(a - 1)*(1 - x^a)^(b - 1)
  # Ensure bounds
  p[x<0|x>1] <- 0
  # Log-probability
  if (log) {
    p <- base::log(p)
  }
  return(p)
}

#' @export
#' @rdname dkumaraswamy
pkumaraswamy <- function(q, a = 1, b = 1, log.p = FALSE, lower.tail = TRUE){
  # Distribution function
  p <- 1 - (1 - q^a)^b
  # Ensure bounds
  p[q<0|q>1] <- 0
  # Lower tail
  if (!lower.tail) {
    p <- 1 - p
  }
  # Log-probability
  if (log.p) {
    p <- base::log(p)
  }
  return(p)
}

#' @export
#' @rdname dkumaraswamy
qkumaraswamy <- function(p, a = 1, b = 1, log.p = FALSE, lower.tail = TRUE){
  # Log-probability
  if (log.p) {
    p <- exp(p)
  }
  # Lower tail
  if (!lower.tail) {
    p <- 1 - p
  }
  # Ensure bounds
  p[p<=0] <- 0
  p[p>=1] <- 1
  # Quantiles
  q <- (1 - (1-p)^(1/b))^(1/a)
  return(q)
}

#' @export
#' @rdname dkumaraswamy
rkumaraswamy <- function(n, a = 1, b = 1){
  # Simulated grades
  u <- runif(n, min = 0, max = 1)
  # Simulated quantiles
  q <- qkumaraswamy(u, a = a, b = b, log.p = FALSE, lower.tail = TRUE)
  return(q)
}

