#' Inverted Gumbel distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the inverted Gumbel distribution.
#'
#' @param x Numeric vector of quantiles.
#' @param q Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param n Number of observations.
#' @param location Numeric location parameter.
#' @param scale Numeric scale parameter.
#' @param log Logical. If `TRUE`, `dinvgumbel()` returns log-densities.
#' @param log.p Logical. If `TRUE`, probabilities are supplied or returned on
#'   the log scale.
#' @param lower.tail Logical. If `TRUE`, probabilities are \eqn{P[X \le x]};
#'   otherwise, \eqn{P[X > x]}.
#'
#' @return
#' - `dinvgumbel()` returns a numeric vector of density values.
#' - `pinvgumbel()` returns a numeric vector of probabilities.
#' - `qinvgumbel()` returns a numeric vector of quantiles.
#' - `rinvgumbel()` returns a numeric vector of random draws.
#'
#' @examples
#' dinvgumbel(c(-1, 0, 1), location = 0, scale = 1)
#' pinvgumbel(c(-1, 0, 1), location = 0, scale = 1)
#' qinvgumbel(c(0.25, 0.75), location = 0, scale = 1)
#'
#' set.seed(1)
#' rinvgumbel(3, location = 0, scale = 1)
#' @name dinvgumbel
#' @rdname dinvgumbel
#' @aliases dinvgumbel
#' @aliases pinvgumbel
#' @aliases qinvgumbel
#' @aliases rinvgumbel
#' @keywords distributions
#' @details Version 1.0.0.
#' @export
dinvgumbel <- function(x, location = 0, scale = 1, log = FALSE){
  # Standardized values
  z <- (x-location)/scale
  # Density
  p <- (1/scale)*exp(z - exp(z))
  # Log probability
  if (log) {
    p <- base::log(p)
  }
  return(p)
}

#' @export
#' @rdname dinvgumbel
pinvgumbel <- function(q, location = 0, scale = 1, log.p = FALSE, lower.tail = TRUE){
  # Standardized values
  z <- (q - location)*scale
  # Distribution
  p <- 1-exp(-exp(z))
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
#' @rdname dinvgumbel
qinvgumbel <- function(p, location = 0, scale = 1, log.p = FALSE, lower.tail = TRUE) {
  # Log probability
  if (log.p) {
    p <- exp(p)
  }
  # Lower tail
  if (!lower.tail) {
    p <- 1 - p
  }
  # Quantiles
  q <- -location + scale*log(-log(1-p))
  return(q)
}

#' @export
#' @rdname dinvgumbel
rinvgumbel <- function(n, location = 0, scale = 1){
  # Simulated grades
  u <- runif(n, min = 0, max = 1)
  # Simulated values
  q <- qinvgumbel(u, location, scale, log.p = FALSE, lower.tail = FALSE)
  return(q)
}

