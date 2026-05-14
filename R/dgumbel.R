#' Gumbel distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the Gumbel distribution.
#'
#' @param x Numeric vector of quantiles.
#' @param q Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param n Number of observations.
#' @param location Numeric location parameter.
#' @param scale Numeric scale parameter.
#' @param log Logical. If `TRUE`, `dgumbel()` returns log-densities.
#' @param log.p Logical. If `TRUE`, probabilities are supplied or returned on
#'   the log scale.
#' @param lower.tail Logical. If `TRUE`, probabilities are \eqn{P[X \le x]};
#'   otherwise, \eqn{P[X > x]}.
#'
#' @return
#' - `dgumbel()` returns a numeric vector of density values.
#' - `pgumbel()` returns a numeric vector of probabilities.
#' - `qgumbel()` returns a numeric vector of quantiles.
#' - `rgumbel()` returns a numeric vector of random draws.
#'
#@references Gumbel distribution [\href{https://en.wikipedia.org/wiki/Gumbel_distribution}{W}].
#'
#' @examples
#' dgumbel(c(-1, 0, 1), location = 0, scale = 1)
#' pgumbel(c(-1, 0, 1), location = 0, scale = 1)
#' qgumbel(c(0.25, 0.75), location = 0, scale = 1)
#'
#' set.seed(1)
#' rgumbel(3, location = 0, scale = 1)
#' @name dgumbel
#' @rdname dgumbel
#' @aliases dgumbel
#' @aliases pgumbel
#' @aliases qgumbel
#' @aliases rgumbel
#' 
#' @keywords distributions
#' @family distributions
#' @note Version 1.0.0.
#' @export
dgumbel <- function(x, location = 0, scale = 1, log = FALSE){
  # Standardized values
  z <- (x-location)/scale
  # Density
  p <- (1/scale)*exp(-(z + exp(-z)))
  # Log probability
  if (log) {
    p <- base::log(p)
  }
  return(p)
}


#' @export
#' @rdname dgumbel
pgumbel <- function(q, location = 0, scale = 1, log.p = FALSE, lower.tail = TRUE){
  # Standardized values
  z <- (q - location)*scale
  # Distribution
  p <- exp(-exp(-z))
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
#' @rdname dgumbel
qgumbel <- function(p, location = 0, scale = 1, log.p = FALSE, lower.tail = TRUE) {
  # Log probability
  if (log.p) {
    p <- exp(p)
  }
  # Lower tail
  if (!lower.tail) {
    p <- 1 - p
  }
  # Quantiles
  q <- location - scale*log(-log(p))
  return(q)
}


#' @export
#' @rdname dgumbel
rgumbel <- function(n, location = 0, scale = 1){
  # Simulated grades
  u <- runif(n, min = 0, max = 1)
  # Simulated values
  q <- qgumbel(u, location, scale, log.p = FALSE, lower.tail = FALSE)
  return(q)
}

