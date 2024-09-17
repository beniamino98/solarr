#' Gumbel random variable
#'
#' Gumbel density, distribution, quantile and random generator.
#'
#' @param x vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If `length(n) > 1`, the length is taken to be the number required.
#' @param location location parameter.
#' @param scale scale parameter.
#' @param log logical; if `TRUE`, probabilities are returned as `log(p)`.
#' @param log.p logical; if `TRUE`, probabilities p are given as `log(p)`.
#' @param lower.tail logical; if TRUE (default), probabilities are `P[X < x]` otherwise, `P[X > x]`.
#'
#' @references Gumbel distribution [\href{https://en.wikipedia.org/wiki/Gumbel_distribution}{W}].
#'
#' @examples
#' # Grid
#' x <- seq(-5, 5, 0.01)
#'
#' # Density function
#' p <- dgumbel(x, location = 0, scale = 1)
#' plot(x, p, type = "l")
#'
#' # Distribution function
#' p <- pgumbel(x, location = 0, scale = 1)
#' plot(x, p, type = "l")
#'
#' # Quantile function
#' qgumbel(0.1)
#' pgumbel(qgumbel(0.1))
#'
#' # Random Numbers
#' rgumbel(1000)
#' plot(rgumbel(1000), type = "l")
#'
#' @name dgumbel
#' @rdname dgumbel
#' @aliases dgumbel
#' @aliases pgumbel
#' @aliases qgumbel
#' @aliases rgumbel
#' @export
dgumbel <- function(x, location = 0, scale = 1, log = FALSE){
  # Standardized values
  z <- (x-location)/scale
  # Density
  probs <- (1/scale)*exp(-(z + exp(-z)))
  # Log probability
  if (log) {
    probs <- base::log(probs)
  }
  return(probs)
}


#' @export
#' @rdname dgumbel
pgumbel <- function(x, location = 0, scale = 1, log.p = FALSE, lower.tail = TRUE){
  # Standardized values
  z <- (x - location)*scale
  # Distribution
  probs <- exp(-exp(-z))
  # Lower tail
  if (!lower.tail) {
    probs <- 1 - probs
  }
  # Log-probability
  if (log.p) {
    probs <- base::log(probs)
  }
  return(probs)
}


#' @export
#' @rdname dgumbel
qgumbel <- function(p, location = 0, scale = 1, log.p = FALSE, lower.tail = TRUE) {
  probs <- p
  # Log probability
  if (log.p) {
    probs <- exp(probs)
  }
  # Lower tail
  if (!lower.tail) {
    probs <- 1 - probs
  }
  # Quantiles
  x <- location + scale*log(-log(probs))
  return(x)
}


#' @export
#' @rdname dgumbel
rgumbel <- function(n, location = 0, scale = 1){
  # Control
  if (length(n) > 1){
    n <- length(n)
  }
  # Simulated grades
  u <- runif(n, min = 0, max = 1)
  # Simulated values
  qgumbel(u, location, scale, log.p = FALSE, lower.tail = TRUE)
}

