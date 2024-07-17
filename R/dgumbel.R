#' Gumbel Random Variable
#'
#' Probability density function for a gumbel random variable
#'
#' @param x vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If `length(n) > 1`, the length is taken to be the number required.
#' @param mean vector of means.
#' @param scale vector of scale parameter.
#' @param log.p logical; if `TRUE`, probabilities p are given as `log(p)`.
#' @param lower.tail logical; if TRUE (default), probabilities are `P[X < x]` otherwise, `P[X > x]`.
#' @param invert logical, use the inverted Gumbel distribution
#'
#' @examples
#' x <- seq(-5, 5, 0.01)
#'
#' # Density function
#' p <- dgumbel(x, mean = 0, scale = 1)
#' plot(x, p, type = "l")
#'
#' # Distribution function
#' p <- pgumbel(x, mean = 0, scale = 1)
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
#' @name gumbel
#' @rdname gumbel
#' @aliases dgumbel
#' @aliases pgumbel
#' @aliases qgumbel
#' @aliases rgumbel
#' @export
dgumbel <- function(x, mean = 0, scale = 1, log.p = FALSE, invert = FALSE){
  z <- (x-mean)/scale
  if (invert) {
    # Inverted Gumbel
    p <- (1/scale)*exp(z - exp(z))
  } else {
    # Gumbel
    p <- (1/scale)*exp(-(z + exp(-z)))
  }
  if (log.p) {
    return(base::log(p))
  }
  return(p)
}


#' @export
#' @rdname gumbel
pgumbel <- function(x, mean = 0, scale = 1, log.p = FALSE, lower.tail = TRUE, invert = FALSE){
  z <- (x - mean)*scale
  if (invert) {
    # Inverted Gumbel
    p <- exp(-exp(z))
  } else {
    # Gumbel
    p <- exp(-exp(-z))
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  if (log.p) {
    p <- base::log(p)
  }
  return(p)
}


#' @export
#' @rdname gumbel
qgumbel <- function(p, mean = 0, scale = 1, log.p = FALSE, lower.tail = TRUE, invert = FALSE) {
  if (log.p) {
    p <- exp(p)
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  x <- mean + ifelse(invert, 1, -1)*scale*log(-log(p))
  return(x)
}


#' @export
#' @rdname gumbel
rgumbel <- function(n, mean = 0, scale = 1, invert = FALSE){
  u <- runif(n, min = 0, max = 1)
  qgumbel(u, mean = mean, scale = scale, log.p = FALSE, lower.tail = TRUE, invert = invert)
}

