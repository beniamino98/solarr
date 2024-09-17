#' Kumaraswamy random variable
#'
#' Kumaraswamy density, distribution, quantile and random generator.
#'
#' @param x vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If `length(n) > 1`, the length is taken to be the number required.
#' @param a parameter `a > 0`.
#' @param b parameter `b > 0`.
#' @param log logical; if `TRUE`, probabilities are returned as `log(p)`.
#' @param log.p logical; if `TRUE`, probabilities p are given as `log(p)`.
#' @param lower.tail logical; if `TRUE`, the default, the computed probabilities are `P[X < x]`. Otherwise, `P[X > x]`.
#'
#' @references Kumaraswamy Distribution [\href{https://en.wikipedia.org/wiki/Kumaraswamy_distribution}{W}].
#'
#' @examples
#' x <- seq(0, 1, 0.01)
#' # Density function
#' plot(x, dkumaraswamy(x, 0.2, 0.3), type = "l")
#' plot(x, dkumaraswamy(x, 2, 1.1), type = "l")
#' # Distribution function
#' plot(x, pkumaraswamy(x, 2, 1.1), type = "l")
#' # Quantile function
#' qkumaraswamy(0.2, 0.4, 1.4)
#' # Random generator
#' rkumaraswamy(20, 0.4, 1.4)
#'
#' @name dkumaraswamy
#' @rdname dkumaraswamy
#' @aliases dkumaraswamy
#' @aliases pkumaraswamy
#' @aliases qkumaraswamy
#' @aliases rkumaraswamy
#' @export
dkumaraswamy <- function(x, a = 1, b = 1, log = FALSE){
  # Density function
  probs <- a*b*x^(a - 1)*(1 - x^a)^(b - 1)
  # Ensure bounds
  probs[x<0|x>1] <- 0
  # Log-probability
  if (log) {
    probs <- base::log(probs)
  }
  return(probs)
}

#' @export
#' @rdname dkumaraswamy
pkumaraswamy <- function(x, a = 1, b = 1, log.p = FALSE, lower.tail = TRUE){
  # Distribution function
  probs <- 1 - (1 - x^a)^b
  # Ensure bounds
  probs[x<0|x>1] <- 0
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
#' @rdname dkumaraswamy
qkumaraswamy <- function(p, a = 1, b = 1, log.p = FALSE, lower.tail = TRUE){
  probs <- p
  # Log-probability
  if (log.p) {
    probs <- exp(probs)
  }
  # Lower tail
  if (!lower.tail){
    probs <- 1 - probs
  }
  # Ensure bounds
  probs[probs<=0] <- 0
  probs[probs>=1] <- 1
  # Quantiles
  x <- (1 - (1-probs)^(1/b))^(1/a)
  return(x)
}

#' @export
#' @rdname dkumaraswamy
rkumaraswamy <- function(n, a = 1, b = 1){
  # Control
  if (length(n) > 1){
    n <- length(n)
  }
  # Simulated grades
  u <- runif(n, min = 0, max = 1)
  # Simulated values
  x <- qkumaraswamy(u, a = a, b = b, log.p = FALSE, lower.tail = TRUE)
  return(x)
}

