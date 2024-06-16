#' Kumaraswamy Random Variable
#'
#' Probability functions for a Kumaraswamy random variable
#'
#' @param x vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If `length(n) > 1`, the length is taken to be the number required.
#' @param a parameter.
#' @param b parameter..
#' @param log.p logical; if `TRUE`, probabilities p are given as `log(p)`.
#' @param lower.tail logical; if TRUE (default), probabilities are `P[X < x]` otherwise, `P[X > x]`.
#'
#' @name kumaraswamy
#' @rdname kumaraswamy
#' @aliases dkumaraswamy
#' @aliases pkumaraswamy
#' @aliases qkumaraswamy
#' @aliases rkumaraswamy
#' @export
dkumaraswamy <- function(x, a = 1, b = 1, log.p = FALSE){
  p <- a*b*x^(a - 1)*(1 - x^a)^(b - 1)
  p[x<0|x>1] <- 0
  if (log.p){
    p <- base::log(p)
  }
  return(p)
}

#' @export
#' @rdname kumaraswamy
pkumaraswamy <- function(x, a = 1, b = 1, log.p = FALSE, lower.tail = TRUE){
  x[x<=0] <- 0
  x[x>=1] <- 1
  p <- 1 - (1 - x^a)^b
  if (!lower.tail){
    p <- 1 - p
  }
  if (log.p){
    p <- base::log(p)
  }
  return(p)
}

#' @export
#' @rdname kumaraswamy
qkumaraswamy <- function(p, a = 1, b = 1, log.p = FALSE, lower.tail = TRUE){
  if (log.p) {
    p <- exp(p)
  }
  if (!lower.tail){
    p <- 1 - p
  }
  p[p<=0] <- 0
  p[p>=1] <- 1
  x <- (1 - (1-p)^(1/b))^(1/a)
  return(x)
}

#' @export
#' @rdname kumaraswamy
rkumaraswamy <- function(n, a = 1, b = 1){
  u <- runif(n, min = 0, max = 1)
  x <- qkumaraswamy(u, a = a, b = b, log.p = FALSE, lower.tail = TRUE)
  return(x)
}





