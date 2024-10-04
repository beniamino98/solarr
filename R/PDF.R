#' Density, distribution and quantile function
#'
#' Return a function of `x` given the specification of a function of `x`.
#'
#' @param .f density function
#' @param cdf cumulative distribution function.
#' @param lower lower bound for integration (domain).
#' @param ... other parameters to be passed to `.f`.
#'
#' @examples
#' # Density
#' pdf <- PDF(dnorm, mean = 0.3, sd = 1.3)
#' pdf(3)
#' dnorm(3, mean = 0.3, sd = 1.3)
#' # Distribution
#' cdf <- CDF(dnorm, mean = 0.3, sd = 1.3)
#' cdf(3)
#' pnorm(3, mean = 0.3, sd = 1.3)
#' # Numeric quantile function
#' pnorm(Quantile(dnorm)(0.9))
#' @name PDF
#' @rdname PDF
#' @aliases PDF
#' @aliases CDF
#' @aliases Quantile
#' @export
PDF <- function(.f, ...){
  function(x, log = FALSE){
    probs <- .f(x, ...)
    # Log-probabilities
    if (log) {
      probs <- base::log(probs)
    }
    return(probs)
  }
}

#' @rdname PDF
#' @export
CDF <- function(.f, lower = -Inf, ...){
  # Density
  pdf <- PDF(.f, ...)
  # Distribution
  cdf <- function(x, pdf) integrate(pdf, lower = lower, upper = x)$value
  # Distribution function
  function(x,  lower.tail = TRUE, log.p = FALSE) {
    probs <- purrr::map_dbl(x, ~cdf(.x, pdf))
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
}

#' @rdname PDF
#' @export
Quantile <- function(cdf, lower = -Inf, x0 = 0){
  # Loss function
  loss_function <- function(x, p) {(cdf(x) - p)^2}
  # Quantile function
  safe_optim <- purrr::quietly(purrr::safely(optim))
  quantile_ <- function(p) purrr::map_dbl(p, ~safe_optim(par = x0, loss_function, p = .x)$result$result$par)

  function(p, log.p = FALSE, lower.tail = TRUE){
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
    x <- quantile_(probs)
    return(x)
  }
}


