#' Construct density, distribution, and quantile functions
#'  
#' @description
#' `r lifecycle::badge("experimental")`
#' Return functions with density-like, distribution-like, or quantile-like
#' interfaces from a supplied density or distribution function.
#'
#' @param .f Function. Density function used by `PDF()` or integrated by
#'   `CDF()`.
#' @param cdf Function. Cumulative distribution function used by `Quantile()`.
#' @param lower Numeric. Lower integration bound used by `CDF()`.
#' @param interval Numeric vector of length two. Search interval used by
#'   `Quantile()` when solving for roots.
#' @param ... Additional arguments passed to `.f`.
#'
#' @return
#' - `PDF()` returns a function with arguments `x` and `log`.
#' - `CDF()` returns a function with arguments `x`, `lower.tail`, and `log.p`.
#' - `Quantile()` returns a function with arguments `p`, `log.p`, and
#'   `lower.tail`.
#'
#' @examples
#' pdf <- PDF(dnorm, mean = 0.3, sd = 1.3)
#' pdf(3)
#'
#' cdf <- CDF(dnorm, mean = 0.3, sd = 1.3)
#' cdf(3)
#'
#' pnorm(Quantile(pnorm)(0.9))
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
Quantile <- function(cdf, interval = c(-100, 100)){

  # Find the quantile numerically
  quantile_root <- function(p, cdf, interval){
    uniroot(function(x) cdf(x) - p,
            interval = interval,
            tol = 10^{-16})$root
  }
  # Quantile function
  quantile_numeric <- function(p) purrr::map_dbl(p, ~quantile_root(.x, cdf, interval))

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
    x <- quantile_numeric(probs)
    return(x)
  }
}




