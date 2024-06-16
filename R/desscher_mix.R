#' Esscher transform of a Gaussian Mixture
#'
#' @param params Gaussian Mixture parameters, mu1, sigma1, mu2, sigma2, p.
#'
#' @aliases desscher_mix
#' @aliases pesscher_mix
#' @export

#' @rdname esscher
#' @export
desscher_mix <- function(params = c(0,1,0,1,0.5)){
  # Mixture Pdf
  dnorm_mixture <- dnorm_mix(params)
  # Esscher Numerator
  esscher_num <- function(x, h = 0) ifelse(is.infinite(exp(h*x)), 0, exp(h*x))
  # Esscher Denominator
  esscher_den <- function(x, h = 0) esscher_num(x, h)*dnorm_mixture(x)
  # Esscher pdf depending on (x; h)
  function(x, h = 0){
    den <- integrate(esscher_den, lower = -Inf, upper = Inf, h = h)$value
    esscher_num(x, h)*dnorm_mixture(x)/den
  }
}

#' @rdname esscher
#' @export
pesscher_mix <- function(params = c(0,1,0,1,0.5)){
  # Density
  pdf <- esscher_mixture_pdf(params)
  function(upper, lower = -Inf, h = 0){
    integrate(function(x) pdf(x, h = h), lower = lower, upper = upper)
  }
}
