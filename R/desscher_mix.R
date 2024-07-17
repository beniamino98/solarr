#' Esscher transform of a Gaussian Mixture
#'
#' @param params Gaussian Mixture parameters, mu1, sigma1, mu2, sigma2, p.
#'
#' @examples
#' library(ggplot2)
#' grid <- seq(-5, 5, 0.01)
#' pdf_1 <- desscher_mix()(grid, h = 0)
#' pdf_2 <- desscher_mix()(grid, h = -0.5)
#' pdf_3 <- desscher_mix()(grid, h = 0.5)
#' ggplot()+
#'  geom_line(aes(grid, pdf_1), color = "black")+
#'  geom_line(aes(grid, pdf_2), color = "green")+
#'  geom_line(aes(grid, pdf_3), color = "red")
#'
#' cdf_1 <- pesscher_mix()(grid, h = 0)
#' cdf_2 <- pesscher_mix()(grid, h = -0.5)
#' cdf_3 <- pesscher_mix()(grid, h = 0.5)
#' ggplot()+
#'   geom_line(aes(grid, cdf_1), color = "black")+
#'   geom_line(aes(grid, cdf_2), color = "green")+
#'   geom_line(aes(grid, cdf_3), color = "red")
#'
#' @aliases desscher_mix
#' @aliases pesscher_mix
#' @export

#' @rdname esscher
#' @export
desscher_mix <- function(params = c(0,0,1,1,0.5)){
  # Gaussian mixture Pdf
  dnorm_mixture <- dnorm_mix(params)
  # Esscher Numerator
  esscher_num <- function(x, h = 0) ifelse(is.infinite(exp(h*x)), 0, exp(h*x))
  # Esscher Denominator
  esscher_den <- function(x, h = 0) esscher_num(x, h)*dnorm_mixture(x)

  # Esscher pdf depending on `x` and `h`
  function(x, h = 0){
    den <- integrate(esscher_den, lower = -Inf, upper = Inf, h = h)$value
    esscher_num(x, h)*dnorm_mixture(x)/den
  }
}

#' @rdname esscher
#'
#' @export
pesscher_mix <- function(params = c(0,0,1,1,0.5)){
  # Esscher pdf
  pdf <- desscher_mix(params)

  # Esscher cdf depending on `x` and `h`
  cdf <- function(x, lower = -Inf, h = 0){
    purrr::map_dbl(x, ~integrate(function(x) pdf(x, h = h), lower = lower, upper = .x)$value)
  }
}
