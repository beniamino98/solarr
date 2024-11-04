#' Esscher transform of a Gaussian Mixture
#'
#' @inheritParams dmixnorm
#' @param theta Esscher parameter, the default is zero.
#'
#' @examples
#' library(ggplot2)
#' grid <- seq(-5, 5, 0.01)
#' # Density
#' pdf_1 <- desscherMixture(mean = c(-3, 3), theta = 0)(grid)
#' pdf_2 <- desscherMixture(mean = c(-3, 3), theta = -0.5)(grid)
#' pdf_3 <- desscherMixture(mean = c(-3, 3), theta = 0.5)(grid)
#' ggplot()+
#'  geom_line(aes(grid, pdf_1), color = "black")+
#'  geom_line(aes(grid, pdf_2), color = "green")+
#'  geom_line(aes(grid, pdf_3), color = "red")
#' # Distribution
#' cdf_1 <- pesscherMixture(mean = c(-3, 3), theta = 0)(grid)
#' cdf_2 <- pesscherMixture(mean = c(-3, 3), theta = -0.2)(grid)
#' cdf_3 <- pesscherMixture(mean = c(-3, 3), theta = 0.2)(grid)
#' ggplot()+
#'   geom_line(aes(grid, cdf_1), color = "black")+
#'   geom_line(aes(grid, cdf_2), color = "green")+
#'   geom_line(aes(grid, cdf_3), color = "red")
#'
#' @rdname desscherMixture
#' @aliases desscherMixture
#' @aliases pesscherMixture
#' @export
desscherMixture <- function(mean = c(0,0), sd = c(1,1), alpha = c(0.5, 0.5), theta = 0){

  # Moment generating function of k-component
  mgf <- function(mean, sd, theta) exp(theta*mean + (theta^2*sd^2)/2)

  num <- c()
  den <- 0
  for(k in 1:length(means)){
    num[k] <- alpha[k]*mgf(mean[k], sd[k], theta)
    den <- den + num[k]
  }
  # Update probabilities
  probs <- num/den
  # Update means parameters
  mean <- mean + theta*sd^2
  # Mixture pdf
  pdf <- function(x) dmixnorm(x, mean, sd, alaph)
  # Esscher pdf
  function(x, log = FALSE){
    probs <- pdf(x)
    # Log-probabilities
    if (log) {
      probs <- base::log(probs)
    }
    return(probs)
  }
}

#' @rdname desscherMixture
#'
#' @export
pesscherMixture <- function(mean = c(0,0), sd = c(1,1), alpha = c(0.5, 0.5), theta = 0){
  # Esscher pdf
  pdf <- desscherMixture(mean, sd, alpha, theta)
  # Esscher cdf depending on `x` and `h`
  cdf <- function(q, lower = -Inf){
    purrr::map_dbl(q, ~integrate(function(x) pdf(x), lower = lower, upper = .x)$value)
  }
}
