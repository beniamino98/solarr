#' Esscher transform of a Gaussian Mixture
#'
#' @inheritParams dmixnorm
#' @param theta Esscher parameter, the default is zero.
#'
#' @examples
#' library(ggplot2)
#' grid <- seq(-5, 5, 0.01)
#' # Density
#' pdf_1 <- desscherMixture(means = c(-3, 3), theta = 0)(grid)
#' pdf_2 <- desscherMixture(means = c(-3, 3), theta = -0.5)(grid)
#' pdf_3 <- desscherMixture(means = c(-3, 3), theta = 0.5)(grid)
#' ggplot()+
#'  geom_line(aes(grid, pdf_1), color = "black")+
#'  geom_line(aes(grid, pdf_2), color = "green")+
#'  geom_line(aes(grid, pdf_3), color = "red")
#' # Distribution
#' cdf_1 <- pesscherMixture(means = c(-3, 3), theta = 0)(grid)
#' cdf_2 <- pesscherMixture(means = c(-3, 3), theta = -0.2)(grid)
#' cdf_3 <- pesscherMixture(means = c(-3, 3), theta = 0.2)(grid)
#' ggplot()+
#'   geom_line(aes(grid, cdf_1), color = "black")+
#'   geom_line(aes(grid, cdf_2), color = "green")+
#'   geom_line(aes(grid, cdf_3), color = "red")
#'
#' @rdname desscherMixture
#' @aliases desscherMixture
#' @aliases pesscherMixture
#' @export
desscherMixture <- function(means = c(0,0), sd = c(1,1), p = c(0.5, 0.5), theta = 0){

  # Moment generating function of k-component
  mgf <- function(mean, sd, theta){
    exp(theta*mean + (theta^2*sd^2)/2)
  }

  num <- c()
  den <- 0
  for(k in 1:length(means)){
    num[k] <- p[k]*mgf(means[k], sd[k], theta)
    den <- den + num[k]
  }
  # Update probabilities
  p <- num/den
  # Update means parameters
  means <- means + theta*sd^2
  # Mixture pdf
  pdf <- function(x) dmixnorm(x, means, sd, p)
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
pesscherMixture <- function(means = c(0,0), sd = c(1,1), p = c(0.5, 0.5), theta = 0){
  # Esscher pdf
  pdf <- desscherMixture(means, sd, p, theta)
  # Esscher cdf depending on `x` and `h`
  cdf <- function(x, lower = -Inf){
    purrr::map_dbl(x, ~integrate(function(x) pdf(x), lower = lower, upper = .x)$value)
  }
}
