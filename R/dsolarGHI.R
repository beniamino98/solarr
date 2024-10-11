#' Solar radiation random variable
#'
#' Solar radiation density, distribution, quantile and random generator.
#'
#' @param x vector of quantiles.
#' @param p vector of probabilities.
#' @param Ct clear sky radiation
#' @param alpha parameter `alpha > 0`.
#' @param beta parameter `beta > 0` and `alpha + beta < 1`.
#' @param pdf_Y density of Y.
#' @param cdf_Y distribution of Y.
#' @param log logical; if `TRUE`, probabilities are returned as `log(p)`.
#' @param log.p logical; if `TRUE`, probabilities p are given as `log(p)`.
#' @param lower.tail logical; if `TRUE`, the default, the computed probabilities are `P[X < x]`. Otherwise, `P[X > x]`.
#' @details Consider a random variable \eqn{Y \in [-\infty, \infty]} with a known density function `pdf_Y`. Then
#' the funtion `dsolarGHI` compute the density function of the following transformed random variable, i.e.
#' \deqn{GHI(Y) = C_t (1-\alpha-\beta \exp(-\exp(Y)))}
#' where \eqn{GHI(Y) \in [C_t(1-\alpha-\beta), C_t(1-\alpha)]}.
#' @examples
#' # Parameters
#' alpha = 0
#' beta = 0.9
#' Ct <- 7
#' # Grid of points
#' grid <- seq(Ct*(1-alpha-beta), Ct*(1-alpha), by = 0.01)
#'
#' # Density
#' dsolarGHI(5, Ct, alpha, beta, function(x) dnorm(x))
#' dsolarGHI(5, Ct, alpha, beta, function(x) dnorm(x, sd=2))
#' plot(grid, dsolarGHI(grid, Ct, alpha, beta, function(x) dnorm(x, mean = -1, sd = 0.9)), type="l")
#'
#' # Distribution
#' psolarGHI(3.993, 7, 0.001, 0.9, function(x) pnorm(x))
#' psolarGHI(3.993, 7, 0.001, 0.9, function(x) pnorm(x, sd=2))
#' plot(grid, psolarGHI(grid, Ct, alpha, beta, function(x) pnorm(x, sd = 0.2)), type="l")
#'
#' # Quantile
#' qsolarGHI(c(0.05, 0.95), 7, 0.001, 0.9, function(x) pnorm(x))
#' qsolarGHI(c(0.05, 0.95), 7, 0.001, 0.9, function(x) pnorm(x, sd=2))
#'
#' # Random generator (I)
#' Ct <- Bologna$seasonal_data$Ct
#' GHI <- purrr::map(Ct, ~rsolarGHI(1, .x, alpha, beta, function(x) pnorm(x, sd=1.4)))
#' plot(1:366, GHI, type="l")
#'
#' # Random generator (II)
#' cdf <- function(x) pmixnorm(x, c(-0.8, 0.5), c(1.2, 0.7), c(0.3, 0.7))
#' GHI <- purrr::map(Ct, ~rsolarGHI(1, .x, 0.001, 0.9, cdf))
#' plot(1:366, GHI, type="l")
#' @rdname dsolarGHI
#' @aliases dsolarGHI
#' @aliases psolarGHI
#' @aliases qsolarGHI
#' @aliases rsolarGHI
#' @export
dsolarGHI  <- function(x, Ct, alpha, beta, pdf_Y, log = FALSE){
  z_x <- (1 - x/Ct - alpha)/beta
  u_x <- suppressWarnings(log(-log(z_x)))
  num <- pdf_Y(u_x)
  den <- Ct*beta*log(z_x^z_x)
  probs <- -num/den
  # Log-probability
  if (log) {
    probs <- base::log(probs)
  }
  return(probs)
}

#' Distribution function for the GHI
#'
#' @rdname dsolarGHI
#' @export
psolarGHI  <- function(x, Ct, alpha, beta, cdf_Y, log.p = FALSE, lower.tail = TRUE){

  z_x <- (1 - x/Ct - alpha)/beta
  u_x <- suppressWarnings(log(-log(z_x)))
  probs <- cdf_Y(u_x)
  probs[x<=Ct*(1-alpha-beta)] <- 0
  probs[x>=Ct*(1-alpha)] <- 1
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


#' Quantile function for the GHI
#'
#' @rdname dsolarGHI
#' @export
qsolarGHI  <- function(p, Ct, alpha, beta, cdf_Y, log.p = FALSE, lower.tail = TRUE){
  # Log-probability
  if (log.p) {
    p <- exp(p)
  }
  # Lower tail
  if (!lower.tail) {
    p <- 1 - p
  }
  # Bounds for GHI
  interval = c(Ct*(1-alpha-beta), Ct*(1-alpha))
  # Density function
  cdf <- function(x) psolarGHI(x, Ct, alpha, beta, cdf_Y)
  # Empirical quantile function
  quantile_numeric <- Quantile(cdf, interval = interval)
  # Quantiles
  x <- quantile_numeric(p)
  return(x)
}

#' Random generator function for the GHI
#'
#' @rdname dsolarGHI
#' @export
rsolarGHI  <- function(n, Ct, alpha, beta, cdf_Y){
  # Simulated grades
  u <- runif(n, 0, 1)
  qsolarGHI(u, Ct, alpha, beta, cdf_Y)
}









