#' Solar risk driver distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the solar risk driver.
#'
#' @param x Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param alpha Numeric scalar. Lower transformation parameter.
#' @param beta Numeric scalar. Scale transformation parameter. Typically
#'   `beta > 0` and `alpha + beta < 1`.
#' @param pdf_Y Function. Density function of the latent variable `Y`.
#' @param cdf_Y Function. Distribution function of the latent variable `Y`.
#' @param log Logical. If `TRUE`, `dsolarX()` returns log-densities.
#' @param log.p Logical. If `TRUE`, probabilities are supplied or returned on
#'   the log scale.
#' @param lower.tail Logical. If `TRUE`, probabilities are \eqn{P[X \le x]};
#'   otherwise, \eqn{P[X > x]}.
#'
#' @return
#' - `dsolarX()` returns a numeric vector of density values.
#' - `psolarX()` returns a numeric vector of probabilities.
#' - `qsolarX()` returns a numeric vector of quantiles.
#' - `rsolarX()` returns a numeric vector of random draws.
#'
#' @details Consider a latent random variable \eqn{Y} with density `pdf_Y` and
#' distribution function `cdf_Y`. The solar risk driver is modeled as
#' \deqn{X(Y) = \alpha+\beta \exp(-\exp(Y))}
#' with support \eqn{[\alpha, \alpha+\beta]}.
#'
#' @examples
#' alpha <- 0.001
#' beta <- 0.9
#' dsolarX(c(0.2, 0.5), alpha, beta, dnorm)
#' psolarX(c(0.2, 0.5), alpha, beta, pnorm)
#' qsolarX(c(0.1, 0.9), alpha, beta, pnorm)
#'
#' set.seed(1)
#' rsolarX(3, alpha, beta, pnorm)
#' 
#' @name dsolarX
#' @rdname dsolarX
#' @aliases dsolarX
#' @aliases psolarX
#' @aliases qsolarX
#' @aliases rsolarX
#' 
#' @keywords distributions
#' @family distributions
#' @note Version 1.0.0.
#' @export
dsolarX  <- function(x, alpha, beta, pdf_Y, log = FALSE){
  z_x <- (x - alpha)/beta
  u_x <- log(-log(z_x))
  probs <- -(pdf_Y(u_x))/(beta*log(z_x^z_x))
  if (log) {
    probs <- base::log(probs)
  }
  return(probs)
}

#' Distribution function for the Solar risk driver
#'
#' @rdname dsolarX
#' @export
psolarX  <- function(x, alpha, beta, cdf_Y, log.p = FALSE, lower.tail = TRUE){

  z_x <- (x - alpha)/beta
  u_x <- log(-log(z_x))
  probs <- 1 - cdf_Y(u_x)
  probs[x<=alpha] <- 0
  probs[x>=(alpha+beta)] <- 1
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

#' Quantile function for the Solar risk driver
#'
#' @rdname dsolarX
#' @export
qsolarX  <- function(p, alpha, beta, cdf_Y, log.p = FALSE, lower.tail = TRUE){
  # Log-probability
  if (log.p) {
    p <- exp(p)
  }
  # Lower tail
  if (!lower.tail) {
    p <- 1 - p
  }
  # Bounds for solar risk driver
  interval <- c(alpha, alpha + beta)
  # Density function
  cdf <- function(x) psolarX(x, alpha, beta, cdf_Y)
  # Empirical quantile function
  quantile_numeric <- Quantile(cdf, interval = interval)
  # Quantiles
  x <- quantile_numeric(p)
  return(x)
}

#' Random generator function for the Solar risk driver
#'
#' @rdname dsolarX
#' @export
rsolarX  <- function(n, alpha, beta, cdf_Y){
  # Simulated grades
  u <- runif(n, 0, 1)
  # Simulated random variable
  qsolarX(u, alpha, beta, cdf_Y)
}
