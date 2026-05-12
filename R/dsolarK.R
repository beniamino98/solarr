#' Clearness index distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the clearness index.
#'
#' @param x Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param alpha Numeric scalar. Lower transformation parameter.
#' @param beta Numeric scalar. Scale transformation parameter. Typically
#'   `beta > 0` and `alpha + beta < 1`.
#' @param pdf_Y Function. Density function of the latent variable `Y`.
#' @param cdf_Y Function. Distribution function of the latent variable `Y`.
#' @param log Logical. If `TRUE`, `dsolarK()` returns log-densities.
#' @param log.p Logical. If `TRUE`, probabilities are supplied or returned on
#'   the log scale.
#' @param lower.tail Logical. If `TRUE`, probabilities are \eqn{P[X \le x]};
#'   otherwise, \eqn{P[X > x]}.
#'
#' @return
#' - `dsolarK()` returns a numeric vector of density values.
#' - `psolarK()` returns a numeric vector of probabilities.
#' - `qsolarK()` returns a numeric vector of quantiles.
#' - `rsolarK()` returns a numeric vector of random draws.
#'
#' @details Consider a latent random variable \eqn{Y} with density `pdf_Y` and
#' distribution function `cdf_Y`. The clearness index is modeled as
#' \deqn{K(Y) = 1-\alpha-\beta \exp(-\exp(Y))}
#' with support \eqn{[1-\alpha-\beta, 1-\alpha]}.
#'
#' @examples
#' alpha <- 0.001
#' beta <- 0.9
#' dsolarK(c(0.2, 0.5), alpha, beta, dnorm)
#' psolarK(c(0.2, 0.5), alpha, beta, pnorm)
#' qsolarK(c(0.1, 0.9), alpha, beta, pnorm)
#'
#' set.seed(1)
#' rsolarK(3, alpha, beta, pnorm)
#' @rdname dsolarK
#' @aliases dsolarK
#' @aliases psolarK
#' @aliases qsolarK
#' @aliases rsolarK
#' @keywords distributions
#' @details Version 1.0.0.
#' @export
dsolarK  <- function(x,  alpha, beta, pdf_Y, log = FALSE){
  z_x <- (1 - x - alpha)/beta
  u_x <- log(-log(z_x))
  probs <- -(pdf_Y(u_x))/(beta*log(z_x^z_x))
  if (log) {
    probs <- base::log(probs)
  }
  return(probs)
}

#' Distribution function for the Clearness index
#'
#' @rdname dsolarK
#' @export
psolarK  <- function(x, alpha, beta, cdf_Y, log.p = FALSE, lower.tail = TRUE){
  z_x <- (1 - x - alpha)/beta
  u_x <- log(-log(z_x))
  probs <- cdf_Y(u_x)
  probs[x<=(1-alpha-beta)] <- 0
  probs[x>=(1-alpha)] <- 1
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

#' Quantile function for the Clearness index
#'
#' @rdname dsolarK
#' @export
qsolarK  <- function(p, alpha, beta, cdf_Y, log.p = FALSE, lower.tail = TRUE){
  # Log-probability
  if (log.p) {
    p <- exp(p)
  }
  # Lower tail
  if (!lower.tail) {
    p <- 1 - p
  }
  # Bounds for Clearness Index
  interval <- c(1-alpha-beta, 1-alpha)
  # Distribution function
  cdf <- function(x) psolarK(x, alpha, beta, cdf_Y)
  # Empirical quantile function
  quantile_numeric <- Quantile(cdf, interval = interval)
  # Quantiles
  x <- quantile_numeric(p)
  return(x)
}

#' Random generator function for the Clearness index
#'
#' @rdname dsolarK
#' @export
rsolarK  <- function(n, alpha, beta, cdf_Y){
  # Simulated grades
  u <- runif(n, 0, 1)
  qsolarK(u, alpha, beta, cdf_Y)
}
