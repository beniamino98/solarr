#' Solar radiation distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' transformed global horizontal irradiance (GHI).
#'
#' @param x Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param Ct Numeric scalar or vector of clear-sky radiation values.
#' @param alpha Numeric scalar. Lower transformation parameter.
#' @param beta Numeric scalar. Scale transformation parameter. Typically
#'   `beta > 0` and `alpha + beta < 1`.
#' @param pdf_Y Function. Density function of the latent variable `Y`.
#' @param cdf_Y Function. Distribution function of the latent variable `Y`.
#' @param log Logical. If `TRUE`, `dsolarGHI()` returns log-densities.
#' @param log.p Logical. If `TRUE`, probabilities are supplied or returned on
#'   the log scale.
#' @param lower.tail Logical. If `TRUE`, probabilities are \eqn{P[X \le x]};
#'   otherwise, \eqn{P[X > x]}.
#' @param link Character string specifying the transformation link. Supported
#'   values are `"invgumbel"`, `"gumbel"`, `"logis"`, and `"norm"`.
#'
#' @return
#' - `dsolarGHI()` returns a numeric vector of density values.
#' - `psolarGHI()` returns a numeric vector of probabilities.
#' - `qsolarGHI()` returns a numeric vector of quantiles.
#' - `rsolarGHI()` returns a numeric vector of random draws.
#'
#' @details Consider a latent random variable \eqn{Y} with density `pdf_Y` and
#' distribution function `cdf_Y`. With the inverse Gumbel link, the transformed
#' solar radiation variable is
#' \deqn{R_t(y) = C(t) (1-\alpha-\beta \exp(-\exp(y)))}
#' with support \eqn{[C(t)(1-\alpha-\beta), C(t)(1-\alpha)]}.
#'
#' @examples
#' alpha <- 0.001
#' beta <- 0.9
#' Ct <- 7
#' dsolarGHI(c(3, 5), Ct, alpha, beta, dnorm)
#' psolarGHI(c(3, 5), Ct, alpha, beta, pnorm)
#' qsolarGHI(c(0.1, 0.9), Ct, alpha, beta, pnorm)
#'
#' set.seed(1)
#' rsolarGHI(3, Ct, alpha, beta, pnorm)
#' @name dsolarGHI
#' @rdname dsolarGHI
#' @aliases dsolarGHI
#' @aliases psolarGHI
#' @aliases qsolarGHI
#' @aliases rsolarGHI
#' 
#' @keywords distributions
#' @family distributions
#' @note Version 1.0.0.
#' @export
dsolarGHI <- function(x, Ct, alpha, beta, pdf_Y, log = FALSE, link = "invgumbel"){
  z_x <- (1 - x/Ct - alpha)/beta
  if (link == "invgumbel") {
    u_x <- suppressWarnings(log(-log(z_x)))
    den <- -Ct*beta*log(z_x^z_x)
  } else if (link == "gumbel") {
    u_x <- suppressWarnings(-log(-log(z_x)))
    den <- -Ct*beta*log(z_x^z_x)
  } else if (link == "logis") {
    u_x <- suppressWarnings(log(z_x/(1-z_x)))
    den <- Ct*beta*z_x*(1-z_x)
  } else if (link == "norm") {
    u_x <- suppressWarnings(qnorm(z_x))
    den <- Ct*beta*dnorm(qnorm(z_x))
  }
  # Numerator
  num <- pdf_Y(u_x)
  # Probability
  probs <- num/den
  # Log-probability
  if (log) {
    probs <- base::log(probs)
  }
  probs[is.nan(probs)] <- 0
  return(probs)
}

#' Distribution function for the GHI
#'
#' @rdname dsolarGHI
#' @export
psolarGHI  <- function(x, Ct, alpha, beta, cdf_Y, log.p = FALSE, lower.tail = TRUE, link = "invgumbel"){
  z_x <- (1 - x/Ct - alpha)/beta
  if (link == "invgumbel") {
    u_x <- suppressWarnings(log(-log(z_x)))
    probs <- cdf_Y(u_x)
  } else if (link == "gumbel") {
    u_x <- suppressWarnings(-log(-log(z_x)))
    probs <- 1-cdf_Y(u_x)
  } else if (link == "logis") {
    u_x <- suppressWarnings(log(z_x/(1-z_x)))
    probs <- 1-cdf_Y(u_x)
  } else if (link == "norm") {
    u_x <- suppressWarnings(qnorm(z_x))
    probs <- 1-cdf_Y(u_x)
  }
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
qsolarGHI  <- function(p, Ct, alpha, beta, cdf_Y, log.p = FALSE, lower.tail = TRUE, link = "invgumbel"){
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
  cdf <- function(x) psolarGHI(x, Ct, alpha, beta, cdf_Y, link = link)
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
rsolarGHI  <- function(n, Ct, alpha, beta, cdf_Y, link = "invgumbel"){
  # Simulated grades
  u <- runif(n, 0, 1)
  qsolarGHI(u, Ct, alpha, beta, cdf_Y, link = link)
}
