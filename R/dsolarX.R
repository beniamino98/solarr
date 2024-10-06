#' Solar risk driver random variable
#'
#' Solar risk driver density, distribution, quantile and random generator.
#'
#' @param x vector of quantiles.
#' @param p vector of probabilities.
#' @param alpha parameter `alpha > 0`.
#' @param beta parameter `beta > 0` and `alpha + beta < 1`.
#' @param pdf_Yt density of Yt.
#' @param log logical; if `TRUE`, probabilities are returned as `log(p)`.
#' @param log.p logical; if `TRUE`, probabilities p are given as `log(p)`.
#' @param lower.tail logical; if `TRUE`, the default, the computed probabilities are `P[X < x]`. Otherwise, `P[X > x]`.
#' @details Consider a random variable \eqn{Y \in [-\infty, \infty]} with a known density function `pdf_Yt`. Then
#' the funtion `dsolarX` compute the density function of the following transformed random variable, i.e.
#' \deqn{X(Y) = \alpha+\beta \exp(-\exp(Y))}
#' where \eqn{X(Y) \in [\alpha, \alpha+\beta]}.
#' @examples
#' # Density
#' dsolarX(0.4, 0.001, 0.9, function(x) dnorm(x))
#' dsolarX(0.4, 0.001, 0.9, function(x) dnorm(x, sd = 2))
#'
#' # Distribution
#' psolarX(0.493, 0.001, 0.9, function(x) dnorm(x))
#' dsolarX(0.493, 0.001, 0.9, function(x) dnorm(x, sd = 2))
#'
#' # Quantile
#' qsolarX(c(0.05, 0.95), 0.001, 0.9, function(x) dnorm(x))
#' qsolarX(c(0.05, 0.95), 0.001, 0.9, function(x) dnorm(x, sd = 1.3))
#'
#' # Random generator (I)
#' Kt <- rsolarX(366, 0.001, 0.9, function(x) dnorm(x, sd = 0.8))
#' plot(1:366, Kt, type="l")
#'
#' # Random generator (II)
#' pdf <- function(x) dmixnorm(x, c(-1.8, 0.9), c(0.5, 0.7), c(0.6, 0.4))
#' Kt <- rsolarX(366, 0.001, 0.9, pdf)
#' plot(1:366, Kt, type="l")
#' @rdname dsolarX
#' @aliases dsolarX
#' @aliases psolarX
#' @aliases qsolarX
#' @aliases rsolarX
#' @export
dsolarX  <- function(x, alpha, beta, pdf_Yt, log = FALSE){
  z_x <- (x - alpha)/beta
  u_x <- log(-log(z_x))
  probs <- -(pdf_Yt(u_x))/(beta*z_x*log(z_x))
  if (log) {
    probs <- base::log(probs)
  }
  probs[is.nan(probs)] <- 0
  return(probs)
}

#' Distribution function for the Solar risk driver
#'
#' @rdname dsolarX
#' @export
psolarX  <- function(x, alpha, beta, pdf_Yt, log.p = FALSE, lower.tail = TRUE){
  probs <- c()
  for(i in 1:length(x)){
    pdf <- function(x) dsolarX(x, alpha, beta, pdf_Yt)
    if (x[i] > alpha+beta) {
      probs[i] <- NA
    } else if (x[i] < alpha) {
      probs[i] <- NA
    } else if (x[i] == alpha) {
      probs[i] <- 0
    } else if (x[i] == alpha+beta) {
      probs[i] <- 1
    } else {
      probs[i] <- integrate(pdf, lower = alpha, upper = x[i])$value
    }
  }
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
qsolarX  <- function(p, alpha, beta, pdf_Yt, log.p = FALSE, lower.tail = TRUE){
  probs <- p
  # Log-probability
  if (log.p) {
    probs <- exp(probs)
  }
  # Lower tail
  if (!lower.tail){
    probs <- 1 - probs
  }

  # Bounds for solar risk driver
  lower_bound <- alpha
  upper_bound <- alpha+beta
  # Initial value for the routine
  init_value <- (upper_bound+lower_bound)/2
  # Density function
  cdf <- function(x) psolarX(x, alpha, beta, pdf_Yt)
  # Empirical quantile function
  quantile_ <- Quantile(cdf, lower = lower_bound, x0 = init_value)
  # Quantiles
  x <- quantile_(p)
  return(x)
}

#' Random generator function for the Solar risk driver
#'
#' @rdname dsolarX
#' @export
rsolarX  <- function(n, alpha, beta, pdf_Yt){
  # Simulated grades
  u <- runif(n, 0, 1)
  qsolarX(u, alpha, beta, pdf_Yt)
}
