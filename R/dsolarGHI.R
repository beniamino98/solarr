#' Solar radiation random variable
#'
#' Solar radiation density, distribution, quantile and random generator.
#'
#' @param x vector of quantiles.
#' @param p vector of probabilities.
#' @param Ct clear sky radiation
#' @param alpha parameter `alpha > 0`.
#' @param beta parameter `beta > 0` and `alpha + beta < 1`.
#' @param pdf_Yt density of Yt.
#' @param log logical; if `TRUE`, probabilities are returned as `log(p)`.
#' @param log.p logical; if `TRUE`, probabilities p are given as `log(p)`.
#' @param lower.tail logical; if `TRUE`, the default, the computed probabilities are `P[X < x]`. Otherwise, `P[X > x]`.
#' @details Consider a random variable \eqn{Y \in [-\infty, \infty]} with a known density function `pdf_Yt`. Then
#' the funtion `dsolarGHI` compute the density function of the following transformed random variable, i.e.
#' \deqn{GHI(Y) = C_t (1-\alpha-\beta \exp(-\exp(Y)))}
#' where \eqn{GHI(Y) \in [Ct(1-\alpha-\beta), Ct(1-\alpha)]}.
#' @examples
#' # Density
#' dsolarGHI(5, 7, 0.001, 0.9, function(x) dnorm(x))
#' dsolarGHI(5, 7, 0.001, 0.9, function(x) dnorm(x, sd=2))
#'
#' # Distribution
#' psolarGHI(3.993, 7, 0.001, 0.9, function(x) dnorm(x))
#' psolarGHI(3.993, 7, 0.001, 0.9, function(x) dnorm(x, sd=2))
#'
#' # Quantile
#' qsolarGHI(c(0.05, 0.95), 7, 0.001, 0.9, function(x) dnorm(x))
#' qsolarGHI(c(0.05, 0.95), 7, 0.001, 0.9, function(x) dnorm(x, sd=2))
#'
#' # Random generator (I)
#' Ct <- Bologna$seasonal_data$Ct
#' GHI <- purrr::map(Ct, ~rsolarGHI(1, .x, 0.001, 0.9, function(x) dnorm(x, sd=0.8)))
#' plot(1:366, GHI, type="l")
#'
#' # Random generator (II)
#' pdf <- function(x) dmixnorm(x, c(-0.8, 0.5), c(1.2, 0.7), c(0.3, 0.7))
#' GHI <- purrr::map(Ct, ~rsolarGHI(1, .x, 0.001, 0.9, pdf))
#' plot(1:366, GHI, type="l")
#' @rdname dsolarGHI
#' @aliases dsolarGHI
#' @aliases psolarGHI
#' @aliases qsolarGHI
#' @aliases rsolarGHI
#' @export
dsolarGHI  <- function(x, Ct, alpha, beta, pdf_Yt, log = FALSE){
  z_x <- (1 - x/Ct - alpha)/beta
  u_x <- log(-log(z_x))
  probs <- -(pdf_Yt(u_x))/(Ct*beta*z_x*log(z_x))
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
psolarGHI  <- function(x, Ct, alpha, beta, pdf_Yt, log.p = FALSE, lower.tail = TRUE){
  probs <- c()
  for(i in 1:length(x)){
    pdf <- function(x) dsolarGHI(x, Ct, alpha, beta, pdf_Yt)
    if (x[i] > Ct*(1-alpha)) {
      probs[i] <- NA
    } else if (x[i] < Ct*(1-alpha-beta)) {
      probs[i] <- NA
    } else if (x[i] == Ct*(1-alpha)) {
      probs[i] <- 1
    } else if (x[i] == Ct*(1-alpha-beta)) {
      probs[i] <- 0
    } else {
      probs[i] <- integrate(pdf, lower = Ct*(1-alpha-beta), upper = x[i])$value
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

#' Quantile function for the GHI
#'
#' @rdname dsolarGHI
#' @export
qsolarGHI  <- function(p, Ct, alpha, beta, pdf_Yt, log.p = FALSE, lower.tail = TRUE){
  probs <- p
  # Log-probability
  if (log.p) {
    probs <- exp(probs)
  }
  # Lower tail
  if (!lower.tail){
    probs <- 1 - probs
  }

  # Bounds for GHI
  lower_bound <- Ct*(1-alpha-beta)
  upper_bound <- Ct*(1-alpha)
  # Initial value for the routine
  init_value <- (lower_bound+upper_bound)/2
  # Density function
  cdf <- function(x) psolarGHI(x, Ct, alpha, beta, pdf_Yt)
  # Empirical quantile function
  quantile_ <- Quantile(cdf, lower = lower_bound, x0 = init_value)
  # Quantiles
  x <- quantile_(p)
  return(x)
}

#' Random generator function for the GHI
#'
#' @rdname dsolarGHI
#' @export
rsolarGHI  <- function(n, Ct, alpha, beta, pdf_Yt){
  # Simulated grades
  u <- runif(n, 0, 1)
  qsolarGHI(u, Ct, alpha, beta, pdf_Yt)
}









