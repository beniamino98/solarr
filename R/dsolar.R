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
#' @rdname dsolar
#' @aliases dsolar
#' @aliases psolarGHI
#' @aliases qsolarGHI
#' @aliases rsolarGHI
#' @export
dsolar  <- function(x, Ct, alpha, beta, pdf_Yt, log = FALSE){
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


#' Clearness index random variable
#'
#' Clearness index density, distribution, quantile and random generator.
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
#' the funtion `dsolarK` compute the density function of the following transformed random variable, i.e.
#' \deqn{K(Y) = 1-\alpha-\beta \exp(-\exp(Y))}
#' where \eqn{K(Y) \in [1-\alpha-\beta, 1-\alpha]}.
#' @examples
#' # Density
#' dsolarK(0.4, 0.001, 0.9, function(x) dnorm(x))
#' dsolarK(0.4, 0.001, 0.9, function(x) dnorm(x, sd = 2))
#'
#' # Distribution
#' psolarK(0.493, 0.001, 0.9, function(x) dnorm(x))
#' psolarK(0.493, 0.001, 0.9, function(x) dnorm(x, sd = 2))
#'
#' # Quantile
#' qsolarK(c(0.05, 0.95), 0.001, 0.9, function(x) dnorm(x))
#' qsolarK(c(0.05, 0.95), 0.001, 0.9, function(x) dnorm(x, sd = 2))
#'
#' # Random generator (I)
#' Kt <- rsolarK(366, 0.001, 0.9, function(x) dnorm(x, sd = 1.3))
#' plot(1:366, Kt, type="l")
#'
#' # Random generator (II)
#' pdf <- function(x) dmixnorm(x, c(-1.8, 0.9), c(0.5, 0.7), c(0.6, 0.4))
#' Kt <- rsolarK(36, 0.001, 0.9, pdf)
#' plot(1:36, Kt, type="l")
#' @rdname dsolarK
#' @aliases dsolarK
#' @aliases psolarK
#' @aliases qsolarK
#' @aliases rsolarK
#' @export
dsolarK  <- function(x,  alpha, beta, pdf_Yt, log = FALSE){
  z_x <- (1 - x - alpha)/beta
  u_x <- log(-log(z_x))
  probs <- -(pdf_Yt(u_x))/(beta*z_x*log(z_x))
  if (log) {
    probs <- base::log(probs)
  }
  probs[is.nan(probs)] <- 0
  return(probs)
}

#' Distribution function for the Clearness index
#'
#' @rdname dsolarK
#' @export
psolarK  <- function(x, alpha, beta, pdf_Yt, log.p = FALSE, lower.tail = TRUE){
  probs <- c()
  for(i in 1:length(x)){
    pdf <- function(x) dsolarK(x, alpha, beta, pdf_Yt)
    if (x[i] > (1-alpha)) {
      probs[i] <- NA
    } else if (x[i] < (1-alpha-beta)) {
      probs[i] <- NA
    } else if (x[i] == (1-alpha)) {
      probs[i] <- 1
    } else if (x[i] == (1-alpha-beta)) {
      probs[i] <- 0
    } else {
      probs[i] <- integrate(pdf, lower = (1-alpha-beta), upper = x[i])$value
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

#' Quantile function for the Clearness index
#'
#' @rdname dsolarK
#' @export
qsolarK  <- function(p, alpha, beta, pdf_Yt, log.p = FALSE, lower.tail = TRUE){
  probs <- p
  # Log-probability
  if (log.p) {
    probs <- exp(probs)
  }
  # Lower tail
  if (!lower.tail){
    probs <- 1 - probs
  }

  # Bounds for Clearness Index
  lower_bound <- (1-alpha-beta)
  upper_bound <- (1-alpha)
  # Initial value for the routine
  init_value <- (upper_bound+lower_bound)/2
  # Density function
  cdf <- function(x) psolarK(x, alpha, beta, pdf_Yt)
  # Empirical quantile function
  quantile_ <- Quantile(cdf, lower = lower_bound, x0 = init_value)
  # Quantiles
  x <- quantile_(p)
  return(x)
}

#' Random generator function for the Clearness index
#'
#' @rdname dsolarK
#' @export
rsolarK  <- function(n, alpha, beta, pdf_Yt){
  # Simulated grades
  u <- runif(n, 0, 1)
  qsolarK(u, alpha, beta, pdf_Yt)
}



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


#' Bivariate PDF GHI
#' @rdname dmvsolarGHI
#' @export
dmvsolarGHI <- function(x, Ct, alpha, beta, joint_pdf_Yt){
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }
  z <- x
  z[,1] <- (1 - x[,1] - alpha[1])/beta[1]
  z[,2] <- (1 - x[,2]/Ct - alpha[2])/beta[2]
  u <- log(-log(z))
  # Compute denominator
  z_prod <- apply(z, 1, prod)
  # Denominator
  den <- Ct*prod(beta)*apply(z, 1, prod)*apply(log(z), 1, prod)
  # Mixture probabilities
  probs <- (1/den)*joint_pdf_Yt(u)
  probs[is.nan(probs)] <- 0
  return(probs)
}
