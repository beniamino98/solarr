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
#'
#' @examples
#' # Density
#' dsolarGHI(5, 7, 0.001, 0.9, function(x) dnorm(x))
#' dsolarGHI(6.993, 7, 0.001, 0.9, function(x) dnorm(x))
#' # Distribution
#' psolarGHI(6.993, 7, 0.001, 0.9, function(x) dnorm(x))
#' # Quantile
#' qsolarGHI(1, 7, 0.001, 0.9, function(x) dnorm(x))
#' qsolarGHI(c(0.05, 0.95), 7, 0.001, 0.9, function(x) dnorm(x))
#' # Random generator
#' rsolarGHI(10, 7, 0.001, 0.9, function(x) dnorm(x))
#' Ct <- Bologna$seasonal_data$Ct
#' GHI <- purrr::map(Ct, ~rsolarGHI(1, .x, 0.001, 0.9, function(x) dsnorm(x, shape = -4)))
#' plot(1:366, GHI, type="l")
#' @rdname dsolarGHI
#' @aliases dsolarGHI
#' @aliases psolarGHI
#' @aliases qsolarGHI
#' @aliases rsolarGHI
#' @export
dsolarGHI  <- function(x, Ct, alpha, beta, pdf_Yt){
  z_x <- (1 - x/Ct - alpha)/beta
  a1 <- 1/log(z_x)
  a2 <- 1/(1-x/Ct - alpha)
  probs <- -(1/Ct)*a1*a2*pdf_Yt(log(-log(z_x)))
  return(probs)
}

#' Distribution function for the GHI
#'
#' @rdname dsolarGHI
#' @export
psolarGHI  <- function(x, Ct, alpha, beta, pdf_Yt){
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
  return(probs)
}

#' Quantile function for the GHI
#'
#' @rdname dsolarGHI
#' @export
qsolarGHI  <- function(p, Ct, alpha, beta, pdf_Yt){
  # Bounds for GHI
  GHI_dw <- Ct*(1-alpha-beta)
  GHI_up <- Ct*(1-alpha)
  # Initial value for the routine
  GHI_init <- (GHI_up+GHI_dw)/2
  # Density function
  cdf <- function(x) psolarGHI(x, Ct, alpha, beta, pdf_Yt)
  # Empirical quantile function
  quantile_ <- numericQuantile(cdf, lower = GHI_dw, x0 = GHI_init)
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




