#' Standard GARCH expected value formula
#'
#' @examples
#' # Forecast horizon
#' h <- 2
#' # GARCH parameters
#' alpha <- 0.08
#' beta <- 0.35
#' omega <- 1*(1 - alpha - beta)
#' # Moments
#' e_x2 = 1
#' e_x4 = 3
#' # Initial values for variance
#' sigma2_t <- 1.2
#' sigma4_t <- sigma2_t^2
#'
#' e_sigma2_h(10, omega, alpha, beta, e_x2, sigma2_t)
#' e_sigma2_h_mix(10, omega, alpha, beta, e_x2[1], sigma2_t)
#' e_sigma4_h_mix(10, omega, alpha, beta, e_x2, e_x4, sigma4_t)
#' v_sigma2_h_mix(10, omega, alpha, beta, e_x2, e_x4, sigma4_t)
#' v_sigma_h_mix(10, omega, alpha, beta, e_x2, e_x4, sigma4_t)
#' e_sigma12_h_mix(10, omega, alpha, beta, e_x2, e_x4, sigma4_t)
#' e_sigma32_h_mix(10, omega, alpha, beta, e_x2, e_x4, sigma4_t)
#' @export
e_sigma2_h <- function(h, omega, alpha, beta, e_x2 = 1, sigma2_t){
  # Derived quantities
  lambda <- alpha * e_x2[1] + beta
  # Long term expectation
  sigma2_inf <- omega / (1 - lambda)
  # Forecasted second moment
  sigma2_h <- sigma2_inf + (alpha + beta)^(1:h) * (sigma2_t - sigma2_inf)
  sigma2_h <- c(sigma2_t, sigma2_h)
  names(sigma2_h) <- paste0("t+", 0:h)
  sigma2_h
}
#' Iterative GARCH expected value formula
#' @export
e_sigma2_h_mix <- function(h, omega, alpha, beta, e_x2 = 1, sigma2_t){
  sigma2_h <- c(sigma2_t)
  # Second moment
  m2 <- e_x2
  if (length(e_x2) == 1 & h > 1){
    m2 <- rep(e_x2, h)
  }
  # Derived quantities
  lambda <- alpha * m2 + beta
  lambda_prod <- lambda
  lambda_prod[1] <- 1
  lambda_prod <- cumprod(lambda_prod)
  sigma2_h <- c(sigma2_h, omega * cumsum(lambda_prod) + sigma2_t * cumprod(lambda))
  names(sigma2_h) <- paste0("t+", 0:h)
  sigma2_h
}
#' Iterative GARCH second moment formula
#' @export
e_sigma4_h_mix <- function(h, omega, alpha, beta, e_x2 = 1, e_x4 = 3, sigma4_t){
  sigma4_h <- c()
  for(i in 0:h){
    if(i == 0){
      sigma4_h <- c(sigma4_t[1])
      next
    }
    # Second moment
    m2 <- e_x2
    if (length(e_x2 == 1) & i > 1){
      m2 <- rep(e_x2, i)
    }
    # Fourth moment
    m4 <- e_x4
    if (length(e_x4 == 1) & i > 1){
      m4 <- rep(e_x4, i)
    }
    gamma_t <- alpha^2 * m4 + 2 * alpha * beta + beta^2
    #gamma_t[i] <- 1
    sigma2_t <- e_sigma2_h_mix(i, omega, alpha, beta, m2, sqrt(sigma4_t))
    lambda_t <- alpha * m2 + beta
    b_t <- omega * (omega + 2 * lambda_t * sigma2_t[1:i])
    sigma4_h <- c(sigma4_h, b_t[i] + sum(b_t[i:1] * cumprod(gamma_t)) + sigma4_t * prod(gamma_t))
  }
  names(sigma4_h) <- paste0("t+", 0:h)
  sigma4_h
}
#' Iterative GARCH variance formula
#' @export
v_sigma2_h_mix <- function(h, omega, alpha, beta, e_x2 = 1, e_x4 = 3, sigma4_t){
  # Second moment GARCH variance
  e_sigma4 <- e_sigma4_h_mix(h, omega, alpha, beta, e_x2, e_x4, sigma4_t)
  # Expectation GARCH variance
  e_sigma2 <- e_sigma2_h_mix(h, omega, alpha, beta, e_x2, sqrt(sigma4_t))
  # Variance
  e_sigma4 - e_sigma2^2
}
#' Iterative GARCH variance formula (approximated)
#' @export
v_sigma_h_mix <- function(h, omega, alpha, beta, e_x2 = 1, e_x4 = 3, sigma4_t){
  # Expectation GARCH variance
  e_sigma2 <- e_sigma2_h_mix(h, omega, alpha, beta, e_x2, sqrt(sigma4_t))
  # Expectation GARCH std. dev
  e_sigma <- e_sigma12_h_mix(h, omega, alpha, beta, e_x2, e_x4, sigma4_t)
  # Variance
  e_sigma2 - e_sigma^2
}
#' Conditional first moment GARCH std. dev (approximated)
#' @export
e_sigma12_h_mix <- function(h, omega, alpha, beta, e_x2 = 1, e_x4 = 3, sigma4_t){
  # Expectation GARCH variance
  e_sigma2 <- e_sigma2_h_mix(h, omega, alpha, beta, e_x2, sqrt(sigma4_t))
  # Variance GARCH variance
  v_sigma2 <- v_sigma2_h_mix(h, omega, alpha, beta, e_x2, e_x4, sigma4_t)
  # Moment to power 1/2 (Approximated)
  e_sigma2^(1/2) - (1/8) * v_sigma2 / sqrt(e_sigma2)^3
}
#' Conditional third moment GARCH std. dev (approximated)
#' @export
e_sigma32_h_mix <- function(h, omega, alpha, beta, e_x2 = 1, e_x4 = 3, sigma4_t){
  # Expectation variance
  e_sigma2 <- e_sigma2_h_mix(h, omega, alpha, beta, e_x2, sqrt(sigma4_t))
  # Variance variance
  v_sigma2 <- v_sigma2_h_mix(h, omega, alpha, beta, e_x2, e_x4, sigma4_t)
  # Moment to power 3/2 (Approximated)
  e_sigma2^(3/2) + (3/8) * v_sigma2 / sqrt(e_sigma2)
}


