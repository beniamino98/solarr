# Compute different types of model errors
modelErrors <- function(x, x_hat, na.rm = FALSE){
  errors <- list()
  errors[["MAX"]] <- max(abs(x-x_hat), na.rm = na.rm)
  errors[["PMAX"]] <- max(abs(x-x_hat)/x, na.rm = na.rm)
  errors[["PMIN"]] <- min(abs(x-x_hat)/x, na.rm = na.rm)
  errors[["MSE"]] <- sd(x-x_hat, na.rm = na.rm)
  errors[["MAE"]] <- mean(abs(x-x_hat), na.rm = na.rm)
  errors[["MAPE"]] <- mean(abs(x-x_hat)/x, na.rm = na.rm)*100
  errors[["MPE"]] <- mean((x-x_hat)/x, na.rm = na.rm)*100
  return(dplyr::bind_rows(errors))
}

AR_GARCH_filter <- function(x, params, x0 = NA, sigma0 = NA, sigma2 = NA, t0 = 0){

  if (missing(x)) {
    stop("Please, provide a vector for `x`!")
  }

  # Length of the time series
  n <- length(x)
  names(params) <- c("phi1", "phi2", "omega0", "omega1", "omega2", "c0", "c1", "c2")

  # AR parameters
  phi <- params[c("phi1", "phi2")]

  # Garch parameters
  omega <- params[c("omega0", "omega1", "omega2")]
  # Long term GARCH variance
  sigma0_2 <- omega[1]/(1 - omega[2] - omega[3])
  if (sigma0_2 < 0 | any(omega > 1) | any(omega < 0)) {
    warning("The GARCH parameters are not valid!")
    return(NA)
  }
  # Unconditional variance
  if (!is.na(sigma2)) {
    omega[1] <- sigma2*(1 - omega[2] - omega[3])
    sigma0_2 <- omega[1]/(1 - omega[2] - omega[3])
  }

  # Seasonal variance
  c_ <- params[c("c0", "c1", "c2")]
  # Positive constraint
  if (c_[1] < sqrt(c_[2]^2 + c_[3]^2)){
    warning("The Seasonal parameters are not valid!")
    return(NA)
  }

  # Initialize std. deviation
  if (is.na(sigma0)) {
    sigma0 <- sqrt(sigma0_2)
  } else {
    sigma0 <- 1
  }
  sigma <- rep(sigma0, n)

  # Default initial residuals
  if (!is.na(x0)) {
    x[1] <- x0
  }
  mu <- x
  loglik <- rep(0, n)
  # AR-GARCH loop
  for(t in 3:n){
    mu[t] <- phi[1]*x[t-1] + phi[2]*x[t-2]
    omega_t <- (2*base::pi*(t0 + t))/365
    sigma_bar <-  sqrt(c_[1] + c_[2]*cos(omega_t) + c_[3]*sin(omega_t))
    sigma[t] <- sqrt(omega[1] + omega[2]*(x[t-1]/sigma_bar)^2 + omega[3]*sigma[t-1]^2)
    loglik[t] <- dnorm(x[t], mean = mu[t], sd = sigma[t]*sigma_bar, log = TRUE)
  }
  message("Log-likelihood: ", sum(loglik), "\r", appendLF = FALSE)

  params <- c(phi, omega, c_)
  # Output
  structure(
    list(
      params = params,
      fitted = dplyr::tibble(
        t = 1:n,
        sigma = sigma,
        sigma_bar = sigma_bar,
        mu_t = mu,
        loglik = loglik,
        u = (x-mu_t)/(sigma*sigma_bar)
      )
    )
  )
}

AR_GARCH_NM_filter <- function(x, params, sigma0 = 1, t0 = 0){

  check_positivity <- function(x){
    if (x[1] < sqrt(x[2]^2 + x[3]^2)){
      warning("The Seasonal parameters are not valid!")
      return(FALSE)
    }
    return(TRUE)
  }

  if (missing(x)) {
    stop("Please, provide a vector for `x`!")
  }
  # Length of the time series
  n <- length(x)

  # Standard names for the parameters
  names(params) <- c("phi1", "phi2", "omega0", "omega1", "omega2", "c0", "c1", "c2",
                     "a_mu_up0", "a_mu_up1", "a_mu_up2", "a_mu_dw0", "a_mu_dw1", "a_mu_dw2",
                     "c_sd_up0", "c_sd_up1", "c_sd_up2", "c_sd_dw0", "c_sd_dw1", "c_sd_dw2")

  # AR parameters
  phi <- params[c("phi1", "phi2")]
  # Garch parameters
  omega <- params[c("omega0", "omega1", "omega2")]
  # Check Long term GARCH variance
  sigma0_2 <- omega[1]/(1 - omega[2] - omega[3])
  if (sigma0_2 < 0 | any(omega > 1) | any(omega < 0)) {
    warning("The GARCH parameters are not valid!")
    return(NA)
  }
  # Seasonal variance
  c_ <- params[c("c0", "c1", "c2")]
  # Positive constraint
  if (!check_positivity(c_)) {
    warning("The Seasonal parameters are not valid!")
    return(NA)
  }
  # Seasonal mean (up)
  a_up <- params[c("a_mu_up0", "a_mu_up1", "a_mu_up2")]
  # Seasonal variance (up)
  c_up <- params[c("c_sd_up0", "c_sd_up1", "c_sd_up2")]
  # Positive constraint
  if (!check_positivity(c_up)) {
    warning("The Seasonal parameters c_up are not valid!")
    return(NA)
  }
  # Seasonal mean (down)
  a_dw <- params[c("a_mu_dw0", "a_mu_dw1", "a_mu_dw2")]
  # Seasonal variance (down)
  c_dw <- params[c("c_sd_dw0", "c_sd_dw1", "c_sd_dw2")]
  # Positive constraint
  if (!check_positivity(c_up)) {
    warning("The Seasonal parameters c_dw are not valid!")
    return(NA)
  }

  t = 1:n
  sigma = rep(sigma[1], n)
  mu_x = rep(NA, n)
  mu_x_up = rep(NA, n)
  mu_x_dw = rep(NA, n)
  sd_x_up = rep(NA, n)
  sd_x_dw = rep(NA, n)
  sigma_bar = rep(NA, n)
  sigma_up = rep(NA, n)
  sigma_dw = rep(NA, n)
  mu_up = rep(NA, n)
  mu_dw = rep(NA, n)
  loglik_up = rep(NA, n)
  loglik_dw = rep(NA, n)
  loglik = rep(NA, n)
  B = rep(NA, n)

  # AR-GARCH loop
  for(t in 3:n){
    # Conditional mean
    mu_x[t] <- phi[1]*x[t-1] + phi[2]*x[t-2]
    # Seasonal variance
    omega_t <- (2*base::pi*(t0 + t))/365
    sigma_bar[t] <-  sqrt(c_[1] + c_[2]*cos(omega_t) + c_[3]*sin(omega_t))
    # GARCH variance
    sigma[t] <- sqrt(omega[1] + omega[2]*(x[t-1]/sigma_bar[t])^2 + omega[3]*sigma[t-1]^2)
    # Mixture components
    # Up
    mu_up[t] <- a_up[1] + a_up[2]*cos(omega_t) + a_up[3]*sin(omega_t)
    sigma_up[t] <-  sqrt(c_up[1] + c_up[2]*cos(omega_t) + c_up[3]*sin(omega_t))
    # Conditional moments under up distribution
    mu_x_up[t] = mu_x[t] + sigma[t]*sigma_bar[t]*mu_up[t]
    sd_x_up[t] = sigma[t]*sigma_bar[t]*sigma_up[t]
    loglik_up[t] <- dnorm(x[t], mean = mu_x_up[t], sd = sd_x_up[t], log = TRUE)
    # Down
    mu_dw[t] <- a_dw[1] + a_dw[2]*cos(omega_t) + a_dw[3]*sin(omega_t)
    sigma_dw[t] <-  sqrt(c_dw[1] + c_dw[2]*cos(omega_t) + c_dw[3]*sin(omega_t))
    # Conditional moments under dw distribution
    mu_x_dw[t] = mu_x[t] + sigma[t]*sigma_bar[t]*mu_dw[t]
    sd_x_dw[t] = sigma[t]*sigma_bar[t]*sigma_dw[t]
    loglik_dw[t] <- dnorm(x[t], mean = mu_x_dw[t], sd = sd_x_dw[t], log = TRUE)

    if(loglik_up[t] > loglik_dw[t]){
      B[t] <- 1
      loglik[t] <- loglik_up[t]
    } else {
      B[t] <- 0
      loglik[t] <- loglik_dw[t]
    }
  }
  message("Log-likelihood: ", sum(loglik, na.rm = TRUE), "\r", appendLF = FALSE)
  params <- c(phi, omega, c_, a_up, a_dw, c_up, c_dw)

  # Output
  structure(
    list(
      params = params,
      fitted = dplyr::tibble(
        t = t,
        sigma = sigma,
        mu_x = mu_x,
        mu_x_up = mu_x_up,
        mu_x_dw = mu_x_dw,
        sd_x_up = sd_x_up,
        sd_x_dw = sd_x_dw,
        sigma_bar = sigma_bar,
        sigma_up = sigma_up,
        sigma_dw = sigma_dw,
        mu_up = mu_up,
        mu_dw = mu_dw,
        loglik_up = loglik_up,
        loglik_dw = loglik_dw,
        loglik = loglik,
        B = B
      ) %>% na.omit()
    )
  )
}

AR_GARCH_logLik <- function(params, x, ...){
  garch_filter <- AR_GARCH_filter(x = x, params = params, ...)
  if (length(garch_filter) == 1 && is.na(garch_filter)){
    return(NA)
  } else {
    sum(garch_filter$fitted$loglik)
  }
}

AR_GARCH_optimize <- function(x, params, sigma2 = NA, ..., control = list(maxit = 1e5, fnscale = -1, ndeps = 1e-5)){

  names(params) <- c("phi1", "phi2", "omega0", "omega1", "omega2", "c0", "c1", "c2")
  opt <- optim(par = params, AR_GARCH_logLik, x = x, sigma2 = sigma2, ..., control = control)

  garch_filter <- AR_GARCH_filter(x, params = opt$par, sigma2 = sigma2, ...)
  structure(
    list(
      init_params = params,
      params = garch_filter$params,
      logLik = sum(garch_filter$fitted$loglik, na.rm = TRUE),
      fitted = garch_filter$fitted
    )
  )
}

#' Separate a named list of parameter into ARCH and GARCH components
#' @param params named vector of GARCH parameters. Tthe ARCH parameters `alpha1` for first lag, `alpha2` for second lag, etc.
#' Finally, the GARCH parameters `beta1` for first lag, `beta2` for second lag and so on. The intercept is named `omega`.
#' @examples
#' init_params <- c(omega=1, beta2=0.2, beta3=0.1, beta1=0.3, alpha1 = 0.01, alpha2 = 0.02)
#' params <- GARCH_parameters(init_params)
#' @keywords internal
#' @noRd
GARCH_parameters <- function(params){
  if (is.list(params)){
    names(params) <- NULL
    params <- unlist(params)
  }
  params_names <- names(params)
  # Extract intercept
  omega <- params[stringr::str_detect(params_names, "omega")]
  # Extract ARCH component
  idx_arch <- which(stringr::str_detect(params_names, "alpha"))
  alpha <- 0
  if (!purrr::is_empty(idx_arch)){
    alpha <- params[idx_arch]
  }
  # Extract GARCH component
  idx_garch <- which(stringr::str_detect(params_names, "beta"))
  beta <- 0
  if (!purrr::is_empty(idx_garch)){
    beta <- params[idx_garch]
    beta <- beta[order(names(beta))]
  }
  list(
    omega = omega,
    alpha = alpha,
    beta = beta
  )
}

#' Filter for a GARCH(1,1)
#'
#' @param x vector, time series to filter.
#' @param params vector of GARCH parameters, in order `omega0`, `omega1`, `omega2`.
#' @param x0 starting value for time series.
#' @param sigma0 starting value for GARCH std. deviation.
#' @param sigma2_bar unconditional variance. If specified will modify `omega0` to match it.
#'
#' @examples
#' set.seed(1)
#' GARCH_filter(x, params = c(0.2, 0.3, 0.5), x0 = NA, sigma0 = 1, sigma2_bar = 1.3)
#' GARCH_filter(x, params = c(0.2, 0.3, 0.5), x0 = 1, sigma0 = 0.8)
#' @rdname GARCH_filter
#' @noRd
GARCH_filter <- function(x, params, x0 = NA, sigma0 = NA, sigma2_bar = NA){

  # Length of the time series
  n <- length(x)
  # Set omega to match unconditional variance
  params[1] <- ifelse(is.na(sigma2_bar), params[1], sigma2_bar*(1-sum(params[2:3])))
  # Unconditional std. deviation
  sigma_bar <- sqrt(params[1]/(1-sum(params[2:3])))
  # Vector of stochastic std. deviations
  sigma0 <- ifelse(is.na(sigma0), sigma_bar, sigma0)
  sigma <- c(sigma0, rep(sigma_bar, n-1))
  # Initialization time series
  u <- x/sigma
  # GARCH(1,1) filter
  for(t in 2:n){
    sigma[t] <- sqrt(params[1] + params[2]*x[t-1]^2 + params[3]*sigma[t-1]^2)
    u[t] <- x[t]/sigma[t]
  }
  # Output
  dplyr::tibble(
    t = 1:n,
    sigma = sigma,
    x = x,
    u = x/sigma
  )
}

#' Simulate a GARCH(1,1) process
#'
#' @param n number of simulations.
#' @inheritParams GARCH_filter
#'
#' @examples
#' GARCH_simulate(10, params = c(0.2, 0.3, 0.5), x0 = NA, sigma0 = 1)
#' GARCH_simulate(10, params = c(0.2, 0.3, 0.5), x0 = 1, sigma0 = 0.8)
#' @rdname GARCH_simulate
#' @noRd
GARCH_simulate <- function(n, params, x0 = NA, sigma0 = NA, seed = 1){
  set.seed(seed)
  # Long-term variance
  sigma_bar <- sqrt(params[1]/(1-sum(params[2:3])))
  # Vector of stochastic std. deviations
  sigma0 <- ifelse(is.na(sigma0), sigma_bar, sigma0)
  sigma <- c(sigma0, rep(sigma_bar, n-1))
  # Simulated residuals
  u <- rnorm(n)
  # Initialization time series
  x <- u*sigma
  x[1] <- ifelse(is.na(x0), x[1], x0)
  # GARCH(1,1) Filter
  for(t in 2:n){
    sigma[t] <- sqrt(params[1] + params[2]*x[t-1]^2 + params[3]*sigma[t-1]^2)
    x[t] <- u[t]*sigma[t]
  }
  # Output
  dplyr::tibble(
    seed = seed,
    t = 1:n,
    sigma = sigma,
    u = u,
    x = u*sigma
  )
}

#' Log-likelihood function GARCH(1,1)
#'
#' @inheritParams GARCH_filter
#' @param weights vector of weights for likelihood.
#' @examples
#' set.seed(1)
#' x <- GARCH_simulate(1000, params = c(0.2, 0.3, 0.5), x0 = NA, sigma0 = 1)
#' GARCH_logLik(x$x, params = c(0.2, 0.32, 0.50), x0 = NA, sigma0 = 1)
#' @rdname GARCH_logLik
#' @keywords internal
#' @noRd
GARCH_logLik <- function(params, x, weights, ...){
  if (any(params < 0)){
    return(NA)
  }
  # Length of the time series
  n <- length(x)
  fitted_garch <- GARCH_filter(x, params, ...)
  loglik <- purrr::map2_dbl(fitted_garch$x, fitted_garch$sigma, ~dnorm(.x, sd = .y, log = TRUE))
  #loglik <- purrr::map_dbl(fitted_garch$u, ~dnorm(.x, log = TRUE))
  # Flexible weights for likelihood
  if (missing(weights)){
    w <- rep(1, n)
  } else {
    w <- ifelse(weights > 0, 1, 0)
  }
  loglik <- sum(loglik*w, na.rm = TRUE)
  return(loglik)
}



#' Fit function GARCH(1,1)
#'
#' @inheritParams GARCH_filter
#'
#' @examples
#' # True parameters
#' coef <- c(0.5, 0.1, 0.3)
#'
#' # Manual optimization
#' x <- GARCH_simulate(10000, params = coef, x0 = NA)
#' # Initial parameters
#' init_coef <- coef * 0.3
#' opt_1 <- GARCH_optimize(init_coef, x$x)
#'
#' # Rugarch optimization
#' spec <- rugarch::ugarchspec(variance.model = list(garchOrder = c(1,1)),
#' mean.model = list(armaOrder = c(0, 0), include.mean = FALSE))
#' opt_2 <- pars <- rugarch::ugarchfit(spec, x$x)
#'
#' # Difference
#' abs(opt_1$params - opt_2@fit$coef)
#'
#' # Log-likelihood
#' GARCH_logLik(opt_1$params, x$x)
#' GARCH_logLik(opt_2@fit$coef, x$x)
#'
#' @param control list of control for `optim`.
#' @rdname GARCH_optimize
#' @keywords internal
#' @noRd
GARCH_optimize <- function(params, x, sigma2 = NA, weights, ..., control = list(maxit = 1e5, fnscale = -1, ndeps = 1e-5)){

  opt <- optim(par = params, GARCH_logLik, x = x, sigma2_bar = sigma2, weights = weights, ..., control = control)
  # Match unconditional variance
  if (!is.na(sigma2)){
    opt$par[1] <- sigma2*(1-opt$par[2]-opt$par[3])
  }

  structure(
    list(
      init_params = params,
      params = opt$par,
      logLik = opt$value,
      fitted = GARCH_filter(x, params = opt$par, ...)
    )
  )
}

