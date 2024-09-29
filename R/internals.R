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

