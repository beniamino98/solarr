#' Separate a named list of parameter into ARCH and GARCH components
#' @param params named vector of GARCH parameters. Tthe ARCH parameters `alpha1` for first lag, `alpha2` for second lag, etc.
#' Finally, the GARCH parameters `beta1` for first lag, `beta2` for second lag and so on. The intercept is named `omega`.
#' @examples
#' init_params <- c(omega=1, beta2=0.2, beta3=0.1, beta1=0.3, alpha1 = 0.01, alpha2 = 0.02)
#' params <- GARCH_parameters(init_params)
#' @keywords internal
#' @noRd
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' params <- c(0.6, 0.1, 0.2)
#' x <- GARCH_simulate(5000, params = params, x0 = NA, sigma0 = 1)
#' GARCH_optimize(params, x$x, sigma2 = 1)
#' @param control list of control for `optim`.
#' @rdname GARCH_optimize
#' @keywords internal
#' @noRd
#' @export
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


#' GARCH models with flexible weights
#'
#' @keywords internal
#' @export
ugarchfit_fp <- function(spec, data, weights, sigma0, ...){

  if (missing(weights)){
    weights <- 1
  } else {
    weights <- ifelse(weights == 0, 0, 1)
  }
  # Initial parameters
  if (is.null(spec@model$start.pars)) {
    # Safe GARCH model
    safe_GARCH <- purrr::safely(rugarch::ugarchfit)
    # Fitted model
    model <- safe_GARCH(data = data, spec = spec, out.sample = 0)$result
    if (model@fit$convergence == 0){
      # Store the parameters inside the specification of the GARCH
      spec@model$fixed.pars <- model@fit$coef
      spec@model$start.pars <- model@fit$coef
    } else {
      fixed.pars <- spec@model$pars[,4][spec@model$pars[,4] == 1]
      fixed.pars <- runif(length(fixed.pars), 0, 0.1)
      fixed.pars <- fixed.pars/sum(fixed.pars)
      fixed.pars[1] <- 1 - sum(fixed.pars[-1])
      names(fixed.pars) <- names(spec@model$pars[,4][spec@model$pars[,4] == 1])
      spec@model$start.pars <- fixed.pars
      spec@model$fixed.pars <- fixed.pars
    }
  }
  # Extract original attributes
  params_attr <- attributes(spec@model$start.pars)
  # Loss function
  log.likelihood <- function(params, spec, sigma0){
    attributes(params) <- params_attr
    if (any(params[stringr::str_detect(names(params), "omega|alpha|beta")] < 0)) return(NA)
    # Unconditional variance
    if (!missing(sigma0)) {
      params[1] <- sigma0*(1-sum(params[stringr::str_detect(names(params), "alpha|beta")]))
    }
    spec@model$fixed.pars <- params
    log.lik <- rugarch::ugarchfilter(spec, data, ...)@filter$log.likelihoods
    sum(-log.lik*weights)
  }
  # Optimization function
  opt <- optim(spec@model$start.pars, log.likelihood, spec = spec, sigma0 = sigma0)
  # Unconditional variance
  if (!missing(sigma0)) {
    opt$par[1] <- sigma0*(1-sum(opt$par[stringr::str_detect(names(opt$par), "alpha|beta")]))
  }
  # Update the parameters
  spec@model$fixed.pars <- opt$par
  return(spec)
}

#' Next step function for a GARCH(p,q)
#'
#' @param omega The intercept
#' @param alpha ARCH parameters
#' @param beta GARCH parameters
#'
#' @keywords internal
#' @export
GARCH_pq_next_step <- function(omega = 1, alpha, beta){

  # ARCH order
  p <- ifelse(missing(alpha), 0, length(alpha))
  # GARCH order
  q <- ifelse(missing(beta), 0, length(beta))

  function(x = 1, sigma = 1){
    # Initialize GARCH variance
    sigma2_next <- omega
    # ARCH(p) component
    if (p > 0) {
      sigma2_next <- sigma2_next + sum(alpha*x^2)
    }
    # GARCH(q) component
    if (q > 0) {
      sigma2_next <- sigma2_next + sum(beta*sigma^2)
    }
    return(sqrt(sigma2_next))
  }
}

AR_next_step <- function(phi_0 = 0, phi){
  # AR order
  p <- ifelse(missing(phi), 0, length(phi))
  function(x = 1){
    # Initialize next step value
    x_next_step <- phi_0
    # AR(p) component
    if (p > 0) {
      x_next_step <- x_next_step + sum(phi*x)
    }
    return(x_next_step)
  }
}




