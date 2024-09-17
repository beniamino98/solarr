#' Gaussian mixture
#'
#' Fit the parameters of a gaussian mixture with k-components.
#'
#' @param x vector
#' @param means vector of initial means parameters.
#' @param sd vector of initial std. deviation parameters.
#' @param p vector of initial probability parameters.
#' @param components number of components.
#' @param match_moments logical. When `TRUE`, the parameters of the second distribution will be estimated such that
#' the empirical first two moments of `x` matches the theoretical Gaussian mixture moments.
#' @param prior_p prior probability for the k-state. If the k-component is not `NA` the probability will be considered as given and
#' the parameter `p[k]` will be equal to `prior_p[k]`.
#' @param weights observations weights, if a weight is equal to zero the observation is excluded, otherwise is included with unitary weight.
#' When `missing` all the available observations will be used.
#' @param maxit maximum number of iterations.
#' @param absotol absolute level for convergence.
#' @param na.rm logical. When `TRUE`, the default, `NA` values will be excluded from the computations.
#'
#' @examples
#' means = c(-3,0,3)
#' sd = rep(1, 3)
#' p = c(0.2, 0.3, 0.5)
#' # Density function
#' pdf <- dmixnorm(means, sd, p)
#' # Distribution function
#' cdf <- pmixnorm(means, sd, p)
#' # Random numbers
#' x <- rgaussianMixture(1000, means, sd, p)
#' gaussianMixture(x$X, means, sd, p, components = 3)
#' gaussianMixture(x$X, means, sd, prior_p = p, components = 3)
#' @return list with clustered components and the optimal parameters.
#'
#' @rdname gaussianMixture
#' @name gaussianMixture
#' @export
gaussianMixture <- function(x, means, sd, p, components = 2, prior_p = rep(NA, components), weights, maxit = 100, abstol = 10e-15,  na.rm = FALSE){

  # Initialization
  n <- length(x)
  # Weights
  if (missing(weights)){
    w <- rep(1, n)
  } else {
    w <- ifelse(weights == 0, 0, 1)
  }

  # Empirical moments
  e_x_hat <- mean(x[w != 0], na.rm = na.rm)
  v_x_hat <- var(x[w != 0], na.rm = na.rm)
  sd_x_hat <- sqrt(v_x_hat)

  # Update n parameter
  n_w <- length(x[w != 0])

  # Default means
  if (missing(means) || any(is.na(means))){
    # means <- rep(e_x_hat, components)
    means <- quantile(x, probs = seq(0.2, 0.8, length.out = components))
  }
  # Default std. deviations
  if (missing(sd) || any(is.na(sd))){
    sd <- rep(sd_x_hat, components)
  }
  # Default probabilities
  if (missing(p) || any(is.na(p))) {
    p <- rep(1/components, components)
  }
  # Fixed prior probabilities
  p <- ifelse(is.na(prior_p), p, prior_p)
  p <- p/sum(p)

  # Routine
  # 0. Initialization
  log_likelihood <- 0
  previous_log_likelihood <- -Inf
  responsibilities <- matrix(0, nrow = n, ncol = components)
  previous_params <- list(mean = means, sd = sd, p = p)
  # EM Algorithm
  for (iteration in 1:maxit) {
    # E-step: posterior probabilities
    for (i in 1:n) {
      for(k in 1:components){
        responsibilities[i, k] <- previous_params$p[k]*dnorm(x[i], previous_params$mean[k], previous_params$sd[k])
      }
      # Normalize the posterior probabilities
      responsibilities[i,] <- (responsibilities[i,])/sum(responsibilities[i,])
    }
    # Optimal parameters
    params <- previous_params
    # M-step: Update the parameters
    for(k in 1:components){
      # Normalizing factor for each group
      n_k <- sum(w*responsibilities[, k], na.rm = na.rm)
      # Mean parameter k-component
      params$mean[k] <- sum(responsibilities[, k]*w*x, na.rm = na.rm)/n_k
      # Std. deviation k-component
      params$sd[k] <- sqrt(sum(responsibilities[, k]*w*(x - params$mean[k])^2, na.rm = na.rm)/n_k)
      # Probability k-component
      params$p[k] <- ifelse(is.na(prior_p[k]), n_k/n_w, prior_p[k])
    }
    # Calculate the log-likelihood
    log_likelihood <- 0
    for(i in 1:n) {
      ll <- 0
      for(k in 1:components){
        ll <- ll + params$p[k]*dnorm(x[i], params$mean[k], params$sd[k])
      }
      log_likelihood <- sum(c(log_likelihood, log(ll)*w[i]), na.rm = TRUE)
    }
    # Check for convergence
    stop_condition <- abs(log_likelihood - previous_log_likelihood) < abstol
    if (stop_condition) {
      break
    } else {
      # Update log-likelihood
      previous_log_likelihood <- log_likelihood
      # Update parameters
      previous_params <- params
    }
  }

  # Final classification of each component
  x_hat <- matrix(0, nrow = n, ncol = components)
  B_hat <- matrix(0, nrow = n, ncol = components)
  for(i in 1:n) {
    ll <- c()
    for(k in 1:components){
      ll[k] <- sum(responsibilities[i,k])
    }
    B_hat[i, which.max(ll)] <- 1
    x_hat[i,] <- B_hat[i,]*x[i]
  }
  colnames(B_hat) <- paste0("B", 1:components)
  colnames(x_hat) <- paste0("x", 1:components)
  B_hat <- dplyr::as_tibble(B_hat)
  x_hat <- dplyr::as_tibble(x_hat)
  z_hat <- dplyr::as_tibble((x_hat-params$mean)/params$sd*B_hat)
  colnames(z_hat) <- paste0("z", 1:components)
  df_fitted <- dplyr::bind_cols(B_hat, x_hat, z_hat)

  # ML-parameters
  params <- previous_params
  # Assign a name to the parameters
  names(params$mean) <- paste0("mu", 1:components)
  names(params$sd) <- paste0("sd", 1:components)
  names(params$p) <- paste0("p", 1:components)
  # Parameters on fitted sub-samples
  params_hat <- params
  for(k in 1:components){
    component_name <- paste0("x", k)
    component <- df_fitted[w == 1,]
    component <- component[component[, component_name] != 0,]
    params_hat$mean[k] <- mean(component[[component_name]], na.rm = na.rm)
    params_hat$sd[k] <- sd(component[[component_name]], na.rm = na.rm)
    params_hat$p[k] <- mean(df_fitted[[paste0("B", k)]], na.rm = na.rm)
  }

  # Reorder the components by decreasing means
  idx_order <- order(params$mean)
  params$mean <- params$mean[idx_order]
  params$sd <- params$sd[idx_order]
  params$p <- params$p[idx_order]
  params_hat$mean <- params_hat$mean[idx_order]
  params_hat$sd <- params_hat$sd[idx_order]
  params_hat$p <- params_hat$p[idx_order]

  colnames(responsibilities) <- paste0("B", 1:components)
  responsibilities <- dplyr::as_tibble(responsibilities)

  # Log-likelihood on fitted parameters
  # Calculate the log-likelihood
  log_likelihood <- 0
  for(i in 1:n) {
    ll <- 0
    for(k in 1:components){
      ll <- ll + params_hat$p[k]*dnorm(x[i], params_hat$mean[k], params_hat$sd[k])
    }
    log_likelihood <- sum(c(log_likelihood, log(ll)*w[i]), na.rm = TRUE)
  }

  structure(
    list(
      iteration = iteration,
      par = params,
      par_hat = params_hat,
      responsibilities = responsibilities,
      log_lik = previous_log_likelihood,
      log_lik_hat = log_likelihood,
      fitted = df_fitted
    ),
    class = c("gaussianMixture")
  )
}

#' Fit a monthly Gaussian Mixture Pdf (??NOT USED)
#'
#' Fit the monthly parameters for the density function of a Gaussian mixture with two components.
#'
#' @param x vector
#' @param date vector of dates
#' @param means matrix of initial means with dimension `12 X components`.
#' @param sd matrix of initial std. deviations with dimension `12 X components`.
#' @param p matrix of initial p with dimension `12 X components`. The rows must sum up to 1.
#' @param prior_p matrix of prior probabilities for the each month. Any element that is different from `NA` will be not optimized and will be considered
#' as given.
#' @param ... other parameters for the optimization function. See \code{\link{gaussianMixture}} for more details.
#'
#' @rdname gaussianMixture_monthly
#' @name gaussianMixture_monthly
#' @export
gaussianMixture_monthly <- function(x, date, means, sd, p, components = 2, prior_p, ...){

  data <- dplyr::tibble(date = date, Month = lubridate::month(date), eps = x)

  # Check matrix of means
  if (missing(means)) {
    means <- matrix(NA, nrow = 12, ncol = components)
  } else {
    if (ncol(mean) != components || nrow(means) != 12) {
      msg <- paste0("The matrix `means` has not dimension: ", "12 X ", components)
      stop(msg)
    }
  }
  # Check matrix of std. deviations
  if (missing(sd)) {
    sd <- matrix(NA, nrow = 12, ncol = components)
  } else {
    if (ncol(sd) != components || nrow(sd) != 12) {
      msg <- paste0("The matrix `sd` has not dimension: ", "12 X ", components)
      stop(msg)
    }
  }
  # Check matrix of probabilities
  if (missing(p)) {
    p <- matrix(NA, nrow = 12, ncol = components)
  } else {
    if (ncol(p) != components || nrow(p) != 12) {
      msg <- paste0("The matrix `p` has not dimension: ", "12 X ", components)
      stop(msg)
    }
    if (any(rowSums(p)!=1)) {
      msg <- paste0("The rows of the matrix `p` do not sum up to 1! ")
      stop(msg)
    }
  }
  # Check matrix of prior probabilities
  if (missing(prior_p)){
    prior_p <- matrix(NA, nrow = 12, ncol = components)
  } else {
    if (ncol(prior_p) != components || nrow(prior_p) != 12) {
      msg <- paste0("The matrix `prior_p` has not dimension: ", "12 X ", components)
      stop(msg)
    }
    if (any(rowSums(prior_p)!=1)) {
      msg <- paste0("The rows of the matrix `prior_p` do not sum up to 1! ")
      stop(msg)
    }
  }

  # Initialization
  GM_model <- list()
  data_months <- list()
  params <- list()
  # Monthly Gaussian Mixture
  for(m in unique(data$Month)){
    # Monthly data
    data_months[[m]] <- dplyr::filter(data, Month == m)
    # Monthly data
    eps <- data_months[[m]]$eps
    # Fitted parameters
    gm <- gaussianMixture(eps, means = means[m,], sd = sd[m,], p = p[m,], prior_p = prior_p[m,])
    # Dataset with fitted parameters
    params[[m]] <- list(means = gm$par$mean,
                        sd = gm$par$sd,
                        p = gm$par$p)
    # Theoretical first moment
    e_u <- sum(gm$par$mean*gm$par$p)
    # Theoretical variance
    v_u <- sum((gm$par$mean^2 + gm$par$sd^2)*gm$par$p) - e_u^2
    # Sample mean
    e_u_hat <- mean(eps, na.rm = TRUE)
    # Sample variance
    v_u_hat <- var(eps, na.rm = TRUE)
    # Store classified series
    data_months[[m]] <- dplyr::bind_cols(data_months[[m]], gm$fitted)
    # Monthly data
    GM_model[[m]] <- dplyr::tibble(Month = m,
                                   loss = gm$log_lik,
                                   nobs = length(eps),
                                   e_x = e_u,
                                   v_x = v_u,
                                   e_x_hat = e_u_hat,
                                   v_x_hat = v_u_hat)
  }

  # Reorder the variables
  GM_model <- dplyr::bind_rows(GM_model)
  names(params) <- lubridate::month(1:12, label = TRUE)
  # Dataset with all the parameters
  l_params <- list(
    means = purrr::map_df(params, ~dplyr::bind_rows(.x$means)),
    sd = purrr::map_df(params, ~dplyr::bind_rows(.x$sd)),
    p = purrr::map_df(params, ~dplyr::bind_rows(.x$p))
  )
  l_params <- purrr::map(l_params, ~dplyr::bind_cols(Month = 1:12, .x))

  structure(
    list(
      fitted = dplyr::bind_rows(data_months),
      params = params,
      l_params = l_params,
      model = GM_model
    )
  )
}








