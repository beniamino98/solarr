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
gaussianMixture <- function(x, means, sd, p, components = 2, weights, maxit = 100, abstol = 10e-15,  na.rm = FALSE){

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
  # Rescale probabilities
  p <- p/sum(p)

  # Routine
  # 0. Initialization
  log_likelihood <- 0
  previous_log_likelihood <- -Inf
  prev_responsibilities <- matrix(0, nrow = n, ncol = components)
  previous_params <- list(mean = means, sd = sd, p = p)
  # EM Algorithm
  for (iteration in 1:maxit) {
    # E-step: posterior probabilities
    responsibilities <- prev_responsibilities
    for (i in 1:n) {
      for(k in 1:components){
        responsibilities[i, k] <- previous_params$p[k]*dnorm(x[i], previous_params$mean[k], previous_params$sd[k])
      }
      # Normalize the posterior probabilities
      responsibilities[i,][is.na(responsibilities[i,])] <- 0
      responsibilities[i,][is.nan(responsibilities[i,])] <- 0
      responsibilities[i,] <- (responsibilities[i,])/sum(responsibilities[i,])
      responsibilities[i,][is.nan(responsibilities[i,])] <- 0
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
      params$p[k] <- n_k/n_w
    }

    if (any(params$p > 0.9)){
      responsibilities <- prev_responsibilities
      break
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


#' Multivariate gaussian mixture
#'
#' @rdname mvgaussianMixture
#' @name mvgaussianMixture
#' @export
mvgaussianMixture <- function(x, means, sd, p, components = 2, maxit = 100, abstol = 10e-15, na.rm = FALSE){

  # Ensure that there are not NAs or NaN observations
  idx_NA <- is.na(x)
  if (any(idx_NA)) {
    x <- na.omit(x)
    wrn <- paste0("Removed ", sum(idx_NA), " NA observations!")
    warning(wrn)
  }
  # Number of observations
  n_w <- nrow(x)
  # Number of variables
  j_w <- ncol(x)
  # Empirical moments
  e_x_hat <- colMeans(x, na.rm = na.rm)
  v_x_hat <- apply(x, 2, var)

  # Default starting means
  if (missing(means) || any(is.na(means))){
    # Initialize a matrix for the components
    means <- matrix(0, nrow = components, ncol = j_w, dimnames = list(1:components, dimnames(x)[[2]]))
    probs <- seq(0.8, 0.2, length.out = components)
    for(k in 1:components){
      means[k,] <- apply(x, 2, quantile, probs = probs[k])
    }
  }

  # Default std. deviations
  if (missing(sd) || any(is.na(sd))){
    sd <- list()
    for(k in 1:components){
      sd[[k]] <- diag(v_x_hat)
    }
  }
  # Default probabilities
  if (missing(p) || any(is.na(p))) {
    p <- rep(1/components, components)
  }

  # Routine
  # 0. Initialization
  log_likelihood <- 0
  previous_log_likelihood <- -Inf
  prev_responsibilities <- matrix(0, nrow = n_w, ncol = components)
  previous_params <- list(mean = means, sd = sd, p = p)
  iteration <- 1
  # EM Algorithm
  for (iteration in 1:maxit) {
    # E-step: posterior probabilities
    responsibilities <- prev_responsibilities
    for (i in 1:n_w) {
      for(k in 1:components){
        responsibilities[i, k] <- previous_params$p[k]*mvtnorm::dmvnorm(x[i,], mean = previous_params$mean[k,], sigma = previous_params$sd[[k]])
      }
      # Normalize the posterior probabilities
      responsibilities[i,] <- (responsibilities[i,])/sum(responsibilities[i,], na.rm = TRUE)
      responsibilities[i,][is.na(responsibilities[i,])] <- 0
    }

    # Optimal parameters
    k <- 1
    params <- previous_params
    # M-step: Update the parameters
    for(k in 1:components){
      # Normalizing factor for each component
      n_k <- sum(responsibilities[, k], na.rm = na.rm)
      # Mean parameters k-component
      params$mean[k,] <- apply(responsibilities[, k]*x, 2, sum)/n_k
      # Covariance matrix k-component
      params$sd[[k]] <- diag(apply(x^2*responsibilities[, k], 2, sum)/n_k - params$mean[k,]^2)
      params$sd[[k]][1,2] <- params$sd[[k]][2,1] <- sum(x[,1]*x[,2]*responsibilities[, k])/n_k - params$mean[k,][1]*params$mean[k,][2]
      # Probability k-component
      params$p[k] <- n_k/n_w
    }

    if(any(params$p > 0.9)){
      warning("Probability greater than 0.9 Break!")
      params <- previous_params
      break
    }

    # Calculate the log-likelihood
    log_likelihood <- 0
    for(i in 1:n_w) {
      ll <- 0
      for(k in 1:components){
        ll <- ll + params$p[k]*mvtnorm::dmvnorm(x[i,], mean = params$mean[,k], sigma = params$sd[[k]])
      }
      log_likelihood <- sum(c(log_likelihood, log(ll)), na.rm = TRUE)
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
    print(log_likelihood)
    if (iteration == maxit) {
      message("Max iteration reached (", iteration, ")")
    }
  }

  # Final classification of each component
  B_hat <- matrix(0, nrow = n_w, ncol = components)
  for(i in 1:n_w) {
    ll <- c()
    for(k in 1:components){
      ll[k] <- sum(responsibilities[i,k])
    }
    B_hat[i, which.max(ll)] <- 1
  }
  colnames(B_hat) <- paste0("B", 1:components)
  B_hat <- dplyr::as_tibble(B_hat)

  # ML-parameters
  params <- previous_params
  # Reorder the components by decreasing means
  colnames(responsibilities) <- paste0("B", 1:components)
  responsibilities <- dplyr::as_tibble(responsibilities)

  # Log-likelihood on fitted parameters
  # Calculate the log-likelihood
  # Calculate the log-likelihood
  log_likelihood <- 0
  for(i in 1:n_w) {
    ll <- 0
    for(k in 1:components){
      ll <- ll + params$p[k]*mvtnorm::dmvnorm(x[i,], mean = params$mean[,k], sigma = params$sd[[k]])
    }
    log_likelihood <- sum(c(log_likelihood, log(ll)), na.rm = TRUE)
  }

  upd_params <- list()
  upd_params$means <- params$mean
  upd_params$sigma2 <- params$mean
  upd_params$rho <- c(0,0)
  upd_params$p <- params$p
  i <- 1
  for(i in 1:length(params$sd)){
    upd_params$sigma2[i,] <- diag(params$sd[[i]])
    upd_params$rho[i] <-  params$sd[[i]][upper.tri(params$sd[[i]])]/prod(sqrt(diag(params$sd[[i]])))
  }

  structure(
    list(
      B_hat = B_hat,
      iteration = iteration,
      params = upd_params,
      responsibilities = responsibilities,
      log_lik = log_likelihood
    ),
    class = c("mvgaussianMixture")
  )
}



