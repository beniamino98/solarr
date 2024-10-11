#' Gaussian mixture
#'
#' Fit the parameters of a gaussian mixture with k-components.
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
#' x <- rmixnorm(5000, means, sd, p)
#' gm <- gaussianMixture$new(components=3)
#' gm$fit(x$X)
#' gm$parameters
#' gm$EM(x$X)
#' gm$parameters
#' gm$fitted
#' @rdname gaussianMixture
#' @name gaussianMixture
#' @export
gaussianMixture <- R6::R6Class("gaussianMixture",
                               public = list(
                                 #' @field maxit Integer, maximum number of iterations.
                                 maxit = 100,
                                 #' @field abstol Numeric, absolute level for convergence.
                                 abstol = 10e-5,
                                 #' @field components Integer, number of components.
                                 components = 2,
                                 #' @field means Numeric vector of means parameters.
                                 means = NA,
                                 #' @field sd Numeric vector of std. deviation parameters.
                                 sd = NA,
                                 #' @field p Numeric vector of probability parameters.
                                 p = NA,
                                 #' @description
                                 #' Initialize a gaussianMixture object
                                 #' @param components Integer, number of components.
                                 #' @param maxit Numeric, maximum number of iterations.
                                 #' @param abstol Numeric, absolute level for convergence.
                                 initialize = function(components = 2, maxit = 500, abstol = 10e-10){
                                   # Control parameters
                                   self$maxit <- maxit
                                   self$abstol <- abstol
                                   self$components <- components
                                 },
                                 #' @description
                                 #' Compute the log-likelihood
                                 #' @param x vector
                                 #' @param params Named list of mixture parameters.
                                 logLik = function(x, params){
                                   if (missing(params)) {
                                     params <- self$parameters
                                   }
                                   # Calculate the log-likelihood
                                   log_likelihood <- 0
                                   for(i in 1:length(x)) {
                                     ll <- 0
                                     for(k in 1:self$components){
                                       ll <- ll + sum(params$p[k]*dnorm(x[i], params$means[k], params$sd[k]))
                                     }
                                     log_likelihood <- sum(c(log_likelihood, log(ll)), na.rm = TRUE)
                                   }
                                   return(log_likelihood)
                                 },
                                 #' @description
                                 #' Compute the posterior probabilities (E-step)
                                 #' @param x vector
                                 #' @param params a list of mixture parameters
                                 E_step = function(x, params){
                                   if (missing(params)) {
                                     params <- self$parameters
                                   }
                                   responsabilities <- matrix(0, nrow = length(x), ncol = self$components)
                                   for(k in 1:self$components){
                                     responsabilities[, k] <- params$p[k]*dnorm(x, params$means[k], params$sd[k])
                                   }
                                   # Normalize the posterior probabilities
                                   #responsabilities <- apply(responsabilities, 2, function(x) ifelse(is.na(x)|is.nan(x), 0, x))
                                   responsabilities <- responsabilities/rowSums(responsabilities)
                                   colnames(responsabilities) <- paste0("B", 1:self$components)
                                   return(dplyr::as_tibble(responsabilities))
                                 },
                                 #' @description
                                 #' Classify the time series in its components
                                 #' @param x vector
                                 classify = function(x){
                                   # Optimal parameters
                                   params <- self$parameters
                                   # Number of observations
                                   n <- length(x)
                                   # E-step: posterior probabilities
                                   responsabilities <- self$E_step(x)
                                   # Classification of each component
                                   x_hat <- matrix(0, nrow = n, ncol = self$components)
                                   B_hat <- matrix(0, nrow = n, ncol = self$components)
                                   classification <- c()
                                   for(i in 1:n) {
                                     classification[i] <- which.max(responsabilities[i,])
                                     B_hat[i, classification[i]] <- 1
                                     x_hat[i,] <- B_hat[i,]*x[i]
                                   }
                                   colnames(B_hat) <- paste0("B", 1:self$components)
                                   colnames(x_hat) <- paste0("x", 1:self$components)
                                   B_hat <- dplyr::as_tibble(B_hat)
                                   x_hat <- dplyr::as_tibble(x_hat)
                                   z_hat <- dplyr::as_tibble((x_hat-params$means)/params$sd*B_hat)
                                   colnames(z_hat) <- paste0("z", 1:self$components)
                                   # Output
                                   dplyr::bind_cols(classification = classification, B_hat, x_hat, z_hat)
                                 },
                                 #' @description
                                 #' Fit the parameters with mclust package
                                 #' @param x vector
                                 #' @param weights observations weights, if a weight is equal to zero the observation is excluded, otherwise is included with unitary weight.
                                 #' When `missing` all the available observations will be used.
                                 fit = function(x, weights){
                                   # Initialization
                                   n <- length(x)
                                   # Weights
                                   if (missing(weights)){
                                     w <- rep(1, n)
                                   } else {
                                     w <- ifelse(weights == 0, 0, 1)
                                   }
                                   w[is.na(x)] <- 0
                                   # Fitted parameters
                                   clust <- mclust::Mclust(x[w!=0], G = self$components, modelNames = c("V"), verbose = FALSE)
                                   # Update the parameters
                                   idx_order <- order(clust$parameters$mean)
                                   self$means <- clust$parameters$mean[idx_order]
                                   self$sd <- clust$parameters$variance$scale[idx_order]
                                   self$p <- clust$parameters$pro[idx_order]
                                   # Assign a name to the parameters
                                   names(self$means) <- paste0("mu", 1:self$components)
                                   names(self$sd) <- paste0("sd", 1:self$components)
                                   names(self$p) <- paste0("p", 1:self$components)

                                   # E-step: posterior probabilities
                                   private$responsabilities <- self$E_step(x)
                                   # Calculate the log-likelihood
                                   private$..loglik <- self$logLik(x[w!=0])
                                   # Final classification of each component
                                   private$..fitted <- self$classify(x)
                                   # Add number of observations
                                   private$..moments$nobs <- length(x[w!=0])
                                   # Add time series
                                   private$x <- x
                                   private$w <- w
                                 },
                                 #' @description
                                 #' Fit the parameters with EM-algorithm
                                 #' @param x vector
                                 #' @param weights observations weights, if a weight is equal to zero the observation is excluded, otherwise is included with unitary weight.
                                 #' When `missing` all the available observations will be used.
                                 EM = function(x, weights){
                                   if (missing(x)){
                                     x <- private$x
                                   }
                                   if (missing(weights)){
                                     w <- private$w
                                   } else {
                                     w <- ifelse(weights == 0, 0, 1)
                                   }
                                   w[is.na(x)] <- 0
                                   # Initialization
                                   n <- length(x)
                                   # Update n parameter
                                   n_w <- length(x[w != 0])
                                   # Initialization
                                   log_likelihood <- 0
                                   previous_log_likelihood <- -Inf
                                   prev_responsabilities <- matrix(0, nrow = n, ncol = self$components)
                                   previous_params <- list(means = self$means, sd = self$sd, p = self$p)
                                   # EM Algorithm
                                   for (iteration in 1:self$maxit){
                                     # E-step: posterior probabilities
                                     responsabilities <- self$E_step(x, previous_params)
                                     # Optimal parameters
                                     params <- previous_params
                                     # M-step: Update the parameters
                                     for(k in 1:self$components){
                                       # Normalizing factor for each group
                                       n_k <- sum(w*responsabilities[, k], na.rm = TRUE)
                                       # Mean parameter k-component
                                       params$means[k] <- sum(responsabilities[, k]*w*x, na.rm = TRUE)/n_k
                                       # Std. deviation k-component
                                       params$sd[k] <- sqrt(sum(responsabilities[, k]*w*(x - params$means[k])^2, na.rm = TRUE)/n_k)
                                       # Probability k-component
                                       params$p[k] <- n_k/n_w
                                     }
                                     # Check divergence
                                     if (any(params$p > 0.9)) {
                                       warning("Probs > 0.9 break")
                                       params <- previous_params
                                       break
                                     }
                                     # Calculate the log-likelihood
                                     log_likelihood <- self$logLik(x[w!=0], params)
                                     # Check for convergence
                                     stop_condition <- abs(log_likelihood - previous_log_likelihood) < self$abstol
                                     if (stop_condition) {
                                       break
                                     } else {
                                       # Update log-likelihood
                                       previous_log_likelihood <- log_likelihood
                                       # Update parameters
                                       previous_params <- params
                                     }
                                   }
                                   # Update the parameters
                                   idx_order <- order(params$means)
                                   self$means <- params$means[idx_order]
                                   self$sd <- params$sd[idx_order]
                                   self$p <- params$p[idx_order]
                                   # Assign a name to the parameters
                                   names(self$means) <- paste0("mu", 1:self$components)
                                   names(self$sd) <- paste0("sd", 1:self$components)
                                   names(self$p) <- paste0("p", 1:self$components)

                                   # E-step: posterior probabilities
                                   private$responsabilities <- self$E_step(x)
                                   # Calculate the log-likelihood
                                   private$..loglik <- self$logLik(x[w!=0])
                                   # Final classification of each component
                                   private$..fitted <- self$classify(x)
                                   # Add number of observations
                                   private$..moments$nobs <- n_w
                                   # Add time series
                                   private$x <- x
                                   private$w <- w
                                 },
                                 #' @description
                                 #' Update the responsabilities, means, sd, p and recompute log-likelihood and fitted data.
                                 #' @param x vector
                                 #' @param weights observations weights, if a weight is equal to zero the observation is excluded, otherwise is included with unitary weight.
                                 #' When `missing` all the available observations will be used.
                                 #' @param means Numeric vector of means parameters.
                                 #' @param sd Numeric vector of std. deviation parameters.
                                 #' @param p Numeric vector of probability parameters.
                                 update = function(x, weights, means, sd, p){
                                   if (missing(x)){
                                     x <- private$x
                                   }
                                   if (missing(weights)){
                                     w <- private$w
                                   } else {
                                     w <- ifelse(weights == 0, 0, 1)
                                   }
                                   w[is.na(x)] <- 0
                                   # Update mean parameters
                                   if (!missing(means)) {
                                     self$means <- means
                                   }
                                   # Update Std. deviations parameters
                                   if (!missing(sd)) {
                                     self$sd <- sd
                                   }
                                   # Update probability parameters
                                   if (!missing(p)) {
                                     self$p <- p
                                   }
                                   # Reorder the parameters
                                   idx_order <- order(self$means)
                                   self$means <- self$means[idx_order]
                                   self$sd <- self$sd[idx_order]
                                   self$p <- self$p[idx_order]
                                   # Assign a name to the parameters
                                   names(self$means) <- paste0("mu", 1:self$components)
                                   names(self$sd) <- paste0("sd", 1:self$components)
                                   names(self$p) <- paste0("p", 1:self$components)

                                   # Update posterior probabilities
                                   private$responsabilities <- self$E_step(x)
                                   # Update log-likelihood
                                   private$..loglik <- self$logLik(x[w!=0])
                                   # Final classification of each component
                                   private$..fitted <- self$classify(x)
                                   # Add number of observations
                                   private$..moments$nobs <- length(x[w!=0])
                                   # Add time series
                                   private$x <- x
                                   private$w <- w
                                 }
                               ),
                               private  = list(
                                 x = NA,
                                 w = NA,
                                 ..loglik = NA,
                                 ..fitted = NA,
                                 responsabilities = NA,
                                 ..moments = list(mean = NA, variance = NA, skewness = NA, kurtosis = NA, nobs = NA)
                               ),
                               active = list(
                                 #' @field parameters named list with mixture parameters.
                                 parameters = function(){
                                   # Output data
                                   list(means = self$means, sd = self$sd, p = self$p)
                                 },
                                 #' @field model Tibble with mixture parameters, in order means, sd, p.
                                 model = function(){
                                   dplyr::bind_cols(as.list(self$means), as.list(self$sd), as.list(self$p))
                                 },
                                 #' @field loglik log-likelihood of the fitted series.
                                 loglik = function(){
                                   private$..loglik
                                 },
                                 #' @field fitted fitted series
                                 fitted = function(){
                                   private$..fitted
                                 },
                                 #' @field moments Tibble with the theoric moments and the number of observations used for fit.
                                 moments = function(){
                                   # Compute moments
                                   private$..moments$mean <- sum(self$means*self$p)
                                   private$..moments$variance <- sum((self$means^2 + self$sd^2)*self$p) - private$..moments$mean^2
                                   # private$..moments$skewness
                                   # private$..moments$kurtosis
                                   return(dplyr::bind_cols(private$..moments))
                                 }
                               ))



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



