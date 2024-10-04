#' Gaussian mixture random variable
#'
#' Gaussian mixture density, distribution, quantile and random generator.
#'
#' @param x vector of quantiles or probabilities.
#' @param n number of observations. If `length(n) > 1`, the length is taken to be the number required.
#' @param means vector of means parameters.
#' @param sd vector of std. deviation parameters.
#' @param p vector of probability parameters.
#' @param log.p logical; if `TRUE`, probabilities p are given as `log(p)`.
#' @param log logical; if `TRUE`, probabilities are returned as `log(p)`.
#' @param lower.tail logical; if TRUE (default), probabilities are `P[X < x]` otherwise, `P[X > x]`.
#'
#' @references Mixture Models [\href{https://en.wikipedia.org/wiki/Mixture_model}{W}].
#'
#' @examples
#' # Parameters
#' means = c(-3,0,3)
#' sd = rep(1, 3)
#' p = c(0.2, 0.3, 0.5)
#' # Density function
#' dmixnorm(3, means, sd, p)
#' # Distribution function
#' dmixnorm(c(1.2, -3), means, sd, p)
#' # Quantile function
#' qmixnorm(0.2, means, sd, p)
#' # Random generator
#' rmixnorm(1000, means, sd, p)
#'
#' @rdname dmixnorm
#' @name dmixnorm
#' @aliases dmixnorm
#' @aliases pmixnorm
#' @aliases qmixnorm
#' @aliases rmixnorm
#' @export
dmixnorm <- function(x, means = rep(0, 2), sd = rep(1, 2), p = rep(1/2, 2),
                     log = FALSE){

  # Number of components
  k <- length(means)
  # List of parameters
  params <- list(mean = means, sd = sd, p = p)
  # Density
  probs <- c()
  for(i in 1:length(x)){
    probs[i] <- 0
    for(s in 1:k){
      probs[i] <- probs[i] + p[s]*dnorm(x[i], mean = means[s], sd = sd[s])
    }
  }
  # Log-probability
  if (log) {
    probs <- base::log(probs)
  }
  return(probs)
}

#' @rdname dmixnorm
#' @export
pmixnorm <- function(x, means = rep(0, 2), sd = rep(1, 2), p = rep(1/2, 2),
                     lower.tail = TRUE, log.p = FALSE){

  # Number of components
  k <- length(means)
  # List of parameters
  params <- list(mean = means, sd = sd, p = p)
  # Distribution
  probs <- c()
  for(i in 1:length(x)){
    probs[i] <- 0
    for(s in 1:k){
      probs[i] <- probs[i] + p[s]*pnorm(x[i], mean = means[s], sd = sd[s], lower.tail = lower.tail)
    }
  }
  # Log-probability
  if (log.p) {
    probs <- base::log(probs)
  }
  return(probs)
}

#' @rdname dmixnorm
#' @export
qmixnorm <- function(x, means = rep(0, 2), sd = rep(1, 2), p = rep(1/2, 2),
                     lower.tail = TRUE, log.p = FALSE) {
  # Log-probabilities
  probs <- x
  if (log.p) {
    probs <- base::exp(probs)
  }
  # Distribution function
  cdf <- function(x) pmixnorm(x, means = means, sd = sd, p = p, lower.tail = lower.tail)
  # Empirical quantile optimizer
  quantile <- function(x, p){(cdf(x) - p)^2}
  # Quantiles
  x <- c()
  for(i in 1:length(probs)){
    x[i] <- suppressWarnings(optim(par = 0, quantile, p = probs[i])$par)
  }
  return(x)
}

#' @rdname dmixnorm
#' @export
rmixnorm <- function(n, means = rep(0, 3), sd = rep(1, 3), p = rep(1/3, 3)){

  # Number of components
  k <- length(means)

  X <- matrix(NA, nrow = n, ncol = k)
  B <- matrix(0, nrow = n, ncol = k)
  index <- 1:n
  for(s in 1:k){
    # Simulated bernoulli
    if (s == k) {
      B[index,][,s] <- 1
    } else {
      B[index,][,s] <- rbinom(n, 1, p[s]/sum(p[s:k]))
    }
    # Simulated component
    X[B[,s] == 1, s] <- rnorm(sum(B[, s]), mean = means[s], sd = sd[s])
    # Update number of remaining elements
    n <- n - sum(B[,s])
    # Update the remaining indexes
    index <- index[!(index %in% which(!is.na(X[,s])) )]
    # Substitue NA values with 0
    X[,s] <- ifelse(is.na(X[,s]), 0, X[,s])
  }
  colnames(X) <- paste0("X", 1:k)
  colnames(B) <- paste0("B", 1:k)
  sim <- dplyr::bind_cols(t = 1:nrow(X), X, B, X = rowSums(X))
  return(sim)
}


#' Multivariate Gaussian mixture random variable
#'
#' Multivariate Gaussian mixture density, distribution, quantile and random generator.
#'
#' @examples
#' # Means components
#' mean_1 = c(-1.8,-0.4)
#' mean_2 = c(0.6, 0.5)
#' # Dimension of the random variable
#' j = length(mean_1)
#' # Matrix of means
#' means = matrix(c(mean_1, mean_2), j,j, byrow = TRUE)
#'
#' # Variance components
#' var_1 = c(1,1.4)
#' var_2 = c(1.3, 1.2)
#' # Matrix of variances
#' sigma2 = matrix(c(var_1, var_2), j,j, byrow = TRUE)
#'
#' # Correlations
#' rho <- c(rho_1 = 0.2, rho_2 = 0.3)
#'
#' # Probability for each component
#' p <- c(0.4, 0.6)
#'
#' x <- matrix(c(0.1,-0.1), nrow = 1)
#' dmvmixnorm(x, means, sigma2, p, rho)
#' pmvmixnorm(x, means, sigma2, p, rho)
#' qmvmixnorm(0.35, means, sigma2, p, rho)
#' @name dmvmixnorm
#' @rdname dmvmixnorm
#' @aliases pmvmixnorm
#' @aliases qmvmixnorm
#' @export
dmvmixnorm <- function(x, means = matrix(0, 2, 2), sigma2 = matrix(1, 2, 2), p = rep(1/2, 2), rho = c(0,0), log = FALSE){

  if(is.vector(x)){
    x <- matrix(x, nrow = 1)
  }
  # Number of observations
  n <- nrow(x)
  # Number of components
  k <- nrow(means)
  # Number of variables
  j <- ncol(means)
  # Covariance matrix
  cv_ <- list()
  for(s in 1:k){
    cv_k <- diag(sigma2[s,])
    cv_k[upper.tri(cv_k)] <- cv_k[lower.tri(cv_k)] <- rho[1]*prod(sigma2[s,])
    cv_ [[s]] <- cv_k
  }
  # Density
  probs <- rep(0, n)
  for(i in 1:n){
    for(s in 1:k){
      probs[i] <- probs[i] + p[s]*mvtnorm::dmvnorm(x[i,], mean = means[s,], sigma  = cv_[[s]])
    }
  }
  # Log-probability
  if (log) {
    probs <- base::log(probs)
  }
  return(probs)
}

#' @rdname dmvmixnorm
#' @export
pmvmixnorm <- function(x, means = matrix(0, 2, 2), sigma2 = matrix(1, 2, 2), p = rep(1/2, 2), rho = c(0,0), lower = -Inf, log.p = FALSE){
  if(is.vector(x)){
    x <- matrix(x, nrow = 1)
  }
  # Number of observations
  n <- nrow(x)
  # Number of components
  k <- nrow(means)
  # Number of variables
  j <- ncol(means)
  # Covariance matrix
  cv_ <- list()
  for(s in 1:k){
    cv_k <- diag(sigma2[s,])
    cv_k[upper.tri(cv_k)] <- cv_k[lower.tri(cv_k)] <- rho[1]*prod(sigma2[s,])
    cv_ [[s]] <- cv_k
  }
  # Density
  probs <- rep(0, n)
  for(i in 1:n){
    for(s in 1:k){
      probs[i] <- probs[i] + p[s]*mvtnorm::pmvnorm(lower = lower, upper = x[i,], mean = means[s,], sigma = cv_[[s]])
    }
  }
  # Log-probability
  if (log.p) {
    probs <- base::log(probs)
  }
  return(probs)
}

#' @rdname dmvmixnorm
#' @export
qmvmixnorm <- function(x, means = matrix(0, 2, 2), sigma2 = matrix(1, 2, 2), p = rep(1/2, 2), rho = c(0,0), log.p = FALSE){
  # Log-probabilities
  probs <- x
  if (log.p) {
    probs <- base::exp(probs)
  }
  # Distribution function
  cdf <- function(x) pmvmixnorm(x, means = means, sigma2 = sigma2, p = p, rho = rho)
  # Empirical quantile optimizer
  quantile_ <- function(x, prob){(cdf(x) - prob)^2}
  # Quantiles
  q <- matrix(0, nrow = length(probs), ncol = ncol(means))
  for(i in 1:length(probs)){
    q[i,] <- suppressWarnings(optim(par = rep(0, ncol(means)), quantile_, prob = probs[i])$par)
  }
  return(q)
}
