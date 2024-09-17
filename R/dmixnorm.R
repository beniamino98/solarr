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
