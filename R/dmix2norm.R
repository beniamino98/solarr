#' Density of a bivariate Gaussian mixture with 2 components.
#'
#' @param x matrix, on the rows the observations, on the columns the locations.
#' @param means matrix of means for each state.
#' @param Sigma List of variance-covariance matrices for each state.
#' @param probs vector of joint probabilities for each state.
#' @param log Logical, when `TRUE`, probabilities are returned as `log(p)`.
#' @param log.p Logical, when `TRUE`, probabilities p are given as `log(p)`.
#' @param lower.tail Logical, when `TRUE`, the default, the computed probabilities are \eqn{\mathbb{P}(X < x)}, otherwise \eqn{\mathbb{P}(X \ge x)}.
#' @examples
#' library(tidyverse)
#' x <- matrix(runif(200, -3, 3), ncol = 2)
#' dmix2norm(x)
#' pmix2norm(x)
#' qmix2norm(x)
#'
#' @rdname d2mixnorm
#' @aliases d2mixnorm
#' @aliases p2mixnorm
#' @aliases q2mixnorm
#' @keywords distributions
#' @details Version 1.0.0.
#' @export
dmix2norm <- function(x, means = matrix(0, 4, 2), Sigma = rep(list(diag(1, 2)), 4), probs = rep(1/4, 4), log = FALSE){
  # Number of joint pairs
  n.obs <- nrow(x)
  # Number of locations
  stopifnot(ncol(x) == 2)
  n.loc <- 2
  # Number of components
  n.com <- 4
  # Components
  x1 <- x[,1]
  x2 <- x[,2]
  # Compute density for each state
  f_s <- matrix(0, nrow = n.obs, ncol = n.com)
  for(s in 1:n.com){
    m1 <- means[s,1]
    m2 <- means[s,2]
    s1 <- sqrt(Sigma[[s]][1,1])
    s2 <- sqrt(Sigma[[s]][2,2])
    # Standardized components
    z1 <- (x1 - m1) / s1
    z2 <- (x2 - m2) / s2
    iD <- diag(1/c(s1, s2))
    # Compute correlation
    rho12 <- Sigma[[s]][1,2] / (s1 * s2)
    f_s[, s] <- probs[s] * exp(- (z1^2 - 2*rho12 * z1 * z2 + z2^2) / (2*(1-rho12^2))) / (2 * base::pi * s1 * s2 * sqrt(1-rho12^2))
  }
  # Total density
  p <- rowSums(f_s)

  # Log-probability
  if (log) {
    p <- base::log(p)
  }
  return(p)
}

#' Distribution of a bivariate Gaussian mixture with 2 components each
#'
#' @rdname dmix2norm
#' @export
pmix2norm <- function(x, means = matrix(0, 4, 2), Sigma = rep(list(diag(1, 2)), 4), probs = rep(1/4, 4), log = FALSE){
  # Number of joint pairs
  n.obs <- nrow(x)
  # Number of locations
  stopifnot(ncol(x) == 2)
  n.loc <- 2
  # Number of components
  n.com <- 4
  # Components
  x1 <- x[,1]
  x2 <- x[,2]
  # Compute distribution for each state
  F_s <- matrix(0, nrow = n.obs, ncol = n.com)
  for(s in 1:n.com){
    m1 <- means[s,1]
    m2 <- means[s,2]
    s1 <- sqrt(Sigma[[s]][1,1])
    s2 <- sqrt(Sigma[[s]][2,2])
    # Standardized components
    z1 <- (x1 - m1) / s1
    z2 <- (x2 - m2) / s2
    iD <- diag(1/c(s1, s2))
    # Compute correlation
    rho12 <- Sigma[[s]][1,2] / (s1 * s2)
    F_s[, s] <- probs[s] * pbivnorm::pbivnorm(z1, z2, rho = rho12)
    F_s[z1 == z2, s] <- 0
  }
  # Total distribution
  p <- rowSums(F_s)
  # Log-probability
  if (log) {
    p <- base::log(p)
  }
  return(p)
}

#' Quantile of a bivariate Gaussian mixture with 2 components each
#'
#' @rdname dmix2norm
#' @export
qmix2norm <- function(p, means = matrix(0, 4, 2), Sigma = rep(list(diag(1, 2)), 4), probs = rep(1/4, 4), v = NULL, x0 = c(-1, 1),
                      range_s = c(-100, 100), log.p = FALSE, lower.tail = TRUE){
  # Log-probability
  if (log.p) {
    p <- exp(p)
  }
  # Lower tail
  if (!lower.tail) {
    p <- 1 - p
  }
  # Number of joint pairs
  n.obs <- length(p)
  # Safety check
  stopifnot(ncol(means) == 2 | length(probs) == 4)
  #all(nrow(c(Sigma, diag(2, 2, 2)) == 3)
  # Number of locations
  n.loc <- 2
  # default direction: all ones
  if (is.null(v)) v <- rep(1, n.loc)
  v <- as.numeric(v)
  v <- v / sqrt(sum(v^2))
  # Reference point (mid-value)
  if (is.null(x0)){
    x0 <- apply(means, 2, sum)/2
  }
  # scalar CDF along the ray
  G <- function(s, target) {
    x <- matrix(x0 + s * v, nrow = 1, ncol = n.loc)
    pmix2norm(x, means, Sigma, probs) - target
  }
  # Initialization
  q_x <- matrix(NA_real_, nrow = n.obs, ncol = n.loc)
  for (t in seq_len(n.obs)) {
    # Root
    root <-  uniroot(
      f = G,
      target = p[t],
      lower = range_s[1],
      upper = range_s[2],
      tol = 1e-8)$root
    # Quantile
    q_x[t, ] <- x0 + root * v
  }
  q_x
}


#' Conditional means component l given Yv = yv
#' @param y_v conditioning value Yv = yv
#' @param M_Yl mean vector for component l (state 1, state0)
#' @param M_Yv mean vector for component v (state 1, state0)
#' @param S_Yv std. deviation vector for component v (state 1, state0)
#' @param Cv_Yl_Yv covariance between Yl and Yv.
#' @keywords distributions
#' @export
mix2norm_M_Yl_given_Yv <- function(y_v, M_Yl, M_Yv, S_Yv, Cv_Yl_Yv){
  c(mu11 = M_Yl[1] + Cv_Yl_Yv[1]/S_Yv[1]^2 * (y_v - M_Yv[1]),
    mu10 = M_Yl[1] + Cv_Yl_Yv[2]/S_Yv[2]^2 * (y_v - M_Yv[2]),
    mu01 = M_Yl[2] + Cv_Yl_Yv[3]/S_Yv[1]^2 * (y_v - M_Yv[1]),
    mu00 = M_Yl[2] + Cv_Yl_Yv[4]/S_Yv[2]^2 * (y_v - M_Yv[2]))
}

#' Conditional std. deviation component l given Yv = yv
#' @param S_Yl std. deviation vector for component l (state 1, state0)
#' @param S_Yv std. deviation vector for component v (state 1, state0)
#' @param Cv_Yl_Yv covariance between Yl and Yv.
#' @keywords distributions
#' @export
mix2norm_S_Yl_given_Yv <- function(S_Yl, S_Yv, Cv_Yl_Yv){
  res <- c(sd11 = c(S_Yl[1]^2 - (Cv_Yl_Yv[1]/S_Yv[1])^2),
           sd10 = c(S_Yl[1]^2 - (Cv_Yl_Yv[2]/S_Yv[2])^2),
           sd01 = c(S_Yl[2]^2 - (Cv_Yl_Yv[3]/S_Yv[1])^2),
           sd00 = c(S_Yl[2]^2 - (Cv_Yl_Yv[4]/S_Yv[2])^2))
  sqrt(res)
}

#' Conditional probabilities component l given Yv = yv
#'
#' @param y_v conditioning value Yv = yv
#' @param M_Yv mean vector for component v (state 1, state0)
#' @param S_Yv std. deviation vector for component v (state 1, state0)
#' @param probs Joint probabilities
#' @keywords distributions
#' @export
mix2norm_probs_given_Yv <- function(y_v, M_Yv, S_Yv, probs){
  # Unconditional pdf of Yv
  pdf_Yv_1 <- dnorm(y_v, M_Yv[1], S_Yv[1])
  pdf_Yv_0 <- dnorm(y_v, M_Yv[2], S_Yv[2])
  # Numerator
  num <- probs * c(pdf_Yv_1, pdf_Yv_0, pdf_Yv_1, pdf_Yv_0)
  # Denominator
  den <- sum(num)
  num/den
}

