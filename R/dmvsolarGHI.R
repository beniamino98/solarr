#' Multivariate PDF GHI
#'
#' @param x vector or quantiles.
#' @param Ct clear sky radiation
#' @param st list of solarTransform.
#' @param pdf_Yt joint density of Y_t.
#'
#' @name dmvsolarGHI
#' @rdname dmvsolarGHI
#' @aliases dmvsolarGHI
#' @aliases pmvsolarGHI
#' @aliases qmvsolarGHI
#'
#' @keywords distributions
#' @details Version 1.0.0.
#' @export
dmvsolarGHI <- function(x, Ct, pdf_Yt, st){
  # Number of joint pairs
  n.obs <- nrow(x)
  # Number of locations
  n.loc <- ncol(x)
  # Normalized risk drivers
  g_prime_eta <- matrix(0, nrow = n.obs, ncol = n.loc)
  # Transformed variable
  Yt <-  matrix(0, nrow = n.obs, ncol = n.loc)
  for(l in 1:n.loc) {
    st_l <- st[[l]]
    g_prime_eta[,l] <- st_l$g_prime(st_l$eta(x[,l], Ct[l])) / (st_l$beta * Ct[l])
    Yt[,l] <- st_l$RY(x[,l], Ct[l])
  }
  # Total jacobian
  jac <- abs(apply(g_prime_eta, 1, prod))
  # Probabilities
  probs <- pdf_Yt(Yt) * jac
  probs[is.nan(probs)] <- 0
  return(probs)
}

#' Distribution function for the GHI
#'
#' @rdname dmvsolarGHI
#' @export
pmvsolarGHI <- function(x, Ct, cdf_Yt, st){
  # Number of joint pairs
  n.obs <- nrow(x)
  # Number of locations
  n.loc <- ncol(x)
  # Convert R into y-coordinates
  Yt <- matrix(0, nrow = n.obs, ncol = n.loc)
  # Bounds
  bounds <- purrr::map2_df(st, Ct, ~setNames(.y*.x$bounds("Kt"), c("R_min", "R_max")))
  for(l in 1:n.loc){
    Yt[,l] <- st[[l]]$RY(x[,l], Ct[l])
    Yt[,l][x[,l] <= bounds$R_min[l]] <- -Inf
    Yt[,l][x[,l] >= bounds$R_max[l]] <- Inf
  }
  cdf_Yt(Yt)
}

#' Directional quantile
#'
#' @rdname dmvsolarGHI
#' @export
qmvsolarGHI <- function(p, Ct, cdf_Yt, st, v = NULL, tol = 1e-8, range_s = c(-100, 100)) {
  # Number of quantiles
  n.obs <- length(p)
  # Reorder range
  range_s <- range(range_s)
  # Number of locations
  n.loc <- length(Ct)
  # Reference point (mid-value)
  x0 <- purrr::map2_dbl(st, Ct, ~ mean(.y * .x$bounds("Kt")))
  # default direction: all ones
  if (is.null(v)) v <- rep(1, n.loc)
  v <- as.numeric(v)
  v <- v / sqrt(sum(v^2))
  # scalar CDF along the ray
  G <- function(s, target) {
    x <- matrix(x0 + s * v, nrow = 1, ncol = n.loc)
    pmvsolarGHI(x, Ct, cdf_Yt, st) - target
  }
  # Initialization
  q_Rt <- matrix(NA_real_, nrow = n.obs, ncol = n.loc)
  for (t in seq_len(n.obs)) {
    # Root
    root <-  uniroot(
      f = G,
      target = p[t],
      lower = range_s[1],
      upper = range_s[2],
      tol = tol)$root
    # Quantile
    q_Rt[t, ] <- x0 + root * v
  }
  q_Rt
}
