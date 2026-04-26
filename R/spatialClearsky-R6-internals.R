#' Control parameters for a `spatialClearsky` object
#'
#'
#' @examples
#' control = control_spatialClearsky()
#' @return Named list of control parameters.
#' @keywords spatialClearsky
#' @note Version 1.0.1
#' @rdname control_spatialClearsky
#' @export
control_spatialClearsky <- function(lambda = 0,
                                    eps = 0.1,
                                    gamma = NA,
                                    normalize_loss = FALSE,  # if TRUE, interpret gamma as normalized weight and set lambda accordingly
                                    penalty = c("none", "ridge"),
                                    mask_time = NULL,        # length p: weights for time-coef blocks (0=no penalty, 1=penalize)
                                    osqp_settings = list(max_iter = 30000L)){
  # Match the penalty
  penalty = match.arg(penalty, choices = c("none", "ridge"))
  # Settings for osqp
  settings <- modifyList(list(verbose = TRUE, eps_abs = 1e-5, eps_rel = 1e-5, polish = TRUE), osqp_settings)

  structure(
    list(
      lambda = lambda,
      gamma = gamma,
      eps = eps,
      normalize_loss = normalize_loss,  # if TRUE, interpret gamma as normalized weight and set lambda accordingly
      penalty = penalty,
      mask_time = mask_time,        # length p: weights for time-coef blocks (0=no penalty, 1=penalize)
      settings = settings
    ),
    class = c("control", "list")
  )
}

#' Specification for `spatialClearsky` model
#'
#' @param formula_params Character, formula for spatial parameters.
#' @param C_tilde_col Integer scalar, number of combinations of sines and cosines.
#' @param R_col Name for GHI.
#' @inheritParams seasonalClearsky_spec
#' @examples
#' spec_Ct = spatialClearsky_spec()
#' @return Named list of control parameters.
#' @keywords spatialClearsky
#' @note Version 1.0.0
#' @rdname spatialClearsky_spec
#' @export
spatialClearsky_spec <- function(formula_params, C_tilde_col = "clearsky", R_col = "GHI",
                                 order = 1, order_H0 = 1, method_H0 = "spencer", period = 365,
                                 include.intercept = TRUE, include.trend = FALSE){
  # Method to compute H0
  method_H0 <- match.arg(method_H0, choices = c("spencer", "cooper"))
  # Build the base formula
  formula_Ct <- paste0(C_tilde_col, " ~ H0")
  # Add orders for H0
  if (order_H0 > 1) {
    for(i in 2:order_H0){
      formula_Ct <- paste0(formula_Ct, " + ", paste0("H0_", i))
    }
  }
  formula_Ct <- ifelse(include.trend, paste0(formula_Ct, " + n"), formula_Ct)
  formula_Ct <- ifelse(include.intercept, formula_Ct, paste0(formula_Ct, "-1"))
  # Convert in formula
  formula_Ct <- seasonalModel_formula(formula_Ct, order = order, t_idx = "n")$formula
  # Number of regressors
  p <- length(formula.tools::rhs.vars(formula_Ct)) + include.intercept
  # Parameters formula
  params.include.intercept <- !stringr::str_detect(formula_params, "-1")
  # Convert in formula
  formula_params <- as.formula(formula_params)
  # Number of base
  K <- length(formula.tools::rhs.vars(formula_params)) + params.include.intercept

  structure(
    list(
      formula_Ct = formula_Ct,
      formula_params = formula_params,
      dim = list(p = p, K = K),
      C_tilde_col = C_tilde_col,
      R_col = R_col,
      order = order,
      order_H0 = order_H0,
      method_H0 = method_H0,
      period = period,
      include.intercept = c(Ct = include.intercept, Base = params.include.intercept),
      include.trend = include.trend
    ),
    class = c("seasonalClearsky_spec", "list")
  )
}

#' Fit the spatial model for clear sky radiation.
#'
#' @rdname spatialClearsky_optimizer_CLS
#' @name spatialClearsky_optimizer_CLS
#' @keywords spatialClearsky
#' @note Version 1.0.0
#' @export
spatialClearsky_optimizer_CLS <- function(formula_Ct,
                                          formula_params,
                                          eps = 0.01,
                                          data_list,
                                          coords,
                                          C_tilde_col = "clearsky",
                                          R_col = "GHI",
                                          lambda = 0,
                                          gamma = NA,
                                          normalize_loss = FALSE,
                                          penalty = c("none", "ridge"),
                                          mask_time = NULL,
                                          settings = list(verbose = TRUE, eps_abs = 1e-5, eps_rel = 1e-5,
                                                          polish = TRUE, max_iter = 4000L)){
  # Required package
  if (!requireNamespace("osqp", quietly = TRUE)) stop("Package 'osqp' required.")
  # Dimension check
  stopifnot(is.list(data_list), nrow(coords) == length(data_list))
  # Number of locations
  L <- length(data_list)
  # Clear-sky
  stopifnot(all(vapply(data_list, function(d) C_tilde_col %in% names(d), logical(1))))
  C_tilde_list <- lapply(data_list, function(d) d[[C_tilde_col]])
  # GHI
  stopifnot(all(vapply(data_list, function(d) R_col %in% names(d), logical(1))))
  R_list <- lapply(data_list, function(d) d[[R_col]])
  stopifnot(length(C_tilde_list) == L, length(R_list) == L)

  # --- build spatial basis B_i for each location (K columns)
  B_mat <- model.matrix(formula_params, data = mutate(coords, x = 2)) # L x K
  # Number of parameters
  K <- ncol(B_mat)
  # --- build time/physics design X_i for each location (T_i x p)
  # formula_cs should be RHS-only (~ ...)
  X_list <- lapply(data_list, function(d) {
    model.matrix(as.formula(formula_Ct), data = d)
  })
  p <- ncol(X_list[[1]])
  stopifnot(all(vapply(X_list, ncol, integer(1)) == p))

  # --- optional time-mask for ridge penalty
  if (is.null(mask_time)) {
    mask_time <- rep(1, p)
  } else {
    stopifnot(length(mask_time) == p, all(mask_time >= 0))
    mask_time <- as.numeric(mask_time)
  }

  # --- build stacked y and b
  y <- unlist(C_tilde_list, use.names = FALSE)
  b <- unlist(Map(function(r) r + eps, R_list), use.names = FALSE)
  N <- length(y)

  # --- build Z as sparse by stacking Z_i = [B_{i1} X_i, ..., B_{iK} X_i]
  # This avoids explicitly forming kron(B_i, X_i) and is faster/sparser.
  # template column names (same for all i)
  cnX <- colnames(X_list[[1]]); if (is.null(cnX)) cnX <- paste0("x", seq_len(p))
  cnB <- colnames(B_mat);       if (is.null(cnB)) cnB <- paste0("B", seq_len(K))
  colnames_Z <- as.vector(outer(cnB, cnX, function(a, b) paste0(a, ":", b)))  # length K*p

  # --- build Z as sparse by stacking Z_i = [B_{i1} X_i, ..., B_{iK} X_i]
  Z_list <- vector("list", L)
  for (i in seq_len(L)) {
    Xi <- X_list[[i]]
    stopifnot(ncol(Xi) == p)
    Xi_sp <- if (inherits(Xi, "Matrix")) Xi else Matrix::Matrix(Xi, sparse = TRUE)
    Xi_sp <- methods::as(Xi_sp, "dgCMatrix")
    Bi <- as.numeric(B_mat[i, ])  # length K
    # blocks: Xi * Bi[j]
    Zi_blocks <- lapply(seq_len(K), function(j) Xi_sp * Bi[j])
    # robust multi-cbind without deprecated cBind
    Zi <- Reduce(Matrix::cbind2, Zi_blocks)
    Zi <- methods::as(Zi, "dgCMatrix")
    colnames(Zi) <- colnames_Z
    Z_list[[i]] <- Zi
  }
  # robust multi-rbind without deprecated rBind
  Z <- Reduce(Matrix::rbind2, Z_list)
  Z <- methods::as(Z, "dgCMatrix")
  stopifnot(nrow(Z) == N, ncol(Z) == p * K)
  # --- penalty matrix P (ridge / none)
  if (penalty == "none" || lambda == 0) {
    Ppen <- Matrix::Diagonal(p * K, x = rep(0, p * K))
  } else {
    # ridge on beta with coefficient-specific weights:
    # beta is stacked as [B1:*Xcols, B2:*Xcols, ..., BK:*Xcols]
    # i.e. K blocks, each of length p (time coefficients).
    # Apply mask_time inside each block.
    P_block <- Matrix::Diagonal(p, x = mask_time)
    Ppen <- Matrix::kronecker(Matrix::Diagonal(K, 1), P_block) # (pK x pK)
  }
  # --- Option A normalization (dimensionless gamma)
  # Normalize objective as:
  #   (1/N)||y - Z beta||^2 + gamma * beta'Pbeta
  # Equivalent to unnormalized:
  #   ||y - Z beta||^2 + (gamma * N) * beta'Pbeta
  if (isTRUE(normalize_loss)) {
    if (is.na(gamma) || is.null(gamma)) stop("normalize_loss=TRUE requires gamma.")
    lambda <- as.numeric(gamma) * N
    message(sprintf("normalize_loss=TRUE: using lambda = gamma * N = %.6g", lambda))
  }

  # --- QP in OSQP form: (1/2) beta' P beta + q'beta  s.t.  Z beta >= b
  # Data term: ||y - Zb||^2 = b' Z'Z b - 2 y'Z b + const
  H_data <- Matrix::crossprod(Z)              # (pK x pK)
  P <- 2 * (H_data + lambda * Ppen)
  P <- Matrix::forceSymmetric(P)
  q <- -2 * Matrix::crossprod(Z, y)
  # Constraints: Z beta >= b  <=>  b <= Z beta <= +Inf
  A <- Z
  l <- b
  u <- rep(Inf, length(b))
  # Solver
  # Settings for osqp
  solver <- osqp::osqp(P = P, q = as.numeric(q), A = A, l = l, u = u, pars = settings)
  res <- solver$Solve()
  # Check for convergence
  if (res$info$status_val < 1 || res$info$status_val > 2) {
    warning(sprintf("OSQP status: %s", res$info$status))
  }
  # Fitted parameters
  beta_hat <- res$x
  # reshape to p x K (columns are spatial basis)
  Beta_hat <- matrix(beta_hat, nrow = p, ncol = K, byrow = FALSE)
  rownames(Beta_hat) <- colnames(X_list[[1]])
  colnames(Beta_hat) <- colnames(B_mat)

  # fitted + slack per location
  fitted_list <- vector("list", L)
  slack_list  <- vector("list", L)
  idx <- 0L
  for (i in seq_len(L)) {
    Ti <- length(C_tilde_list[[i]])
    Zi <- Z[(idx + 1L):(idx + Ti), , drop = FALSE]
    fit_i <- drop(Zi %*% beta_hat)
    fitted_list[[i]] <- fit_i
    slack_list[[i]] <- fit_i - (R_list[[i]] + eps)
    idx <- idx + Ti
  }

  list(
    Beta_hat = Beta_hat,
    coords_basis = B_mat,
    osqp = res,
    lambda = lambda,
    gamma = gamma,
    fitted_list = fitted_list,
    slack_list = slack_list
  )
}



