#' Control parameters for a `seasonalClearsky` object
#'
#' @param delta0 Numeric scalar, initial value for the optimization. The estimated clear sky is increased by delta0.
#' @param quiet Logical, when `FALSE`, the default, the functions displays warning or messages.
#' @inheritParams clearsky_delta
#' @inheritParams clearsky_optimizer_RLS
#' @inheritParams clearsky_optimizer_CLS
#'
#' @details The parameters `ntol`, `lower`, `upper` and `by` are used exclusively in \code{\link{clearsky_delta}}.
#' @examples
#' control = control_seasonalClearsky()
#' @return Named list of control parameters.
#' @keywords clearsky
#' @note Version 1.0.1
#' @rdname control_seasonalClearsky
#' @export
control_seasonalClearsky <- function(delta0 = 1.4, lower = 0, upper = 3, by = 0.001, ntol = 0, eps = 0.01, lambda = 100, quiet = FALSE){

  structure(
    list(
      delta0 = delta0,
      lower = lower,
      upper = upper,
      by = by,
      ntol = ntol,
      eps = eps,
      lambda = lambda,
      quiet = quiet
    ),
    class = c("control", "list")
  )
}

#' Specification for `seasonalClearsky` model
#'
#' @param method Character, method of fitting: can be:  `OLS`, `RLS`, `CLS`
#' @param order Integer scalar, number of combinations of sines and cosines.
#' @param order_H0 Integer scalar, order of the polynomial for extra-terrestrial radiation.
#' @param method_H0 Character, method of computation of the extra-terrestrial radiation: can be:  `spencer`, `cooper`.
#' @param include.intercept Logical, when `TRUE`, the default, the intercept will be included in the clear sky model.
#' @param period Integer scalar, seasonality period. The default is 365.
#'
#' @examples
#' spec_Ct = spec_seasonalClearsky()
#' @return Named list of control parameters.
#' @keywords clearsky
#' @note Version 1.0.0
#' @rdname seasonalClearsky_spec
#' @export
seasonalClearsky_spec <- function(method = "OLS", order = 1, order_H0 = 1, method_H0 = "spencer", period = 365, include.intercept = TRUE, include.trend = FALSE){

  method <- match.arg(method, choices = c("CLS", "RLS", "OLS"))
  method_H0 <- match.arg(method_H0, choices = c("spencer", "cooper"))

  structure(
    list(
      method = method,
      order = order,
      order_H0 = order_H0,
      method_H0 = method_H0,
      period = period,
      include.intercept = include.intercept,
      include.trend = include.trend
    ),
    class = c("seasonalClearsky_spec", "list")
  )
}

#' Fit the seasonal model for clear sky radiation.
#'
#' @param x Numeric vector, time series of solar radiation.
#' @param dates Character or Date vector, time series of dates.
#' @param clearsky Numeric vector, time series of target clear sky radiation.
#' @param lat Numeric scalar, reference latitude. Overwritten when `seasonal_model_Ct` is specified.
#' @param spec See the function \code{\link{seasonalClearsky_spec}} for more details.
#' Overwritten when `seasonal_model_Ct` is specified.
#' @param control See the function \code{\link{control_seasonalClearsky}} for more details.
#' Overwritten when `seasonal_model_Ct` is specified.
#' @param seasonal_model_Ct Optional object of the class \code{\link{seasonalClearsky}}.
#' When submitted `lat`, `spec` and `control` parameters used are recovered from the specified model.
#' @examples
#' x  = solarr::spec$data$GHI
#' dates  = solarr::spec$data$date
#' lat = solarr::spec$coords$lat
#' clearsky = solarr::spec$data$clearsky
#' # Without a clearsky model
#' seasonal_model_Ct <- clearsky_optimizer(x, dates, clearsky, lat)
#' # With a clearsky model
#' clearsky_optimizer(x, dates, clearsky, seasonal_model_Ct = seasonal_model_Ct)
#' @rdname clearsky_optimizer
#' @name clearsky_optimizer
#' @keywords clearsky
#' @note Version 1.0.0
#' @export
clearsky_optimizer <- function(x, dates, clearsky, lat, spec = seasonalClearsky_spec(), control = control_seasonalClearsky(), seasonal_model_Ct){

  if (missing(seasonal_model_Ct)) {
    # Initialize a clear-sky model
    seasonal_model_Ct <- seasonalClearsky$new(spec = spec, control = control)
    # Initialize a fit with OLS estimate
    seasonal_model_Ct$fit(x = x,
                          dates = dates,
                          lat = lat,
                          clearsky = clearsky)
  } else {
    control <- seasonal_model_Ct$control
    spec <- seasonal_model_Ct$spec
    lat <- seasonal_model_Ct$lat
  }

  if (spec$method == "OLS") {
    return(seasonal_model_Ct)
  }

  if (missing(x) || missing(dates) || missing(clearsky)) stop("`x`, `dates` and `cleasky` should be specified when `seasonal_model_Ct` is supplied!")
  # Build a dataset for RLS or CLS refinements
  newdata <- dplyr::tibble(GHI = x,
                           dates = dates,
                           n = number_of_day(dates),
                           clearsky = clearsky,
                           H0 = seasonal_model_Ct$ssf$Hon(n, lat))

  # Refine the fit depending on the method
  if (spec$method == "CLS") {
    seasonal_model_Ct <- clearsky_optimizer_CLS(seasonal_model_Ct, newdata, eps = control$eps)
  } else if (spec$method == "RLS") {
    seasonal_model_Ct <- clearsky_optimizer_RLS(seasonal_model_Ct, newdata, lambda = control$lambda)
  }
  return(seasonal_model_Ct)
}

#' Penalized Least-Squares
#'
#' @param seasonal_model_Ct An object of the class `seasonalClearsky`. See \code{\link{seasonalClearsky}} for more details.
#' @param newdata dataset
#' @param lambda penalization parameter.
#'
#' @rdname clearsky_optimizer_RLS
#' @name clearsky_optimizer_RLS
#' @keywords clearsky
#' @note Version 1.0.0
#' @export
clearsky_optimizer_RLS <- function(seasonal_model_Ct, newdata, lambda = 0, delta_optimize = TRUE){
  # Penality function
  penality_fun <- function(x) pmin(x, 0)^2
  # Target clear-sky
  Ct_target <- newdata$clearsky
  # Target Radiation
  Rt_target <- newdata$GHI
  # Clone the model
  sm <- seasonal_model_Ct$clone(TRUE)
  # Number of observations
  n <- nrow(newdata)
  # Number of parameters
  p <- length(sm$coefficients)
  # Loss function to minimize
  Q_fun <- function(params, sm, data, lambda){
    sm$update(params)
    # Prediction
    C_hat <- sm$predict(newdata = data)
    # Full loss function
    sum((Ct_target - C_hat)^2) + lambda * sum(penality_fun(C_hat - Rt_target))
  }
  # Optimal parameters
  opt <- optim(sm$coefficients, Q_fun, sm = sm, data = newdata, lambda = lambda, hessian = TRUE)
  # Update the parameters
  sm$update(opt$par)
  # Fitted clear-sky
  C_hat <- sm$predict(newdata$n)
  # Residuals
  eps <- newdata$clearsky - C_hat
  # Approx covariance: sigma^2 * H^{-1}
  sigma2 <- sum(eps^2) / (n - p)
  # Variance-covariance matrix
  Vcov <- tryCatch(sigma2 * solve(opt$hessian), error = function(e) NULL)
  # Standard errors
  std.errors <- if (is.null(Vcov)) rep(NA_real_, p) else sqrt(pmax(diag(Vcov), 0))
  std.errors <- setNames(std.errors, names(sm$std.errors))
  # Update the std. errors
  sm$update_std.errors(std.errors)
  if (delta_optimize) {
    sm <- clearsky_optimizer_delta(sm, newdata)
  }
  return(sm)
}

#' Fully constrained optimization (Constrained Least Squares)
#'
#' @inheritParams clearsky_optimizer_RLS
#' @param eps threshold parameter
#'
#' @rdname clearsky_optimizer_CLS
#' @name clearsky_optimizer_CLS
#' @keywords clearsky
#' @note Version 1.0.0
#' @export
clearsky_optimizer_CLS <- function(seasonal_model_Ct, newdata, eps = 0.01, delta_optimize = TRUE){
  # Target clear-sky
  Ct_target <- newdata$clearsky
  # Target Radiation
  Rt_target <- newdata$GHI
  # Clone the model
  sm <- seasonal_model_Ct$clone(TRUE)
  newdata$H0 <- seasonal_model_Ct$ssf$Hon(newdata$n, seasonal_model_Ct$lat)
  # Number of observations
  n <- nrow(newdata)
  # Number of parameters
  p <- length(sm$coefficients)
  # Custom formula
  formula <- sm$.__enclos_env__$private$formula

  # Add Extraterrestrial orders
  if (sm$spec$order_H0 > 1){
    for(i in 2:sm$spec$order_H0){
      newdata[[paste0("H0_", i)]] <- newdata$H0^i
    }
  }

  # Matrix of regressors
  X <- model.matrix(formula, data = newdata)
  # QP matrices: Constrained least squares
  Dmat <- 2 * crossprod(X)                    # 2 X'X
  dvec <- 2 * crossprod(X, Ct_target)  # 2 X'R  (quadprog uses -dvec, but see solve.QP form below)
  # quadprog solves: min 1/2 b'Db - d'b  s.t. A^T b >= b0
  Amat <- t(X)                    # 4 x T  -> A^T b >= b0 is X b >= b0
  bvec <- Rt_target + eps
  # Safe optimization
  safe_solve.QP <- purrr::safely(quadprog::solve.QP)
  opt <- safe_solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 0)
  if (is.null(opt$result)) {
    warning("Errors in the solution of solve.QP!")
    return(NULL)
  }
  # Optimized parameters (constraint)
  c_hat <- setNames(opt$result$solution, names(sm$coefficients))
  # Update parameters
  sm$update(c_hat)
  # Fitted clear-sky
  C_hat <- drop(X %*% c_hat)
  # Active set
  idx_A <- drop(C_hat - bvec == 0)
  # Number of constraints
  k <- sum(idx_A)
  # Residuals
  eps <- Ct_target - C_hat
  # sigma^2
  sigma2 <- sum(eps^2) / (n - p + k)

  H <- crossprod(X)
  if (k == 0) {
    Vcov <- tryCatch(sigma2 * solve(H),  error = function(e) NULL)
  } else {
    A <- X[idx_A,,drop = FALSE]
    Vcov <- tryCatch(sigma2 * (solve(H) - solve(H) %*% t(A) %*% solve(A %*% solve(H) %*% t(A)) %*% A %*% solve(H)), error = function(e) NULL)
  }

  # Standard errors
  std.errors <- if (is.null(Vcov)) rep(NA_real_, p) else sqrt(pmax(diag(Vcov), 0))
  std.errors <- setNames(std.errors, names(sm$std.errors))
  # Update the std. errors
  sm$update_std.errors(std.errors)
  if (delta_optimize) {
    sm <- clearsky_optimizer_delta(sm, newdata)
  }
  return(sm)
}

#' Optimizer for Solar Clear sky
#'
#' Find the best parameter delta for fitting clear sky radiation.
#'
#' @param x Numeric vector, time series of solar radiation.
#' @param Ct Numeric vector, time series of  clear sky radiation.
#' @param lower Numeric scalar, lower bound for grid of delta parameters used for optimization. Default is `0`.
#' @param upper Numeric scalar, upper bound for grid of delta parameters used for optimization. Default is `3`.
#' @param by Numeric scalar, step for the grid. Default is `0.01`.
#' @param ntol Integer scalar, Tolerance for the maximum number of violations admitted of the condition `Ct > x`. The default is `0`.
#' @details
#' Detect the best parameter `delta` such that the condition `x > delta * Ct` is satisfied for all the observations except for `ntol` observations.
#'
#' @name clearsky_delta
#' @rdname clearsky_delta
#' @keywords clearsky
#' @note Version 1.0.1
#' @export
clearsky_delta <- function(x, Ct, lower = 0, upper = 3, by = 0.01, ntol = 0){
  # Grid of points
  grid <- seq(lower, upper, by = by)
  # Loss
  opt <- data.frame(delta = grid, loss = purrr::map_dbl(grid, ~sum(.x*Ct - x < 0)))
  # Return minimum delta satisfying the constraint
  delta <- opt[which(opt$loss <= ntol)[1],]$delta
  return(delta)
}

#' Delta optimizer on a clear-sky model
#'
#' @inheritParams clearsky_optimizer_RLS
#' @rdname clearsky_optimizer_delta
#' @name clearsky_optimizer_delta
#' @keywords clearsky
#' @note Version 1.0.0
#' @export
clearsky_optimizer_delta <- function(seasonal_model_Ct, newdata){
  # Clone the seasonal model
  sm <- seasonal_model_Ct$clone(TRUE)
  # Control
  control <- sm$control
  # Original std.errors
  std.errors_1 <- sm$std.errors
  # Initial fit average clear sky
  newdata$Ct_hat_1 <- sm$predict(newdata$n)
  # SSE
  sse_1 <- sum((newdata$clearsky - newdata$Ct_hat_1)^2)

  # 2. Optimize the fit
  # Optimization: delta_init*Ct_max ~ delta*(a_0 + a_1 H0 + a_2 cos(.) + a_3 sin(.) + a_4 cos(2*.) + a_5 sin(2*.) + ...)
  delta <- clearsky_delta(newdata$GHI, newdata$Ct_hat_1 * control$delta0,
                                    control$lower, control$upper, control$by, control$ntol)
  # Optimized parameters
  c_hat_2 <- sm$coefficients * delta * control$delta0
  # Original std.errors
  std.errors_2 <- std.errors_1 * delta * control$delta0
  # Store original coefficients
  sm$extra_params[["coefficients_orig"]] <- sm$coefficients
  sm$extra_params[["std.errors_orig"]] <- sm$std.errors
  # Store delta parameter
  sm$extra_params[["delta"]] <- delta * control$delta0
  # Update model
  sm$update(c_hat_2)
  sm$update_std.errors(std.errors_2)
  return(sm)
}

#' Impute clear sky outliers
#'
#' Detect and impute outliers with respect to a maximum level of radiation (Ct)
#'
#' @param x Numeric vector, time series of solar radiation.
#' @param Ct Numeric vector, time series of fitted clear sky radiation.
#' @param date Character or Date vector, time series of dates used to precisely impute solar radiation according to the realized values in the same day of the year.
#' @param threshold Numeric, scalar, threshold value used for imputation. Default is `0.0001`.
#' @param quiet Logical.
#' @examples
#' clearsky_outliers(c(1,2,3), 2)
#' @details
#' The function will detect the observations for which `x > Ct`, `x < 0` or `is.na(x)`. Then, if
#' \describe{
#' \item{`x < 0`}{If a value is below 0 for a day it will be imputed to be equal to `min(x)` for that day}.
#' \item{`x > Ct`}{If a value is above the maximum clear sky Ct it will be imputed to be `Ct*(1-threshold)`}.
#' \item{`is.na(x)`}{If a value is NA it will be imputed to be the average `mean(x)` for that day.}.
#' }
#'
#' @rdname clearsky_outliers
#' @name clearsky_outliers
#' @keywords clearsky
#' @note Version 1.0.0
#' @export
clearsky_outliers <- function(x, Ct, date, threshold = 0.0001, quiet = FALSE){
  # Initialize a dataset
  data <- dplyr::tibble(Ct = Ct, x = x)
  # Eventually add a date for non-stationary data
  if (!missing(date)){
    data <- dplyr::mutate(data,
                          date = as.Date(date),
                          Month = lubridate::month(date),
                          Day = lubridate::day(date))
  }
  # Number of observations
  nobs <- nrow(data)
  # Detect problems and violations
  outliers_na <- which(is.na(data$x))
  outliers_lo <- which(data$x <= 0 & !is.na(data$x))
  outliers_hi <- which(data$x >= data$Ct & !is.na(data$x))
  # Complete outliers index
  idx_outliers <- unique(c(outliers_na, outliers_lo, outliers_hi))
  # Check presence of outliers and impute them
  if (purrr::is_empty(idx_outliers)) {
    if (!quiet) message("No outliers!")
    data_clean <- data
  } else {
    # Verbose message
    if (!quiet) message("Outliers: ", length(idx_outliers), " (", format(length(idx_outliers)/nobs*100, digits = 2), " %)")
    # Dataset without outliers
    data_no_outliers <- data[-c(idx_outliers),]
    # Imputed dataset
    data_clean <- data
    # Impute outliers
    for (i in idx_outliers) {
      df_n <- data[i,]
      if (!missing(date)) {
        # Data for the same day and month (without outliers)
        df_day <- dplyr::filter(data_no_outliers, Month == df_n$Month & Day == df_n$Day & date != df_n$date)
      } else {
        df_day <- data_no_outliers
      }
      # Imputed data depending on outliers type
      if (i %in% outliers_na) {
        # Outlier is an NA
        data_clean[i,]$x <- mean(df_day$x, na.rm = TRUE)
      } else if (i %in% outliers_lo) {
        # Outlier is under minimum value for the day
        data_clean[i,]$x <- min(df_day$x, na.rm = TRUE)
      } else if (i %in% outliers_hi) {
        # Outlier is above maximum value for the day
        data_clean[i,]$x <- data_clean$Ct[i] * (1-threshold)
      }
    }
  }

  if (missing(date)) {
    date <- NA
  } else {
    date <- data$date[idx_outliers]
  }

  # Structure output data
  structure(
    list(
      # Complete imputed series
      x = data_clean$x,
      # Values before imputation
      original = data$x[idx_outliers],
      # Values after imputation
      imputed = data_clean$x[idx_outliers],
      # Index of imputed observations
      index = idx_outliers,
      # Index of type of outlier
      index_type = list(na = outliers_na, lo = outliers_lo, hi = outliers_hi),
      # Dates of imputed observations
      date = date,
      # Number of outliers
      n = length(idx_outliers),
      # Error of imputed series with respect to original data
      MAPE = mean(abs((data$x - data_clean$x)/data$x))*100,
      MSE = sd(data$x - data_clean$x),
      # Threshold for imputation
      threshold = threshold
    )
  )
}


