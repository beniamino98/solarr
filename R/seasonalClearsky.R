#' Control parameters for a `seasonalClearsky` object
#'
#' @param include.intercept Logical, when `TRUE`, the default, the intercept will be included in the clear sky model.
#' @param order Integer scalar, number of combinations of sines and cosines.
#' @param period Integer scalar, seasonality period. The default is 365.
#' @param delta0 Numeric scalar, initial value for delta.
#' @param quiet Logical, when `FALSE`, the default, the functions displays warning or messages.
#' @inheritParams clearsky_optimizer
#' @details The parameters `ntol`, `lower`, `upper` and `by` are used exclusively in \code{\link{clearsky_optimizer}}.
#' @note Version 1.0.0
#' @examples
#' control = control_seasonalClearsky()
#' @return Named list of control parameters.
#' @note Version 1.0.0
#' @rdname control_seasonalClearsky
#' @export
control_seasonalClearsky <- function(order = 1, order_H0 = 1, period = 365, include.intercept = TRUE, include.trend = FALSE,
                                     delta0 = 1.4, lower = 0, upper = 3, by = 0.001, ntol = 0, quiet = FALSE){
  structure(
    list(
      order = order,
      order_H0 = order_H0,
      period = period,
      include.intercept = include.intercept,
      include.trend = include.trend,
      delta0 = delta0,
      lower = lower,
      upper = upper,
      by = by,
      ntol = ntol,
      quiet = quiet
    ),
    class = c("control", "list")
  )
}

#' R6 implementation for a clear sky seasonal model
#'
#' @rdname seasonalClearsky
#' @name seasonalClearsky
#' @export
seasonalClearsky <- R6::R6Class("seasonalClearsky",
                                inherit = seasonalModel,
                                # ====================================================================================================== #
                                #                                             Public slots
                                # ====================================================================================================== #
                                public = list(
                                  #' @field control Named list. Control parameters. See the function \code{\link{control_seasonalClearsky}} for details.
                                  control = list(),
                                  #' @field lat Numeric, scalar. Latitude of the location considered.
                                  lat = NA_integer_,
                                  #' @description
                                  #' Initialize a `seasonalClearsky` object.
                                  #' @param control Named list. Control parameters. See the function \code{\link{control_seasonalClearsky}} for details.
                                  initialize = function(control = control_seasonalClearsky()){
                                    # Store control parameters
                                    self$control <- control
                                    # Update order
                                    private$..order <- control$order
                                    # Update period
                                    private$..period <- control$period
                                  },
                                  #' @description
                                  #' Fit the seasonal model for clear sky radiation.
                                  #' @param x Numeric vector. Time series of solar radiation.
                                  #' @param date Character or Date vector. Time series of dates.
                                  #' @param lat Numeric scalar. Reference latitude.
                                  #' @param clearsky Numeric vector. Time series of CAMS clear sky radiation.
                                  fit = function(x, date, lat, clearsky){
                                    # Self arguments
                                    control = self$control
                                    # Control parameters
                                    include.intercept = control$include.intercept
                                    include.trend = control$include.trend
                                    # Ensure clearsky is specified
                                    if (missing(clearsky)) {
                                      stop('`clearsky` time series must be specified.')
                                    }
                                    # Add the function to compute extraterrestrial radiation
                                    private$..ssf <- seasonalSolarFunctions$new("spencer")
                                    # Store reference latitude
                                    self$lat <- lat[1]
                                    # Initialize the dataset
                                    data <- dplyr::tibble(date = as.Date(date))
                                    data <- dplyr::mutate(data,
                                                          Year = lubridate::year(date),
                                                          Month = lubridate::month(date),
                                                          Day = lubridate::day(date),
                                                          t = Year - max(Year),
                                                          n = number_of_day(date),
                                                          Rt = x,
                                                          H0 = self$H0(n),
                                                          clearsky = clearsky)
                                    # Method: Estimate with Extraterrestrial and clear sky radiation
                                    # ========================================================================
                                    # 1. Daily maximum clearsky: Ct_max ~ a_0 + a_1 H0 + a_2 cos(.) + a_3 sin(.) + a_4 cos(2*.) + a_5 sin(2*.) + ...
                                    # ========================================================================
                                    # Compute maximum clear sky for each day
                                    base_formula <- "clearsky ~ H0"
                                    if (control$order_H0 > 1){
                                      for(i in 2:control$order_H0){
                                        data[[paste0("H0_", i)]] <- data$H0^i
                                        base_formula <- paste0(base_formula, " + ", paste0("H0_", i))
                                      }
                                    }
                                    base_formula <- ifelse(include.trend, paste0(base_formula, " + t"), base_formula)
                                    base_formula <- ifelse(include.intercept, base_formula, paste0(base_formula, "-1"))
                                    # Fit the coefficients of the clear sky max model
                                    super$fit(formula = base_formula, data = data)
                                    # Initial fit average clear sky
                                    data$Ct_hat <- self$predict(newdata = data)
                                    # Optimize the fit
                                    data <- dplyr::select(data, n, t, H0, Rt, Ct_hat)
                                    # ========================================================================
                                    # 2. Optimization: delta_init*Ct_max ~ delta*(a_0 + a_1 H0 + a_2 cos(.) + a_3 sin(.) + a_4 cos(2*.) + a_5 sin(2*.) + ...)
                                    # ========================================================================
                                    # Optimize the fit
                                    delta <- clearsky_delta_optimizer(data$Rt, data$Ct_hat*control$delta0, control$lower, control$upper, control$by, control$ntol)
                                    # Standard names for coefficients
                                    coefs_names <- c()
                                    orig_names <- super$.__enclos_env__$private$..model$coefficients_names
                                    if (include.intercept) {
                                      coefs_names[1] <- "delta_0"
                                      orig_names <- orig_names[-c(1)]
                                      for(i in 1:control$order_H0){
                                        coefs_names[i+1] <- paste0("delta_extra", i)
                                        orig_names <- orig_names[-c(1)]
                                      }
                                    } else {
                                      for(i in 1:control$order_H0){
                                        coefs_names[i] <- paste0("delta_extra", i)
                                        orig_names <- orig_names[-c(1)]
                                      }
                                    }
                                    if (include.trend) {
                                      coefs_names <- c(coefs_names, "t")
                                      orig_names <- orig_names[-c(1)]
                                    }
                                    if (self$order > 0) {
                                      coefs_names <- c(coefs_names, paste0("delta_", orig_names))
                                    }
                                    # Store original coefficients
                                    private$coefficients_orig <- self$coefficients
                                    # Store delta parameter
                                    private$delta <- delta * control$delta0
                                    # Update coefficients values and names
                                    super$.__enclos_env__$private$..model$coefficients <- super$.__enclos_env__$private$..model$coefficients * private$delta
                                    super$.__enclos_env__$private$..model$coefficients_names <- coefs_names
                                    # Update std errors values and names
                                    super$.__enclos_env__$private$..std.errors <- super$.__enclos_env__$private$..std.errors * private$delta
                                    names(super$.__enclos_env__$private$..std.errors) <- coefs_names
                                  },
                                  #' @description
                                  #' Compute the extraterrestrial radiation at a given location.
                                  #' @param n Integer, scalar or vector. Number of day of the year.
                                  H0 = function(n){
                                    private$..ssf$H0(n, self$lat)
                                  },
                                  #' @description
                                  #' Predict method for `seasonalClearsky` object.
                                  #' @param n Integer, scalar or vector. number of day of the year.
                                  predict = function(n, newdata){
                                    if (missing(newdata)) {
                                      if (missing(n)) {
                                        predict.lm(private$..model)
                                      } else {
                                        H0 <- self$H0(n)
                                        newdata <- data.frame(n = n, H0 = H0)
                                        if (self$control$order_H0 > 1){
                                          for(i in 2:self$control$order_H0){
                                            newdata[[paste0("H0_", i)]] <- newdata$H0^i
                                          }
                                        }
                                        predict.lm(private$..model, newdata = newdata)
                                      }
                                    } else {
                                      newdata$H0 <- self$H0(newdata$n)
                                      if (self$control$order_H0 > 1){
                                        for(i in 2:self$control$order_H0){
                                          newdata[[paste0("H0_", i)]] <- newdata$H0^i
                                        }
                                      }
                                      predict.lm(private$..model, newdata = newdata)
                                    }
                                  },
                                  #' @description
                                  #' Print method for `seasonalClearsky` object.
                                  print = function(){
                                    cat(paste0("----------------------- seasonalClearsky ----------------------- \n"))
                                    msg_1 <- paste0(" - Order: ", self$order, "\n - Period: ", self$period, "\n")
                                    msg_2 <- paste0("- External regressors: 1 (H0) \n")
                                    msg_3 <- paste0("- Version: ", private$version, "\n")
                                    cat(c(msg_1, msg_2, msg_3))
                                    cat(paste0("--------------------------------------------------------------\n"))
                                    print(self$model)
                                  }
                                ),
                                private = list(
                                  version = "1.0.0",
                                  coefficients_orig = NA,
                                  delta = NA,
                                  ..ssf = NA
                                )
                              )

#' Optimizer for Solar Clear sky
#'
#' Find the best parameter delta for fitting clear sky radiation.
#'
#' @param x Numeric vector. Time series of solar radiation.
#' @param Ct Numeric vector. Time series of  clear sky radiation.
#' @param ntol Integer scalar. Tolerance for the maximum number of violations admitted of the condition `clearsky > GHI`. Default is `0`.
#' @param lower Numeric scalar. Lower bound for grid of delta parameters used for optimization. Default is `0`.
#' @param upper Numeric scalar. Upper bound for grid of delta parameters used for optimization. Default is `3`.
#' @param by Numeric scalar, step for the grid. Default is `0.01`.
#'
#' @return Numeric, scalar the optimal delta parameter.
#'
#' @name clearsky_delta_optimizer
#' @rdname clearsky_delta_optimizer
#' @keywords internal
#' @export
clearsky_delta_optimizer <- function(x, Ct, lower = 0, upper = 3, by = 0.01, ntol = 0){
  # Grid of points
  grid <- seq(lower, upper, by = by)
  # Loss
  opt <- data.frame(delta = grid, loss = purrr::map_dbl(grid, ~sum(.x*Ct - x < 0)))
  # Return minimum delta satisfying the constraint
  delta <- opt[which(opt$loss <= ntol)[1],]$delta
  return(delta)
}
#' Optimizer for Solar Clear sky
#'
#' @name clearsky_optimizer
#' @rdname clearsky_optimizer
#' @keywords internal
#' @export
clearsky_optimizer <- function(seasonal_model_Ct, data, ntol = 0){
  # Clone the model
  sm <- seasonal_model_Ct$clone(TRUE)
  # Loss function
  loss <- function(params, sm, data){
    sm$update(params)
    # Prediction
    pred <- sm$predict(newdata = data)
    # Violations
    violations <- ifelse(data$GHI - pred > 0, 1, 0)
    # Check number of violations lower than ntol
    mse <- sum((data$GHI - pred)^2) + 1000000*(sum(violations) - ntol)
    return(mse)
  }
  # Optimal parameters
  opt <- optim(sm$coefficients, loss, sm = sm, data = data)
  # Update the parameters
  sm$update(opt$par)
  return(sm)
}


#' Impute clear sky outliers
#'
#' Detect and impute outliers with respect to a maximum level of radiation (Ct)
#' \describe{
#' \item{`GHI < 0`}{If a value is below 0 for a day it will be imputed to be equal to min(GHI) for that day}.
#' \item{`GHI > Ct`}{If a value is above the maximum clear sky Ct it will be imputed to be Ct*(1-threshold)}.
#' \item{`is.na(GHI)`}{If a value is NA it will be imputed to be the average for that day mean(GHI)}.
#' }
#'
#' @param x Numeric vector. Time series of solar radiation.
#' @param Ct Numeric vector. Time series of clear sky radiation.
#' @param date Character or Date vector, optional. Time series of dates.
#' Used to precisely impute solar radiation according to the realized values in the same day of the year.
#' @param threshold Numeric, scalar. Threshold value used for imputation. Default is `0.0001`.
#' @examples
#' clearsky_outliers(c(1,2,3), 2)
#'
#' @rdname clearsky_outliers
#' @name clearsky_outliers
#' @keywords internal
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



