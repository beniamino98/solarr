#' Control parameters for a `clearskySeasonal` object
#'
#' @param method character, method for clearsky estimate, can be `I` or `II`.
#' @param include.intercept logical. When `TRUE`, the default, the intercept will be included in the model.
#' @param order numeric, of sine and cosine elements.
#' @param period numeric, periodicity. The default is `365`.
#' @param seed numeric, random seed for reproducibility. It is used to random impute the violations.
#' @param delta_init Value for delta init in the clear sky model.
#' @param quiet logical. When `FALSE`, the default, the functions displays warning or messages.
#' @inheritParams clearsky_optimize
#'
#' @details The parametes `tol`, `lower`, `upper` and `by` are used exclusively in \code{\link{clearskyModel_optimize}}.
#'
#' @rdname control_clearskyModel
#' @export
control_clearskyModel <- function(method = "II", include.intercept = TRUE, order = 1, period = 365, seed = 1,
                                  delta_init = 1.4, ntol = 30, lower = 0, upper = 2, by = 0.001, quiet = FALSE){
  structure(
    list(
      method = method,
      include.intercept = include.intercept,
      order = order,
      period = period,
      seed = seed,
      delta_init = delta_init,
      ntol = ntol,
      delta0 = c(lower = lower, upper = upper, by = by),
      quiet = quiet
    ),
    class = c("control", "list")
  )
}

#' Optimizer for Solar Clear sky
#'
#' Find the best parameter delta for fitting clear sky radiation.
#'
#' @param x vector of solar radiation
#' @param Ct vector of extraterrestrial radiation
#' @param ntol integer, tolerance for `clearsky > GHI` condition. Maximum number of violations admitted.
#' @param lower numeric, lower bound for delta grid.
#' @param upper numeric, upper bound for delta grid.
#' @param by numeric, step for delta grid.
#'
#' @return the optimal delta
#'
#' @name clearsky_optimize
#' @rdname clearsky_optimize
#' @export
clearsky_optimize <- function(x, Ct, lower = 0, upper = 3, by = 0.01, ntol = 30){

  # Grid of parameters
  grid <- seq(lower, upper, by = by)
  # Loss
  opt <- data.frame(delta = grid, loss = purrr::map_dbl(grid, ~sum(.x*Ct - x < 0)))
  # Optimal parameter
  delta <- opt[which(opt$loss <= ntol)[1],]$delta
  return(delta)
}

#' Detect outliers with respect to clearsky radiation
#'
#' @param x time series of solar radiation
#' @param Ct time series of clearsky radiation
#' @param date time series of dates
#'
#' @rdname clearsky_outliers
#' @name clearsky_outliers
#' @export
clearsky_outliers <- function(x, Ct, date, quiet = FALSE){

  # Initialize a dataset
  data <- dplyr::tibble(date = as.Date(date),
                        Month = lubridate::month(date),
                        Day = lubridate::day(date),
                        Ct = Ct, x = x)
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
      # Dataset without outlier
      df_day <- dplyr::filter(data_no_outliers, Month == df_n$Month & Day == df_n$Day & date != df_n$date)
      if (i %in% outliers_na) {
        # Imputed data
        data_clean[i,]$x <- mean(df_day$x)
      } else if (i %in% outliers_lo) {
        # Imputed data
        data_clean[i,]$x <- min(df_day$x)*0.999
      } else {
        # Imputed data
        data_clean[i,]$x <- min(c((max(df_day$x) + df_n$Ct)/2, df_n$Ct*0.999))
      }
    }
  }

  structure(
    list(
      x = data_clean$x,
      # Original series
      original = data$x[idx_outliers],
      # Imputed series
      imputed = data_clean$x[idx_outliers],
      # Additional information
      index = idx_outliers,
      date = data$date[idx_outliers],
      n = length(idx_outliers)
    )
  )
}


#' Seasonal model for clear sky radiation
#'
#' Fit a seasonal model for clear sky radiation in a location.
#'
#' @param data dataset
#' @param seasonal_data dataset with two columns: `n` with the number of the day in 1:365
#' and `H0` with the extraterrestrial radiation.
#' @param method character, method for clearsky estimate, can be `I` or `II`.
#' @param control list of control parameters. See \code{\link{control_clearskyModel}} for details.
#'
#' @rdname clearskySeasonal
#' @name clearskySeasonal
#' @export
clearskySeasonal <- function(data, seasonal_data, method = "II", control = control_clearskyModel()){

  # Control parameters
  include.intercept = control$include.intercept
  order = control$order
  period = control$period
  delta_init = control$delta_init
  ntol = control$ntol
  delta0 = control$delta0
  # data <- dplyr::left_join(data, seasonal_data, by = c("n"))
  # Method I: Estimate clear sky as function of Extraterrestrial radiation
  # 1. LM: Ct ~ a_0 + a_1 H0 + a_2 cos(.) + a_3 sin(.) + a_4 cos(2*.) + a_5 sin(2*.) + ...
  # 2. Optimization: delta_init*Ct ~ delta*(a_0 + a_1 H0 + a_2 cos(.) + a_3 sin(.) + a_4 cos(2*.) + a_5 sin(2*.) + ...)
  df_fit <- data
  if (method == "I") {
    formula <- paste0("GHI ~ H0", ifelse(include.intercept, " + 1", " - 1"))
    # Fit the coefficients of the clear sky max model
    seasonal_model_Ct <- seasonalModel$new(formula = formula, order = order, period = period, data = df_fit)
    # Fitted clear sky max
    df_fit$Ct <- seasonal_model_Ct$predict()
    # Optimize the fit
    delta <- clearsky_optimize(df_fit$GHI, df_fit$Ct*delta_init, lower = delta0[1], upper = delta0[2], by = delta0[3], ntol = ntol)
    # Store old coefficients inside the `params` slot
    seasonal_model_Ct$.__enclos_env__$private$params$coefficients_orig <- seasonal_model_Ct$coefficients
    seasonal_model_Ct$.__enclos_env__$private$params$delta <- delta*delta_init
    # Update coefficients
    seasonal_model_Ct$update(seasonal_model_Ct$coefficients*delta*delta_init)
    # Method II: Estimate with Extraterrestrial and clear sky radiation
    # 1. LM: Ct_max ~ a_0 + a_1 H0 + a_2 cos(.) + a_3 sin(.) + a_4 cos(2*.) + a_5 sin(2*.) + ...
    # 2. Optimization: delta_init*Ct_max ~ delta*(a_0 + a_1 H0 + a_2 cos(.) + a_3 sin(.) + a_4 cos(2*.) + a_5 sin(2*.) + ...)
  } else if (method == "II") {
    # Compute maximum clear sky for each day
    df_fit <- dplyr::group_by(data, Month, Day)
    df_fit <- dplyr::reframe(df_fit, Ct_max = max(clearsky), H0 = mean(H0))
    df_fit$n <- 1:nrow(df_fit)
    formula <- paste0("Ct_max ~ H0", ifelse(include.intercept, " + 1", " - 1"))
    # Fit the coefficients of the clear sky max model
    seasonal_model_Ct <- seasonalModel$new(formula = formula, order = order, period = period, data = df_fit)
    # Fitted clear sky max
    df_fit$Ct_max_fit <- seasonal_model_Ct$predict()
    # Optimize the fit
    df_fit <- dplyr::select(df_fit, Month, Day, Ct_max, Ct_max_fit)
    df_fit <- dplyr::left_join(data, df_fit, by = c("Month", "Day"))
    # Store extraterrestrial for predictions
    seasonal_model_Ct$.__enclos_env__$private$seasonal_data <- seasonal_data
    df_fit$Ct <- seasonal_model_Ct$predict(df_fit$n)
    # Optimize the fit
    delta <- clearsky_optimize(df_fit$GHI, df_fit$Ct*delta_init, lower = delta0[1], upper = delta0[2], by = delta0[3], ntol = ntol)
    # Store old coefficients inside the `lm` object
    seasonal_model_Ct$.__enclos_env__$private$params$coefficients_orig <- seasonal_model_Ct$coefficients
    seasonal_model_Ct$.__enclos_env__$private$params$delta <- delta*delta_init
    # Update coefficients
    seasonal_model_Ct$update(seasonal_model_Ct$coefficients*delta*delta_init)
  }
  # Store period, extraterrestrial and  control inside the `lm` object
  seasonal_model_Ct$.__enclos_env__$private$params$period <- control$period
  seasonal_model_Ct$.__enclos_env__$private$params$control <- control
  seasonal_model_Ct$.__enclos_env__$private$seasonal_data <- seasonal_data
  # Class
  class(seasonal_model_Ct) <- c("clearskySeasonal", class(seasonal_model_Ct))
  return(seasonal_model_Ct)
}



