#' Control parameters for a `seasonalClearsky` object
#'
#' @param method Character, method used for estimate clear sky radiation. Can be `I` or `II`.
#' @param include.intercept Logical, when `TRUE`, the default, the intercept will be included in the clear sky model.
#' @param order Integer scalar, number of combinations of sines and cosines.
#' @param period Integer scalar, seasonality period. The default is 365.
#' @param delta0 Numeric scalar, initial value for delta.
#' @param quiet Logical, when `FALSE`, the default, the functions displays warning or messages.
#' @inheritParams clearsky_optimizer
#' @details The parameters `ntol`, `lower`, `upper` and `by` are used exclusively in \code{\link{clearsky_optimizer}}.
#' @examples
#' control = control_seasonalClearsky()
#' @return Named list of control parameters.
#'
#' @rdname control_seasonalClearsky
#' @export
control_seasonalClearsky <- function(method = "II", include.intercept = TRUE, order = 1, period = 365, delta0 = 1.4,
                                     lower = 0, upper = 3, by = 0.001, ntol = 30, quiet = FALSE){
  structure(
    list(
      method = method,
      include.intercept = include.intercept,
      order = order,
      period = period,
      delta0 = delta0,
      ntol = ntol,
      lower = lower,
      upper = upper,
      by = by,
      quiet = quiet
    ),
    class = c("control", "list")
  )
}

#' Clear sky seasonal model
#'
#' @examples
#' library(ggplot2)
#' # Arguments
#' place <- "Palermo"
#' # solarModel specification
#' spec <- solarModel_spec(place, target = "GHI")
#' # Extract the required elements
#' x <- spec$data$GHI
#' date <- spec$data$date
#' lat <- spec$coords$lat
#' clearsky <- spec$data$clearsky
#'
#' # Initialize the model
#' model <- seasonalClearsky$new()
#' # Fit the model
#' model$fit(x, date, lat, clearsky)
#' # Predict the seasonal values
#' spec$data$Ct <- model$predict(spec$data$n)
#'
#' @rdname seasonalClearsky
#' @name seasonalClearsky
#' @export
seasonalClearsky <- R6::R6Class("seasonalClearsky",
                                inherit = seasonalModel,
                                public = list(
                                  #' @field control Named list of control parameters. See the function \code{\link{control_seasonalClearsky}} for details.
                                  control = list(),
                                  #' @field lat latitude of the place considered.
                                  lat = NA_integer_,
                                  #' @method initialize seasonalClearsky
                                  #' @description
                                  #' Initialize a `seasonalClearsky` model
                                  #' @param control Named list of control parameters. See the function \code{\link{control_seasonalClearsky}} for details.
                                  initialize = function(control = control_seasonalClearsky()){
                                    self$control <- control
                                    private$..order <- control$order
                                    private$..period <- control$period
                                  },
                                  #' @method fit seasonalClearsky
                                  #' @description
                                  #' Fit a seasonal model for clear sky radiation.
                                  #' @param x time series of solar radiation.
                                  #' @param date time series of dates.
                                  #' @param lat Numeric, reference latitude.
                                  #' @param clearsky Numeric, optional time series of observed clerasky radiation.
                                  fit = function(x, date, lat, clearsky){
                                    # Self arguments
                                    control = self$control
                                    # Control parameters
                                    include.intercept = control$include.intercept
                                    order = control$order
                                    period = control$period
                                    delta0 = control$delta0
                                    ntol = control$ntol
                                    method = control$method
                                    # Initialize the dataset
                                    data <- dplyr::tibble(date = as.Date(date), n = number_of_day(date),
                                                          Month = lubridate::month(date), Day = lubridate::day(date))
                                    # Add target variable
                                    data$x <- x
                                    # Add extraterrestrial radiation
                                    data$H0 <- seasonalSolarFunctions$new("spencer")$H0(data$n, lat)
                                    # Add clearsky
                                    if (method == "II") {
                                      if (missing(clearsky)) {
                                        stop('`clearsky` time series must be specified when `method = "II"`')
                                      } else {
                                        data$clearsky <- clearsky
                                      }
                                    }
                                    # Fit dataset
                                    df_fit <- data
                                    # Method I: Estimate directly clear sky as function of extraterrestrial radiation
                                    # 1. LM: Ct ~ a_0 + a_1 H0 + a_2 cos(.) + a_3 sin(.) + a_4 cos(2*.) + a_5 sin(2*.) + ...
                                    # 2. Optimization: delta0*Ct ~ delta*(a_0 + a_1 H0 + a_2 cos(.) + a_3 sin(.) + a_4 cos(2*.) + a_5 sin(2*.) + ...)
                                    if (method == "I") {
                                      formula <- paste0("x ~ H0", ifelse(include.intercept, " + 1", " - 1"))
                                      # Fit the coefficients of the clear sky max model
                                      super$fit(formula = formula, data = df_fit)
                                      # Method II: Estimate with Extraterrestrial and clear sky radiation
                                      # 1. LM: Ct_max ~ a_0 + a_1 H0 + a_2 cos(.) + a_3 sin(.) + a_4 cos(2*.) + a_5 sin(2*.) + ...
                                      # 2. Optimization: delta_init*Ct_max ~ delta*(a_0 + a_1 H0 + a_2 cos(.) + a_3 sin(.) + a_4 cos(2*.) + a_5 sin(2*.) + ...)
                                    } else {
                                      # Compute maximum clear sky for each day
                                      df_fit <- dplyr::group_by(df_fit, n)
                                      df_fit <- dplyr::summarise(df_fit, Ct_max = max(clearsky, na.rm = TRUE), H0 = mean(H0))
                                      formula <- paste0("Ct_max ~ H0", ifelse(include.intercept, " + 1", " - 1"))
                                      # Fit the coefficients of the clear sky max model
                                      super$fit(formula = formula, data = df_fit)
                                      # Fitted clear sky max
                                      df_fit$Ct_max_fit <- self$predict(df_fit$n)
                                      # Optimize the fit
                                      df_fit <- dplyr::select(df_fit, n, Ct_max, Ct_max_fit)
                                      df_fit <- dplyr::left_join(data, df_fit, by = c("n"))
                                    }
                                    # Fitted clear sky max
                                    df_fit$Ct <- self$predict(df_fit$n)
                                    # Optimize the fit
                                    delta <- clearsky_optimizer(df_fit$x, df_fit$Ct*delta0, control$lower, control$upper, control$by, control$ntol)

                                    # Control parameters
                                    self$control <- control
                                    # Latitude
                                    self$lat <- lat
                                    # Original coefficients
                                    private$coefficients_orig <- self$coefficients
                                    # Delta parameter
                                    private$delta <- delta*delta0

                                    # Standard names for coefficients
                                    coefs_names <- c()
                                    params <- self$coefficients*delta*delta0
                                    if (include.intercept) {
                                      coefs_names[1] <- "delta_0"
                                      coefs_names[2] <- "delta_extra"
                                    } else {
                                      coefs_names[1] <- "delta_extra"
                                    }
                                    if (order > 0) {
                                      base_names <- paste0("delta_", rep(c("sin_", "cos_")))
                                      coefs_names <- c(coefs_names, unlist(purrr::map(1:order, ~paste0(base_names, .x))))
                                    }
                                    names(params) <- coefs_names
                                    # Update coefficients
                                    self$update(params)
                                    },
                                  #' @method updateH0 seasonalClearsky
                                  #' @description
                                  #' Update the time series of extraterrestrial radiation for a given latitude.
                                  #' @param lat reference latitude
                                  updateH0 = function(lat){
                                    self$lat <- lat
                                    self$seasonal_data$H0 <- seasonalSolarFunctions$new("spencer")$H0(self$seasonal_data$n, self$lat)
                                  },
                                  #' @description
                                  #' Print method for the class `seasonalClearsky`
                                  print = function(){
                                    cat(paste0("----------------------- seasonalClearsky ----------------------- \n"))
                                    msg_1 <- paste0(" - Order: ", self$order, "\n - Period: ", self$period, "\n")
                                    n_external_reg <- ncol(self$seasonal_data)-1
                                    if (n_external_reg == 0) {
                                      msg_2 <- paste0("- External regressors: ", n_external_reg, "\n")
                                    } else {
                                      msg_2 <- paste0("- External regressors: ", n_external_reg, " (", names(self$seasonal_data)[-1], ")\n")
                                    }
                                    cat(c(msg_1, msg_2))
                                    cat(paste0("--------------------------------------------------------------\n"))
                                    print(self$model)
                                  }
                                ),
                                private = list(
                                  coefficients_orig = NA,
                                  delta = NA
                                )
                              )


#' Optimizer for Solar Clear sky
#'
#' Find the best parameter delta for fitting clear sky radiation.
#'
#' @param x Numeric vector, realized solar radiation.
#' @param Ct Numeric vector, clear sky radiation.
#' @param ntol Integer, tolerance for `clearsky > GHI` condition. Maximum number of violations admitted.
#' @param lower Numeric scalar, lower bound for delta grid.
#' @param upper Numeric scalar, upper bound for delta grid.
#' @param by Numeric scalar, step for delta grid.
#'
#' @return Numeric, scalar the optimal delta parameter.
#'
#' @name clearsky_optimizer
#' @rdname clearsky_optimizer
#' @keywords internal
#' @export
clearsky_optimizer <- function(x, Ct, lower = 0, upper = 3, by = 0.01, ntol = 30){
  # Grid of points
  grid <- seq(lower, upper, by = by)
  # Loss
  opt <- data.frame(delta = grid, loss = purrr::map_dbl(grid, ~sum(.x*Ct - x < 0)))
  # Return minimum delta satisfying the constraint
  delta <- opt[which(opt$loss <= ntol)[1],]$delta
  return(delta)
}

#' Impute clear sky outliers
#'
#' Detect and impute outliers with respect to a maximum level of radiation (Ct)
#'
#' @param x Numeric vector, realized solar radiation.
#' @param Ct Numeric vector, clear sky radiation.
#' @param date Vector, optional time series of dates.
#' It will be used for a more precise imputation when a solar radiation value is `NA`.
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
    data <- dplyr::mutate(data, date = as.Date(date),
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
      if (!missing(date)){
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
        # Outlier is under minum value for the day
        data_clean[i,]$x <- min(df_day$x, na.rm = TRUE)
      } else {
        # Outlier is above maximum value for the day
        data_clean[i,]$x <- data_clean$Ct[i]*(1-threshold)
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


