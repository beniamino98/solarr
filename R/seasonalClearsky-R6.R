#' Seasonal clear-sky model
#'
#' @description R6 implementation of a clear-sky seasonal model for solar radiation,
#' using extraterrestrial radiation and seasonal harmonics.
#'
#' @docType class
#' @rdname seasonalClearsky
#' @name seasonalClearsky
#' @keywords clearsky
#' @note Version 1.0.3
#' @export
seasonalClearsky <- R6::R6Class("seasonalClearsky",
                                inherit = seasonalModel,
                                public = list(
                                  #' @field lat Numeric, scalar, latitude of the location considered.
                                  lat = NA_integer_,
                                  #' @field spec List, model specification
                                  spec = seasonalClearsky_spec(),
                                  #' @description
                                  #' Initialize a `seasonalClearsky` object.
                                  #' @param spec Named list, model's specification. See the function \code{\link{seasonalClearsky_spec}} for more details.
                                  #' @param control Named list, control parameters. See the function \code{\link{control_seasonalClearsky}} for more details.
                                  initialize = function(spec = seasonalClearsky_spec(), control = control_seasonalClearsky()){
                                    # Store specification
                                    self$spec <- spec
                                    # Store control
                                    self$control <- control

                                    # Build the base formula
                                    base_formula <- "clearsky ~ H0"
                                    if (spec$order_H0 > 1){
                                      for(i in 2:spec$order_H0){
                                        base_formula <- paste0(base_formula, " + ", paste0("H0_", i))
                                      }
                                    }
                                    base_formula <- ifelse(spec$include.trend, paste0(base_formula, " + n"), base_formula)
                                    base_formula <- ifelse(spec$include.intercept, base_formula, paste0(base_formula, "-1"))
                                    # Add the function to compute extraterrestrial radiation
                                    private$..ssf <- seasonalSolarFunctions$new(spec$method_H0)
                                    # Initialization from `seasonalModel`
                                    super$initialize(base_formula, order = spec$order, period = spec$period)
                                    # ========================================================================
                                    # 2. Standardize parameters names
                                    # ========================================================================
                                    new_coefs_names <- self$coefs_names
                                    if (spec$include.intercept) {
                                      new_coefs_names["intercept"] <- "delta_0"
                                    }
                                    if (spec$include.trend) {
                                      new_coefs_names["n"] <- "delta_n"
                                    }

                                    new_coefs_names["H0"] <- "delta_extra_1"

                                    if (spec$order_H0 > 1){
                                      for(i in 2:spec$order_H0){
                                        new_coefs_names[paste0("H0_", i)] <- paste0("delta_extra", i)
                                      }
                                    }
                                    # Add delta on standard coefficients
                                    idx_sincos <- stringr::str_detect(new_coefs_names, "sin|cos")
                                    new_coefs_names[idx_sincos] <- paste0("delta_", new_coefs_names[idx_sincos])
                                    self$update_coefs_names(new_coefs_names)
                                  },
                                  #' @description
                                  #' Fit the seasonal model for clear sky radiation.
                                  #' @param x Numeric vector, time series of solar radiation.
                                  #' @param dates Character or Date vector, time series of dates.
                                  #' @param lat Numeric scalar, reference latitude.
                                  #' @param clearsky Numeric vector, time series of target clear sky radiation.
                                  fit = function(x, dates, lat, clearsky){
                                    # Self arguments
                                    spec = self$spec
                                    control = self$control
                                    # Ensure clearsky is specified
                                    if (missing(clearsky)) {
                                      stop('`clearsky` time series must be specified.')
                                    }
                                    # Store reference latitude
                                    self$lat <- lat[1]
                                    # Initialize the dataset
                                    data <- dplyr::tibble(date = as.Date(dates),
                                                          n = number_of_day(date),
                                                          Rt = x,
                                                          H0 = self$ssf$Hon(n, self$lat),
                                                          clearsky = clearsky)
                                    # Method: Estimate with Extraterrestrial and clear sky radiation
                                    # ========================================================================
                                    # 1. Daily maximum clearsky: Ct_max ~ a_0 + a_1 H0 + a_2 cos(.) + a_3 sin(.) + a_4 cos(2*.) + a_5 sin(2*.) + ...
                                    # ========================================================================
                                    # Add Extraterrestrial orders
                                    if (spec$order_H0 > 1){
                                      for(i in 2:spec$order_H0){
                                        data[[paste0("H0_", i)]] <- data$H0^i
                                      }
                                    }
                                    # Fit the coefficients of the clear sky max model
                                    super$fit(data = data)
                                    # Initial fit average clear sky
                                    data$Ct_hat <- super$predict(newdata = data)
                                    # 2. Optimize the fit
                                    # Optimization: delta_init*Ct_max ~ delta*(a_0 + a_1 H0 + a_2 cos(.) + a_3 sin(.) + a_4 cos(2*.) + a_5 sin(2*.) + ...)
                                    delta <- clearsky_delta(data$Rt, data$Ct_hat * control$delta0,
                                                            control$lower, control$upper, control$by, control$ntol)

                                    # ========================================================================
                                    # 2. Standardize parameters names
                                    # ========================================================================
                                    # Store original coefficients
                                    self$extra_params[["coefficients_orig"]] <- self$coefficients
                                    self$extra_params[["std.errors_orig"]] <- self$std.errors
                                    # Store delta parameter
                                    self$extra_params[["delta"]] <- delta * control$delta0
                                    # Update std errors
                                    std.errors <- self$std.errors * self$extra_params[["delta"]]
                                    # Update coefficients values and names
                                    self$update(self$extra_params[["coefficients_orig"]] * delta * control$delta0)
                                    # Update std errors values and names
                                    self$update_std.errors(std.errors)
                                  },
                                  #' @description
                                  #' Predict method for `seasonalClearsky` object.
                                  #' @param n Integer, scalar or vector. number of day of the year.
                                  #' @param newdata An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
                                  predict = function(n, newdata){
                                    if (!missing(n)) {
                                      H0 <- self$ssf$Hon(n, self$lat)
                                      newdata <- data.frame(n = n, H0 = H0)
                                    } else {
                                      newdata[["H0"]] <- self$ssf$Hon(newdata$n, self$lat)
                                    }
                                    # Add orders of H0
                                    if (self$spec$order_H0 > 1){
                                      for(i in 2:self$spec$order_H0){
                                        newdata[[paste0("H0_", i)]] <- newdata$H0^i
                                      }
                                    }
                                    super$predict(newdata = newdata)
                                  },
                                  #' @description
                                  #' Differential method for `seasonalClearsky` object.
                                  #' @param n Integer, scalar or vector. number of day of the year.
                                  #' @param newdata An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
                                  differential = function(n, newdata){
                                    if (!missing(n)) {
                                      H0 <- self$ssf$Hon(n, self$lat, deriv = FALSE)
                                      dH0 <- self$ssf$Hon(n, self$lat, deriv = TRUE)
                                      newdata <- data.frame(n = n, H0 = dH0)
                                    } else {
                                      H0 <- self$ssf$Hon(newdata$n, self$lat, deriv = FALSE)
                                      dH0 <- self$ssf$Hon(newdata$n, self$lat, deriv = TRUE)
                                      newdata[["H0"]] <- H0
                                    }

                                    # Add derivatives of the powerts of H0
                                    if (self$spec$order_H0 > 1){
                                      for(i in 2:self$spec$order_H0){
                                        newdata[[paste0("H0_", i)]] <- i * H0^(i-1) * dH0
                                      }
                                    }
                                    # Differential of the trend
                                    if (self$spec$include.trend){
                                      super$differential(newdata = newdata) - self$coefficients["delta_n"]*newdata[["n"]] + self$coefficients["delta_n"]
                                    } else {
                                      super$differential(newdata = newdata)
                                    }
                                  },
                                  #' @description
                                  #' Print method for `seasonalClearsky` object.
                                  print = function(){
                                    msg_0 <- "-------------------------- seasonalClearsky --------------------------"
                                    cat(paste0(msg_0, "\n"))
                                    msg_1 <- paste0("- Order: ", self$order, "\n",
                                                    "- Order (H0): ", self$spec$order_H0, "\n",
                                                    "- Period: ", self$period, "\n",
                                                    "- External regressors: 1 (H0)", "\n",
                                                    "- Linear trend: ", self$spec$include.trend, "\n",
                                                    "- Method: ", self$spec$method, "\n",
                                                    "- Version: ", private$version, "\n")
                                    msg_line <- paste0(rep("-", length(str_split(msg_0, "")[[1]])), collapse = "")
                                    cat(msg_1)
                                    cat(paste0(msg_line, "\n"))
                                    print(self$coefficients)
                                  }
                                ),
                                private = list(
                                  version = "1.0.3",
                                  ..ssf = NA
                                ),
                                active = list(
                                  #' @field ssf See \code{\link{seasonalSolarFunctions}} for more details.
                                  ssf = function(){
                                    private$..ssf
                                  }
                                )
                              )



