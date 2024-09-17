#' Seasonal model object
#'
#' @rdname seasonalModel
#' @name seasonalModel
#' @export
seasonalModel <- R6::R6Class("seasonalModel",
                             public = list(
                               #' @description
                               #' Fit a seasonal model as a linear combination of sine and cosine functions.
                               #'
                               #' @param formula formula, an object of class `formula` (or one that can be coerced to that class).
                               #' It is a symbolic description of the model to be fitted and can be used to include or exclude the intercept or external regressors in `data`.
                               #' @param order numeric, of sine and cosine elements.
                               #' @param period numeric, periodicity. The default is `2*base::pi/365`.
                               #' @param data 	an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
                               #' If not found in data, the variables are taken from environment(formula), typically the environment from which `lm` is called.
                               #' @param ... other parameters to be passed to the function lm.
                               initialize = function(formula = "Yt ~ 1", order = 1, period = 365, data, ...){
                                 # Model formula
                                 if (order > 0) {
                                   for (i in 1:order){
                                     formula <- paste0(formula, " + ", "cos((2*base::pi)/", eval(period), "*n*", i, ")", " + ",
                                                       "sin((2*base::pi)/", eval(period), "*n*", i, ")")
                                   }
                                 }
                                 # Fit seasonal model
                                 private$model <- lm(as.formula(formula), data = data, ...)
                                 # Store model, period, order and formula
                                 private$period = period
                                 private$order = order
                                 private$formula = as.formula(formula)
                               },
                               #' @description
                               #' Convert a date into the respective number of day of the year.
                               #' @param x character or Date object.
                               from_date_to_n = function(x){
                                 purrr::map_int(x, ~ifelse(is.character(.x) | lubridate::is.Date(.x), number_of_day(.x), .x))
                               },
                               #' @description
                               #' Predict method for a `seasonalModel`
                               #' @param n integer, number of day of the year.
                               predict = function(n){
                                 if (missing(n)){
                                   predict.lm(private$model)
                                 } else {
                                   n <- self$from_date_to_n(n)
                                   newdata <- data.frame(n = n %% private$period)
                                   newdata$n <- ifelse(newdata$n == 0, 1, newdata$n)
                                   newdata <- dplyr::left_join(newdata, private$seasonal_data, by = "n")
                                   predict.lm(private$model, newdata = newdata)
                                 }
                               },
                               #' @description
                               #' Update the parameters of a `seasonalModel`.
                               #' @param n integer, number of day of the year.
                               update = function(coefficients){
                                 old_coef <- private$model$coefficients
                                 # Check lenght
                                 if (length(old_coef) != length(coefficients)){
                                   warning("The lenght of `coefficients` do not match the length of the old coefficients.")
                                   return(invisible(NULL))
                                 }
                                 # Update parameters
                                 private$model$coefficients <- coefficients
                               }
                             ),
                             private = list(
                               model = NA,
                               period = NA,
                               order = NA,
                               formula = NA,
                               params = list(),
                               seasonal_data = dplyr::tibble(n = 1:366)
                             ),
                             active = list(
                               coefficients = function(){
                                 private$model$coefficients
                               }
                             )

)

