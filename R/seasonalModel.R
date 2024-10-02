#' Seasonal Model Object
#'
#' @description
#' The `seasonalModel` class implements a seasonal regression model as a linear combination of sine and cosine functions.
#' This model is designed to capture periodic effects in time series data, particularly for applications involving seasonal trends.
#'
#' @details
#' The seasonal model is fitted using a specified formula, which allows for the inclusion of external regressors along with sine and cosine terms
#' to model seasonal variations. The periodicity can be customized, and the model can be updated with new coefficients after fitting.
#'
#' @rdname seasonalModel
#' @name seasonalModel
#' @export
seasonalModel <- R6::R6Class("seasonalModel",
                             public = list(
                               #' @field seasonal_data Slot that contains eventual externals seasonal regressors used for fitting.
                               seasonal_data = dplyr::tibble(n = c(seq(1, 59), 59.5, seq(60, 365))),
                               #' @field extra_params Slot used for containing eventual extra parameters.
                               extra_params = list(),
                               #' @method initialize seasonalModel
                               #' @description
                               #' Initialize an object of the class `seasonalModel`.
                               #' @param order numeric, number of sine and cosine used in fitting.
                               #' @param period numeric, seasonal periodicity. The default is \eqn{\frac{2\pi}{365}}.
                               initialize = function(order = 1, period = 365){
                                 # Store period and order
                                 private$..period = period
                                 private$..order = order
                               },
                               #' @method fit seasonalModel
                               #' @description
                               #' Fit a seasonal model as a linear combination of sine and cosine functions and
                               #' eventual external regressors specified in the formula. The external regressors used should
                               #' have the same periodicity, i.e. not stochastic regressors are allowed.
                               #' @param formula formula, an object of class `formula` (or one that can be coerced to that class).
                               #' It is a symbolic description of the model to be fitted and can be used to include or exclude the intercept or external regressors in `data`.
                               #' @param data 	an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
                               #' If not found in data, the variables are taken from environment(formula), typically the environment from which `lm` is called.
                               #' @param ... other parameters to be passed to the function lm.
                               fit = function(formula, data, ...){
                                 # Model formula
                                 if (self$order > 0) {
                                   for (i in 1:self$order){
                                     formula <- paste0(formula, " + ", "cos((2*base::pi)/", eval(self$period), "*n*", i, ")", " + ",
                                                       "sin((2*base::pi)/", eval(self$period), "*n*", i, ")")
                                   }
                                 }
                                 # Store the formula
                                 private$formula <- as.formula(formula)
                                 # Fit seasonal model
                                 private$..model <- lm(private$formula, data = data, ...)
                                 # Extract external seasonal regressors
                                 external_regressors <- formula.tools::get.vars(formula(private$..model))[-c(1)]
                                 external_regressors <- external_regressors[!stringr::str_detect(external_regressors, "cos|sin")]
                                 # Store external seasonal regressors inside `seasonal_data` slot
                                 if (!purrr::is_empty(external_regressors)) {
                                   external_regressors <- data[, c("n", external_regressors)]
                                   external_regressors <- external_regressors[!duplicated(external_regressors$n),]
                                   self$seasonal_data <- dplyr::left_join(self$seasonal_data, external_regressors, by = c("n"))
                                 }
                               },
                               #' @method predict seasonalModel
                               #' @description
                               #' Predict method for the class `seasonalModel`.
                               #' @param n integer, number of day of the year.
                               predict = function(n){
                                 if (missing(n)){
                                   predict.lm(private$..model)
                                 } else {
                                   n <- number_of_day(n)
                                   newdata <- data.frame(n = n %% self$period)
                                   newdata$n <- ifelse(newdata$n == 0, 1, newdata$n)
                                   # Add the external regressors only when specified
                                   if (ncol(self$seasonal_data) > 1) {
                                     newdata <- dplyr::left_join(newdata, self$seasonal_data, by = "n")
                                   }
                                   predict.lm(private$..model, newdata = newdata)
                                 }
                               },
                               #' @method update seasonalModel
                               #' @description
                               #' Update the parameters inside the model.
                               #' @param coefficients vector of parameters.
                               update = function(coefficients){
                                 old_coef <- self$coefficients
                                 # Check length
                                 if (length(old_coef) != length(coefficients)) {
                                   warning("The lenght of `coefficients` do not match the length of the old coefficients.")
                                   return(invisible(NULL))
                                 }
                                 # Update parameters
                                 private$..model$coefficients <- coefficients
                               }
                             ),
                             private = list(
                               ..model = NA,
                               formula = NA,
                               ..period = 1,
                               ..order = 365
                             ),
                             active = list(
                               #' @field coefficients Get a vector with the fitted coefficients.
                               coefficients = function(){
                                  private$..model$coefficients
                                  },
                               #' @field model Get the fitted `lm` object.
                               model = function(){
                                  private$..model
                               },
                               #' @field period Get the seasonality in days.
                               period = function(){
                                 private$..period
                               },
                               #' @field order Get the number of combinations of sines and cosines used.
                               order = function(){
                                 private$..order
                               }
                             )
)


