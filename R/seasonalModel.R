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
                               initialize = function(formula, order = 1, period = 365, data, ...){
                                 # Model formula
                                 if (order > 0) {
                                   for (i in 1:order){
                                     formula <- paste0(formula, " + ", "cos((2*base::pi)/", eval(period), "*n*", i, ")", " + ",
                                                       "sin((2*base::pi)/", eval(period), "*n*", i, ")")
                                   }
                                 }
                                 # Store the formula
                                 private$formula <- as.formula(formula)
                                 # Fit seasonal model
                                 private$..model <- lm(private$formula, data = data, ...)
                                 # Extract external seasonal regressors
                                 external_regressors <- formula.tools::get.vars(formula(private$..model))[-c(1)]
                                 external_regressors <- external_regressors[!stringr::str_detect(external_regressors, "cos|sin")]
                                 if (!purrr::is_empty(external_regressors)){
                                   external_regressors <- data[, c("n", external_regressors)]
                                   external_regressors <- external_regressors[!duplicated(external_regressors$n),]
                                   private$external_regressors <- dplyr::left_join(private$external_regressors, external_regressors, by = c("n"))
                                 }
                                 # Store model, period and order
                                 private$period = period
                                 private$order = order
                               },
                               #' @description
                               #' Predict method for a `seasonalModel`,
                               #' @param n integer, number of day of the year.
                               predict = function(n){
                                 if (missing(n)){
                                   predict.lm(private$..model)
                                 } else {
                                   n <- number_of_day(n)
                                   newdata <- data.frame(n = n %% private$period)
                                   newdata$n <- ifelse(newdata$n == 0, 1, newdata$n)
                                   newdata <- dplyr::left_join(newdata, private$external_regressors, by = "n")
                                   predict.lm(private$..model, newdata = newdata)
                                 }
                               },
                               #' @description
                               #' Update the parameters of a `seasonalModel`.
                               #' @param coefficients vector of parameters.
                               update = function(coefficients){
                                 old_coef <- private$..model$coefficients
                                 # Check lenght
                                 if (length(old_coef) != length(coefficients)){
                                   warning("The lenght of `coefficients` do not match the length of the old coefficients.")
                                   return(invisible(NULL))
                                 }
                                 # Update parameters
                                 private$..model$coefficients <- coefficients
                               }
                             ),
                             private = list(
                               ..model = NA,
                               period = NA,
                               order = NA,
                               formula = NA,
                               params = list(),
                               external_regressors = dplyr::tibble(n = c(seq(1, 59), 59.5, seq(60, 365)))
                             ),
                             active = list(
                                coefficients = function(){
                                  private$..model$coefficients
                                  },
                                model = function(){
                                  private$..model
                               }
                             )
)

