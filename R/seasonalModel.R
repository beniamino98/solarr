#' Seasonal Model
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
                               #' @field extra_params List customizable to contain eventual extra parameters.
                               extra_params = list(),
                               #' @method initialize seasonalModel
                               #' @description
                               #' Initialize a `seasonalModel` object.
                               #' @param order Integer, number of combinations of sines and cosines.
                               #' @param period Integer, seasonality period. The default is 365.
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
                               #' @param ... other parameters to be passed to the function `lm`.
                               fit = function(formula, data, ...){
                                 base_formula <- formula
                                 # Model formula
                                 if (self$order > 0) {
                                   for (i in 1:self$order){
                                     formula <- paste0(formula, " + ", "sin((2*base::pi)/", eval(self$period), "*n*", i, ")", " + ",
                                                       "cos((2*base::pi)/", eval(self$period), "*n*", i, ")")
                                   }
                                 }
                                 # Store the formula
                                 private$mformula <- as.formula(formula)
                                 # Fit seasonal model
                                 private$..model <- lm(private$mformula, data = data, ...)
                                 # Extract external regressors
                                 regressors <- formula.tools::get.vars(formula(private$..model))[-c(1)]
                                 external_regressors <- which(!stringr::str_detect(regressors, "sin|cos"))
                                 private$external_regressors <- regressors[external_regressors]
                                 # Reparametrize the coefficients
                                 coefs <- self$coefficients2
                                 coef_phi <- coefs[stringr::str_detect(names(coefs), "phi")]
                                 coef_B <- coefs[stringr::str_detect(names(coefs), "B")]
                                 formula <- ""
                                 if (self$order > 0) {
                                   for (i in 1:self$order){
                                     formula <- paste0(formula, " +", "I((2*base::pi)/", eval(self$period),
                                                       "*cos((2*base::pi)/", eval(self$period), "*n*", i, " + ", eval(coef_phi[i]), "))")
                                   }
                                 }
                                 # Model formula
                                 formula <- paste0(formula.tools::lhs(private$mformula), " ~ ", formula, " - 1")
                                 # Store the formula
                                 private$dformula <- as.formula(formula)
                                 # Fit seasonal model
                                 private$..dmodel <- lm(private$dformula, data = data)
                                 private$..dmodel$coefficients <- coef_B
                               },
                               #' @method predict seasonalModel
                               #' @description
                               #' Predict method for the class `seasonalModel`.
                               #' @param n integer, number of day of the year.
                               #' @param dt Numeric, time step.
                               predict = function(n, newdata, dt = 1){
                                 if (missing(newdata)) {
                                   if (missing(n)) {
                                     predict.lm(private$..model)
                                   } else {
                                     n <- number_of_day(n) %% self$period
                                     predict.lm(private$..model, newdata = data.frame(n = n))
                                   }
                                 } else {
                                   predict.lm(private$..model, newdata = newdata)
                                 }
                               },
                               #' @method differential seasonalModel
                               #' @description
                               #' Compute the differential of the sinusoidal function. It do not consider the differential of
                               #' eventual external regressors.
                               #' @param n Integer, number of day of the year.
                               #' @param dt Numeric, time step.
                               differential = function(n, newdata, dt = 1){
                                 if (missing(newdata)) {
                                   if (missing(n)) {
                                     predict.lm(private$..dmodel)
                                   } else {
                                     n <- (number_of_day(n) + dt) %% self$period
                                     predict.lm(private$..dmodel, newdata = data.frame(n = n))
                                   }
                                 } else {
                                   predict.lm(private$..dmodel, newdata = newdata)
                                 }
                               },
                               #' @method update seasonalModel
                               #' @description
                               #' Update the parameters inside the model.
                               #' @param coefficients A named vector with coefficients.
                               update = function(coefficients){
                                 old_coef <- self$coefficients
                                 # Check length
                                 if (length(old_coef) != length(coefficients)) {
                                   warning("The lenght of `coefficients` do not match the length of the old coefficients.")
                                   return(invisible(NULL))
                                 }
                                 # Update parameters
                                 private$..model$coefficients <- coefficients
                               },
                               #' @description
                               #' Print method for the class `seasonalModel`.
                               print = function(){
                                 cat(paste0("----------------------- seasonalModel ----------------------- \n"))
                                 msg_1 <- paste0(" - Order: ", self$order, "\n - Period: ", self$period, "\n")
                                 if (any(is.na(private$external_regressors))) {
                                   msg_2 <- paste0("- External regressors: 0 \n")
                                 } else {
                                   n_external_regressors <- length(private$external_regressors)
                                   msg_2 <- paste0("- External regressors: ", n_external_regressors, " (", private$external_regressors, ")\n")
                                 }
                                 cat(c(msg_1, msg_2))
                                 cat(paste0("--------------------------------------------------------------\n"))
                                 print(self$model)
                               }
                             ),
                             private = list(
                               ..model = NA,
                               ..dmodel = NA,
                               mformula = NA,
                               dformula = NA,
                               ..period = 1,
                               ..order = 365,
                               external_regressors = NA
                             ),
                             active = list(
                               #' @field coefficients A named vector with the fitted coefficients.
                               coefficients = function(){
                                 private$..model$coefficients
                               },
                               #' @field coefficients2 A named vector with the coefficients reparametrized
                               #' to obtain a linear combination of only shifted sine functions.
                               coefficients2 = function(){
                                 # Extract seasonal regressors
                                 regressors <- formula.tools::get.vars(formula(private$..model))[-c(1)]
                                 seasonal_regressors <- which(stringr::str_detect(regressors, "sin|cos")) + 1
                                 external_regressors <- which(!stringr::str_detect(regressors, "sin|cos")) + 1
                                 # Reparametrized parameters
                                 coefs <- private$..model$coefficients[c(1, seasonal_regressors)]
                                 index <- seq(1, length(coefs)-1, by = 2)
                                 param <- c(private$..model$coefficients[c(1, external_regressors)])
                                 names(param) <- c("A", names(private$..model$coefficients[external_regressors]))
                                 k <- 1
                                 i <- 1
                                 for(i in index){
                                   coefs_names <- names(param)
                                   param <- c(param, sqrt(coefs[i + 1]^2 + coefs[i + 2]^2))
                                   param <- c(param, atan(coefs[i + 2]/coefs[i + 1]) + ifelse(coefs[i + 1] < 0, -base::pi, 0))
                                   names(param) <- c(coefs_names, paste0(c("B", "phi"), k))
                                   k <- k + 1
                                 }
                                 return(param)
                               },
                               #' @field model A slot with the fitted `lm` object.
                               model = function(){
                                 private$..model
                               },
                               #' @field period Integer, the seasonality period.
                               period = function(){
                                 private$..period
                               },
                               #' @field order Integer, number of combinations of sines and cosines.
                               order = function(){
                                 private$..order
                               }
                             )
)
