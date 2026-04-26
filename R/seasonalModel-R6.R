#' Seasonal Model
#'
#' @description
#' The `seasonalModel` class implements a seasonal regression model as a linear combination of sine and cosine functions.
#' This model is designed to capture periodic effects in time series data, particularly for applications involving seasonal trends.
#'
#' @details
#' The seasonal model is fitted using a specified formula, which allows for the inclusion of external regressors along with sine and cosine terms
#' to model seasonal variations. The periodicity can be customized, and the model can be updated with new coefficients after the initial fit.
#'
#' @examples
#' # Sample time series
#' n <- 1000
#' eps <- rnorm(n, sd = 0.2)
#' Yt <- 0.5 + 0.2 * sin(2*base::pi/365*c(1:n)) + 0.3 * cos(2*base::pi/365*c(1:n)) + eps
#'
#' # Initialize the model
#' sm <- seasonalModel$new("Yt ~ 1", 1, 365)
#' # Fit the model
#' data = data.frame(Yt = Yt, n = 1:1000)
#' sm$fit(data = data)
#' data$Yt_hat <- sm$predict(n = data$n)
#' # Fitted parameters
#' sm$coefficients
#'
#' # Plot the fitted data
#' library(ggplot2)
#' ggplot()+
#' geom_line(data = data, aes(n, Yt))+
#' geom_line(data = data, aes(n, Yt_hat), color = "red")
#'
#' @rdname seasonalModel
#' @name seasonalModel
#' @keywords seasonalModel
#' @note Version 1.0.3
#' @export
seasonalModel <- R6::R6Class("seasonalModel",
                             public = list(
                               #' @field extra_params List to contain custom extra parameters.
                               extra_params = list(),
                               #' @field control List to contain custom control parameters.
                               control = list(include.intercept = TRUE),
                               #' @description
                               #' Initialize a `seasonalModel` object.
                               #' @param order Integer, number of combinations of sines and cosines.
                               #' @param period Integer, seasonality period. The default is 365.
                               #' @param formula formula, an object of class `formula` (or one that can be coerced to that class).
                               #' It is a symbolic description of the model to be fitted and can be used to include or exclude the intercept or external regressors in `data`.
                               initialize = function(formula, order = 1, period = 365){
                                 # Store period and order
                                 private[["..period"]] = period
                                 private[["..order"]] = order
                                 # Formula with standard names
                                 formulas <- seasonalModel_formula(formula, order = order, period = period)

                                 # Main model
                                 private[["formula"]] <- formulas$formula
                                 # Differential model w.r.t. t
                                 private[["formula_dt"]] <- as.formula(formulas$formula_dt)

                                 # Extract regressors from the formula excluding target variable
                                 regressors <- formula.tools::rhs.vars(private$formula)
                                 private$target <- formula.tools::lhs.vars(private$formula)
                                 # Detect number of regressors
                                 n_regressors <- length(regressors)
                                 # Detect intercept from the formula
                                 self$control$include.intercept <- !stringr::str_detect(as.character(private$formula), "-1") & any(!stringr::str_detect(regressors, "-1"))
                                 # Seasonal regressors
                                 idx_seasonal_regressors <- which(stringr::str_detect(regressors, "sin|cos"))
                                 n_seasonal_reg <- length(idx_seasonal_regressors)
                                 # External regressors
                                 idx_external_regressors <- which(!stringr::str_detect(regressors, "sin|cos"))
                                 private[["external_regressors"]] <- regressors[idx_external_regressors]
                                 n_external_regressors <- length(idx_external_regressors)
                                 # Standard names
                                 coefs_names <- regressors
                                 coefs_names[idx_seasonal_regressors] <- formulas$coefs_names
                                 coefs_names[idx_external_regressors] <- regressors[idx_external_regressors]
                                 # Intercept
                                 if (self$control$include.intercept) {
                                   coefs_names <- c("intercept", coefs_names)
                                 }
                                 # Initialize a vector of coefficients
                                 private[["..coefs_names"]] <- setNames(coefs_names, coefs_names)
                                 private[["..coefficients"]] <- setNames(rep(0, length(coefs_names)), coefs_names)
                                 private[["..dcoefficients"]] <- setNames(rep(0, length(coefs_names)), coefs_names)
                                 private[["..dcoefficients"]][stringr::str_detect(coefs_names, "intercept")] <- 0
                                 private[["..std.errors"]] <- setNames(rep(0, length(coefs_names)), coefs_names)
                               },
                               #' @description
                               #' Fit a seasonal model as a linear combination of sine and cosine functions and
                               #' eventual external regressors specified in the formula. The external regressors used should
                               #' have the same periodicity, i.e. not stochastic regressors are allowed.
                               #' @param data A data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
                               #' @param ... other parameters to be passed to the function `lm`.
                               fit = function(data, ...){
                                 # Fit seasonal model
                                 lm_model <- lm(private$formula, data = data, ...)
                                 # Standard coefficients names
                                 coefs_names <- self$coefs_names
                                 # Extract the coefficients
                                 private[["..coefficients"]] <- setNames(lm_model$coefficients, coefs_names)
                                 # Update differential parameters
                                 dcoefs <- self$coefficients
                                 dcoefs[stringr::str_detect(coefs_names, "cos")] <- -dcoefs[stringr::str_detect(coefs_names, "cos")]
                                 if (self$control$include.intercept) dcoefs[1] <- 0
                                 private[["..dcoefficients"]] <- dcoefs
                                 # Update std.errors
                                 private[["..std.errors"]] <- setNames(broom::tidy(lm_model)$std.error, coefs_names)
                               },
                               #' @description
                               #' Predict method for the class `seasonalModel`.
                               #' @param n Integer vector, numbers of day of the year.
                               #' @param newdata an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
                               predict = function(n, newdata){
                                 if (missing(newdata)) {
                                   newdata <- data.frame(x = 0, n = n)
                                   colnames(newdata) <- c(private$target, "n")
                                   if (!is.null(private$external_regressors[1])){
                                     for(regressor in private$external_regressors){
                                       newdata[[regressor]] <- 0
                                     }
                                   }
                                   drop(model.matrix(private$formula, data = newdata) %*% self$coefficients)
                                 } else {
                                   newdata[[private$target]] <- 0
                                   drop(model.matrix(private$formula, data = newdata) %*% self$coefficients)
                                 }
                               },
                               #' @description
                               #' Compute the differential with respect to time (`n`).
                               #' @param n Integer, number of day of the year.
                               #' @param newdata an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
                               differential = function(n, newdata){
                                 if (missing(newdata)) {
                                   newdata <- data.frame(x = 0, n = n)
                                   newdata <- data.frame(x = 0, n = n)
                                   colnames(newdata) <- c(private$target, "n")
                                   if (!is.null(private$external_regressors[1])){
                                     for(regressor in private$external_regressors){
                                       newdata[[regressor]] <- 0
                                     }
                                   }
                                   drop(model.matrix(private$formula_dt, data = newdata) %*% self$dcoefficients)
                                 } else {
                                   newdata[[private$target]] <- 0
                                   drop(model.matrix(private$formula_dt, data = newdata) %*% self$dcoefficients)
                                 }
                               },
                               #' @description
                               #' Update the parameters.
                               #' @param coefficients Named vector, with new parameters, the names should match the names of the coefficients.
                               update = function(coefficients){
                                 # Extract old coefficients
                                 new_coefs <- self$coefficients
                                 # Extract names
                                 names_old <- self$coefs_names
                                 names_new <- names(coefficients)
                                 # Update only if they are present
                                 for(i in 1:length(coefficients)){
                                   condition <- names_new[i] %in% names_old
                                   if (condition) {
                                     old_coef <- new_coefs[names_new[i]]
                                     if (old_coef != coefficients[i]){
                                       private[["..coefficients"]][names_new[i]] <- coefficients[i]
                                       private[["..std.errors"]][names_new[i]]   <- NA_integer_
                                     }
                                   }
                                 }
                                 # ************************************************************
                                 # Update coefficients of the differential
                                 new_dcoefs <- self$coefficients
                                 new_dcoefs[stringr::str_detect(names_old, "cos")] <- -new_dcoefs[stringr::str_detect(names_old, "cos")]
                                 if (self$control$include.intercept) new_dcoefs[1] <- 0
                                 private[["..dcoefficients"]] <- new_dcoefs
                               },
                               #' @description
                               #' Update the names of the parameters.
                               #' @param coefs_names Named vector, new parameters names, the names should match the names of the old coefficients.
                               update_coefs_names = function(coefs_names){
                                 # Extract old coefficients names
                                 new_coefs_names <- self$coefs_names
                                 # Extract names
                                 names_old <- names(new_coefs_names)
                                 names_new <- names(coefs_names)
                                 # Update only if they are present
                                 for(i in 1:length(coefs_names)){
                                   condition <- names_new[i] %in% names_old
                                   if (condition) {
                                     private[["..coefs_names"]][names_new[i]] <- coefs_names[i]
                                   }
                                 }
                                 # Store updated coefficients names
                                 private[["..coefs_names"]] <- setNames(private[["..coefs_names"]], private[["..coefs_names"]])
                                 # Update all the names
                                 coefs_names <- self$coefs_names
                                 private[["..coefficients"]] <- setNames(private[["..coefficients"]], coefs_names)
                                 private[["..dcoefficients"]] <- setNames(private[["..dcoefficients"]], coefs_names)
                                 private[["..std.errors"]] <- setNames(private[["..std.errors"]], coefs_names)
                               },
                               #' @description
                               #' Update the standard errors of the parameter.
                               #' @param std.errors Named vector with new standard errors, the names should match the names of the coefficients.
                               update_std.errors = function(std.errors){
                                 # Extract old coefficients
                                 new_std.errors <- self$std.errors
                                 # Extract names
                                 names_old <- names(new_std.errors)
                                 names_new <- names(std.errors)
                                 # Update only if they are present
                                 for(i in 1:length(std.errors)){
                                   condition <- names_new[i] %in% names_old
                                   if (condition) {
                                     private[["..std.errors"]][names_new[i]] <- std.errors[i]
                                   }
                                 }
                               },
                               #' @description
                               #' Print method for the class `seasonalModel`.
                               print = function(){
                                 cat(paste0("----------------------- seasonalModel ----------------------- \n"))
                                 msg_1 <- paste0("-  Order: ", self$order, "\n - Period: ", self$period, "\n")
                                 if (any(is.na(private$external_regressors))) {
                                   msg_2 <- paste0("- External regressors: 0 \n")
                                 } else {
                                   n_external_regressors <- length(private$external_regressors)
                                   msg_2 <- paste0("- External regressors: ", n_external_regressors, " (", private$external_regressors, ")\n")
                                 }
                                 msg_3 <- paste0("- Version: ", private$version, "\n")
                                 cat(c(msg_1, msg_2, msg_3))
                                 cat("-------------------------------------------------------------\n")
                                 print(self$coefficients)
                               }
                             ),
                             private = list(
                               version = "1.0.3",
                               target = "",
                               formula = NA,
                               formula_dt = NA,
                               ..coefs_names = NULL,
                               ..coefficients = c(),
                               ..dcoefficients = c(),
                               ..period = 1,
                               ..order = 365,
                               ..std.errors = c(),
                               external_regressors = NULL
                             ),
                             active = list(
                               #' @field period Integer vector, periodicities of the seasonality.
                               period = function(){
                                 private$..period
                               },
                               #' @field order Integer vector, number of combinations of sines and cosines.
                               order = function(){
                                 private$..order
                               },
                               #' @field omega Integer, periodicity computed as \eqn{\frac{2 \pi}{\text{period}}}.
                               omega = function(){
                                 2 * base::pi / private$..period
                               },
                               #' @field coefficients Named vector, model's parameters.
                               coefficients = function(){
                                 private$..coefficients
                               },
                               #' @field dcoefficients Named vector, model's parameters for the differential with respect to `n`.
                               dcoefficients = function(){
                                 private$..dcoefficients
                               },
                               #' @field std.errors Named vector, standard errors of the model's parameters.
                               std.errors = function(){
                                 private$..std.errors
                               },
                               #' @field coefs_names Named vector, names of the model's parameters.
                               coefs_names = function(){
                                 private$..coefs_names
                               },
                               #' @field tidy A tibble with estimated parameters and std. errors.
                               tidy = function(){
                                 dplyr::tibble(
                                   term = private$coefs_names,
                                   estimate = self$coefficients,
                                   std.error = private$..std.errors
                                 )
                               }
                             )
)
