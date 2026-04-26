#' ARMA model
#'
#' @description
#' R6 class for ARMA(p, q) model
#'
#' @rdname ARMA_modelR6
#' @name ARMA_modelR6
#' @keywords ARMA
#' @note Version 1.0.2
#' @seealso [stats::arima()] which is wrapped in the method `fit`.
#' @export
ARMA_modelR6 <- R6::R6Class("ARMA_modelR6",
                          public = list(
                            #' @field control list, to contain custom control parameters.
                            control = list(include.intercept = FALSE),
                            #' @field model_name character(1), standard model name.
                            model_name = "ARMA(0, 0)",
                            #' @description
                            #' Initialize an ARMA model
                            #' @param arOrder integer(1), order of the Autoregressive component.
                            #' @param maOrder integer(1), order of the Moving-Average component.
                            #' @param include.intercept logical(1), the default is `FALSE`. When `TRUE` the intercept will be included.
                            initialize = function(arOrder = 1, maOrder = 1, include.intercept = FALSE){
                              # Logical value for the intercept parameter
                              self$control[["include.intercept"]] <- include.intercept
                              # Store AR order
                              private[["..arOrder"]] <- arOrder
                              if (arOrder > 0){
                                private[["..phi"]] <- setNames(rep(0, arOrder), paste0("phi_", 1:arOrder))
                              }
                              # Store MA order
                              private[["..maOrder"]] <- maOrder
                              if (maOrder > 0){
                                private[["..theta"]] <- setNames(rep(0, maOrder), paste0("theta_", 1:maOrder))
                              }
                              # Standard errors
                              private[["..std.errors"]] <- c(intercept = NA, self$phi, self$theta)
                              # Pre-compute vector b
                              private[["..b"]] <- ARMA_vector_b(arOrder, maOrder)
                              # Model name
                              self$model_name <- paste0("ARMA", " (", arOrder, ", ", maOrder, ")")
                            },
                            #' @description
                            #' Fit the ARMA model with `arima` function.
                            #' @param x numeric vector, time series to fit.
                            fit = function(x){
                              # Fitted model
                              ARMA_model <- arima(x, order = c(self$order[1], 0, self$order[2]),
                                                  include.mean = self$control$include.intercept, method = "CSS")
                              # Standardize parameters names
                              params <- ARMA_model$coef
                              # Extract the std.errors
                              std.errors <- sqrt(diag(ARMA_model$var.coef))
                              names(std.errors) <- names(params)
                              # Extract intercept
                              intercept <- c(intercept = 0)
                              std.errors_intercept <- c(intercept = NA)
                              if (self$control$include.intercept) {
                                index <- which(names(params) == "intercept")
                                intercept <- c(intercept = params[[index]])
                                std.errors_intercept <- c(intercept = std.errors[[index]])
                                std.errors <- std.errors[-c(index)]
                                params <- params[-c(index)]
                              }
                              # AR coefficients without intercept
                              phi <- c(phi_1 = 0)
                              std.errors_ar <- c(phi_1 = NA)
                              if (self$order[1] > 0){
                                index <- stringr::str_detect(names(params), "ar")
                                phi <- params[index]
                                std.errors_ar <- std.errors[index]
                                names(phi) <- names(std.errors_ar) <- paste0("phi_", 1:self$order[1])
                              }
                              # MA coefficients without intercept
                              theta <- c(theta_1 = 0)
                              std.errors_ma <- c(theta_1 = NA)
                              if (self$order[2] > 0){
                                index <- stringr::str_detect(names(params), "ma")
                                theta <- params[index]
                                std.errors_ma <- std.errors[index]
                                names(theta) <- names(std.errors_ma) <- paste0("theta_", 1:self$order[2])
                              }
                              # *********************************************
                              # Store coefficients
                              private$..intercept <- intercept
                              private$..phi <- phi
                              private$..theta <- theta
                              # Compute companion matrix
                              private$..A <- ARMA_companion_matrix(phi, theta)
                              # Store the std. errors
                              private$..std.errors <- c(std.errors_intercept, std.errors_ar, std.errors_ma)
                              # Update fitted variance
                              private$..sigma2 <- sqrt(ARMA_model$sigma2)
                            },
                            #' @description
                            #' Filter the time-series and compute fitted values and residuals. See the function \code{\link{ARMA_filter}} for more details.
                            #' @param x numeric vector, time series to filter.
                            filter = function(x){
                              ARMA_filter(x, self$A, self$b, self$intercept)
                            },
                            #' @description
                            #' Next step function. See the function \code{\link{ARMA_next_step}} for more details.
                            #' @param x Numeric vector, state vector with past observations and residuals.
                            #' @param n.ahead integer(1), number of steps ahead.
                            #' @param eps optional numeric vector of length n.ahead, next step realized residuals.
                            next_step = function(x, n.ahead = 1, eps = 0){
                              ARMA_next_step(n.ahead, x, self$A, self$b, self$intercept, eps)
                            },
                            #' @description
                            #' Forecast expected value. See the function \code{\link{ARMA_expectation}} for more details.
                            #' @param h integer(1), number of steps ahead.
                            #' @param X0 numeric vector of size p + q. State vector of past values.
                            expectation = function(h = 1, X0){
                              if (missing(X0)){
                                X0 <- rep(0, sum(self$order))
                              }
                              ARMA_expectation(h, X0, self$A, self$b, self$intercept)
                            },
                            #' @description
                            #' Forecast variance. See the function \code{\link{ARMA_variance}} for more details.
                            #' @param h integer(1), number of steps ahead.
                            #' @param sigma2 integer(1), standard deviation of the residuals.
                            variance = function(h = 1, sigma2 = 1){
                              ARMA_variance(h, self$A, self$b, sigma2)
                            },
                            #' @description
                            #' Update the model's parameters
                            #' @param coefficients Numeric named vector, model's coefficients. If missing nothing will be updated.
                            update = function(coefficients){
                              if (missing(coefficients)) {
                                return(invisible(NULL))
                              }
                              # Extract old coefficients
                              new_coefs <- self$coefficients
                              # Extract names
                              names_old <- names(new_coefs)
                              names_new <- names(coefficients)
                              # Update only if they are present
                              for(i in 1:length(coefficients)){
                                # Check if new parameters matches the parameters in the model
                                condition <- names_new[i] %in% names_old
                                # Update the parameter
                                if (condition) {
                                  old_coef <- new_coefs[names_new[i]]
                                  if (old_coef != coefficients[i]){
                                    new_coefs[names_new[i]] <- coefficients[i]
                                    private[["..std.errors"]][names_new[i]] <- NA_integer_
                                  }
                                }
                              }
                              # Update intercept
                              if (self$control$include.intercept){
                                private[["..intercept"]] <- new_coefs["intercept"]
                              }
                              # Update AR parameters
                              if (self$order[1] > 0) {
                                private[["..phi"]] <- new_coefs[stringr::str_detect(names_old, "phi")]
                              }
                              # Update MA parameters
                              if (self$order[2] > 0) {
                                private[["..theta"]] <- new_coefs[stringr::str_detect(names_old, "theta")]
                              }
                              if (!self$control$include.intercept){
                                new_coefs <- new_coefs[-which(names_old == "intercept")]
                              } else {
                                index <- which(names_old == "intercept")
                                new_coefs <- c(new_coefs[-index], new_coefs[index])
                              }
                              # Update companion matrix
                              private[["..A"]] <- ARMA_companion_matrix(self$phi, self$theta)
                            },
                            #' @description
                            #' Update the standard errors of the parameters.
                            #' @param std.errors Numeric named vector, parameters' standard errors. If missing nothing will be updated.
                            update_std.errors = function(std.errors){
                              if (missing(std.errors) || purrr::is_empty(std.errors)) {
                                return(invisible(NULL))
                              }
                              # Extract old coefficients
                              new_std.errors <- private$..std.errors
                              # Extract names
                              names_old <- names(new_std.errors)
                              names_new <- names(std.errors)
                              # Update only if they are present
                              for(i in 1:length(std.errors)){
                                if (names_new[i] %in% names_old) {
                                  new_std.errors[names_new[i]] <- std.errors[i]
                                }
                              }
                              # Update private std. errors
                              private[["..std.errors"]] <- new_std.errors
                            },
                            #' @description
                            #' Update the variance of the residuals.
                            #' @param sigma2 integer(1), variance of the residuals.
                            update_sigma2 = function(sigma2){
                              if (!missing(sigma2)){
                                private[["..sigma2"]] <- sigma2
                              }
                            },
                            #' @description
                            #' Print method for `ARMA_modelR6` class.
                            print = function(){
                              green <- function(x) paste0("\033[1;32m", x, "\033[0m")
                              red <- function(x) paste0("\033[1;31m", x, "\033[0m")
                              msg_col <- function(x) ifelse(x, green(x), red(x))
                              # Intercept
                              par.intercept <- format(self$intercept, digits = 4)
                              std.error.intercept <- private$..std.errors[names(private$..std.errors) == names(par.intercept)]
                              format.intercept <- paste0(par.intercept, " (", format(std.error.intercept, digits = 3), ")")
                              # AR parameters
                              par.ar <- format(self$phi, digits = 4)
                              std.error.ar <- private$..std.errors[names(private$..std.errors) %in% names(par.ar)]
                              format.ar <- paste0(par.ar, " (", format(std.error.ar, digits = 3), ")")
                              # MA parameters
                              par.ma <- format(self$theta, digits = 4)
                              std.error.ma <- private$..std.errors[names(private$..std.errors) %in% names(par.ma)]
                              format.ma <- paste0(par.ma, " (", format(std.error.ma, digits = 3), ")")
                              # ********************************************************
                              msg0 <- paste0("--------------------- ", self$model_name, " ---------------------")
                              cat(msg0,  "\n")
                              cat(" Include Intercept: ", msg_col(self$control$include.intercept), "\n",
                                  "AR: ", msg_col(self$order[1] != 0), "\n",
                                  "MA: ", msg_col(self$order[2] != 0), "\n",
                                  "Intercept: ", format.intercept, " \n",
                                  "AR parameters: ", format.ar, " \n",
                                  "MA parameters: ", format.ma, " \n",
                                  "Version: ", private$version, "\n",
                                  paste0(rep("-", length(strsplit(msg0, "")[[1]])-1), collapse = ""), "\n")
                            }
                          ),
                          private = list(
                            version = "1.0.2",
                            ..arOrder = 0,
                            ..maOrder = 0,
                            ..intercept = c(intercept = 0),
                            ..phi = c(phi_1 = 0),
                            ..theta = c(theta_1 = 0),
                            ..b = c(),
                            ..A = c(),
                            ..sigma2 = 1,
                            ..std.errors = c()
                          ),
                          active = list(
                            #' @field arOrder integer(1), Autoregressive order.
                            arOrder = function(){
                              private$..arOrder
                            },
                            #' @field maOrder integer(1), Moving-Average order.
                            maOrder = function(){
                              private$..maOrder
                            },
                            #' @field order named numeric vector of size 2. Orders of the ARMA model: the first element is the `arOrder`, while the second the `maOrder`.
                            order = function(){
                              c(AR = private$..arOrder, MA = private$..maOrder)
                            },
                            #' @field intercept named integer(1), intercept of the model.
                            intercept = function(){
                              private$..intercept
                            },
                            #' @field phi named numeric vector of size arOrder, AR parameters. If `arOrder = 0`, the parameter is zero.
                            phi = function(){
                              private$..phi
                            },
                            #' @field theta named numeric vector of size maOrder, MA parameters. If `maOrder = 0`, the parameter is zero.
                            theta = function(){
                              private$..theta
                            },
                            #' @field coefficients named numeric vector of size arOrder + maOrder + 1.
                            #' The first element ie sthe intercept, then the ARMA parameters.
                            coefficients = function(){
                              c(private$..intercept, private$..phi, private$..theta)
                            },
                            #' @field std.errors Numeric named vector, standard errors of the intercept and ARMA parameters.
                            std.errors = function(){
                              private$..std.errors
                            },
                            #' @field sigma2 integer(1), std.errors of the residuals.
                            sigma2 = function(){
                              private$..sigma2
                            },
                            #' @field A numeric matrix of size (arOrder + maOrder) x (arOrder + maOrder). See the function \code{\link{ARMA_companion_matrix}} for more details.
                            A = function(){
                              private$..A
                            },
                            #' @field b numeric vector of size arOrder + maOrder. See the function \code{\link{ARMA_vector_b}} for more details.
                            b = function(){
                              private$..b
                            },
                            #' @field tidy Tibble with estimated parameters and standard errors.
                            tidy = function(){
                              dplyr::tibble(
                                term = names(self$coefficients),
                                estimate = self$coefficients,
                                std.error = self$std.errors
                              )
                            }
                          ))


