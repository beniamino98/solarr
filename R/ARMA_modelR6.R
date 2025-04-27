#' ARMA(p, q) model implementation in a class R6
#'
#'
#' @examples
#' model <- solarModel$new(spec)
#' model$fit()
#' # Reference time series
#' x <- model$data$Yt_tilde
#' # Initialize the model without intercept
#' arma <- ARMA_modelR6$new(arOrder = 2, maOrder = 2)
#' # Fit the model
#' arma$fit(x)
#' arma
#' arma$coefficients
#' arma$variance
#' # Next step prediction
#' arma$next_step(c(0.4, 0.2, 0.2, -0.2))
#' arma$next_step(c(0.4, 0.2, 0.2, -0.2), n.ahead = 10)
#' # Update coefficients
#' params <- arma$coefficients*1.01
#' arma$update(params)
#' arma
#' # All sample prediction
#' arma$filter(x, arma$model$residuals[c(1,2)])
#' @rdname ARMA_modelR6
#' @name ARMA_modelR6
#' @note Version 1.0.0
#' @export
#self <- arma$.__enclos_env__$self
#private <- arma$.__enclos_env__$private
ARMA_modelR6 <- R6::R6Class("ARMA_modelR6",
                          public = list(
                            #' @description
                            #' Initialize an ARMA model
                            #' @param arOrder Numeric, scalar. Order for Autoregressive component.
                            #' @param maOrder Numeric, scalar. Order for Moving-Average component.
                            #' @param include.intercept Logical. When `TRUE` the intercept will be included. The default is `FALSE`.
                            initialize = function(arOrder = 1, maOrder = 1, include.intercept = FALSE){
                              private$include.intercept <- include.intercept
                              private$..arOrder <- arOrder
                              private$..maOrder <- maOrder
                            },
                            #' @description
                            #' Fit the ARMA model with `arima` function.
                            #' @param x Numeric, vector. Time series to fit.
                            fit = function(x){
                              # Fitted model
                              ARMA_model <- arima(x, order = c(self$order[1], 0, self$order[2]),
                                                  include.mean = private$include.intercept, method = "CSS")
                              # Standardize parameters names
                              params <- ARMA_model$coef
                              # Extract intercept
                              intercept <- c(intercept = 0)
                              if (private$include.intercept) {
                                index <- which(names(params) == "intercept")
                                intercept <- c(intercept = params[[index]])
                                params <- params[-c(index)]
                              }
                              # AR coefficients without intercept
                              phi <- c()
                              if (self$order[1] > 0){
                                phi <- params[stringr::str_detect(names(params), "ar")]
                                names(phi) <- paste0("phi_", 1:self$order[1])
                              }
                              # MA coefficients without intercept
                              theta <- c()
                              if (self$order[2] > 0){
                                theta <- params[stringr::str_detect(names(params), "ma")]
                                names(theta) <- paste0("theta_", 1:self$order[2])
                              }
                              # Extract the std.errors
                              std.errors <- broom::tidy(ARMA_model)$std.error
                              names(std.errors) <- names(params)
                              # Store coefficients
                              private$..intercept <- intercept
                              private$..phi <- phi
                              private$..theta <- theta
                              # Store the fitted model
                              private$..model <- ARMA_model
                              # Store the std. errors
                              private$..std.errors <- std.errors
                            },
                            #' @description
                            #' Filter the time-series and compute fitted values and residuals.
                            #' @param x Numeric, vector. Time series to filter.
                            #' @param eps0 Numeric vector. Initial residuals of the same length of the MA order.
                            filter = function(x, eps0){
                              # ARMA order
                              k <- max(self$order)
                              # Length of the time series
                              n <- length(x)
                              # Vector to store the fitted residuals
                              e_hat <- eps0
                              # Vector to store the fitted time series
                              x_hat <- x
                              i <- k+1
                              # Initialize the state vector
                              x_t <- c(x[k:(k-self$order[1]+1)], eps0[self$order[2]:1])
                              i <- k+1
                              for(i in (k+1):n){
                                x_t <- self$next_step(x_t)
                                # Fitted series
                                x_hat[i] <- x_t[1,1]
                                # Fitted residuals
                                e_hat[i] <- x[i] - x_hat[i]
                                # Update state vector
                                x_t <- x_t + self$b * e_hat[i]
                              }
                              dplyr::tibble(fitted = x_hat, residuals = e_hat)
                            },
                            #' @description
                            #' Next step function
                            #' @param x Numeric, vector. State vector with past observations and residuals.
                            #' @param n.ahead Numeric, scalar. Number of steps ahead.
                            #' @param eps Numeric vector. Optional realized residuals.
                            next_step = function(x, n.ahead = 1, eps = 0){
                              if (n.ahead > 1 & (length(eps) == 1 && eps == 0)) {
                                eps <- rep(0, n.ahead)
                              } else if (length(eps) != n.ahead){
                                stop("eps must be of the same length of n.ahead when specified!")
                              }
                              # Forecasting loop
                              x_t <- x
                              for(step in 1:n.ahead){
                                # Next step state space vector
                                x_t <-  self$Phi %*% x_t + self$b * eps[step]
                                x_t[1] <- x_t[1] + self$intercept
                              }
                              return(x_t)
                            },
                            #' @description
                            #' Update the model's parameters
                            #' @param coefficients Numeric, named vector. Model's coefficients. If missing nothing is updated.
                            update = function(coefficients){
                              # 1) Update the parameters
                              if (!missing(coefficients)) {
                                names(coefficients) <- names(self$coefficients)
                                # Intercept
                                if(private$include.intercept){
                                  # Update intercept
                                  private$..intercept <- c(intercept = coefficients[[1]])
                                } else {
                                  # Set intercept equal to zero
                                  private$..intercept <- c(intercept = 0)
                                  # Remove intercept from coefficients
                                  coefficients <- coefficients[-c(1)]
                                }
                                # Update the parameters inside the ARMA model
                                private$..model$coef <- coefficients
                                # Update AR parameters
                                if (self$order[1] > 0){
                                  private$..phi <- coefficients[stringr::str_detect(names(coefficients), "phi")]
                                }
                                # Update MA parameters
                                if (self$order[2] > 0){
                                  private$..theta <- coefficients[stringr::str_detect(names(coefficients), "theta")]
                                }
                                # Set std. errors to NA
                                private$..std.errors <- rep(NA, length(coefficients))
                              }
                            },
                            #' @description
                            #' Update the std. errors of the parameters.
                            #' @param std.errors Numeric, named vector. Parameters std. errors.
                            update_std.errors = function(std.errors){
                              # Update the vector of std. errors
                              std.errors_updated <- self$coefficients
                              std.errors_updated[!(names(std.errors_updated) %in% names(std.errors))] <- NA_integer_
                              std.errors_updated[names(std.errors_updated) %in% names(std.errors)] <- std.errors
                              # Update private std. errors
                              private$..std.errors <- std.errors_updated
                            },
                            #' @description
                            #' Print method for `AR_modelR6` class.
                            print = function(){
                              cat(paste0("--------------------- ", "ARMA", "(", self$order[1], ", ", self$order[2], ")", "--------------------- \n"))
                              cat(paste0("intercept: ",  paste0(format(self$intercept, digits = 4), collapse = " "), "\n"))
                              cat(paste0("phi: ", paste0(format(self$phi, digits = 4), collapse = " "), "\n"))
                              cat(paste0("theta: ", paste0(format(self$theta, digits = 4), collapse = " "), "\n"))
                              cat(paste0("------------------------------------------------\n"))
                              cat(paste0("Intercept: ", private$include.intercept, "\n"))
                              cat(paste0("Version: ", private$version, "\n"))
                            }
                          ),
                          private = list(
                            version = "1.0.0",
                            ..model = NA,
                            ..arOrder = 0,
                            ..maOrder = 0,
                            ..intercept = 0,
                            ..phi = c(),
                            ..theta = c(),
                            ..std.errors = NA,
                            include.intercept = FALSE
                          ),
                          active = list(
                            #' @field intercept Numeric named scalar. Intercept.
                            intercept = function(){
                              private$..intercept
                            },
                            #' @field phi Numeric named vector. AR parameters.
                            phi = function(){
                              private$..phi
                            },
                            #' @field theta Numeric named vector. MA parameters.
                            theta = function(){
                              private$..theta
                            },
                            #' @field coefficients Numeric named vector. Intercept and ARMA parameters.
                            coefficients = function(){
                              c(self$intercept, self$phi, self$theta)
                            },
                            #' @field order Numeric named vector. ARMA order.
                            order = function(){
                              c(AR = private$..arOrder, MA = private$..maOrder)
                            },
                            #' @field mean Numeric scalar. Long term expectation.
                            mean = function(){
                              self$intercept / (1 - sum(self$phi))
                            },
                            #' @field variance Numeric scalar. Long term variance.
                            variance = function(){
                              var <- NA
                              phi <- self$phi
                              sigma2 <- 1
                              if (self$order[2] == 0 & self$order[1] <= 4) {
                                if (length(phi) == 1) {
                                  var <- sigma2/(1-phi[1]^2)
                                } else if (length(phi) == 2) {
                                  var <- sigma2*(1-phi[2])/((1 - phi[2])*(1 - phi[1]^2 - phi[2]^2) - 2*phi[1]^2*phi[2])
                                } else if (length(phi) == 3) {
                                  phi_tilde_1 <- (phi[1] + phi[2]*phi[3])/(1 - phi[2] - phi[3]^2 - phi[1]*phi[3])
                                  phi_tilde_2 <- (phi[1] + phi[3])*phi_tilde_1 + phi[2]
                                  phi_tilde_3 <- (phi[1]*phi_tilde_2 + phi[2]*phi_tilde_1 + phi[3])
                                  phi_tilde_0 <- 1/(1 - phi[1]*phi_tilde_1 - phi[2]*phi_tilde_2 - phi[3]*phi_tilde_3)
                                  var <- phi_tilde_0*sigma2
                                } else if (length(phi) == 4) {
                                  phi_1 <- (phi[1] + phi[3])/(1 - phi[4])
                                  phi_0 <- phi[2]/(1 - phi[4])
                                  psi_1 <- (phi_1*phi[3] + phi[1]*phi_0*phi[4]+ phi[1] + phi[3]*phi[4])
                                  psi_1 <- psi_1/(1-phi_1*(phi[3] + phi[1]*phi[4]) - phi[2]*(1 + phi[4]) - phi[4]^2)
                                  psi_2 <- phi_1*psi_1 + phi_0
                                  psi_3 <- phi[1]*psi_2 + phi[2]*psi_1 + phi[4]*psi_1 + phi[3]
                                  psi_4 <- phi[1]*psi_3 + phi[2]*psi_2 + phi[3]*psi_1 + phi[4]
                                  var <- sigma2/(1 - phi[1]*psi_1 - phi[2]*psi_2 - phi[3]*psi_3 - phi[4]*psi_4)
                                }
                              } else {
                                # Numerical computation
                                A <- self$Phi
                                pow_sum <- sigma2
                                for(i in 1:100){
                                  A_j <- pow_matrix(A, i)
                                  pow_sum <- pow_sum + A_j %*% self$b %*% t(self$b) %*% t(A_j)
                                }
                                var <- pow_sum[1,1] * sigma2
                              }
                              return(var)
                            },
                            #' @field model Fitted ARMA model from the function `arima`.
                            model = function(){
                              private$..model
                            },
                            #' @field Phi Numeric, matrix. Companion matrix.
                            Phi = function(){
                              par <- self$coefficients[-1]
                              # Max order
                              k <- max(self$order)
                              # Vector for shifts
                              v_shift <- rep(0, sum(self$order))

                              j <- 1
                              # Companion matrix
                              A <- matrix(par, nrow = 1, ncol = sum(self$order))
                              if (ncol(A) > 1) {
                                while(j < ncol(A)){
                                  m <- v_shift
                                  if (j < self$order[1] & self$order[1] > 0){
                                    m[j] <- 1
                                  } else if (j > self$order[1]){
                                    m[j] <- 1
                                  } else {
                                    m[j] <- 0
                                  }
                                  A <- rbind(A, m)
                                  j <- j + 1
                                }
                              }
                              rownames(A) <- NULL
                              colnames(A) <- names(par)
                              return(A)
                            },
                            #' @field b Numeric, vector. Vector for matrix form of the residuals.
                            b = function(){
                              # Companion vector AR
                              if (self$order[1] > 1){
                                b_ar <- c(1, rep(0, self$order[1]-1))
                              } else {
                                b_ar <- 1
                              }
                              # Companion vector MA
                              if (self$order[2] == 1) {
                                b_ma <- 1
                              } else if (self$order[2] > 1) {
                                b_ma <- c(1, rep(0, self$order[2]-1))
                              } else {
                                b_ma <- c()
                              }
                              # Combine the two parts
                              b <- c(b_ar, b_ma)
                              return(b)
                            },
                            #' @field tidy Method tidy for the estimated parameters
                            tidy = function(){
                              dplyr::tibble(
                                term = names(private$..model$coef),
                                estimate = private$..model$coef,
                                std.error = private$..std.errors
                              )
                            }
                          ))
