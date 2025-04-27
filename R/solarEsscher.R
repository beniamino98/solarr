#' Function to Calibrate Esscher Bounds and parameters
#'
#' @rdname solarEsscher
#' @name solarEsscher
#' @export
solarEsscher <- R6::R6Class("solarEsscher",
                            public = list(
                              #' @field control list containing the control parameters
                              control = list(),
                              #' @field grid list containing the grids
                              grid = NA,
                              #' @field theta_up list containing the grids
                              theta_up = NA,
                              #' @field theta_bar list containing the grids
                              theta_bar = NA,
                              #' @field theta_dw list containing the grids
                              theta_dw = NA,
                              #' @description
                              #' Initialize the settings for calibration of Esscher parameter.
                              #'
                              #' @param n_key_points integer, number of key points used to create the grid for interpolation.
                              #' @param init_lambda numeric, initial value for the Esscher parameter.
                              #' @param lower_lambda numeric, lower value for the Esscher parameter.
                              #' @param upper_lambda numeric, upper value for the Esscher parameter.
                              #' @param quiet logical
                              #' @param put logical, when `TRUE` will be considered a put contract otherwise a call contract.
                              #' @param control_options control function. See \code{\link{control_solarOption}} for details.
                              initialize = function(n_key_points = 15, init_lambda = 0, lower_lambda = -1, upper_lambda = 1,
                                                    put = TRUE, quiet = FALSE, control_options = control_solarOption()) {
                                self$control$n_key_points = n_key_points
                                self$control$init_lambda = init_lambda
                                self$control$lower_lambda = lower_lambda
                                self$control$upper_lambda = upper_lambda
                                self$control$quiet = quiet
                                self$control$put = put
                                self$control$control_options = control_options
                              },
                              #' @description
                              #' Calibrate the optimal Esscher parameter given a target price
                              #'
                              #' @param model solar model
                              #' @param target_price the `target_price` represent the model price under the target Q-measure.
                              #' @param nmonths month or months
                              calibrate_theta = function(model, mom, target_price){
                                # Esscher control parameters
                                put = self$control$put
                                control_options = self$control$control_options
                                lower_lambda = self$control$lower_lambda
                                upper_lambda = self$control$upper_lambda
                                init_lambda = self$control$init_lambda
                                quiet = self$control$quiet
                                # Generator of loss functions for Esscher theta
                                loss_generator <- function(model, mom, target_price, put, control_options, quiet){
                                  function(theta) {
                                    # Model premium
                                    model_price <- solarOption_model(model, mom, theta = theta, put = put, control_options = control_options)
                                    # Compute the difference from target price
                                    price_difference <- abs(model_price$payoff_year$premium - target_price)
                                    # Expected return in absolute value
                                    expected_return <- price_difference/target_price
                                    if (!quiet) message("Price difference: ", price_difference,
                                                        " Expected return: |", format(expected_return*100, digits = 4),
                                                        "| Lambda: ", format(theta, digits = 5), "\r", appendLF = FALSE)
                                    flush.console()
                                    return(price_difference^2)
                                  }
                                }
                                # Calibrator function
                                calibrator <- function(loss){
                                  # Optimal Esscher parameter
                                  opt <- optim(par = init_lambda, loss, method = "Brent", lower = lower_lambda, upper = upper_lambda)
                                  return(opt)
                                }
                                # Specify the loss for a month
                                loss <- loss_generator(model, mom, target_price, put, control_options, quiet)
                                if (!quiet) cat(paste0("\033[1;35m---------------\033[0m",
                                                       " Calibrating Yearly Esscher parameter ",
                                                       "\033[1;32m",
                                                       "\033[1;35m---------------\033[0m \n"))
                                # Optimal Esscher parameter
                                par <- calibrator(loss)$par
                                return(par)
                              },
                              #' @description
                              #' Create a grid of optimal theta and expected returns with respect of the benchmark price.
                              #' Fit the model to predict the optimal Esscher parameters given the grid.
                              #' @param model object with the class `solarModel`. See the function \code{\link{solarModel}} for details.
                              #' @param benchmark_price benchmark price for an expected return equal to zero.
                              #' @param lower_price lower price in the grid.
                              #' @param upper_price upper price in the grid.
                              create_grid = function(model, mom, benchmark_price, lower_price, upper_price){
                                # Esscher control parameters
                                option_side = ifelse(self$control$put, "put", "call")
                                n_key_points = self$control$n_key_points
                                # Grid of prices inside the Esscher corridor (decreasing)
                                grid_prices <- seq(upper_price, lower_price, length.out = n_key_points)
                                # Grid of expected returns (decreasing)
                                grid_rates <- benchmark_price/grid_prices - 1
                                # Implied grid of Esscher lambda
                                grid <- dplyr::tibble(prices = grid_prices, rates = grid_rates, theta = NA)
                                # Add the theta correspondent to lower and upper price if specified
                                key_points <- 1:n_key_points
                                if (!purrr::is_empty(key_points)) {
                                  for(i in key_points){
                                    if (!self$control$quiet) message(rep("-", 20), " ", i, "/", n_key_points, " ", rep("-", 20))
                                    grid$theta[i] <- self$calibrate_theta(model, mom, grid$prices[i])
                                    if (!self$control$quiet) message("r: ", grid$rates[i], " h: ", grid$theta[i])
                                  }
                                }
                                # Add manually r = 0
                                idx_pos <- which(grid_rates > 0)
                                idx_neg <- which(grid_rates < 0)
                                grid <- dplyr::bind_rows(grid[idx_neg,],
                                                         dplyr::tibble(prices = benchmark_price, rates = 0, theta = 0),
                                                         grid[idx_pos,])
                                self$grid <- dplyr::bind_cols(side = option_side, grid)
                                # Fit and store the model
                                private$..model <- lm(theta ~ rates + I(rates^2), data = self$grid)
                              },
                              #' @description
                              #' Predict the optimal Esscher parameters given a certain level of expected return.
                              #' @param r expected return
                              theta = function(r){
                                stats::predict.lm(self$model, newdata = data.frame(rates = r))
                              },
                              #' @description
                              #' Print method for `solarEsscher` class.
                              print = function(){
                                # Complete data specifications
                                data <- self$dates$data
                                cat(paste0("--------------------- ", "esscherBounds", "--------------------- \n"))
                                cat(paste0("Bounds: ", self$theta_dw, " - ", self$theta_up, "\n"))
                                cat(paste0("Q-measure: ", self$theta_bar, "\n"))
                                cat(paste0("Version: ", private$version, "\n"))
                              }
                            ),
                            private = list(
                              version = "1.0.0",
                              ..model = NA
                            ),
                            active = list(
                              #' @field model Models to predict the optimal theta given the expected return.
                              model = function(){
                                private$..model
                              }
                            )
)
