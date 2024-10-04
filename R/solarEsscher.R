#' Calibrate Esscher Bounds and parameters
#'
#' @rdname solarEsscher
#' @name solarEsscher
#' @export
solarEsscher <- R6::R6Class("solarEsscher",
                            public = list(
                              #' @field control list containing the control parameters
                              control = list(),
                              #' @field grid list containing the grids
                              grid = list(Yt = NA, Rt = NA),
                              #' @description
                              #' Initialize the settings for calibration of Esscher parameter.
                              #'
                              #' @param n_key_points integer, number of key points for interpolation.
                              #' @param init_lambda numeric, initial value for the Esscher parameter.
                              #' @param lower_lambda numeric, lower value for the Esscher parameter.
                              #' @param upper_lambda numeric, upper value for the Esscher parameter.
                              #' @param quiet logical
                              #' @param put logical, when `TRUE` will be considered a put contract otherwise a call contract.
                              #' @param target.Yt logical, when `TRUE` will be distorted with esscher parameter the pdf of Yt otherwise the pdf of the GHI.
                              #' @param control_options control function. See \code{\link{control_solarOption}} for details.
                              initialize = function(n_key_points = 15, init_lambda = 0, lower_lambda = -1, upper_lambda = 1,
                                                    put = TRUE, target.Yt = TRUE, quiet = FALSE, control_options = control_solarOption()) {
                                self$control$n_key_points = n_key_points
                                self$control$init_lambda = init_lambda
                                self$control$lower_lambda = lower_lambda
                                self$control$upper_lambda = upper_lambda
                                self$control$quiet = quiet
                                self$control$put = put
                                self$control$target.Yt = target.Yt
                                self$control$control_options = control_options
                              },
                              #' @description
                              #' Calibrate the optimal Esscher parameter given a target price
                              #'
                              #' @param model solar model
                              #' @param target_price the `target_price` represent the model price under the target Q-measure.
                              #' @param nmonths month or months
                              #' @param target.Yt logical, when `TRUE` will be distorted with esscher parameter the pdf of Yt otherwise the pdf of the GHI.
                              calibrator = function(model, target_price, nmonths = 1:12, target.Yt){
                                # Esscher control parameters
                                put = self$control$put
                                target.Yt = ifelse(missing(target.Yt), self$control$target.Yt, target.Yt)
                                control_options = self$control$control_options
                                lower_lambda = self$control$lower_lambda
                                upper_lambda = self$control$upper_lambda
                                init_lambda = self$control$init_lambda
                                quiet = self$control$quiet
                                # Generator of loss functions for Esscher theta
                                loss_generator <- function(model, nmonth, target_price, target.Yt, put, control_options, quiet){
                                  function(theta) {
                                    # Model premium for months "nmonth"
                                    model_price <- solarOption_model(model, theta = theta, put = put, target.Yt = target.Yt, nmonths = nmonth, control_options = control_options)
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

                                if (length(nmonths) > 1 & length(target_price) == 1) {
                                  # Specify the loss for a month
                                  loss <- loss_generator(model, nmonths, target_price, target.Yt, put, control_options, quiet)
                                  if (!quiet) cat(paste0("\033[1;35m---------------\033[0m", " Calibrating Yearly Esscher parameter ", "\033[1;32m", "\033[1;35m---------------\033[0m \n"))
                                  # Optimal Esscher parameter
                                  par <- calibrator(loss)$par
                                } else {
                                  par <- c()
                                  for(nmonth in nmonths){
                                    flush.console()
                                    message("", appendLF = TRUE)
                                    # Test tolerance parameter
                                    if (!quiet) cat(paste0("\033[1;35m---------------\033[0m", " Calibrating Monthly Esscher parameter (", "\033[1;32m",
                                                           lubridate::month(nmonth, label = TRUE, abbr = FALSE), "\033[0m", ")",  " \033[1;35m---------------\033[0m \n"))
                                    # Specify the loss for a month
                                    loss <- loss_generator(model, nmonth, target_price[nmonth], target.Yt, put, control_options, quiet)
                                    # Optimal Esscher parameter
                                    opt_par <- calibrator(loss)$par
                                    names(opt_par) <- lubridate::month(nmonth, label = TRUE)
                                    par <- c(par, opt_par)
                                  }
                                }
                                return(par)
                              },
                              #' @description
                              #' Calibrate Esscher upper and lower bounds
                              #'
                              #' @param model object with the class `solarModel`. See the function \code{\link{solarModel}} for details.
                              #' @param payoffs object with the class `solarOptionPayoffs`. See the function \code{\link{solarOptionPayoffs}} for details.
                              #' @param target.Yt logical, when `TRUE` will be distorted with esscher parameter the pdf of Yt otherwise the pdf of the GHI.
                              calibrate_bounds = function(model, payoffs, target.Yt){
                                # Match the type of contract
                                option_type <- ifelse(self$control$put, "put", "call")
                                target.Yt <- self$control$target.Yt
                                # Extract payoffs
                                # Bootstrapped payoff (historical)
                                payoff_boot <- payoffs[[option_type]]$model$boot
                                # Historical payoff (realized)
                                payoff_hist <- payoffs[[option_type]]$historical
                                # Model payoff
                                payoff_model_P <- payoffs[[option_type]]$model$P

                                # ********* Calibrate the optimal Esscher theta for the upper and lower bounds *********
                                # Default benchmark for worst case price (buyer point of view)
                                benchmark_Qup <- payoff_boot$payoff_year$premium_up
                                # Optimal upper Esscher parameter (up, worse case)
                                opt_theta_up <- self$calibrator(model, benchmark_Qup, 1:12, target.Yt)
                                # Optimal up and down parameters
                                theta <- solarEsscher_theta_bounds(opt_theta_up)[[1]]
                                # Add Bounds
                                if (target.Yt) {
                                  private$..bounds$Yt <- list(Qdw = theta$dw, Q = 0.5*(theta$up + theta$dw), Qup = theta$up)
                                } else {
                                  private$..bounds$Rt <- list(Qdw = theta$dw, Q = 0.5*(theta$up + theta$dw), Qup = theta$up)
                                }
                              },
                              #' @description
                              #' Create a grid of optimal theta and expected returns with respect of the benchmark price.
                              #'
                              #' @param model object with the class `solarModel`. See the function \code{\link{solarModel}} for details.
                              #' @param benchmark_price benchmark price for an expected return equal to zero.
                              #' @param lower_price lower price in the grid.
                              #' @param upper_price upper price in the grid.
                              #' @param target.Yt logical, when `TRUE` will be distorted with esscher parameter the pdf of Yt otherwise the pdf of the GHI.
                              create_grid = function(model, benchmark_price, lower_price, upper_price, target.Yt){
                                # Esscher control parameters
                                option_side = ifelse(self$control$put, "put", "call")
                                n_key_points = self$control$n_key_points
                                target.Yt = ifelse(missing(target.Yt), self$control$target.Yt, target.Yt)

                                # Grid of prices inside the Esscher corridor (decreasing)
                                grid_prices <- seq(upper_price, lower_price, length.out = n_key_points)
                                # Grid of expected returns (decreasing)
                                grid_rates <- benchmark_price/grid_prices - 1
                                # Implied grid of Esscher lambda
                                df_grid <- dplyr::tibble(prices = grid_prices, rates = grid_rates, theta = NA)
                                # Add the theta correspondent to lower and upper price if specified
                                key_points <- 1:n_key_points
                                if (!purrr::is_empty(key_points)) {
                                  for(i in key_points){
                                    if (!self$control$quiet) message(rep("-", 20), " ", i, "/", n_key_points, " ", rep("-", 20))
                                    df_grid$theta[i] <- self$calibrator(model, df_grid$prices[i], 1:12, target.Yt)
                                    if (!self$control$quiet) message("r: ", df_grid$rates[i], " h: ", df_grid$theta[i])
                                  }
                                }
                                # Add manually r = 0
                                idx_pos <- which(grid_rates > 0)
                                idx_neg <- which(grid_rates < 0)
                                df_grid <- dplyr::bind_rows(df_grid[idx_neg,],
                                                            dplyr::tibble(prices = benchmark_price, rates = 0, theta = 0),
                                                            df_grid[idx_pos,])
                                df_grid <- dplyr::bind_cols(side = option_side, target.Yt = target.Yt, df_grid)
                                if (target.Yt) {
                                  self$grid$Yt <- df_grid
                                } else {
                                  self$grid$Rt <- df_grid
                                }
                              },
                              #' @description
                              #' Fit the models to predict the optimal Esscher parameters given the grid.
                              fit_theta = function(){
                                if (is.na(self$grid[[1]]) || is.na(self$grid[[2]])){
                                  stop("The slots `grid are currently empty`! Create a grid first!")
                                }
                                # Extract the grid
                                grid <- self$grid
                                # Create a model for both types of theta
                                theta_esscher_Yt <- lm(theta ~ rates + I(rates^2), data = grid$Yt)
                                theta_esscher_Rt <- lm(theta ~ rates + I(rates^2), data = grid$Rt)
                                # Create a dataset with both theta to connect the models
                                grid_theta <- dplyr::bind_cols(theta_Yt = grid$Yt$theta, theta_Rt = grid$Rt$theta)
                                theta_esscher_Yt_Rt <- lm(theta_Yt ~ theta_Rt + I(theta_Rt^2), data = grid_theta)
                                theta_esscher_Rt_Yt <- lm(theta_Rt ~ theta_Yt + I(theta_Yt^2), data = grid_theta)
                                # Output
                                private$..models$theta_Yt <- theta_esscher_Yt
                                private$..models$theta_Rt <- theta_esscher_Rt
                                private$..models$theta_Yt_Rt <- theta_esscher_Yt_Rt
                                private$..models$theta_Rt_Yt <- theta_esscher_Rt_Yt

                                if (self$control$target.Yt) {
                                  theta_Rt <- stats::predict.lm(theta_esscher_Rt_Yt, newdata = data.frame(theta_Yt = unlist(self$bounds$Yt)))
                                  private$..bounds$Rt <- list(Qdw = theta_Rt[1], Q = theta_Rt[2], Qup = theta_Rt[3])
                                } else {
                                  theta_Yt <- stats::predict.lm(theta_esscher_Yt_Rt, newdata = data.frame(theta_Rt = unlist(self$bounds$Rt)))
                                  private$..bounds$Yt <- list(Qdw = theta_Yt[1], Q = theta_Yt[2], Qup = theta_Yt[3])
                                }
                              },
                              #' @description
                              #' Predict the optimal Esscher parameters given a certain level of expected return.
                              #' @param r expected return
                              #' @param target.Yt logical, when `TRUE` will be distorted with esscher parameter the pdf of Yt otherwise the pdf of the GHI.
                              predict = function(r, target.Yt = FALSE){
                                theta <- stats::predict.lm(private$..models$theta_Rt, newdata = data.frame(rates = r))
                                if (target.Yt) {
                                  theta <- stats::predict.lm(private$..models$theta_Yt_Rt, newdata = data.frame(theta_Rt = theta))
                                }
                                return(theta)
                              }
                            ),
                            private = list(
                              ..bounds = list(Rt = list(Qdw = NA, Q = NA, Qup = NA), Yt = list(Qdw = NA, Q = NA, Qup = NA)),
                              ..models = list(theta_Yt = NA, theta_Rt = NA, theta_Yt_Rt = NA, theta_Rt_Yt = NA)
                            ),
                            active = list(
                              #' @field bounds calibrated bounds with respect to bootstrapped payoff.
                              bounds = function(){
                                private$..bounds
                              },
                              #' @field models models to predict the optimal theta given the expected return.
                              models = function(){
                                private$..models
                              }
                            )
                          )






