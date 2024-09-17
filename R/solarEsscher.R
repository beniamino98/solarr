#' Control for Esscher calibration.
#'
#' Control parameters for calibration of Esscher parameters.
#'
#' @param nsim integer, number of simulations used to bootstrap the premium bounds.
#' @param ci numeric, confidence interval for bootstrapping. See `solar_option_payoff_bootstrap()`.
#' @param seed integer, random seed for reproducibility.
#' @param n_key_points integer, number of key points for interpolation.
#' @param init_lambda numeric, initial value for the Esscher parameter.
#' @param lower_lambda numeric, lower value for the Esscher parameter.
#' @param upper_lambda numeric, upper value for the Esscher parameter.
#' @param quiet logical
#'
#' @rdname control_solarEsscher
#' @name control_solarEsscher
#' @export
control_solarEsscher <- function(nsim = 200, ci = 0.05, seed = 1, n_key_points = 15,
                                 init_lambda = 0, lower_lambda = -1, upper_lambda = 1, quiet = FALSE){
  stopifnot(n_key_points >= 2)
  list(
    nsim = nsim,
    ci = ci,
    seed = seed,
    n_key_points = n_key_points,
    init_lambda = init_lambda,
    lower_lambda = lower_lambda,
    upper_lambda = upper_lambda,
    quiet = quiet
  )
}
#' Calibrate monthly Esscher parameter given the expected return
#'
#' Calibrator function for the monthly Esscher parameter of a solarOption given a
#' desired level of expected return at maturity.
#'
#' @param model solar model
#' @param nmonth month
#' @param expected_return expected return at maturity. The benchmark for the `target_price` to match will be the
#' mean cumulated net payoff on the last day of the month plus the model price paid under the Esscher measure.
#' The return of the `target_price` with respect to the model price will match the parameter `expected_return`. For example `0.01`, `0.02`, ecc.
#' @param target_price alternative to the `expected_return` parameter. Submitting a `target_price` will imply that the `expected_return = 0` so that
#' the model price under the Esscher measure matches the `target_price`
#' @param control_esscher control
#' @param control_options control
#'
#' @rdname solarEsscher_calibrator_month
#' @name solarEsscher_calibrator_month
#' @export
solarEsscher_calibrator_month <- function(model, nmonth = 1, expected_return = 0, target_price = NA, control_esscher = control_solarEsscher(), control_options = control_solarOption()){

  # Esscher control parameters
  lower_lambda = control_esscher$lower_lambda
  upper_lambda = control_esscher$upper_lambda
  init_lambda = control_esscher$init_lambda
  quiet = control_esscher$quiet

  # Compute realized historical payoffs
  payoff_hist <- solarOption_historical(model, nmonths = 1:12, control_options = control_options)$payoff

  # Calibrator function for Esscher theta
  calibrator <- function(theta, nmonth, expected_return = 0, target_price = NA) {
    # Model premium for month "nmonth"
    model_price <- solarOption_model(model, theta = theta, nmonths = nmonth, control_options = control_options)
    # Extract model daily premium
    payoff_model <- dplyr::select(model_price$payoff_month_day, Month, Day, premium)
    # Extract model monthly premium
    model_price <- sum(model_price$payoff_month$premium)
    # Default target price
    if (is.na(target_price)) {
      # Compute the expected cumulated net payoff on historical data with model premium
      net_payoff <- dplyr::left_join(payoff_hist, payoff_model, by = c("Month", "Day")) %>%
        dplyr::filter(Month == nmonth) %>%
        dplyr::group_by(Year) %>%
        dplyr::mutate(net_payoff = cumsum(payoff - premium)) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(Month, Day) %>%
        dplyr::mutate(e_net_payoff = mean(net_payoff, na.rm = TRUE))
      # Extract the expected net cumulated payoff from the last day of the month
      day_max <- lubridate::days_in_month(net_payoff$Month[1])
      e_net_payoff <- dplyr::filter(net_payoff, Day == day_max)$e_net_payoff[1]
      # Compute the target price
      target_price <- model_price + e_net_payoff
    } else {
      expected_return <- 0
    }
    # Compute expected return with respect to the target price
    e_return <- target_price/model_price
    # Loss to match the expected return at maturity
    loss <- abs((e_return - expected_return)^2 - 1)
    if (!quiet) message("Loss: ", loss, " Expected return: ", format((e_return-1)*100, digits = 4), " Lambda: ", format(theta, digits = 5))
    return(loss)
  }

  # Optimal Esscher parameter
  opt <- optim(par = init_lambda, calibrator, nmonth = nmonth,
               expected_return = expected_return, target_price = target_price,
               method = "Brent", lower = lower_lambda, upper = upper_lambda)

  return(opt)
}

#' Calibrator for Esscher parameter
#'
#' @param model object with the class `solarModel`. See the function \code{\link{solarModel}} for details.
#' @param target_price target price for the calibration.
#' @param nmonths months used in the model computation.
#' @param control_options control function. See \code{\link{control_solarOption}} for details.
#' @param control_esscher control function. See \code{\link{control_solarEsscher}} for details.
#'
#' @rdname solarEsscher_calibrator
#' @name solarEsscher_calibrator
#'
#' @export
solarEsscher_calibrator <- function(model, target_price, nmonths = 1:12, control_esscher = control_solarEsscher(), control_options = control_solarOption()){

  # Bounds for parameters
  lower_lambda = control_esscher$lower_lambda
  upper_lambda = control_esscher$upper_lambda
  # Initial theta parameter
  init_lambda = control_esscher$init_lambda
  # Verbose parameter
  quiet = control_esscher$quiet

  # Loss function
  loss_function <- function(theta, target_price) {
    # Price under model
    model_price <- solarOption_model(model, theta = theta, nmonths = nmonths, control_options = control_options)$payoff_year$premium
    # Loss function
    loss <- (model_price - target_price)^2
    if (!quiet) message("Loss: ", loss, " Theta: ", theta)
    return(loss)
  }
  opt <- optim(par = init_lambda, loss_function, target_price = target_price,
               method = "Brent", lower = lower_lambda, upper = upper_lambda)
  return(opt)
}

#' Function to establish up and down parameters of the Esscher transform
#' A positive theta identify a `down` parameter, a negative theta identify an `up` parameter.
#'
#' @param theta Esscher parameter.
#' @return A `list` with first element named `up` with the positive parameter and second element named `dw` with the negative one.
#'
#' @keywords internals
#' @noRd
#' @export
solarEsscherTheta <- function(theta){
  # dw parameter = positive theta
  # up parameter = negative theta
  if (theta >= 0) {
    par_dw <-  theta
    par_up <- -theta/(1+theta)
  } else {
    par_dw <- -theta/(1-theta)
    par_up <-  theta
  }
  params <- list(up = par_up, dw = par_dw)
  return(params)
}

#' solarEsscherFunction
#'
#' @keywords internals
#' @noRd
#' @export
solarEsscherFunction <- function(par = c(h1 = 1, h2 = 1, h3 = 1), imp_r = NA){

  is_r_implied <- !is.na(imp_r[1])
  if (is_r_implied) {
    stopifnot(length(imp_r) == 12)
  }

  # Parametric functions for the h1 e h2
  h1_tilde <- function(r_imp){par[1] + par[2]*r_imp + par[3]*r_imp^2}
  h2_tilde <- function(r_imp){2*par[3]*r_imp + par[2]}

  if (is_r_implied) {
    # Parametric function for monthly Esscher theta
    function(r){
      function(nmonth) {
        h1_tilde(imp_r[nmonth]) + h2_tilde(imp_r[nmonth])*r + par[3]*r^2
      }
    }
  } else {
    # Parametric function for Esscher theta
    function(r){
      h1_tilde(0) + h2_tilde(0)*r + par[3]*r^2
    }

  }
}


#' Calibrate Esscher Bounds and parameters
#'
#' @param model object with the class `solarModel`. See the function \code{\link{solarModel}} for details.
#' @param control_options control function. See \code{\link{control_solarOption}} for details.
#' @param control_esscher control function. See \code{\link{control_solarEsscher}} for details.
#'
#' @rdname solarEsscher_bounds
#' @name solarEsscher_bounds
#' @export
solarEsscher_bounds <- function(model, control_options = control_solarOption(), control_esscher = control_solarEsscher()){

  # Bootstrap parameters
  nsim = control_esscher$nsim
  ci = control_esscher$ci
  seed = control_esscher$seed
  # Grid controls
  n_key_points = control_esscher$n_key_points
  # Verbose parameter
  quiet = control_esscher$quiet

  # Compute payoffs
  # Fair Payoff bootstrapped (historical)
  payoff_boot <- solarOption_bootstrap(model, nsim = nsim, ci = ci, seed = seed, control_options = control_options)
  # Historical payoff computed on data (historical)
  payoff_hist <- solarOption_historical(model, control_options = control_options)
  # P-Payoff computed with model
  payoff_model_P <- solarOption_model(model, theta = 0, control_options = control_options)
  # Discrepance between historical and model P-price
  error_hist_P <- payoff_hist$payoff_year$premium - payoff_model_P$payoff_year$premium
  # 1) ------- Calibrate the optimal theta for the bounds -------
  # Default benchmark for worste price
  benchmark_Qup <- payoff_boot$payoff_year$premium_up + error_hist_P
  # Optimal upper Esscher parameter (up, worse case)
  opt_lambda_up <- solarEsscher_calibrator(model, target_price = benchmark_Qup, control_esscher = control_esscher, control_options = control_options)
  # Optimal up and down parameters
  lambda_up <- solarEsscherTheta(opt_lambda_up$par)$up
  lambda_dw <- solarEsscherTheta(opt_lambda_up$par)$dw
  lambda_bar <- (lambda_up+lambda_dw)/2
  # Compute payoffs
  # Qdw-payoff: best case scenario (model)
  payoff_model_Qdw <- solarOption_model(model, theta = lambda_dw, control_options = control_options)
  # Qup-payoff: worste case scenario (model)
  payoff_model_Qup <- solarOption_model(model, theta = lambda_up, control_options = control_options)
  # Q-payoff as average between lambda up and down (model)
  payoff_model_Q <- solarOption_model(model, theta = lambda_bar, control_options = control_options)

  # 3) ------- Calibrate the Esscher parameter for different levels of r -------
  # Benchmark price with r = 0
  benchmark_Q <- payoff_model_Q$payoff_year$premium
  # Benchmark for best case prices with r > 0
  lower <- payoff_model_Qdw$payoff_year$premium
  # Benchmark for worste case prices with r < 0
  upper <- payoff_model_Qup$payoff_year$premium
  # Grid of prices in the Esscher corridor (decreasing)
  grid_prices <- seq(upper, lower, length.out = n_key_points)
  # Implied grid of risk-free rates (decreasing)
  grid_rates <- benchmark_Q/grid_prices - 1
  # Implied grid of Esscher lambda
  df_grid <- dplyr::tibble(prices = grid_prices, rates = grid_rates, theta = NA)
  df_grid$theta[1] <- lambda_up
  df_grid$theta[n_key_points] <- lambda_dw

  key_points <- c(1:n_key_points)[-c(1, n_key_points)]
  if(!purrr::is_empty(key_points)){
    for(i in key_points){
      if (!quiet) message(rep("-", 20), " ", i, "/", n_key_points, " ", rep("-", 20))
      df_grid$theta[i] <- solarEsscher_calibrator(model, target_price = df_grid$prices[i],
                                                  control_esscher = control_esscher, control_options = control_options)$par
      if (!quiet) message("r: ", df_grid$rates[i], " h: ", df_grid$theta[i])
    }
  }

  # Add manually r = 0
  idx_pos <- which(grid_rates > 0)
  idx_neg <- which(grid_rates < 0)
  df_grid <- dplyr::bind_rows(df_grid[idx_neg,],
                              dplyr::tibble(prices = benchmark_Q, rates = 0, theta = lambda_bar),
                              df_grid[idx_pos,])

  # 4) ------- Fit a parametric model for the Esscher parameter -------
  # theta = h0 + h1*r + h2*r^2 model in terms of risk-free rates
  model_param_esscher <- lm(theta ~ rates + I(rates^2), data = df_grid)
  # Fitted values
  df_grid$lambda_fit <- predict(model_param_esscher, newdata = df_grid)
  # Fitted parameters
  esscher_params <- as.list(model_param_esscher$coefficients)
  names(esscher_params) <- c("h0","h1","h2")
  # Add Esscher parametric model
  model$esscher <- list(
    grid = df_grid,
    params = list(Qdw = lambda_dw, Q = lambda_bar, Qup = lambda_up),
    model = model_param_esscher,
    coefficients = esscher_params,
    theta = solarEsscherFunction(par = unlist(esscher_params), imp_r = NA),
    theta_m = NA,
    implied_r = NA,
    control = control_esscher
  )
  # Add Historical payoffs to the model
  model$payoffs$hist <- payoff_hist
  # Compute implied monthly returns
  implied_expected_returns <- solarOption_implied_return(model, nmonths = 1:12, control_options = control_options)$implied_r
  # Function for Optimal Esscher parameters for monthly case
  model$esscher$theta_m <- solarEsscherFunction(par = unlist(esscher_params), imp_r = implied_expected_returns)
  model$esscher$implied_r <- implied_expected_returns

  # Add payoffs to the model
  # Bootstrapped
  model$payoffs$boot <- payoff_boot
  # Model
  model$payoffs$model$P <- payoff_model_P
  model$payoffs$model$Q <- payoff_model_Q
  model$payoffs$model$Qdw <- payoff_model_Qdw
  model$payoffs$model$Qup <- payoff_model_Qup
  return(model)
}

