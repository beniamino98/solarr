#' Control for Esscher calibration.
#'
#' Control parameters for calibration of Esscher parameters.
#'
#' @param nsim integer, number of simulations used to bootstrap the premium bounds.
#' @param ci numeric, confidence interval for bootstrapping. See \code{\link{solarOption_bootstrap}}.
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
  structure(
    list(
      nsim = nsim,
      ci = ci,
      seed = seed,
      n_key_points = n_key_points,
      init_lambda = init_lambda,
      lower_lambda = lower_lambda,
      upper_lambda = upper_lambda,
      quiet = quiet
    ),
    class = c("control", "list")
  )
}


#' Calibrate an Esscher parameter given a target price
#'
#' Calibrator function for the monthly Esscher parameter of a solarOption
#'
#' @param model solar model
#' @param nmonths month or months
#' @param target_price the `target_price` represent the model price under the target Q-measure.
#' @param control_esscher control
#' @param control_options control
#'
#' @examples
#' model <- Bologna
#' # Compute realized historical payoffs
#' payoff_hist <- solarOption_historical(model, nmonths = 1:12)
#' # Monthly calibration
#' solarEsscher_calibrator(model, 1:3, payoff_hist$payoff_month$premium[1:3])
#' # Yearly calibration
#' solarEsscher_calibrator(model, 1:12, payoff_hist$payoff_year$premium)
#'
#' @rdname solarEsscher_calibrator
#' @name solarEsscher_calibrator
#' @export
solarEsscher_calibrator <- function(model, nmonths = 1, target_price, control_esscher = control_solarEsscher(), control_options = control_solarOption()){

  # Esscher control parameters
  lower_lambda = control_esscher$lower_lambda
  upper_lambda = control_esscher$upper_lambda
  init_lambda = control_esscher$init_lambda
  quiet = control_esscher$quiet

  # Loss function for Esscher theta
  loss_function <- function(model, nmonth, target_price, control_options, quiet) {
    function(theta) {
      # Model premium for months "nmonth"
      model_price <- solarOption_model(model, theta = theta, nmonths = nmonth, control_options = control_options)
      # Compute the difference from target price
      price_difference <- abs(model_price$payoff_year$premium - target_price)
      # Expected return in absolute value
      expected_return <- price_difference/target_price
      if (!quiet) message("Error: ", price_difference, " Expected return: |", format(expected_return*100, digits = 4), "| Lambda: ", format(theta, digits = 5),
                          "\r", appendLF = FALSE)
      flush.console()
      return(price_difference^2)
    }
  }
  # Calibrator function
  calibrator <- function(loss) {
    # Optimal Esscher parameter
    opt <- optim(par = init_lambda, loss, method = "Brent", lower = lower_lambda, upper = upper_lambda)
    return(opt)
  }
  if (length(nmonths) > 1 & length(target_price) == 1){
    # Specify the loss for a month
    loss <- loss_function(model, nmonths, target_price, control_options, quiet)
    if (!quiet) cat(paste0("\033[1;35m---------------\033[0m", " Calibrating Yearly Esscher parameter ", "\033[1;32m", "\033[1;35m---------------\033[0m \n"))
    # Optimal Esscher parameter
    par <- calibrator(loss)$par
  }  else {
    par <- c()
    for(nmonth in nmonths) {
      flush.console()
      message("", appendLF = TRUE)

      # Test tolerance parameter
      if (!quiet) cat(paste0("\033[1;35m---------------\033[0m", " Calibrating Monthly Esscher parameter (", "\033[1;32m",
                  lubridate::month(nmonth, label = TRUE, abbr = FALSE), "\033[0m", ")",  " \033[1;35m---------------\033[0m \n"))
      # Specify the loss for a month
      loss <- loss_function(model, nmonth, target_price[nmonth], control_options, quiet)
      # Optimal Esscher parameter
      opt_par <- calibrator(loss)$par
      names(opt_par) <- lubridate::month(nmonth, label = TRUE)
      par <- c(par, opt_par)
    }
  }
  return(par)
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
  n_key_points = control_esscher$n_key_points
  quiet = control_esscher$quiet
  # Set target.Yt always on TRUE
  target.Yt <- control_options$target.Yt

  # Compute payoffs
  # Fair Payoff bootstrapped (historical)
  payoff_boot <- solarOption_bootstrap(model, nsim = nsim, ci = ci, seed = seed, control_options = control_options)$payoff_year
  # Historical payoff computed on data (historical)
  payoff_hist <- solarOption_historical(model, control_options = control_options)
  # P-Payoff computed with model
  payoff_model_P <- solarOption_model(model, theta = 0, control_options = control_options)$payoff_year
  # Difference between bootstrapped and model P-price
  error_boot_P <- payoff_boot$premium - payoff_model_P$premium
  # Percentage error between bootstrapped and model P-price
  perc_error_boot_P <- error_boot_P/payoff_boot$premium*100

  # 1) ------- Calibrate the optimal theta for the bounds -------
  # Default benchmark for worste price
  benchmark_Qup <- payoff_boot$premium_up + error_boot_P
  # Optimal upper Esscher parameter (up, worse case)
  opt_theta_up <- solarEsscher_calibrator(model, nmonths = 1:12, target_price = benchmark_Qup,
                                          control_esscher = control_esscher, control_options = control_options)
  # Optimal up and down parameters
  esscher_theta <- solarEsscher_theta_bounds(opt_theta_up)[[1]]
  theta_up <- esscher_theta$up
  theta_dw <- esscher_theta$dw
  theta_bar <- 0.5*(theta_up + theta_dw)
  # Compute payoffs
  # Qdw-payoff: best case scenario (model)
  payoff_model_Qdw <- solarOption_model(model, theta = theta_dw, control_options = control_options)
  # Qup-payoff: worste case scenario (model)
  payoff_model_Qup <- solarOption_model(model, theta = theta_up, control_options = control_options)
  # Q-payoff as average between lambda up and down (model)
  payoff_model_Q <- solarOption_model(model, theta = theta_bar, control_options = control_options)

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
  df_grid$theta[1] <- theta_up
  df_grid$theta[n_key_points] <- theta_dw

  key_points <- c(1:n_key_points)[-c(1, n_key_points)]
  if(!purrr::is_empty(key_points)){
    for(i in key_points){
      if (!quiet) message(rep("-", 20), " ", i, "/", n_key_points, " ", rep("-", 20))
      df_grid$theta[i] <- solarEsscher_calibrator(model, nmonths = 1:12, target_price = df_grid$prices[i],
                                                  control_esscher = control_esscher, control_options = control_options)$par
      if (!quiet) message("r: ", df_grid$rates[i], " h: ", df_grid$theta[i])
    }
  }

  # Add manually r = 0
  idx_pos <- which(grid_rates > 0)
  idx_neg <- which(grid_rates < 0)
  df_grid <- dplyr::bind_rows(df_grid[idx_neg,],
                              dplyr::tibble(prices = benchmark_Q, rates = 0, theta = theta_bar),
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
    params = list(Qdw = theta_dw, Q = theta_bar, Qup = theta_up),
    model = model_param_esscher,
    coefficients = esscher_params,
    theta = solarEsscherTheta(esscher_params$h0, esscher_params$h1, esscher_params$h2),
    control = control_esscher
  )
  # Add Historical payoffs to the model
  #model$payoffs$hist <- payoff_hist
  # Compute implied monthly returns
  #implied_expected_returns <- solarOption_implied_return(model, nmonths = 1:12, control_options = control_options)$implied_r
  # Function for Optimal Esscher parameters for monthly case
  #model$esscher$theta_m <- solarEsscherFunction(par = unlist(esscher_params), imp_r = implied_expected_returns)
  #model$esscher$implied_r <- implied_expected_returns

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


#' Function to establish up and down parameters of the Esscher transform
#' A positive theta identify a `down` parameter, a negative theta identify an `up` parameter.
#'
#' @param theta Esscher parameter.
#' @return A `list` with first element named `up` with the positive parameter and second element named `dw` with the negative one.
#'
#' @keywords internals
#' @noRd
#' @export
solarEsscher_theta_bounds <- function(theta){
  # Identify the up and down parameters of Esscher transform
  # dw parameter (positive): theta > 0
  # up parameter (negative): theta < 0
  i <- 1
  par <- list()
  for(i in 1:length(theta)){
    par[[i]] <- list(up = 0, dw = 0)
    if (theta[i] >= 0) {
      par[[i]][["dw"]] <- theta[i]
      par[[i]][["up"]] <- invertEsscherTheta(theta[i])
    } else {
      par[[i]][["dw"]] <- invertEsscherTheta(theta[i])
      par[[i]][["up"]] <- theta[i]
    }
  }
  return(par)
}

#' invertEsscherTheta
#'
#' @keywords internals
#' @noRd
#' @export
invertEsscherTheta <- function(theta){
  # dw parameter = positive theta
  # up parameter = negative theta
  if (theta >= 0) {
    -theta/(1+theta)
  } else {
    -theta/(1-theta)
  }
}


#' solarEsscherFunction
#'
#' @keywords internals
#' @noRd
#' @export
solarEsscherTheta <- function(h1 = 1, h2 = 1, h3 = 1){
  # Parametric functions for the h1 e h2
  h1_tilde <- function(r_imp){par[1] + par[2]*r_imp + par[3]*r_imp^2}
  h2_tilde <- function(r_imp){2*par[3]*r_imp + par[2]}
  # Parametric function for esscher Thetta
  function(r, r_imp = 0){
    h1_tilde(r_imp) + h2_tilde(r_imp)*r + par[3]*r^2
  }
}



