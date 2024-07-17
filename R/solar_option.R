#' Payoff on Historical Data
#'
#' @param data slot `data` from `solarModel` object.
#' @param nmonth index for the months.
#' @param control control function, see `control_solarOption`.
#'
#' @rdname solar_option_payoff_historical
#' @name solar_option_payoff_historical
#' @export
solar_option_payoff_historical <- function(data, nmonth = 1:12, control = control_solarOption()){

  nyears <- control$nyears
  K <- control$K
  put <- control$put

  # Historical Daily Payoffs
  df_payoff <- dplyr::filter(data, Year >= nyears[1] & Year <= nyears[2])
  df_payoff <- dplyr::filter(df_payoff, Month %in% nmonth)
  # Option strike
  df_payoff$strike <- df_payoff$GHI_bar*exp(K)
  # Payoff
  if (put) {
    df_payoff$exercise <- ifelse(df_payoff$GHI < df_payoff$strike, 1, 0)
    df_payoff$payoff <- (df_payoff$strike - df_payoff$GHI)*df_payoff$exercise # Put
    df_payoff$side <- "put"
  } else {
    df_payoff$exercise <- ifelse(df_payoff$GHI > df_payoff$strike, 1, 0)
    df_payoff$payoff <- (df_payoff$GHI - df_payoff$strike)*df_payoff$exercise # Call
    df_payoff$side <- "call"
  }
  df_payoff <- dplyr::select(df_payoff, side, date, Year, Month, Day, n, GHI, strike, payoff, exercise)

  # Historical Payoffs aggregated by Year and Month
  df_year_month <- df_payoff %>%
    dplyr::group_by(Year, Month, side) %>%
    dplyr::reframe(ndays = dplyr::n(),
                   payoff = sum(payoff),
                   exercise = mean(exercise),
                   GHI = sum(GHI),
                   strike = sum(strike))

  # Fair Premium (Historical)
  # Historical Premium by Month and Day (fair prices)
  df_month_day <- df_payoff %>%
    dplyr::group_by(Month, Day, side) %>%
    dplyr::reframe(premium = mean(payoff),
                   exercise = mean(exercise),
                   GHI = mean(GHI),
                   strike = mean(strike))
  df_month_day$n <- 1:nrow(df_month_day)

  # Historical Premium by Month (fair prices)
  df_month <- df_month_day %>%
    dplyr::group_by(Month, side)  %>%
    dplyr::reframe(ndays = dplyr::n(),
                   premium = sum(premium),
                   daily_premium = premium/ndays,
                   exercise = mean(exercise),
                   GHI = sum(GHI),
                   strike = sum(strike))
  # Historical Premium by Year (fair prices)
  df_year <- df_month %>%
    dplyr::group_by(side)  %>%
    dplyr::reframe(ndays = sum(ndays),
                   premium = sum(premium),
                   exercise = mean(exercise),
                   GHI = sum(GHI),
                   strike = sum(strike))
  structure(
    list(
      aggregate = aggregate_payoff(df_month_day, control),
      payoff = df_payoff,
      payoff_year_month = df_year_month,
      payoff_month_day = df_month_day,
      payoff_month = df_month,
      payoff_year = df_year
    ),
    class = c("solarOptionPayoff", "list")
  )
}

#' Payoff on Simulated Data
#'
#' @param sim slot `sim` from `solarModel` object.
#' @param nmonth index for the months.
#' @param nsim number of simulation to use.
#' @param control control function, see `control_solarOption`.
#'
#' @rdname solar_option_payoff_scenarios
#' @name solar_option_payoff_scenarios
#' @export
solar_option_payoff_scenarios <- function(sim, nmonth = 1:12, nsim = NULL, control = control_solarOption()){

  # Control parameters
  nyears <- control$nyears
  K <- control$K
  put <- control$put

  df_payoff <- dplyr::filter(sim$sim, Year >= nyears[1] & Year <= nyears[2])
  df_payoff <- dplyr::filter(df_payoff, Month %in% nmonth)
  if (!is.null(nsim)) {
    nsim <- min(c(nsim, nrow(df_payoff$data[[1]])))
    df_payoff <- dplyr::mutate(df_payoff, data = purrr::map(data, ~.x[1:nsim,]))
  }

  # Simulated Payoffs
  df_payoff <- df_payoff %>%
    dplyr::mutate(strike = purrr::map_dbl(data, ~.x$GHI_bar[1]*exp(K)),
                  side = ifelse(put, "put", "call")) %>%
    dplyr::mutate(
      data = purrr::map(data, ~dplyr::mutate(.x, strike = GHI_bar*exp(K))),
      data = purrr::map(data, ~dplyr::mutate(.x, side = ifelse(put, "put", "call"))),
      data = purrr::map(data,
                        ~dplyr::mutate(.x,
                                       exercise = dplyr::case_when(
                                         side == "put" & GHI <= strike ~ 1,
                                         side == "put" & GHI > strike ~ 0,
                                         side == "call" & GHI <= strike ~ 0,
                                         side == "call" & GHI > strike ~ 1))),
      data = purrr::map(data,
                        ~dplyr::mutate(.x,
                                       payoff = dplyr::case_when(
                                         side == "put" ~ (strike - GHI)*exercise,
                                         side == "call" ~ (GHI - strike)*exercise)))
      ) %>%
    dplyr::mutate(
      payoff = purrr::map_dbl(data, ~mean(.x$payoff)),
      GHI = purrr::map_dbl(data, ~mean(.x$GHI)),
      exercise = purrr::map_dbl(data, ~mean(.x$exercise))) %>%
    dplyr::select(side, date, Year, Month, Day, n, GHI, strike, payoff, exercise)

  # Simulated Payoffs aggregated by Year and Month
  df_year_month <- df_payoff %>%
    dplyr::group_by(Year, Month, side) %>%
    dplyr::reframe(premium = sum(payoff),
                   exercise = mean(exercise),
                   GHI = sum(GHI),
                   strike = sum(strike),
                   ndays = dplyr::n())

  # Simulated Payoffs aggregated by Month and Day
  df_month_day <- df_payoff %>%
    dplyr::group_by(Month, Day, side) %>%
    dplyr::reframe(premium = mean(payoff),
                   exercise = mean(exercise),
                   GHI = mean(GHI),
                   strike = mean(strike))
  df_month_day$n <- 1:nrow(df_month_day)

  # Simulated Premium aggregated by Month (fair price)
  df_month <- df_month_day %>%
    dplyr::group_by(Month, side) %>%
    dplyr::reframe(ndays = dplyr::n(),
                   premium = sum(premium),
                   daily_premium = premium/ndays,
                   exercise = mean(exercise),
                   GHI = sum(GHI),
                   strike = sum(strike))

  # Simulated Premium aggregated by Year (fair price)
  df_year <- df_month %>%
    dplyr::group_by(side) %>%
    dplyr::reframe(ndays = sum(ndays),
                   premium = sum(premium),
                   exercise = mean(exercise),
                   GHI = sum(GHI),
                   strike = sum(strike))
  structure(
    list(
      aggregate = aggregate_payoff(df_month_day, control),
      payoff = df_payoff,
      payoff_year_month = df_year_month,
      payoff_month_day = df_month_day,
      payoff_month = df_month,
      payoff_year = df_year
    ),
    class = c("solarOptionPayoff", "list")
  )
}

#' Bootstrap a fair price from historical data
#'
#' @param model an object of the class `solarModel`.
#' @param nsim number of simulation to bootstrap.
#' @param ci confidence interval for quantile
#' @param seed random seed.
#' @param control control function, see `control_solarOption`.
#'
#' @rdname solar_option_payoff_bootstrap
#' @name solar_option_payoff_bootstrap
#' @export
solar_option_payoff_bootstrap <- function(model, nsim = 500, ci = 0.05, seed = 1, control = control_solarOption()){

  # Historical Payoff
  payoff_hist <- solar_option_payoff_historical(model$data, nmonth = 1:12, control = control)
  # Dataset with empirical quantiles for each Month and day
  df_quantiles <- payoff_hist$payoff %>%
    dplyr::group_by(Month, Day) %>%
    dplyr::select(Month, Day, payoff) %>%
    tidyr::nest() %>%
    dplyr::mutate(quant = purrr::map(data, ~function(u) quantile(.x$payoff, probs = u))) %>%
    dplyr::select(-data)

  ndays <- nrow(df_quantiles)
  df_sim <- df_quantiles
  df_sim$u <- runif(ndays, 0, 1)
  df_sim <- dplyr::mutate(df_sim, payoff = purrr::map2_dbl(quant, u, ~.x(.y)))

  ndays <- nrow(df_quantiles)
  df_sim <- df_quantiles
  boot <- list()
  set.seed(seed)
  for(i in 1:nsim){
    # Simulate grades
    df_sim$u <- runif(ndays, 0, 1)
    # Bootstrapped historical payoff
    df_sim$i <- i
    df_sim <- dplyr::mutate(df_sim, payoff = purrr::map2_dbl(quant, u, ~.x(.y)))
    df_sim$exercise <- ifelse(df_sim$payoff > 0, 1, 0)
    # Store simulation
    boot[[i]] <- dplyr::select(df_sim, -u, -quant)
  }

  df_month_day <- dplyr::bind_rows(boot) %>%
    dplyr::group_by(Month, Day) %>%
    dplyr::reframe(premium_boot = mean(payoff)) %>%
    dplyr::right_join(payoff_hist$payoff_month_day, by = c("Month", "Day")) %>%
    dplyr::mutate(premium = premium*0.5 + premium_boot*0.5) %>%
    dplyr::select(-premium_boot) %>%
    dplyr::select(Month, Day, side, premium, exercise, GHI, strike, n)

  df_month <- df_month_day %>%
    dplyr::group_by(Month, side) %>%
    dplyr::reframe(ndays = dplyr::n(),
                   premium = sum(premium),
                   daily_premium = premium/ndays,
                   exercise = mean(exercise),
                   GHI = sum(GHI),
                   strike = sum(strike))

  df_year <- df_month %>%
    dplyr::group_by(side) %>%
    dplyr::reframe(ndays = sum(ndays),
                   premium = sum(premium),
                   exercise = mean(exercise),
                   GHI = sum(GHI),
                   strike = sum(strike))

  df_year_boot <- dplyr::bind_rows(boot) %>%
    dplyr::mutate(premium_boot = payoff) %>%
    dplyr::select(-payoff, -exercise)%>%
    dplyr::right_join(payoff_hist$payoff_month_day, by = c("Month", "Day")) %>%
    dplyr::mutate(premium = premium*0.5 + premium_boot*0.5) %>%
    dplyr::group_by(i, side) %>%
    dplyr::reframe(premium = sum(premium),
                   exercise = mean(exercise)) %>%
    dplyr::group_by(side) %>%
    dplyr::reframe(ndays = df_year$ndays,
                   premium_dw = quantile(premium, probs = ci),
                   premium_up = quantile(premium, probs = 1-ci),
                   premium = mean(premium),
                   exercise = mean(exercise),
                   nsim = nsim)
  structure(
    list(
      payoff_month_day = df_month_day,
      payoff_month = df_month,
      payoff_year = df_year,
      payoff_boot = df_year_boot
    ),
    class = c("solarOptionPayoff", "list")
  )
}

#' Pricing function for a solar model (for all the year)
#'
#' @param model an object of the class `solarModel`.
#' @param lambda Esscher parameter
#' @param vol unconditional GARCH variance, when `NA` will be used the fitted one,
#' @param nmonths index for the months.
#' @param control control function, see `control_solarOption`.
#'
#' @rdname solar_option_payoff_model
#' @name solar_option_payoff_model
#' @export
solar_option_payoff_model <- function(model, lambda = 0, vol = NA, nmonths = 1:12, control = control_solarOption()){

  K <- control$K
  put <- control$put

  # Default volatility
  if (is.na(vol)) {
    vol <- model$GARCH$vol
  }

  df <- dplyr::filter(model$seasonal_data, Month %in% nmonths)
  nmonths <- df$Month
  ndays <- df$Day

  #' @examples
  #' lambda = 0
  #' vol = 1
  #' nmonth = 3
  #' nday = 1

  # Pricing function
  pricing_month_day <- function(model, nmonth = 1, nday = 1, lambda = 0, vol = 1){
    #message("Day: ", nday, " Month: ", nmonth)
    # Seasonal Data
    df_n <- dplyr::filter(model$data, Month == nmonth & Day == nday)[1,]
    # Normal Mixture Data
    df_nm <- dplyr::filter(model$NM_model, Month == df_n$Month)
    # Option Strike
    df_n$strike <- df_n$GHI_bar*exp(K)
    # Normal Mixture parameters
    params <- list(
      mu_up = df_nm$mu_up*df_n$sigma_bar*vol + df_n$Yt_bar,
      mu_dw = df_nm$mu_dw*df_n$sigma_bar*vol + df_n$Yt_bar,
      sd_up = df_nm$sd_up*df_n$sigma_bar*vol,
      sd_dw = df_nm$sd_dw*df_n$sigma_bar*vol,
      p_up = df_nm$p_up
    )
    # Esscher Mixture Pdf
    esscher_mixture_pdf <- desscher_mix(unlist(params))
    # Expected value function (GHI)
    e_GHI <- function(x, h = 0){
      df_n$Ct*(1 - model$Xt(x))*esscher_mixture_pdf(x, h = h)
    }
    # Integration Point
    z_tk <- (1/model$params$beta)*(1 - model$params$alpha - df_n$strike/df_n$Ct)
    z <- log(-log(z_tk))
    # Option pricing
    if (put) {
      # Expected value (Put)
      df_n$GHI <- integrate(e_GHI, lower = -Inf, upper = z, h = lambda)$value
      # Probability of exercise
      df_n$exercise <- integrate(esscher_mixture_pdf, lower = -Inf, upper = z, h = lambda)$value
      # Option expected value
      df_n$premium <- df_n$strike*df_n$exercise - df_n$GHI
      # Option type
      df_n$side <- "put"
    } else {
      # Expected value (Call)
      df_n$GHI <- integrate(e_GHI, lower = z, upper = Inf)$value
      # Probability of exercise
      df_n$exercise <- integrate(esscher_mixture_pdf, lower = z, upper = Inf, h = lambda)$value
      # Option expected value
      df_n$premium <- df_n$GHI - df_n$strike*df_n$exercise
      # Option type
      df_n$side <- "call"
    }
    # Expected value
    df_n$GHI <- integrate(e_GHI, lower = -Inf, upper = Inf)$value
    # Select only relevant variables
    df_n <- dplyr::select(df_n, Month, Day, side, premium, exercise, GHI, strike, n)
    return(df_n)
  }

  # Model premium for each day and month
  df_month_day <- purrr::map2_df(nmonths, ndays, ~pricing_month_day(model, nmonth = .x, nday = .y, lambda = lambda, vol = vol))
  df_month_day$n <- 1:nrow(df_month_day)

  # Model premium aggregated for each Month
  df_month <- df_month_day %>%
    dplyr::group_by(Month, side) %>%
    dplyr::reframe(ndays = dplyr::n(),
                   premium = sum(premium),
                   daily_premium = premium/ndays,
                   exercise = mean(exercise),
                   GHI = sum(GHI),
                   strike = sum(strike))

  # Model premium aggregated by Year
  df_year <- df_month %>%
    dplyr::group_by(side) %>%
    dplyr::reframe(ndays = sum(ndays),
                   premium = sum(premium),
                   exercise = mean(exercise),
                   GHI = sum(GHI),
                   strike = sum(strike))

  structure(
    list(
      payoff_month_day = df_month_day,
      payoff_month = df_month,
      payoff_year = df_year
    ),
    class = c("solarOptionPayoff", "list")
  )
}


#' Calibrate Esscher Bounds and parameters
#'
#' @param model an object of the class `solarModel`.
#' @param sim simulations object.
#' @param control_options control function, see `control_solarOption`.
#' @param control control function, see `control_solarEsscher`.
#'
#' @rdname solar_option_esscher_calibrator
#' @name solar_option_esscher_calibrator
#' @export
solar_option_esscher_calibrator <- function(model, sim, control_options = control_solarOption(), control = control_solarEsscher()){

  # Esscher controls
  nsim = control$nsim
  ci = control$ci
  seed = control$seed
  quiet = control$quiet
  n_key_points = control$n_key_points
  # Bounds for parameters
  lower_lambda = control$lower_lambda
  upper_lambda = control$upper_lambda
  # Initial theta parameter
  init_lambda = control$init_lambda

  # Calibrator function for Esscher parameter
  loss_esscher_lambda <- function(lambda, benchmark_price){
    # Price under measure Q
    Q_price <- solar_option_payoff_model(model, lambda = lambda, vol = NA, control = control_options)$payoff_year$premium
    loss <- (Q_price - benchmark_price)^2
    if (!quiet) message("Loss: ", loss, " Lambda: ", lambda)
    return(loss)
  }

  # 1) ------- Calibrate the optimal theta bounds -------
  # Fair Payoff bootstrapped (historical)
  payoff_boot <- solar_option_payoff_bootstrap(model = model, nsim = nsim, ci = ci, seed = seed, control = control_options)
  # Benchmark for worste case price
  benchmark_price <- payoff_boot$payoff_boot$premium_up
  # Optimal upper Esscher parameter (up, worse case)
  opt_lambda_up <- optim(par = init_lambda, loss_esscher_lambda, benchmark_price = benchmark_price,
                         method = "Brent", lower = lower_lambda, upper = upper_lambda)
  # Function to establish up and down parameters
  esscherLambda <- function(lambda){
    # dw parameter = positive lambda
    # up parameter = negative lambda
    if (lambda >= 0) {
      param_dw <- lambda
      param_up <- -lambda/(1+lambda)
    } else {
      param_dw <- -lambda/(1-lambda)
      param_up <- lambda
    }
    return(list(up = param_up, dw = param_dw))
  }

  # Optimal up and down parameters
  lambda_up <- esscherLambda(opt_lambda_up$par)$up
  lambda_dw <- esscherLambda(opt_lambda_up$par)$dw
  # Qdw-payoff: best case scenario (model)
  payoff_model_Qdw <- solar_option_payoff_model(model, vol = NA, lambda = lambda_dw, control = control_options)
  # Qup-payoff: worste case scenario (model)
  payoff_model_Qup <- solar_option_payoff_model(model, vol = NA, lambda = lambda_up, control = control_options)

  # 2) ------- Calibrator function for optimal theta parameters -------
  # Historical payoff computed on data (historical)
  payoff_hist <- solar_option_payoff_historical(model$data, nmonth = 1:12, control = control_options)
  # P-Payoff computed by means of scenarios (simulated)
  payoff_sim <- solar_option_payoff_scenarios(sim, nmonth = 1:12, nsim = NULL, control = control_options)
  # Q-Payoff computed with model
  payoff_model_Q <- solar_option_payoff_model(model, vol = NA, lambda = 0, control = control_options)
  # Benchmark for price under P-measure
  benchmark_price <- payoff_sim$payoff_year$premium
  # Optimal Esscher parameter (from Q to P)
  opt_lambda_Q <- optim(par = init_lambda, fn = loss_esscher_lambda, benchmark_price = benchmark_price,
                        method = "Brent", lower = lower_lambda, upper = upper_lambda)
  # Benchmark for market price (depends on r)
  benchmark_price <- payoff_hist$payoff_year$premium*(control_options$B(1)^(1))
  # Optimal Esscher parameter  (from Q to Qr)
  opt_lambda_Qr <- optim(par = init_lambda, fn = loss_esscher_lambda, benchmark_price = benchmark_price,
                         method = "Brent", lower = lower_lambda, upper = upper_lambda)
  # Qr payoff calibrated on market price (model)
  payoff_model_Qr <- solar_option_payoff_model(model, vol = NA, lambda = opt_lambda_Qr$par, control = control_options)
  # P-payoff calibrated on P-price (model)
  payoff_model_P <- solar_option_payoff_model(model, vol = NA, lambda = opt_lambda_Q$par, control = control_options)

  # 3) Calibrate the Esscher parameter for different levels of r
  # Benchmark for market prices with r = 0
  fair <- payoff_hist$payoff_year$premium
  # Benchmark for best case prices with r > 0
  lower <- payoff_model_Qdw$payoff_year$premium
  # Benchmark for worste case prices with r < 0
  upper <- payoff_model_Qup$payoff_year$premium

  # Grid of prices in the Esscher corridor (decreasing)
  grid_prices <- seq(upper, lower, length.out = n_key_points)
  # Implied grid of risk-free rates (decreasing)
  grid_rates <- fair/grid_prices - 1
  # Implied grid of Esscher lambda
  df_grid <- dplyr::tibble(prices = grid_prices, rates = grid_rates, theta = NA)
  for(i in 1:n_key_points){
    if (!quiet) message(rep("-", 20), " ", i, "/", n_key_points, " ", rep("-", 20))
    df_grid$theta[i] <- optim(init_lambda, loss_esscher_lambda,  benchmark_price = df_grid$prices[i],
                               method = "Brent", lower = lower_lambda, upper = upper_lambda)$par
    if (!quiet) message("r: ", df_grid$rates[i], " h: ", df_grid$theta[i])
  }

  # 4) Fit a parametric model for the Esscher parameter:
  #   theta = h0 + h1 r model in terms of risk-free rates
  model_param_esscher <- lm(theta ~ rates, data = df_grid)
  # Fitted values
  df_grid$lambda_fit <- predict(model_param_esscher, newdata = df_grid)
  # Fitted parameters
  esscher_params <- list(h0 = model_param_esscher$coefficients[1], h1 = model_param_esscher$coefficients[2])

  # - Add Esscher parametric model
  model$esscher$model <- list(
    grid = df_grid,
    model = model_param_esscher,
    params = esscher_params,
    theta = function(x) esscher_params$h0 + x*esscher_params$h1,
    r = function(x) (x - esscher_params$h0)/esscher_params$h1)
  # - Add Esscher parameter to the model
  model$esscher$params <- list(Qdw = lambda_dw,
                               P_to_Q  = opt_lambda_Q$par,
                               Q_to_Qr = opt_lambda_Qr$par,
                               P_to_Qr = -(-opt_lambda_Q$par + opt_lambda_Qr$par),
                               P_to_Qdw = opt_lambda_Q$par - lambda_dw,
                               P_to_Qup = opt_lambda_Q$par - lambda_up,
                               Qup = lambda_up)
  # - Add Esscher control setups to the model
  model$esscher$control <- control
  # - Add payoffs to the model
  # Historical
  model$payoffs$hist <- payoff_hist
  # Bootstrapped
  model$payoffs$boot <- payoff_boot
  # Simulated
  model$payoffs$sim$P <- payoff_sim
  # Model
  model$payoffs$model$P <- payoff_model_P
  model$payoffs$model$Q <- payoff_model_Q
  model$payoffs$model$Qdw <- payoff_model_Qdw
  model$payoffs$model$Qup <- payoff_model_Qup
  model$payoffs$model$Qr <- payoff_model_Qr

  return(model)
}

#' Structure payoffs
#'
#' @param model an object of the class `solarModel`.
#' @param type can be `sim` or `model`.
#' @param exact_daily_premium when `TRUE` the historical premium is computed as daily average among all the years.
#' Otherwise the monthly premium is computed and then divided by the number of days of the month.
#'
#' @rdname solar_option_payoff_structure
#' @name solar_option_payoff_structure
#' @export
solar_option_payoff_structure <- function(model, type = "sim", exact_daily_premium = TRUE){

  type <- match.arg(type, choices = c("sim", "model"))
  payoff <- model$payoffs[c("hist", type)]
  # Yearly Premium
  df_year <- dplyr::tibble(
    side =  payoff$hist$payoff_year$side,
    ndays =  payoff$hist$payoff_year$ndays,
    # Premiums for the option
    premium = payoff$hist$payoff_year$premium,
    premium_P = payoff[[type]]$P$payoff_year$premium,
    premium_Q = payoff[[type]]$Q$payoff_year$premium,
    premium_Qdw = payoff[[type]]$Qdw$payoff_year$premium,
    premium_Qup = payoff[[type]]$Qup$payoff_year$premium,
    premium_Qr = payoff[[type]]$Qr$payoff_year$premium,
    # Probabilities of exercise the option
    exercise = payoff$hist$payoff_year$exercise,
    exercise_P = payoff[[type]]$P$payoff_year$exercise,
    exercise_Q = payoff[[type]]$Q$payoff_year$exercise,
    exercise_Qup = payoff[[type]]$Qup$payoff_year$exercise,
    exercise_Qdw = payoff[[type]]$Qdw$payoff_year$exercise,
    exercise_Qr = payoff[[type]]$Qr$payoff_year$exercise
  )

  # Monthly Premium
  df_month <- dplyr::tibble(
    Month = payoff$hist$payoff_month$Month,
    side = payoff$hist$payoff_month$side,
    n = payoff$hist$payoff_month$ndays,
    # Premiums for the option
    premium = payoff$hist$payoff_month$premium,
    premium_P = payoff[[type]]$P$payoff_month$premium,
    premium_Q = payoff[[type]]$Q$payoff_month$premium,
    premium_Qup = payoff[[type]]$Qup$payoff_month$premium,
    premium_Qdw = payoff[[type]]$Qdw$payoff_month$premium,
    premium_Qr = payoff[[type]]$Qr$payoff_month$premium,
    # Probabilities of exercise the option
    exercise = payoff$hist$payoff_month$exercise,
    exercise_P = payoff[[type]]$P$payoff_month$exercise,
    exercise_Q = payoff[[type]]$Q$payoff_month$exercise,
    exercise_Qup = payoff[[type]]$Qup$payoff_month$exercise,
    exercise_Qdw = payoff[[type]]$Qdw$payoff_month$exercise,
    exercise_Qr = payoff[[type]]$Qr$payoff_month$exercise
  )

  # Monthly Daily Premium
  df_month_day_mean <- dplyr::tibble(
    Month = payoff$hist$payoff_month$Month,
    side = payoff$hist$payoff_month$side,
    n = payoff$hist$payoff_month$ndays,
    # Premiums for the option
    premium = payoff$hist$payoff_month$daily_premium,
    premium_P = payoff[[type]]$P$payoff_month$daily_premium,
    premium_Q = payoff[[type]]$Q$payoff_month$daily_premium,
    premium_Qup = payoff[[type]]$Qup$payoff_month$daily_premium,
    premium_Qdw = payoff[[type]]$Qdw$payoff_month$daily_premium,
    premium_Qr = payoff[[type]]$Qr$payoff_month$daily_premium,
  )

  # Monthly-Day Premium
  df_month_day <- dplyr::tibble(
    Month = payoff$hist$payoff_month_day$Month,
    Day = payoff$hist$payoff_month_day$Day,
    side = payoff$hist$payoff_month_day$side,
    # Premiums for the option
    premium = payoff$hist$payoff_month_day$premium,
    premium_P = payoff[[type]]$P$payoff_month_day$premium,
    premium_Q = payoff[[type]]$Q$payoff_month_day$premium,
    premium_Qup = payoff[[type]]$Qup$payoff_month_day$premium,
    premium_Qdw = payoff[[type]]$Qdw$payoff_month_day$premium,
    premium_Qr = payoff[[type]]$Qr$payoff_month_day$premium,
    # Probabilities of exercise the option
    exercise = payoff$hist$payoff_month_day$exercise,
    exercise_P = payoff[[type]]$P$payoff_month_day$exercise,
    exercise_Q = payoff[[type]]$Q$payoff_month_day$exercise,
    exercise_Qup = payoff[[type]]$Qup$payoff_month_day$exercise,
    exercise_Qdw = payoff[[type]]$Qdw$payoff_month_day$exercise,
    exercise_Qr = payoff[[type]]$Qr$payoff_month_day$exercise,
    n = payoff$hist$payoff_month_day$n
  )

  # Historical Daily Premiums
  df_payoff <- dplyr::select(payoff$hist$payoff, -exercise, -GHI, -strike)
  if (!exact_daily_premium) {
    df_payoff <- dplyr::left_join(df_payoff, dplyr::select(df_month_day_mean, Month, premium:premium_Qr), by = c("Month"))
  } else {
    df_payoff <- dplyr::left_join(df_payoff, dplyr::select(df_month_day, Month, Day, premium:premium_Qr), by = c("Month", "Day"))
  }
  df_payoff <- dplyr::mutate(df_payoff, net_payoff = payoff - premium)

  # Cumulated Historical Payoff minus (daily) Premiums
  j <- 1
  cumulated_payoff <- list()
  seq_years <- seq(min(df_payoff$Year), max(df_payoff$Year), 1)
  for (j in 1:length(seq_years)){
    df_cum <- dplyr::filter(df_payoff, Year == seq_years[j])
    # Remove the 29-02 for graphic purposes
    df_cum <- df_cum[paste0(df_cum$Month, "-", df_cum$Day) != "2-29",]
    df_cum <- dplyr::mutate(df_cum,
                            cum_net_payoff = NA,
                            cum_net_payoff_Qdw = NA,
                            cum_net_payoff_P = NA,
                            cum_net_payoff_Q = NA,
                            cum_net_payoff_Qr = NA,
                            cum_net_payoff_Qup = NA)
    ndays <- nrow(df_cum)
    for(i in 1:nrow(df_cum)){
      df_cum$cum_net_payoff[i] <- (df_year$premium + sum(df_cum$payoff[1:i]) - sum(df_cum$premium[1:i]))
      df_cum$cum_net_payoff_disc[i] <- (df_year$premium*model$payoffs$control$B(1) + sum(df_cum$payoff[1:i]) - sum(df_cum$premium[1:i])*model$payoffs$control$B(1))
      df_cum$cum_net_payoff_P[i] <- df_year$premium_P + sum(df_cum$payoff[1:i]) - sum(df_cum$premium_P[1:i])
      df_cum$cum_net_payoff_Q[i] <- df_year$premium_Q + sum(df_cum$payoff[1:i]) - sum(df_cum$premium_Q[1:i])
      df_cum$cum_net_payoff_Qdw[i] <- df_year$premium_Qdw + sum(df_cum$payoff[1:i]) - sum(df_cum$premium_Qdw[1:i])
      df_cum$cum_net_payoff_Qup[i] <- df_year$premium_Qup + sum(df_cum$payoff[1:i]) - sum(df_cum$premium_Qup[1:i])
      df_cum$cum_net_payoff_Qr[i] <- df_year$premium_Qr + sum(df_cum$payoff[1:i]) - sum(df_cum$premium_Qr[1:i])
    }
    cumulated_payoff[[j]] <- df_cum
  }

  # Compute the expected trajectory for each day
  df_cum <- dplyr::bind_rows(cumulated_payoff) %>%
    dplyr::group_by(Month, Day) %>%
    dplyr::mutate(
      e_cum_net_payoff = mean(cum_net_payoff),
      e_cum_net_payoff_disc = mean(cum_net_payoff_disc),
      e_cum_net_payoff_P = mean(cum_net_payoff_P),
      e_cum_net_payoff_Q = mean(cum_net_payoff_Q),
      e_cum_net_payoff_Qdw = mean(cum_net_payoff_Qdw),
      e_cum_net_payoff_Qup = mean(cum_net_payoff_Qup),
      e_cum_net_payoff_Qr = mean(cum_net_payoff_Qr)
    ) %>%
    dplyr::ungroup()
  model$payoffs[[type]]$structured <- list(payoff = df_payoff,
                                           payoff_year = df_year,
                                           payoff_month = df_month,
                                           payoff_month_day = df_month_day,
                                           payoff_cum = df_cum)
  return(model)
}

