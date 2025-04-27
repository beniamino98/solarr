#' Control parameters for a solar option
#'
#' @param nyears numeric vector. Interval of years considered. The first element will be the minimum and the second the maximum years used in
#' the computation of the fair payoff.
#' @param K numeric, level for the strike with respect to the seasonal mean. The seasonal mean is multiplied by `exp(K)`.
#' @param leap_year logical, when `FALSE`, the default, the year will be considered of 365 days, otherwise 366.
#' @param nsim integer, number of simulations used to bootstrap the premium's bounds. See \code{\link{solarOption_historical_bootstrap}}.
#' @param ci numeric, confidence interval for bootstrapping. See \code{\link{solarOption_historical_bootstrap}}.
#' @param seed integer, random seed for reproducibility. See \code{\link{solarOption_historical_bootstrap}}.
#' @param B function. Discount factor function. Should take as input a number (in years) and return a discount factor.
#'
#' @examples
#' control_options <- control_solarOption()
#'
#' @rdname control_solarOption
#' @name control_solarOption
#' @export
control_solarOption <- function(nyears = c(2005, 2024), K = 0, leap_year = FALSE, nsim = 200, ci = 0.05, seed = 1, B = discountFactor()){

  structure(
    list(
      nyears = nyears,
      from = as.Date(paste0(nyears[1], "-01-01")),
      to = as.Date(paste0(nyears[2], "-01-01")),
      K = K,
      leap_year = leap_year,
      nsim = nsim,
      ci = ci,
      seed = seed,
      B = B
    ),
    class = c("controlSolarOption", "list")
  )
}

#' Payoff on Historical Data
#'
#' @param model An object with the class `solarModel`. See the function \code{\link{solarModel}} for details.
#' @param nmonths Numeric vector. Months in which the payoff should be computed. Can vary from 1 (January) to 12 (December).
#' @param put Logical. When `TRUE`, the default, will be computed the price for a `put` contract, otherwise for a `call` contract.
#' @param control_options Named list. Control parameters, see \code{\link{control_solarOption}} for more details.
#'
#' @return An object of the class `solarOptionPayoff`.
#' @examples
#' model <- Bologna
#' solarOption_historical(model, put=TRUE)
#' solarOption_historical(model, put=FALSE)
#'
#' @rdname solarOption_historical
#' @name solarOption_historical
#' @export
solarOption_historical <- function(model, nmonths = 1:12, put = TRUE, control_options = control_solarOption()){

  # Options control parameters
  K <- control_options$K
  # Target and seasonal mean
  target <- model$target
  target_bar <- paste0(model$target, "_bar")
  # Complete data
  data <- dplyr::select(model$data, date, n, Year, Month, Day, tidyr::any_of(c(target, target_bar)))

  # Filter for control years
  data <- dplyr::filter(data, date >= control_options$from & date <= control_options$to)
  # Filter for selected months
  data <- dplyr::filter(data, Month %in% nmonths)
  # Rename target variable
  data$Rt <- data[[target]]
  # Compute strike price
  data$strike <- data[[target_bar]]*exp(K)
  # Select only relevant variables
  df_payoff <- dplyr::select(data, date, n, Year, Month, Day, Rt, strike)

  # Daily option payoffs
  if (put) {
    # Payoff for a put
    df_payoff$exercise <- ifelse(df_payoff$Rt < df_payoff$strike, 1, 0)
    df_payoff$payoff <- (df_payoff$strike - df_payoff$Rt) * df_payoff$exercise
    df_payoff$side <- "put"
  } else {
    # Payoff for a call
    df_payoff$exercise <- ifelse(df_payoff$Rt > df_payoff$strike, 1, 0)
    df_payoff$payoff <- (df_payoff$Rt - df_payoff$strike) * df_payoff$exercise
    df_payoff$side <- "call"
  }

  # Aggregation of Historical Payoffs
  df_payoff$premium <- df_payoff$payoff
  # Reorder variables
  df_payoff <- dplyr::select(df_payoff, side, date, Year, Month, Day, Rt, strike, premium, payoff, exercise)
  # Output structure
  solarOptionPayoff(df_payoff, control_options$leap_year)
}

#' Bootstrap a fair premium from historical data
#'
#' @inheritParams solarOption_historical
#'
#' @return An object of the class `solarOptionBoot`.
#' @examples
#' model <- Bologna
#' solarOption_historical_bootstrap(model, control_options = control_solarOption(ci = 0.4, nsim = 1000))
#'
#' @rdname solarOption_historical_bootstrap
#' @name solarOption_historical_bootstrap
#' @export
solarOption_historical_bootstrap <- function(model, put = TRUE, control_options = control_solarOption()){

  # Control parameters
  set.seed(control_options$seed)
  leap_year <- control_options$leap_year
  ci = control_options$ci
  nsim = control_options$nsim
  option_type = ifelse(put, "put", "call")
  # Historical monthly payoffs
  payoff_hist <- solarOption_historical(model, nmonths = 1:12, put = put, control_options = control_options)

  # Group by month and day and create an empirical quantile function for each group
  df_quantiles <- payoff_hist$payoff %>%
    dplyr::group_by(Month) %>%
    dplyr::select(Month, Day, payoff) %>%
    tidyr::nest() %>%
    dplyr::mutate(quant = purrr::map(data, ~function(u) quantile(.x$payoff, probs = u)),
                  ndays = lubridate::days_in_month(Month)) %>%
    dplyr::select(-data)

  if (leap_year) {
    df_quantiles$ndays[2] <- 29
  }

  # 1) Initialize the dataset
  df_sim <- df_quantiles
  # Compute the simulated daily payoffs by applying the quantile function to the simulated grades
  df_sim <- dplyr::mutate(df_sim, u = purrr::map2(quant, ndays, ~.x(runif(.y))))
  # Evaluate the exercise set based on the simulated payoff
  df_sim <- dplyr::mutate(df_sim, exercise = purrr::map_dbl(u, ~mean(ifelse(.x > 0, 1, 0))))
  # Compute the monthly payoffs
  df_sim <- dplyr::mutate(df_sim, payoff = purrr::map_dbl(u, ~sum(.x)))
  # 2) Bootstrap the payoff
  boot <- dplyr::tibble()
  for(j in 1:nsim){
    df_sim$j <- j
    # Compute the simulated daily payoffs by applying the quantile function to the simulated grades
    df_sim <- dplyr::mutate(df_sim, u = purrr::map2(quant, ndays, ~.x(runif(.y))))
    # Evaluate the exercise set based on the simulated payoff
    df_sim <- dplyr::mutate(df_sim, exercise = purrr::map_dbl(u, ~mean(ifelse(.x > 0, 1, 0))))
    # Compute the monthly payoffs
    df_sim <- dplyr::mutate(df_sim, payoff = purrr::map_dbl(u, ~sum(.x)))
    # Store the simulation
    boot <- dplyr::bind_rows(boot, dplyr::select(dplyr::select(df_sim, -u, -quant), j, dplyr::everything()))
  }
  # 3) Aggregate bootstraps
  # Aggregate by Month
  df_month <- boot %>%
    dplyr::group_by(Month) %>%
    dplyr::summarise(
      side = option_type,
      premium_dw = quantile(payoff, probs = ci),
      premium = mean(payoff),
      premium_up = quantile(payoff, probs = 1-ci),
      exercise_dw = quantile(exercise, probs = ci),
      exercise = mean(exercise),
      exercise_up = quantile(exercise, probs = 1-ci),
      ndays = mean(ndays)) %>%
    dplyr::mutate(daily_premium = premium/ndays)

  # Aggregate by a simulation of 1 Year
  df_year <- boot %>%
    dplyr::group_by(j) %>%
    dplyr::summarise(ndays = sum(ndays),
                     exercise = mean(exercise),
                     payoff = sum(payoff)) %>%
    dplyr::summarise(
      side = option_type,
      premium_dw = quantile(payoff, probs = ci),
      premium = mean(payoff),
      premium_up = quantile(payoff, probs = 1-ci),
      exercise_dw = quantile(exercise, probs = ci),
      exercise = mean(exercise),
      exercise_up = quantile(exercise, probs = 1-ci),
      ndays = mean(ndays),
      nsim = nsim)

  structure(
    list(
      payoff = boot,
      payoff_month = df_month,
      payoff_year = df_year
    ),
    class = c("solarOptionBoot", "list")
  )
}

#' Payoff on simulated Data
#'
#' @param scenario object with the class `solarModelScenario`. See the function \code{\link{solarModel_scenarios}} for details.
#' @param nsim number of simulation to use for computation.
#' @inheritParams solarOption_historical
#'
#' @param control_options control function, see \code{\link{control_solarOption}} for details.
#'
#' @return An object of the class `solarOptionPayoff`.
#' @examples
#' model <- Bologna
#' scenario <- solarScenario(model, from = "2011-01-01", to = "2012-01-01", by = "1 month", nsim = 10, seed = 3)
#' solarOption_scenario(model, scenario)
#' solarOption_historical(model)
#'
#' @rdname solarOption_scenario
#' @name solarOption_scenario
#' @export
solarOption_scenario <- function(model, scenario, nmonths = 1:12, put = TRUE, nsim, control_options = control_solarOption()){

  # Control parameters
  nyears <- control_options$nyears
  K <- control_options$K
  leap_year = control_options$leap_year
  # Target and seasonal mean
  target <- scenario$target
  target_bar <- paste0(target, "_bar")

  # Simulated daily Payoffs
  sim <- scenario$sim
  # Filter for control years
  df_payoff <- dplyr::filter(sim, date >= control_options$from & date <= control_options$to)
  # Filter for selected months
  df_payoff <- dplyr::filter(df_payoff, Month %in% nmonths)
  # Eventually reduce the number of simulations used
  if (!missing(nsim)) {
    nsim <- min(c(nsim, nrow(df_payoff$data[[1]])))
    df_payoff <- dplyr::mutate(df_payoff, data = purrr::map(data, ~.x[1:nsim,]))
  }
  # Add seasonal mean (strike)
  data <- dplyr::left_join(df_payoff, model$data[,c("date", target_bar)], by = "date")
  df_payoff$strike <- data[[target_bar]]*exp(K)
  df_payoff <- tidyr::unnest(df_payoff, cols = "data")
  df_payoff$Rt <- df_payoff[[target]]
  # Compute simulated payoffs
  df_payoff <- df_payoff %>%
    dplyr::mutate(side = ifelse(put, "put", "call"),
                  exercise = dplyr::case_when(side == "put" & Rt <= strike ~ 1,
                                              side == "put" & Rt > strike ~ 0,
                                              side == "call" & Rt <= strike ~ 0,
                                              side == "call" & Rt > strike ~ 1),
                  payoff_sim = dplyr::case_when(side == "put" ~ (strike - Rt)*exercise,
                                                side == "call" ~ (Rt - strike)*exercise)) %>%
    dplyr::group_by(date) %>%
    dplyr::group_by(side, date, Year, Month, Day) %>%
    dplyr::reframe(
      Rt = mean(Rt),
      strike = mean(strike),
      premium = mean(payoff_sim),
      exercise = mean(exercise)
    ) %>%
    dplyr::ungroup()
  # Add and compute realized payoff
  data <- dplyr::left_join(scenario$emp, model$data[,c("date", target_bar)], by = "date")
  scenario$emp$Rt <- scenario$emp[[target]]
  scenario$emp$strike <- data[[target_bar]]*exp(K)
  scenario$emp$side <- ifelse(put, "put", "call")
  scenario$emp <- dplyr::mutate(scenario$emp,
                                payoff = dplyr::case_when(
                                  side == "put" ~ (strike - Rt)*ifelse(Rt < strike, 1, 0),
                                  TRUE ~ (Rt - strike)*ifelse(Rt > strike, 1, 0)))

  df_payoff <- dplyr::left_join(df_payoff, scenario$emp[,c("date", "payoff")], by = "date")
  # Reorder variables
  df_payoff <- dplyr::select(df_payoff, side, date, Year, Month, Day, Rt, strike, premium, payoff, exercise)
  # Output structure
  solarOptionPayoff(df_payoff, control_options$leap_year)
}

#' Compute the price of a `solarOption`
#'
#' @param moments description
#' @param sorad An object of the class `solarOption`.
#' @inheritParams solarOption_historical
#'
#' @examples
#' model <- Bologna$clone(TRUE)
#' moments <- filter(model$moments$conditional, Year == 2022)
#' # Pricing without contracts
#' solarOption_pricing(moments[1,])
#' # Pricing with contracts specification
#' sorad <- solarOption$new()
#' sorad$set_contract("2021-12-31", "2022-01-01", "2022-04-20", moments$GHI_bar[1])
#' solarOption_pricing(moments[1,], sorad)
#' solarOption_pricing(moments[1,], sorad, theta = 0.02)
#' solarOption_pricing(moments[1,], sorad, theta = -0.02)
#'
#' @rdname solarOption_pricing
#' @name solarOption_pricing
#'
#' @export
solarOption_pricing <- function(moments, sorad, theta = 0, put = TRUE, control_options = control_solarOption()){
  # Control parameters
  K <- control_options$K

  df_n <- moments
  # Create datasets for the moments
  comb <- dplyr::tibble(mean = c(df_n$M_Y1, df_n$M_Y0), sd = c(df_n$S_Y1, df_n$S_Y0), probs = c(df_n$p1, 1-df_n$p1))

  # Esscher parameter
  if (theta == 0) {
    # Mixture Pdf
    pdf_Yt <- function(x) dmixnorm(x, comb$mean, comb$sd, comb$probs)
    # Distribution Yt
    cdf_Yt <- function(x) pmixnorm(x, comb$mean, comb$sd, comb$probs)
  } else {
    # Esscher Mixture Pdf
    pdf_Yt <- desscherMixture(comb$mean, comb$sd, comb$probs, theta)
    # Esscher Distribution Yt
    cdf_Yt <- pesscherMixture(comb$mean, comb$sd, comb$probs, theta)
  }

  # G-function
  G <- function(lower, upper) integrate(function(y) exp(-exp(y)) * pdf_Yt(y), lower = lower, upper = upper)$value
  # Strike
  if (missing(sorad)){
    strike <- df_n$GHI_bar*exp(K)
    t_hor <- df_n$date
  } else {
    strike <- sorad$strike*exp(K)
    t_hor <- sorad$t_hor
  }
  # Strike price in terms of Y
  K_Y <- log(log(df_n$beta) -log(1 - strike / df_n$Ct - df_n$alpha))

  # Option pricing
  if (put) {
    # Probability of exercise
    exercise <- cdf_Yt(K_Y)
    # Value first component
    V1 <- (strike - df_n$Ct * (1 - df_n$alpha)) * exercise
    # Value second component
    V2 <- df_n$Ct * df_n$beta * G(-Inf, K_Y)
    # Option expected value
    premium <- V1 + V2
  } else {
    # Probability of exercise
    exercise <- 1 - cdf_Yt(K_Y)
    # Value first component
    V1 <- (df_n$Ct * (1 - df_n$alpha) - strike) * exercise
    # Value second component
    V2 <- - df_n$Ct * df_n$beta * G(K_Y, Inf)
    # Option expected value
    premium <- V1 + V2
  }

  # Select only relevant variables
  df_n <- dplyr::tibble(
    date = t_hor,
    Year = lubridate::year(date),
    Month = lubridate::month(date),
    Day = lubridate::day(date),
    side = ifelse(put, "put", "call"),
    premium = V1 + V2,
    exercise = exercise,
    strike = strike)
  return(df_n)
}

#' Compute the price of a `solarOptionPortfolio`
#'
#' @param moments description
#' @param portfolio A list of objects of the class `solarOptionPortfolio`.
#' @inheritParams solarOption_historical
#'
#' @examples
#' # Model
#' model <- Bologna$clone(TRUE)
#' # Pricing without portfolio
#' moments <- model$moments$unconditional
#' # Premium
#' premium_Vt <- solarOption_model(model, moments, theta = 0.0, put = TRUE)
#' # Pricing date
#' t_now <- as.Date("2021-12-31")
#' # Inception date
#' t_init <- as.Date("2022-01-01")
#' # Maturity date
#' t_hor <- as.Date("2022-12-31")
#' # SoRad portfolio
#' portfolio <- SoRadPorfolio(model, t_now, t_init, t_hor)
#' # Moments
#' moments <- purrr::map_df(portfolio, ~model$Moments(t_now, .x$t_hor))
#' # Premium
#' premium_Vt <- solarOption_model(model, moments, portfolio, theta = 0.0, put = TRUE)
#' premium_Vt$payoff_year$premium
#' @rdname solarOption_model
#' @name solarOption_model
#'
#' @export
solarOption_model <- function(model, moments, portfolio, nmonths = 1:12, theta = 0, implvol = 1, put = TRUE, control_options = control_solarOption()){

  # Options control parameters
  K <- control_options$K
  leap_year <- control_options$leap_year

  if (missing(portfolio)) {
    moments <- dplyr::filter(moments, date >= control_options$from & date <= control_options$to)
    moments <- dplyr::filter(moments, Month %in% nmonths)
    moments$S_Y0 <- moments$S_Y0 * implvol
    moments$S_Y1 <- moments$S_Y1 * implvol
    data <- purrr::map_df(1:nrow(moments), ~solarOption_pricing(moments[.x,], theta = theta, put = put, control_options = control_options))
  } else {
    moments$S_Y0 <- moments$S_Y0 * implvol
    moments$S_Y1 <- moments$S_Y1 * implvol
    data <- purrr::map_df(portfolio, ~solarOption_pricing(dplyr::filter(moments, date == .x$t_hor), .x,
                                                          theta = theta, put = put, control_options = control_options))
  }
  # Add realized GHI
  data <- dplyr::left_join(data, dplyr::select(model$data, date, Rt = "GHI"), by = "date")
  # Compute realized payoff
  data <- dplyr::mutate(data,
                        payoff = dplyr::case_when(
                          side == "put" ~ (strike - Rt)*ifelse(Rt < strike, 1, 0),
                          TRUE ~ (Rt - strike)*ifelse(Rt > strike, 1, 0)))
  # Reorder variables
  df_payoff <- dplyr::select(data, side, date, Year, Month, Day, Rt, strike, premium, payoff, exercise)
  solarOptionPayoff(df_payoff, leap_year)
}

#' Calibrator function for solarOptions
#' Recalibrate and adjust the mixture parameters such that the model premium
#' matches exactly the historical premium for all the months.
#'
#' @inheritParams solarOption_historical
#' @param abstol The absolute convergence tolerance. Only useful for non-negative functions, as a tolerance for reaching zero.
#' @param reltol Relative convergence tolerance. The algorithm stops if it is unable to reduce the value by a factor of reltol * (abs(val) + reltol) at a step.
#' Defaults to `sqrt(.Machine$double.eps)`, typically about 1e-8.
#' @param conditional Logical. When `TRUE` the target will be the option price computed with conditional moments.
#' @examples
#' model <- Bologna
#' nmonth <- 5
#' model_cal <- solarOption_calibrator(model, nmonths = nmonth, reltol=1e-3)
#' # Compare log-likelihoods
#' model$loglik
#' model_cal$loglik
#' # Compare parameters
#' model$NM_model$coefficients[nmonth,]
#' model_cal$NM_model$coefficients[nmonth,]
#' # Compare moments
#' model$NM_model$moments[nmonth,]
#' model_cal$NM_model$moments[nmonth,]
#' @export
solarOption_calibrator <- function(model, nmonths = 1:12, conditional = TRUE, abstol = 1e-4, reltol = 1e-4, control_options = control_solarOption()){
  # Clone the model
  model_cal <- model$clone(deep = TRUE)
  # Historical monthly payoff
  payoff_call <- solarOption_historical(model_cal, put = FALSE, control_options = control_options)$payoff_month$premium
  payoff_put <- solarOption_historical(model_cal, put = TRUE, control_options = control_options)$payoff_month$premium
  # Loss function
  loss_function <- function(params, model_cal, nmonth, init_params, conditional){
    # Update the parameters
    model_cal$NM_model$model[[nmonth]]$update(means = params[1:2], sd = params[3:4])
    model_cal$update_moments()
    moments <- model_cal$moments$conditional
    if (!conditional) {
      nyear <- lubridate::year(model_cal$dates$train$to)
      moments <- dplyr::filter(model_cal$moments$unconditional, Year == nyear)
    }
    # Put pricing function
    df_put <- solarOption_model(model_cal, moments, put = TRUE, nmonth = nmonth, control_options = control_options)
    # Put price
    price_put <- df_put$payoff_year$premium
    # Call pricing function
    df_call <- solarOption_model(model_cal, moments, put = FALSE, nmonth = nmonth, control_options = control_options)
    # Call price
    price_call <- df_call$payoff_year$premium
    # Loss function
    loss <- abs(price_call - payoff_call[nmonth]) + abs(price_put - payoff_put[nmonth])

    params <- format(params, digits = 3)
    params <- purrr::map2_chr(params, init_params, ~paste0(.x, " (", format(.y, digits = 3), ")"))
    message("Loss: ", loss, " Params: ", format(params, digits = 3), "\r", appendLF = FALSE)
    return(loss^2)
  }
  # Monthly calibration
  for(nmonth in nmonths){
    message("------------------------------------ Month: ", nmonth, " ------------------------------------ ")
    params <- unlist(model_cal$NM_model$coefficients[nmonth, c(2:6)])[-5]
    opt <- optim(params, loss_function, model_cal = model_cal, nmonth = nmonth,
                 init_params = params, conditional = conditional,
                 control = list(abstol = abstol, reltol = reltol))
    params <- purrr::map2_chr(opt$par, params, ~paste0(format(.x, digits = 3), " (", format(.y, digits = 3), ")"))
    message("Loss: ", opt$value, " Params: ", paste0(params, collapse = " "))
    params <- opt$par
    # Update the parameters
    model_cal$NM_model$model[[nmonth]]$update(means = params[1:2], sd = params[3:4])
  }
  # Update log-likelihood and conditional moments
  model_cal$update_moments()
  model_cal$update_logLik()
  return(model_cal)
}

