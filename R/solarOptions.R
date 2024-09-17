#' Option payoff function
#'
#' Compute the payoffs of an option at maturity.
#'
#' @param x numeric, vector of values at maturity.
#' @param strike numeric, option strike.
#' @param v0 numeric, price of the option.
#' @param put logical, when `TRUE`, the default, the payoff function is a put othewise a call.
#'
#' @examples
#' optionPayoff(10, 9, 1, put = TRUE)
#' mean(optionPayoff(seq(0, 20), 9, 1, put = TRUE))
#'
#' @rdname optionPayoff
#' @name optionPayoff
#' @export
optionPayoff <- function(x, strike = 0, c0 = 0, put = TRUE){
  payoff <- c()
  put <- ifelse(put, -1, 1)
  payoff <- (x - strike)*put
  payoff[payoff<0] <- 0
  payoff <- payoff - c0
  return(payoff)
}

#' Discount factor function
#'
#' @param r level of yearly constant risk-free rate
#' @param discrete logical, when `TRUE`, the default, discrete compounding will be used. Otherwise continuous compounding.
#'
#' @rdname discountFactor
#' @name discountFactor
#' @export
discountFactor <- function(r = 0.03, discrete = TRUE) {
  risk_free <- r
  if (discrete) {
    function(tau){
      (1 + risk_free/365)^(-tau)
    }
  } else {
    function(tau){
      exp(-risk_free/365*tau)
    }
  }
}

#' Control parameters for a solar option
#'
#' @param nyears numeric vector. Interval of years considered. The first element will be the minimum and the second the maximum years used in
#' the computation of the fair payoff.
#' @param K numeric, level for the strike with respect to the seasonal mean. The seasonal mean is multiplied by `exp(K)`.
#' @param put logical, when `TRUE`, the default, the computations will consider a `put` contract. Otherwise a `call`.
#' @param leap_year logical, when `FALSE`, the default, the year will be considered of 365 days, otherwise 366.
#' @param B function. Discount factor function. Should take as input a number (in years) and return a discount factor.
#'
#' @rdname control_solarOption
#' @name control_solarOption
#' @export
control_solarOption <- function(nyears = c(2005, 2023), K = 0, put = TRUE, leap_year = FALSE, B = discountFactor()){

  structure(
    list(
      nyears = nyears,
      from = as.Date(paste0(nyears[1], "-01-01")),
      to = as.Date(paste0(nyears[2], "-01-01")),
      K = K,
      put = put,
      leap_year = leap_year,
      B = B
    )
  )
}


#' Optimal number of contracts
#'
#' Compute the optimal number of contracts given a particular setup.
#'
#' @param model object with the class `solarModel`. See the function \code{\link{solarModel}} for details.
#' @param type character, method used for computing the premium. Can be `model` (Model with integral) or `sim` (Monte Carlo).
#' @param premium character, premium used. Can be `P`, `Qdw`, `Qup`, or `Q`.
#' @param nyear integer, actual year. The optimization will be performed excluding the year `nyear` and the following.
#' @param tick numeric, conversion tick for the monetary payoff of a contract.
#' @param efficiency numeric, mean efficiency of the solar panels.
#' @param n_panels numeric, number of meters squared of solar panels.
#' @param pun numeric, reference electricity price at which the energy is sold for computing the cash-flows.
#'
#' @rdname solarOption_contracts
#' @name solarOption_contracts
#'
#' @export
solarOption_contracts  <- function(model, type = "model", premium = "Q", nyear = 2021, tick = 0.06, efficiency = 0.2, n_panels = 2000, pun = 0.06){
  # All the payoffs
  payoffs <- model$payoffs
  # Extract historical payoff
  payoff_hist <- model$payoffs$hist$payoff
  # Match the type of computation
  type <- match.arg(type, names(payoffs))
  payoffs <- payoffs[[type]]
  # Match the type of premium
  type <- match.arg(premium, names(payoffs))
  # Extract daily premium
  payoffs <- payoffs[[premium]]$payoff_month_day
  # Select only relevant column
  payoffs <- dplyr::select(payoffs, Month, Day, premium)
  # Merge realized GHI and premium
  df_hedged <- dplyr::left_join(payoff_hist, payoffs, by = c("Month", "Day"))

  # Loss function depending on the number of contracts
  loss_function <- function(n_contracts, df_hedged){
    # Compute hedged cash-flows
    df_ <- dplyr::mutate(df_hedged, hedged = pun*n_panels*efficiency*GHI + tick*n_contracts*(payoff - premium))
    # Exclude nyear and the following from loss computation
    loss <- sd(dplyr::filter(df_, Year < nyear)$hedged)
    return(loss)
  }
  # Optimize the number of contracts
  opt <- optim(par = 10, fn = loss_function, method = "Brent",
               lower = 1, upper = n_panels*efficiency*10, df_hedged = df_hedged)

  list(
    tick = tick,
    efficiency = efficiency,
    n_panels = n_panels,
    pun = pun,
    nyear = nyear + 1,
    n_contracts = trunc(opt$par),
    type = type,
    premium = premium,
    sd = opt$value
  )
}


#' Payoff on Historical Data
#'
#' @param model object with the class `solarModel`. See the function \code{\link{solarModel}} for details.
#' @param nmonths numeric, months of which the payoff will be computed.
#' @param control_options control list, see \code{\link{control_solarOption}} for more details.
#'
#' @rdname solarOption_historical
#' @name solarOption_historical
#' @export
solarOption_historical <- function(model, nmonths = 1:12, control_options = control_solarOption()){

  # Control parameters
  nyears <- control_options$nyears
  K <- control_options$K
  put <- control_options$put
  leap_year <- control_options$leap_year
  # Complete data
  data <- model$data
  # Seasonal data
  seasonal_data <- dplyr::select(model$seasonal_data, Month, Day, GHI_bar)

  # Historical Daily Payoffs
  data <- dplyr::left_join(data, seasonal_data, by = c("Month", "Day"))
  # Select only relevant variables
  df_payoff <- dplyr::select(data, date, n, Year, Month, Day, GHI, GHI_bar)
  # Filter for control years
  df_payoff <- dplyr::filter(df_payoff, Year >= nyears[1] & Year <= nyears[2])
  # Filter for selected months
  df_payoff <- dplyr::filter(df_payoff, Month %in% nmonths)
  # Daily option strike
  df_payoff$strike <- df_payoff$GHI_bar*exp(K)
  # Daily option payoffs
  if (put) {
    # Payoff for a put
    df_payoff$exercise <- ifelse(df_payoff$GHI < df_payoff$strike, 1, 0)
    df_payoff$payoff <- (df_payoff$strike - df_payoff$GHI)*df_payoff$exercise
    df_payoff$side <- "put"
  } else {
    # Payoff for a call
    df_payoff$exercise <- ifelse(df_payoff$GHI > df_payoff$strike, 1, 0)
    df_payoff$payoff <- (df_payoff$GHI - df_payoff$strike)*df_payoff$exercise
    df_payoff$side <- "call"
  }

  # Aggregation of Historical Payoffs
  # Reorder variables
  df_payoff <- dplyr::select(df_payoff, side, date, Year, Month, Day, n, GHI, strike, payoff, exercise)

  # Aggregation by Year and Month
  df_year_month <- df_payoff %>%
    dplyr::group_by(Year, Month, side) %>%
    dplyr::reframe(ndays = dplyr::n(),
                   payoff = sum(payoff),
                   exercise = mean(exercise),
                   GHI = sum(GHI),
                   strike = sum(strike))

  # Include or not 29-th of February from computation
  if (leap_year) {
    df_payoff_ <- df_payoff
  } else {
    df_payoff_ <- dplyr::filter(df_payoff, !(Month == 2 & Day == 29))
  }

  # Aggregation by Month and Day
  df_month_day <- df_payoff_ %>%
    dplyr::group_by(Month, Day, side) %>%
    dplyr::reframe(premium = mean(payoff),
                   exercise = mean(exercise),
                   GHI = mean(GHI),
                   strike = mean(strike))
  # Add number of the day
  df_month_day$n <- 1:nrow(df_month_day)

  # Aggregation by Month
  df_month <- df_month_day %>%
    dplyr::group_by(Month, side)  %>%
    dplyr::reframe(ndays = dplyr::n(),
                   premium = sum(premium),
                   daily_premium = premium/ndays,
                   exercise = mean(exercise),
                   GHI = sum(GHI),
                   strike = sum(strike))

  # Aggregation by Year
  df_year <- df_month %>%
    dplyr::group_by(side)  %>%
    dplyr::reframe(ndays = sum(ndays),
                   premium = sum(premium),
                   exercise = mean(exercise),
                   GHI = sum(GHI),
                   strike = sum(strike))

  # Structured output
  structure(
    list(
      payoff = df_payoff,
      payoff_year_month = df_year_month,
      payoff_month_day = df_month_day,
      payoff_month = df_month,
      payoff_year = df_year
    ),
    class = c("solarOption", "list")
  )
}


#' Payoff on Simulated Data
#'
#' @param sim simulated scenarios with the function \code{\link{solarModel_scenarios}}.
#' @param nmonths numeric, months of which the payoff will be computed.
#' @param nsim number of simulation to use for computation.
#' @param control_options control function, see \code{\link{control_solarOption}} for details.
#'
#' @rdname solarOption_scenario
#' @name solarOption_scenario
#' @export
solarOption_scenario <- function(sim, nmonths = 1:12, nsim, control_options = control_solarOption()){

  # Control parameters
  nyears <- control_options$nyears
  K <- control_options$K
  put <- control_options$put
  leap_year = control_options$leap_year

  # Simulated Daily Payoffs
  # Filter for control years
  df_payoff <- dplyr::filter(sim, Year >= nyears[1] & Year <= nyears[2])
  # Filter for selected months
  df_payoff <- dplyr::filter(df_payoff, Month %in% nmonths)
  # Eventually reduce the number of simulations used
  if (!missing(nsim)) {
    nsim <- min(c(nsim, nrow(df_payoff$data[[1]])))
    df_payoff <- dplyr::mutate(df_payoff, data = purrr::map(data, ~.x[1:nsim,]))
  }
  # Compute simulated payoffs
  df_payoff <- df_payoff %>%
    dplyr::mutate(strike = purrr::map_dbl(data, ~.x$GHI_bar[1]*exp(K)),
                  side = ifelse(put, "put", "call")) %>%
    dplyr::mutate(
      data = purrr::map(data, ~dplyr::mutate(.x, strike = GHI_bar*exp(K))),
      data = purrr::map(data, ~dplyr::mutate(.x, side = ifelse(put, "put", "call"))),
      data = purrr::map(data, ~dplyr::mutate(.x, exercise =
                                               dplyr::case_when(side == "put" & GHI <= strike ~ 1,
                                                                side == "put" & GHI > strike ~ 0,
                                                                side == "call" & GHI <= strike ~ 0,
                                                                side == "call" & GHI > strike ~ 1))),
      data = purrr::map(data, ~dplyr::mutate(.x, payoff =
                                               dplyr::case_when(side == "put" ~ (strike - GHI)*exercise,
                                                                side == "call" ~ (GHI - strike)*exercise)))) %>%
    dplyr::mutate(payoff = purrr::map_dbl(data, ~mean(.x$payoff)),
                  GHI = purrr::map_dbl(data, ~mean(.x$GHI)),
                  exercise = purrr::map_dbl(data, ~mean(.x$exercise))) %>%
    dplyr::select(side, date, Year, Month, Day, n, GHI, strike, payoff, exercise)

  # Aggregation of Simulated Payoffs
  # Aggregation by Year and Month
  df_year_month <- df_payoff %>%
    dplyr::group_by(Year, Month, side) %>%
    dplyr::reframe(premium = sum(payoff),
                   exercise = mean(exercise),
                   GHI = sum(GHI),
                   strike = sum(strike),
                   ndays = dplyr::n())

  if (leap_year) {
    df_payoff_ <- df_payoff
  } else {
    df_payoff_ <- dplyr::filter(df_payoff, !(Month == 2 & Day == 29))
  }

  # Aggregation by Month and Day
  df_month_day <- df_payoff_ %>%
    dplyr::group_by(Month, Day, side) %>%
    dplyr::reframe(premium = mean(payoff),
                   exercise = mean(exercise),
                   GHI = mean(GHI),
                   strike = mean(strike))
  df_month_day$n <- 1:nrow(df_month_day)
  # Aggregation by Month
  df_month <- df_month_day %>%
    dplyr::group_by(Month, side) %>%
    dplyr::reframe(ndays = dplyr::n(),
                   premium = sum(premium),
                   daily_premium = premium/ndays,
                   exercise = mean(exercise),
                   GHI = sum(GHI),
                   strike = sum(strike))
  # Aggregation by Year
  df_year <- df_month %>%
    dplyr::group_by(side) %>%
    dplyr::reframe(ndays = sum(ndays),
                   premium = sum(premium),
                   exercise = mean(exercise),
                   GHI = sum(GHI),
                   strike = sum(strike))
  structure(
    list(
      payoff = df_payoff,
      payoff_year_month = df_year_month,
      payoff_month_day = df_month_day,
      payoff_month = df_month,
      payoff_year = df_year
    ),
    class = c("solarOption", "list")
  )
}


#' Bootstrap a fair premium from historical data
#'
#' @param model object with the class `solarModel`. See the function \code{\link{solarModel}} for details.
#' @param nsim number of simulation to bootstrap.
#' @param ci confidence interval for quantile
#' @param seed random seed.
#' @param control_options control function, see \code{\link{control_solarOption}} for details.
#'
#' @return An object of the class `solarOptionPayoffBoot`.
#'
#' @rdname solarOption_bootstrap
#' @name solarOption_bootstrap
#' @export
solarOption_bootstrap <- function(model, nsim = 500, ci = 0.05, seed = 1, control_options = control_solarOption()){

  leap_year <- control_options$leap_year

  # 0) Historical monthly quantiles
  # Compute historical payoffs
  payoff_hist <- solarOption_historical(model, nmonths = 1:12, control_options = control_options)
  # Group by month and day and create an empirical quantile function for each group
  df_quantiles <- payoff_hist$payoff %>%
    dplyr::group_by(Month) %>%
    dplyr::select(Month, Day, payoff) %>%
    tidyr::nest() %>%
    dplyr::mutate(quant = purrr::map(data, ~function(u) quantile(.x$payoff, probs = u)),
                  ndays = lubridate::days_in_month(Month)) %>%
    dplyr::select(-data)

  if (leap_year){
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
  set.seed(seed)
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
    class = c("solarOption", "list")
  )
}


#' Pricing function under the solar model
#'
#' @inheritParams solarOption_historical
#' @param theta Esscher parameter
#' @param implvol implied unconditional GARCH variance, the default is `1`.
#' @examples
#' control <- control_solarModel(outliers_quantile = 0.0005)
#' spec <- solarModel_spec("Berlino", from="2005-01-01", to="2024-01-01", control_model = control)
#' model <- solarModel(spec)
#' solarOption_model(model)
#' solarOption_historical(model)
#' @rdname solarOption_model
#' @name solarOption_model
#' @export
solarOption_model <- function(model, nmonths = 1:12, theta = 0, implvol = 1, control_options = control_solarOption()){

  # Options control
  K <- control_options$K
  put <- control_options$put
  leap_year = control_options$leap_year
  alpha_ <- model$params$alpha
  beta_ <- model$params$beta

  # AR(2) stationary variance
  par <- model$AR_model_Yt$coefficients
  ar_variance <- (1-par[2])/((1 - par[2])*(1 - par[1]^2 - par[2]^2) - 2*par[1]^2*par[2])
  # Seasonal data
  seasonal_data <- dplyr::filter(model$seasonal_data, Month %in% nmonths)
  if (!leap_year) {
    seasonal_data <- dplyr::filter(seasonal_data, !(Month == 2 & Day == 29))
  }
  # Option Strike
  seasonal_data$strike <- seasonal_data$GHI_bar*exp(K)
  # Add seasonal parameters
  GM_model <- dplyr::select(model$NM_model, Month, mu_up:p_dw)
  seasonal_data <- dplyr::left_join(seasonal_data, GM_model, by = "Month")
  seasonal_data <- dplyr::left_join(seasonal_data, model$monthly_data, by = "Month")

  # Add Esscher parameter
  if (is.function(theta)) {
    seasonal_data$theta <- purrr::map_dbl(seasonal_data$Month, ~theta(.x))
  } else {
    seasonal_data$theta <- theta
  }
  # Compute mixture parameters
  seasonal_data <- dplyr::mutate(seasonal_data,
                                 mu_up = Yt_bar + Yt_tilde_uncond + sigma_bar*sigma_m*mu_up*sqrt(ar_variance),
                                 mu_dw = Yt_bar + Yt_tilde_uncond + sigma_bar*sigma_m*mu_dw*sqrt(ar_variance),
                                 sd_up = implvol*sqrt(ar_variance)*sd_up*sigma_bar*sigma_m,
                                 sd_dw = implvol*sqrt(ar_variance)*sd_dw*sigma_bar*sigma_m)

  # Select only necessary variables
  seasonal_data <- dplyr::select(seasonal_data, Month, Day, Ct, GHI_bar, strike:p_dw, theta)

  # Pricing function
  pricing_month_day <- function(model, seasonal_data, nmonth = 1, nday = 1, put = TRUE){
    #message("Month: ", nmonth, " Day: ", nday)
    # Seasonal Data
    df_n <- dplyr::filter(seasonal_data, Month == nmonth & Day == nday)

    # Mixture Pdf
    pdf_Yt <- desscherMixture(means = c(df_n$mu_up, df_n$mu_dw), sd = c(df_n$sd_up, df_n$sd_dw), p = c(df_n$p_up, 1-df_n$p_up), theta = df_n$theta)
    # Density for GHI
    pdf_GHI <- function(x) dsolarGHI(x, df_n$Ct, alpha_, beta_, pdf_Yt)
    # Esscher density
    pdf_GHI <- desscher(pdf_GHI, theta = df_n$theta, lower = df_n$Ct*(1-alpha_-beta_), upper = df_n$Ct*(1-alpha_))
    # Distribution for GHI
    cdf_GHI <- CDF(pdf_GHI, lower = df_n$Ct*(1-alpha_-beta_))
    # Expected value function (GHI)
    e_GHI <- function(lower, upper){
      if (missing(lower)){
        lower <- df_n$Ct*(1-alpha_-beta_)
      }
      if (missing(upper)){
        upper <- df_n$Ct*(1-alpha_)
      }
      integrate(function(x) x*pdf_GHI(x), lower = lower, upper = upper)$value
    }

    # Option pricing
    if (put) {
      # Expected value (Put)
      df_n$GHI_plus <- e_GHI(upper = df_n$strike)
      # Probability of exercise
      df_n$exercise <- cdf_GHI(df_n$strike)
      # Option expected value
      df_n$premium <- df_n$strike*df_n$exercise - df_n$GHI_plus
      # Option type
      df_n$side <- "put"
    } else {
      # Expected value (Call)
      df_n$GHI_plus <-  e_GHI(lower = df_n$strike)
      # Probability of exercise
      df_n$exercise <- 1 - cdf_GHI(df_n$strike)
      # Option expected value
      df_n$premium <- df_n$GHI_plus - df_n$strike*df_n$exercise
      # Option type
      df_n$side <- "call"
    }
    # Expected value
    df_n$GHI <- e_GHI()
    # Select only relevant variables
    df_n <- dplyr::select(df_n, Month, Day, side, premium, exercise, GHI_plus, strike, GHI)
    return(df_n)
  }

  # Model premium for each day and month
  df_month_day <- purrr::map2_df(seasonal_data$Month, seasonal_data$Day,
                                 ~pricing_month_day(model, seasonal_data, nmonth = .x, nday = .y, put = put))
  day_date <- paste0(ifelse(leap_year, "2020-", "2019-"), df_month_day$Month, "-", df_month_day$Day)
  df_month_day$n <- number_of_day(day_date)

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
    class = c("solarOption", "list")
  )
}


#' Test errors solar Option model
#'
#'@export
solarOption_model_test <- function(model, control_options = control_solarOption()){

  payoff_model <- solarOption_model(model, control_options = control_options)$payoff_month
  payoff_hist <- solarOption_historical(model, control_options = control_options)$payoff_month

  # Monthly payoff
  payoff_month <- dplyr::tibble(
    Month = payoff_model$Month,
    Model = payoff_model$premium,
    Hist = payoff_hist$premium,
    Diff = Hist - Model,
    Error = round(Diff/Hist*100, 3)
  )
  # Yearly payoff
  payoff_year <- dplyr::summarise_all(payoff_month,  sum)
  payoff_year <- dplyr::mutate(payoff_year, Month = "Total", Error = round(Diff/Hist*100, 3))
  # Tranform months in character
  payoff_month$Month <- as.character(lubridate::month(payoff_month$Month, label = TRUE))
  dplyr::bind_rows(payoff_month, payoff_year)
}


#' Pricing function under the solar model
#'
#' @inheritParams solarOption_model
#' @inheritParams spatialModel_neighborhoods
#'
#' @rdname solarOption_model_spatial
#' @name solarOption_model_spatial
#' @export
solarOption_model_spatial <- function(object, lat, lon, nmonths = 1:12, theta = 0, implvol = 1, control_options = control_solarOption()){

  # Options control
  K <- control_options$K
  put <- control_options$put
  leap_year = control_options$leap_year
  # Interpolate the model
  model <- spatialModel_interpolate(object, lat, lon, n = 4, beta = 2)
  # Seasonal data
  seasonal_data <- dplyr::filter(model$seasonal_data, Month %in% nmonths)
  if (!leap_year) {
    seasonal_data <- dplyr::filter(seasonal_data, !(Month == 2 & Day == 29))
  }
  seasonal_data <- dplyr::left_join(seasonal_data, model$monthly_data, by = c("Month"))

  # Compute all possible combinations
  all_combinations <- spatialModel_combinations(object, lat, lon)

  # AR(2) stationary variance
  par <- model$AR_model_Yt$coefficients
  ar_variance <- (1-par[2])/((1 - par[2])*(1 - par[1]^2 - par[2]^2) - 2*par[1]^2*par[2])

  # Option Strike
  seasonal_data$strike <- seasonal_data$GHI_bar*exp(K)

  if (is.function(theta)) {
    seasonal_data$theta <- purrr::map_dbl(seasonal_data$Month, ~theta(.x))
  } else {
    seasonal_data$theta <- theta
  }
  # Compute mixture parameters
  seasonal_data <- dplyr::mutate(seasonal_data,
                                 z_tk = (1 - model$params$alpha - strike/Ct)/model$params$beta,
                                 z = log(-log(z_tk)))

  # Pricing function
  pricing_month_day <- function(model, seasonal_data, nmonth = 1, nday = 1, put = TRUE){

    message("Month: ", nmonth, " Day: ", nday, "\r", appendLF = FALSE)
    combinations <- dplyr::filter(all_combinations, Month == nmonth)
    # Seasonal Data
    df_n <- dplyr::filter(seasonal_data, Month == nmonth & Day == nday)
    # Parameters
    combinations$mean <- df_n$Yt_bar + df_n$Yt_tilde_uncond + sqrt(ar_variance)*df_n$sigma_bar*df_n$sigma_m*combinations$mean
    combinations$sd <- implvol*sqrt(ar_variance)*df_n$sigma_bar*df_n$sigma_m*combinations$sd
    # Mixture Pdf
    mixture_pdf <- desscherMixture(means = combinations$mean, sd = combinations$sd, p = combinations$probs, theta = df_n$theta)

    # Expected value function (GHI)
    e_GHI <- function(x){
      df_n$Ct*(1 - model$transform$Yt(x, inverse = TRUE))*mixture_pdf(x)
    }

    # Option pricing
    if (put) {
      # Expected value (Put)
      df_n$GHI_plus <- integrate(e_GHI, lower = -Inf, upper = df_n$z)$value
      # Probability of exercise
      df_n$exercise <- integrate(mixture_pdf, lower = -Inf, upper = df_n$z)$value
      # Option expected value
      df_n$premium <- df_n$strike*df_n$exercise - df_n$GHI_plus
      # Option type
      df_n$side <- "put"
    } else {
      # Expected value (Call)
      df_n$GHI_plus <- integrate(e_GHI, lower = df_n$z, upper = Inf)$value
      # Probability of exercise
      df_n$exercise <- integrate(mixture_pdf, lower = df_n$z, upper = Inf)$value
      # Option expected value
      df_n$premium <- df_n$GHI_plus - df_n$strike*df_n$exercise
      # Option type
      df_n$side <- "call"
    }
    # Expected value
    df_n$GHI <- integrate(e_GHI, lower = -Inf, upper = Inf)$value
    # Select only relevant variables
    df_n <- dplyr::select(df_n, Month, Day, side, premium, exercise, GHI_plus, strike, GHI)
    return(df_n)
  }

  # Model premium for each day and month
  df_month_day <- purrr::map2_df(seasonal_data$Month, seasonal_data$Day,
                                 ~pricing_month_day(model, seasonal_data, nmonth = .x, nday = .y, put = put))
  day_date <- paste0(ifelse(leap_year, "2020-", "2019-"), df_month_day$Month, "-", df_month_day$Day)
  df_month_day$n <- number_of_day(day_date)

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
    class = c("solarOption", "list")
  )
}

#' Implied expected return at maturity
#'
#' @export
solarOption_implied_return <- function(model, target_prices = NA, nmonths = 1:12, control_options = control_solarOption()){

  # Control parameters
  quiet <- model$control$quiet

  # Default target price is the historical price
  if (is.na(target_prices)) {
    target_prices <- model$payoffs$hist$payoff_month$premium
    stopifnot(!is.null(target_prices))
  }

  # Loss function
  loss_function <- function(r, nmonth = 1, target_price = NA){
    stopifnot(!is.na(target_price))
    premium_model <- solarOption_model(model, nmonths = nmonth, theta = model$esscher$theta(r), control_options = control_options)
    l <- (target_price - premium_model$payoff_month$premium[1])^2
    if(!quiet) message("Loss: ", round(l, 10), ", r (implied): ", format(r*100, digits = 5), " %")
    return(l)
  }

  # Compute implied returns
  implied_r <- c()
  for(m in nmonths){
    if(!quiet) message("--------------------- Month: ", m, "---------------------")
    implied_r[m] <- optim(par = 0, loss_function, lower = -0.5, upper = 0.5, method = "Brent",
                          nmonth = nmonths[m], target_price = target_prices[m])$par
  }

  dplyr::tibble(
    Month = nmonths,
    implied_r = implied_r
  )
}


#' Structure payoffs
#'
#' @param model object with the class `solarModel`. See the function \code{\link{solarModel}} for details.
#' @param type method used for computing the premium. If `model`, the default will be used the analytic model,
#' otherwise with `sim` the monte carlo scenarios stored inside the `model$scenarios$P`.
#' @param exact_daily_premium when `TRUE` the historical premium is computed as daily average among all the years.
#' Otherwise the monthly premium is computed and then divided by the number of days of the month.
#'
#' @rdname solarOption_structure
#' @name solarOption_structure
#' @export
solarOption_structure <- function(model, type = "model", exact_daily_premium = TRUE){

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
    # Probabilities of exercise the option
    exercise = payoff$hist$payoff_year$exercise,
    exercise_P = payoff[[type]]$P$payoff_year$exercise,
    exercise_Q = payoff[[type]]$Q$payoff_year$exercise,
    exercise_Qup = payoff[[type]]$Qup$payoff_year$exercise,
    exercise_Qdw = payoff[[type]]$Qdw$payoff_year$exercise,
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
    # Probabilities of exercise the option
    exercise = payoff$hist$payoff_month$exercise,
    exercise_P = payoff[[type]]$P$payoff_month$exercise,
    exercise_Q = payoff[[type]]$Q$payoff_month$exercise,
    exercise_Qup = payoff[[type]]$Qup$payoff_month$exercise,
    exercise_Qdw = payoff[[type]]$Qdw$payoff_month$exercise,
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
    # Probabilities of exercise the option
    exercise = payoff$hist$payoff_month_day$exercise,
    exercise_P = payoff[[type]]$P$payoff_month_day$exercise,
    exercise_Q = payoff[[type]]$Q$payoff_month_day$exercise,
    exercise_Qup = payoff[[type]]$Qup$payoff_month_day$exercise,
    exercise_Qdw = payoff[[type]]$Qdw$payoff_month_day$exercise,
    n = payoff$hist$payoff_month_day$n
  )

  # Historical Daily Premiums
  df_payoff <- dplyr::select(payoff$hist$payoff, -exercise, -GHI, -strike)
  if (!exact_daily_premium) {
    df_payoff <- dplyr::left_join(df_payoff, dplyr::select(df_month_day_mean, Month, premium:premium_Qdw), by = c("Month"))
  } else {
    df_payoff <- dplyr::left_join(df_payoff, dplyr::select(df_month_day, Month, Day, premium:premium_Qdw), by = c("Month", "Day"))
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
                            cum_net_payoff_Qup = NA)
    ndays <- nrow(df_cum)
    for(i in 1:nrow(df_cum)){
      df_cum$cum_net_payoff[i] <- (df_year$premium + sum(df_cum$payoff[1:i]) - sum(df_cum$premium[1:i]))
      df_cum$cum_net_payoff_P[i] <- df_year$premium_P + sum(df_cum$payoff[1:i] - df_cum$premium_P[1:i])
      df_cum$cum_net_payoff_Q[i] <- df_year$premium_Q + sum(df_cum$payoff[1:i] - df_cum$premium_Q[1:i])
      df_cum$cum_net_payoff_Qdw[i] <- df_year$premium_Qdw + sum(df_cum$payoff[1:i] - df_cum$premium_Qdw[1:i])
      df_cum$cum_net_payoff_Qup[i] <- df_year$premium_Qup + sum(df_cum$payoff[1:i] - df_cum$premium_Qup[1:i])
    }
    cumulated_payoff[[j]] <- df_cum
  }

  # Compute the expected trajectory for each day
  df_cum <- dplyr::bind_rows(cumulated_payoff) %>%
    dplyr::group_by(Month, Day) %>%
    dplyr::mutate(
      e_cum_net_payoff = mean(cum_net_payoff),
      e_cum_net_payoff_P = mean(cum_net_payoff_P),
      e_cum_net_payoff_Q = mean(cum_net_payoff_Q),
      e_cum_net_payoff_Qdw = mean(cum_net_payoff_Qdw),
      e_cum_net_payoff_Qup = mean(cum_net_payoff_Qup),
    ) %>%
    dplyr::ungroup()
  model$payoffs[[type]]$structured <- list(payoff = df_payoff,
                                           payoff_year = df_year,
                                           payoff_month = df_month,
                                           payoff_month_day = df_month_day,
                                           payoff_cum = df_cum)
  return(model)
}


