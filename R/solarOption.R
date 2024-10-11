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
control_solarOption <- function(nyears = c(2005, 2023), K = 0, leap_year = FALSE, nsim = 200, ci = 0.05, seed = 1, B = discountFactor()){

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


#' solarOptionPayoff
#'
#' @param model solarModel
#' @param control_options control list, see \code{\link{control_solarOption}} for more details.
#'
#' @return An object of the class `solarOptionPayoffs`.
#'
#' @rdname solarOptionPayoffs
#' @name solarOptionPayoffs
#' @export
solarOptionPayoffs <- function(model, control_options = control_solarOption()){

  # Initialize a list for call options
  payoff_call = list(
    historical = solarOption_historical(model, put = FALSE, control_options = control_options),
    scenarios = list(P = NA, Q = NA, Qup = NA, Qdw = NA, Qr = NA, structured = NA),
    model = list(P = NA, Q = NA, Qup = NA, Qdw = NA, Qr = NA, boot = NA, structured = NA)
  )
  # Initialize a list for put options
  payoff_put = list(
    historical = solarOption_historical(model,  put = TRUE, control_options = control_options),
    scenarios = list(P = NA, Q = NA, Qup = NA, Qdw = NA, Qr = NA, structured = NA),
    model = list(P = NA, Q = NA, Qup = NA, Qdw = NA, Qr = NA, boot = NA, structured = NA)
  )

  structure(
    list(
      call = payoff_call,
      put = payoff_put,
      esscher = list(),
      control_options = control_options
    ),
    class = c("solarOptionPayoffs", "list")
  )
}


#' Payoff on Historical Data
#'
#' @param model object with the class `solarModel`. See the function \code{\link{solarModel}} for details.
#' @param nmonths numeric vector of months in which the payoff is computed. Range from 1 to 12.
#' @param put logical, when `TRUE`, the default, the computations will consider a `put` contract. Otherwise a `call`.
#' @param control_options control list, see \code{\link{control_solarOption}} for more details.
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
  nyears <- control_options$nyears
  K <- control_options$K
  leap_year <- control_options$leap_year
  # Target and seasonal mean
  target <- model$target
  target_bar <- paste0(model$target, "_bar")
  # Complete data
  data <- dplyr::select(model$data, date, n, Year, Month, Day, tidyr::any_of(c(target, target_bar)))

  # Filter for control years
  data <- dplyr::filter(data, Year >= nyears[1] & Year <= nyears[2])
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
    df_payoff$payoff <- (df_payoff$strike - df_payoff$Rt)*df_payoff$exercise
    df_payoff$side <- "put"
  } else {
    # Payoff for a call
    df_payoff$exercise <- ifelse(df_payoff$Rt > df_payoff$strike, 1, 0)
    df_payoff$payoff <- (df_payoff$Rt - df_payoff$strike)*df_payoff$exercise
    df_payoff$side <- "call"
  }

  # Aggregation of Historical Payoffs
  # Reorder variables
  df_payoff <- dplyr::select(df_payoff, side, date, Year, Month, Day, n, Rt, strike, payoff, exercise)

  # Aggregation by Year and Month
  df_year_month <- df_payoff %>%
    dplyr::group_by(Year, Month, side) %>%
    dplyr::reframe(ndays = dplyr::n(),
                   payoff = sum(payoff),
                   exercise = mean(exercise),
                   Rt = sum(Rt),
                   strike = sum(strike))

  # Include or not 29-th of February from computation
  if (leap_year) {
    df_payoff_ <- df_payoff
    data_month_day <- dplyr::select(model$seasonal_data, Month, Day, n)
  } else {
    df_payoff_ <- dplyr::filter(df_payoff, !(Month == 2 & Day == 29))
    data_month_day <- dplyr::select(model$seasonal_data, Month, Day, n)
    data_month_day <- dplyr::filter(data_month_day, !(Month == 2 & Day == 29))
  }

  # Aggregation by Month and Day
  df_month_day <- df_payoff_ %>%
    dplyr::group_by(Month, Day, side) %>%
    dplyr::reframe(premium = mean(payoff),
                   exercise = mean(exercise),
                   Rt = mean(Rt),
                   strike = mean(strike)) %>%
    dplyr::right_join(data_month_day, by = c("Month", "Day"))

  # Aggregation by Month
  df_month <- df_month_day %>%
    dplyr::group_by(Month, side)  %>%
    dplyr::reframe(ndays = dplyr::n(),
                   premium = sum(premium),
                   daily_premium = premium/ndays,
                   exercise = mean(exercise),
                   Rt = sum(Rt),
                   strike = sum(strike))

  # Aggregation by Year
  df_year <- df_month %>%
    dplyr::group_by(side)  %>%
    dplyr::reframe(ndays = sum(ndays),
                   premium = sum(premium),
                   exercise = mean(exercise),
                   Rt = sum(Rt),
                   strike = sum(strike))

  # Rename columns to match target name
  colnames(df_payoff)[colnames(df_payoff) == "Rt"] <- target
  colnames(df_year_month)[colnames(df_year_month) == "Rt"] <- target
  colnames(df_month_day)[colnames(df_month_day) == "Rt"] <- target
  colnames(df_month)[colnames(df_month) == "Rt"] <- target
  colnames(df_year)[colnames(df_year) == "Rt"] <- target

  # Structured output
  structure(
    list(
      payoff = df_payoff,
      payoff_year_month = df_year_month,
      payoff_month_day = df_month_day,
      payoff_month = df_month,
      payoff_year = df_year
    ),
    class = c("solarOptionPayoff", "list")
  )
}


#' Bootstrap a fair premium from historical data
#'
#' @inheritParams solarOption_historical
#'
#' @return An object of the class `solarOptionBoot`.
#' @examples
#' model <- Bologna
#' solarOption_historical_bootstrap(model)
#'
#' @rdname solarOption_historical_bootstrap
#' @name solarOption_historical_bootstrap
#' @export
solarOption_historical_bootstrap <- function(model, put = TRUE, control_options = control_solarOption()){

  leap_year <- control_options$leap_year
  ci = control_options$ci
  nsim = control_options$nsim
  seed = control_options$seed

  # 0) Historical monthly quantiles
  # Compute historical payoffs
  payoff_hist <- solarOption_historical(model, nmonths = 1:12, put = put, control_options = control_options)
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
      side = ifelse(put, "put", "call"),
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
      side = ifelse(put, "put", "call"),
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
#' scenario <- solarScenario(model, from = "2011-01-01", to = "2012-01-01", by = "1 month", nsim = 1, seed = 3)
#' solarOption_scenario(scenario)
#' solarOption_historical(model)
#'
#' @rdname solarOption_scenario
#' @name solarOption_scenario
#' @export
solarOption_scenario <- function(scenario, nmonths = 1:12, put = TRUE, nsim, control_options = control_solarOption()){

  # Control parameters
  nyears <- control_options$nyears
  K <- control_options$K
  leap_year = control_options$leap_year
  # Target and seasonal mean
  target <- scenario$target
  target_bar <- paste0(target, "_bar")

  # Simulated Daily Payoffs
  sim <- scenario$sim
  # Filter for control years
  df_payoff <- dplyr::filter(sim, Year >= nyears[1] & Year <= nyears[2])
  # Filter for selected months
  df_payoff <- dplyr::filter(df_payoff, Month %in% nmonths)
  # Eventually reduce the number of simulations used
  if (!missing(nsim)) {
    nsim <- min(c(nsim, nrow(df_payoff$data[[1]])))
    df_payoff <- dplyr::mutate(df_payoff, data = purrr::map(data, ~.x[1:nsim,]))
  }

  df_payoff <- tidyr::unnest(df_payoff, cols = "data")
  colnames(df_payoff)[colnames(df_payoff) %in% c(target, target_bar)] <- c("Rt", "Rt_bar")

  # Compute simulated payoffs
  df_payoff <- df_payoff %>%
    dplyr::mutate(strike = Rt_bar*exp(K),
                  side = ifelse(put, "put", "call"),
                  exercise = dplyr::case_when(side == "put" & Rt <= strike ~ 1,
                                              side == "put" & Rt > strike ~ 0,
                                              side == "call" & Rt <= strike ~ 0,
                                              side == "call" & Rt > strike ~ 1),
                  payoff_sim = dplyr::case_when(side == "put" ~ (strike - Rt)*exercise,
                                                side == "call" ~ (Rt - strike)*exercise)) %>%
    dplyr::group_by(side, date, Year, Month, Day) %>%
    dplyr::reframe(
      Rt = mean(Rt),
      strike = mean(strike),
      payoff = mean(payoff_sim),
      exercise = mean(exercise)
    ) %>%
    dplyr::ungroup()

  # Aggregation of Simulated Payoffs
  # Aggregation by Year and Month
  df_year_month <- df_payoff %>%
    dplyr::group_by(Year, Month, side) %>%
    dplyr::reframe(premium = sum(payoff),
                   exercise = mean(exercise),
                   Rt = sum(Rt),
                   strike = sum(strike),
                   ndays = dplyr::n())

  # Include or not 29-th of February from computation
  if (leap_year) {
    df_payoff_ <- df_payoff
  } else {
    df_payoff_ <- dplyr::filter(df_payoff, !(Month == 2 & Day == 29))
  }

  # Aggregation by Month and Day
  df_month_day <- df_payoff_ %>%
    dplyr::group_by(Month, Day, side) %>%
    dplyr::reframe(n = dplyr::n(),
                   premium = mean(payoff),
                   exercise = mean(exercise),
                   Rt = mean(Rt),
                   strike = mean(strike)) %>%
    dplyr::filter(!is.na(premium))
  day_date <- paste0(ifelse(leap_year, "2020-", "2019-"), df_month_day$Month, "-", df_month_day$Day)
  df_month_day$n <- number_of_day(day_date)


  # Aggregation by Month
  df_month <- df_month_day %>%
    dplyr::group_by(Month, side) %>%
    dplyr::reframe(ndays = dplyr::n(),
                   premium = sum(premium),
                   daily_premium = premium/ndays,
                   exercise = mean(exercise),
                   Rt = sum(Rt),
                   strike = sum(strike))

  # Aggregation by Year
  df_year <- df_month %>%
    dplyr::group_by(side) %>%
    dplyr::reframe(ndays = sum(ndays),
                   premium = sum(premium),
                   exercise = mean(exercise),
                   Rt = sum(Rt),
                   strike = sum(strike))

  # Rename columns to match target name
  colnames(df_payoff)[colnames(df_payoff) == "Rt"] <- target
  colnames(df_year_month)[colnames(df_year_month) == "Rt"] <- target
  colnames(df_month_day)[colnames(df_month_day) == "Rt"] <- target
  colnames(df_month)[colnames(df_month) == "Rt"] <- target
  colnames(df_year)[colnames(df_year) == "Rt"] <- target

  structure(
    list(
      payoff = df_payoff,
      payoff_year_month = df_year_month,
      payoff_month_day = df_month_day,
      payoff_month = df_month,
      payoff_year = df_year
    ),
    class = c("solarOptionPayoff", "list")
  )
}


#' Pricing function under the solar model
#'
#' @inheritParams solarOption_historical
#' @param theta Esscher parameter
#' @param implvol implied unconditional GARCH variance, the default is `1`.
#' @param combinations list of 12 elements with gaussian mixture components.
#' @param target.Yt logical, when `TRUE`, the default, the computations will consider the pdf of `Yt` otherwise the pdf of solar radiation.
#'
#' @return An object of the class `solarOptionPayoff`.
#' @examples
#' model <- Bologna
#' solarOption_model(model, put=FALSE)
#' solarOption_model(model, put=TRUE)
#'
#' @rdname solarOption_model
#' @name solarOption_model
#' @export
solarOption_model <- function(model, nmonths = 1:12, theta = 0, combinations = NA, implvol = 1, put = TRUE, target.Yt = TRUE, control_options = control_solarOption()){

  # Options control
  K <- control_options$K
  leap_year = control_options$leap_year
  # Target and seasonal mean
  target <- model$target
  target_bar <- paste0(model$target, "_bar")
  target_plus <- paste0(model$target, "_plus")
  # AR(2) stationary variance
  par <- model$AR_model_Yt$coefficients
  ar_variance <- (1-par[2])/((1 - par[2])*(1 - par[1]^2 - par[2]^2) - 2*par[1]^2*par[2])

  # Extract seasonal data
  seasonal_data <- model$seasonal_data
  # Filter for the selected months
  seasonal_data <- dplyr::filter(seasonal_data, Month %in% nmonths)
  # Remove 29 of February from computations
  if (!leap_year) {
    seasonal_data <- dplyr::filter(seasonal_data, !(Month == 2 & Day == 29))
  }
  # Option Strike
  seasonal_data$strike <- seasonal_data[[target_bar]]*exp(K)
  # Add upper and lower bounds and unconditional moments
  seasonal_data <- dplyr::mutate(seasonal_data,
                                 e_Yt = Yt_bar + Yt_tilde_uncond,
                                 sd_Yt = implvol*sigma_bar*sigma_m*sqrt(ar_variance))
  # Add transform parameters
  seasonal_data$alpha <- model$transform$alpha
  seasonal_data$beta <- model$transform$beta
  if (target.Yt) {
    # Bounds for integration
    seasonal_data <- dplyr::mutate(seasonal_data,
                                   lower = -Inf,
                                   upper = Inf,
                                   z_x = model$transform$Y(1-strike/Ct),
                                   f_x = list(function(x, Ct) model$transform$GHI_y(x, Ct)))
  } else {
    # Bounds for integration
    seasonal_data <- dplyr::mutate(seasonal_data,
                                   lower = Ct*(1-alpha-beta),
                                   upper = Ct*(1-alpha),
                                   z_x = strike,
                                   f_x = list(function(x, Ct) x))
  }

  # Add Esscher parameter
  if (is.function(theta)) {
    seasonal_data$theta <- purrr::map_dbl(seasonal_data$Month, ~theta(.x))
  } else {
    seasonal_data$theta <- theta
  }
  # Select only relevant variables
  seasonal_data <- dplyr::select(seasonal_data, Month, Day, strike, Ct, theta, z_x, e_Yt, sd_Yt, alpha, beta,
                                 mu1, mu2, sd1, sd2, p1, lower, upper, f_x)

  # Pricing function
  pricing_month_day <- function(seasonal_data, combinations = NA, nmonth = 1, nday = 1, put = TRUE){

    # Seasonal Data
    df_n <- dplyr::filter(seasonal_data, Month == nmonth & Day == nday)
    # Mixture combinations
    if (is.na(combinations) && length(combinations) == 1) {
      df_n <- dplyr::mutate(df_n,
                            e_Yt_up = e_Yt + sd_Yt*mu1,
                            e_Yt_dw = e_Yt + sd_Yt*mu2,
                            sd_Yt_up = sd_Yt*sd1,
                            sd_Yt_dw = sd_Yt*sd2)
      # Create combinations table
      comb <- dplyr::tibble(mean = c(df_n$e_Yt_up, df_n$e_Yt_dw), sd = c(df_n$sd_Yt_up, df_n$sd_Yt_dw), probs = c(df_n$p1, 1-df_n$p1))
    } else {
      comb <- dplyr::filter(combinations, Month == nmonth)
      # Update mixture parameters
      comb$mean <- df_n$e_Yt + df_n$sd_Yt*comb$mean
      comb$sd <- df_n$sd_Yt*comb$sd
    }

    if (target.Yt) {
      if (df_n$theta == 0) {
        # Mixture Pdf
        pdf_Yt <- function(x) dmixnorm(x, means = comb$mean, sd = comb$sd, p = comb$probs)
        # Distribution Yt
        cdf_Rt <- function(x) pmixnorm(x, means = comb$mean, sd = comb$sd, p = comb$probs)
      } else {
        # Esscher Mixture Pdf
        pdf_Yt <- desscherMixture(means = comb$mean, sd = comb$sd, p = comb$probs, theta = df_n$theta)
        # Esscher Distribution Yt
        cdf_Rt <- function(x) integrate(pdf_Yt, lower = df_n$lower, upper = x)$value
      }
      # First moment for solar radiation
      e_Rt <- function(lower, upper) integrate(function(x) df_n$f_x[[1]](x, df_n$Ct)*pdf_Yt(x), lower = lower, upper = upper)$value
    } else {
      # Mixture Pdf
      pdf_Yt <- function(x) dmixnorm(x, means = comb$mean, sd = comb$sd, p = comb$probs)
      # Density for solar radiation
      pdf_Rt <- function(x) dsolarGHI(x, df_n$Ct, df_n$alpha, df_n$beta, pdf_Yt)

      # Esscher transform
      if (df_n$theta != 0) {
        pdf_Rt <- desscher(pdf_Rt, theta = df_n$theta, lower = df_n$lower, upper = df_n$upper)
        # Distribution for solar radiation
        cdf_Rt <- function(x) integrate(pdf_Rt, lower = df_n$lower, upper = x)$value
      } else {
        # Distribution for solar radiation
        cdf_Rt <- function(x) psolarGHI(x, df_n$Ct, model$transform$alpha, model$transform$beta, pdf_Yt)
      }
      # First moment for solar radiation
      e_Rt <- function(lower, upper) integrate(function(x) x*pdf_Rt(x), lower = lower, upper = upper)$value
    }

    # Option pricing
    if (put) {
      # Expected value (Put)
      df_n$Rt_plus <- e_Rt(df_n$lower, df_n$z_x)
      # Probability of exercise
      df_n$exercise <- cdf_Rt(df_n$z_x)
      # Option expected value
      df_n$premium <- df_n$strike*df_n$exercise - df_n$Rt_plus
      # Option type
      df_n$side <- "put"
    } else {
      # Expected value (Call)
      df_n$Rt_plus <-  e_Rt(df_n$z_x, df_n$upper)
      # Probability of exercise
      df_n$exercise <- 1 - cdf_Rt(df_n$z_x)
      # Option expected value
      df_n$premium <- df_n$Rt_plus - df_n$strike*df_n$exercise
      # Option type
      df_n$side <- "call"
    }
    # Expected value (GHI)
    df_n$Rt <- e_Rt(df_n$lower, df_n$upper)
    # Select only relevant variables
    df_n <- dplyr::select(df_n, Month, Day, side, premium, exercise, Rt_plus, Rt, strike)
    return(df_n)
  }

  # Model premium for each day and month
  df_month_day <- purrr::map2_df(seasonal_data$Month, seasonal_data$Day,
                                 ~pricing_month_day(seasonal_data, combinations, nmonth = .x, nday = .y, put = put))

  day_date <- paste0(ifelse(leap_year, "2020-", "2019-"), df_month_day$Month, "-", df_month_day$Day)
  df_month_day$n <- number_of_day(day_date)
  # Model premium aggregated for each Month
  df_month <- df_month_day %>%
    dplyr::group_by(Month, side) %>%
    dplyr::reframe(ndays = dplyr::n(),
                   premium = sum(premium),
                   daily_premium = premium/ndays,
                   exercise = mean(exercise),
                   Rt = sum(Rt),
                   Rt_plus = sum(Rt_plus),
                   strike = sum(strike))

  # Model premium aggregated by Year
  df_year <- df_month %>%
    dplyr::group_by(side) %>%
    dplyr::reframe(ndays = sum(ndays),
                   premium = sum(premium),
                   exercise = mean(exercise),
                   Rt = sum(Rt),
                   Rt_plus = sum(Rt_plus),
                   strike = sum(strike))

  # Rename columns to match target name
  rename_columns <- function(data, old_names, new_names){
    col_names <- colnames(data)
    for(i in 1:length(new_names)) {
      idx_old_names <- which(col_names == old_names[i])
      col_names[idx_old_names] <- new_names[i]
    }
    colnames(data) <- col_names
    return(data)
  }
  df_month_day <- rename_columns(df_month_day, c("Rt", "Rt_plus"), c(target, target_plus))
  df_month <- rename_columns(df_month, c("Rt", "Rt_plus"), c(target, target_plus))
  df_year <- rename_columns(df_year, c("Rt", "Rt_plus"), c(target, target_plus))

  structure(
    list(
      payoff_month_day = df_month_day,
      payoff_month = df_month,
      payoff_year = df_year
    ),
    class = c("solarOptionPayoff", "list")
  )
}



solarOption_model <- function(model, nmonths = 1:12, theta = 0, combinations = NA, implvol = 1, put = TRUE, target.Yt = TRUE, control_options = control_solarOption()){

  # Options control
  K <- control_options$K
  leap_year = control_options$leap_year
  # Target and seasonal mean
  target <- model$target
  target_bar <- paste0(model$target, "_bar")
  target_plus <- paste0(model$target, "_plus")
  # AR(2) stationary variance
  par <- model$AR_model_Yt$coefficients
  ar_variance <- (1-par[2])/((1 - par[2])*(1 - par[1]^2 - par[2]^2) - 2*par[1]^2*par[2])

  # Extract seasonal data
  seasonal_data <- model$seasonal_data
  # Filter for the selected months
  seasonal_data <- dplyr::filter(seasonal_data, Month %in% nmonths)
  # Remove 29 of February from computations
  if (!leap_year) {
    seasonal_data <- dplyr::filter(seasonal_data, !(Month == 2 & Day == 29))
  }
  # Option Strike
  seasonal_data$strike <- seasonal_data[[target_bar]]*exp(K)
  # Add upper and lower bounds and unconditional moments
  seasonal_data <- dplyr::mutate(seasonal_data,
                                 e_Yt = Yt_bar + Yt_tilde_uncond,
                                 sd_Yt = implvol*sigma_bar*sigma_m*sqrt(ar_variance))
  # Add transform parameters
  seasonal_data$alpha <- model$transform$alpha
  seasonal_data$beta <- model$transform$beta
  if (target.Yt) {
    # Bounds for integration
    seasonal_data <- dplyr::mutate(seasonal_data,
                                   lower = -Inf,
                                   upper = Inf,
                                   z_x = model$transform$Y(1-strike/Ct),
                                   f_x = list(function(x, Ct) model$transform$GHI_y(x, Ct)))
  } else {
    # Bounds for integration
    seasonal_data <- dplyr::mutate(seasonal_data,
                                   lower = Ct*(1-alpha-beta),
                                   upper = Ct*(1-alpha),
                                   z_x = strike,
                                   f_x = list(function(x, Ct) x))
  }

  # Add Esscher parameter
  if (is.function(theta)) {
    seasonal_data$theta <- purrr::map_dbl(seasonal_data$Month, ~theta(.x))
  } else {
    seasonal_data$theta <- theta
  }
  # Select only relevant variables
  seasonal_data <- dplyr::select(seasonal_data, Month, Day, strike, Ct, theta, z_x, e_Yt, sd_Yt, alpha, beta,
                                 mu1, mu2, sd1, sd2, p1, lower, upper, f_x)

  # Pricing function
  pricing_month_day <- function(seasonal_data, combinations = NA, nmonth = 1, nday = 1, put = TRUE){

    # Seasonal Data
    df_n <- dplyr::filter(seasonal_data, Month == nmonth & Day == nday)
    # Mixture combinations
    if (is.na(combinations) && length(combinations) == 1) {
      df_n <- dplyr::mutate(df_n,
                            e_Yt_up = e_Yt + sd_Yt*mu1,
                            e_Yt_dw = e_Yt + sd_Yt*mu2,
                            sd_Yt_up = sd_Yt*sd1,
                            sd_Yt_dw = sd_Yt*sd2)
      # Create combinations table
      comb <- dplyr::tibble(mean = c(df_n$e_Yt_up, df_n$e_Yt_dw), sd = c(df_n$sd_Yt_up, df_n$sd_Yt_dw), probs = c(df_n$p1, 1-df_n$p1))
    } else {
      comb <- dplyr::filter(combinations, Month == nmonth)
      # Update mixture parameters
      comb$mean <- df_n$e_Yt + df_n$sd_Yt*comb$mean
      comb$sd <- df_n$sd_Yt*comb$sd
    }

    if (target.Yt) {
      if (df_n$theta == 0) {
        # Mixture Pdf
        pdf_Yt <- function(x) dmixnorm(x, means = comb$mean, sd = comb$sd, p = comb$probs)
        # Distribution Yt
        cdf_Rt <- function(x) pmixnorm(x, means = comb$mean, sd = comb$sd, p = comb$probs)
      } else {
        # Esscher Mixture Pdf
        pdf_Yt <- desscherMixture(means = comb$mean, sd = comb$sd, p = comb$probs, theta = df_n$theta)
        # Esscher Distribution Yt
        cdf_Rt <- function(x) integrate(pdf_Yt, lower = df_n$lower, upper = x)$value
      }
      # First moment for solar radiation
      e_Rt <- function(lower, upper) integrate(function(x) df_n$f_x[[1]](x, df_n$Ct)*pdf_Yt(x), lower = lower, upper = upper)$value
    } else {
      # Mixture Pdf
      pdf_Yt <- function(x) dmixnorm(x, means = comb$mean, sd = comb$sd, p = comb$probs)
      # Density for solar radiation
      pdf_Rt <- function(x) dsolarGHI(x, df_n$Ct, df_n$alpha, df_n$beta, pdf_Yt)

      # Esscher transform
      if (df_n$theta != 0) {
        pdf_Rt <- desscher(pdf_Rt, theta = df_n$theta, lower = df_n$lower, upper = df_n$upper)
        # Distribution for solar radiation
        cdf_Rt <- function(x) integrate(pdf_Rt, lower = df_n$lower, upper = x)$value
      } else {
        # Distribution for solar radiation
        cdf_Rt <- function(x) psolarGHI(x, df_n$Ct, model$transform$alpha, model$transform$beta, pdf_Yt)
      }
      # First moment for solar radiation
      e_Rt <- function(lower, upper) integrate(function(x) x*pdf_Rt(x), lower = lower, upper = upper)$value
    }

    # Option pricing
    if (put) {
      # Expected value (Put)
      df_n$Rt_plus <- e_Rt(df_n$lower, df_n$z_x)
      # Probability of exercise
      df_n$exercise <- cdf_Rt(df_n$z_x)
      # Option expected value
      df_n$premium <- df_n$strike*df_n$exercise - df_n$Rt_plus
      # Option type
      df_n$side <- "put"
    } else {
      # Expected value (Call)
      df_n$Rt_plus <-  e_Rt(df_n$z_x, df_n$upper)
      # Probability of exercise
      df_n$exercise <- 1 - cdf_Rt(df_n$z_x)
      # Option expected value
      df_n$premium <- df_n$Rt_plus - df_n$strike*df_n$exercise
      # Option type
      df_n$side <- "call"
    }
    # Expected value (GHI)
    df_n$Rt <- e_Rt(df_n$lower, df_n$upper)
    # Select only relevant variables
    df_n <- dplyr::select(df_n, Month, Day, side, premium, exercise, Rt_plus, Rt, strike)
    return(df_n)
  }

  # Model premium for each day and month
  df_month_day <- purrr::map2_df(seasonal_data$Month, seasonal_data$Day,
                                 ~pricing_month_day(seasonal_data, combinations, nmonth = .x, nday = .y, put = put))

  day_date <- paste0(ifelse(leap_year, "2020-", "2019-"), df_month_day$Month, "-", df_month_day$Day)
  df_month_day$n <- number_of_day(day_date)
  # Model premium aggregated for each Month
  df_month <- df_month_day %>%
    dplyr::group_by(Month, side) %>%
    dplyr::reframe(ndays = dplyr::n(),
                   premium = sum(premium),
                   daily_premium = premium/ndays,
                   exercise = mean(exercise),
                   Rt = sum(Rt),
                   Rt_plus = sum(Rt_plus),
                   strike = sum(strike))

  # Model premium aggregated by Year
  df_year <- df_month %>%
    dplyr::group_by(side) %>%
    dplyr::reframe(ndays = sum(ndays),
                   premium = sum(premium),
                   exercise = mean(exercise),
                   Rt = sum(Rt),
                   Rt_plus = sum(Rt_plus),
                   strike = sum(strike))

  # Rename columns to match target name
  rename_columns <- function(data, old_names, new_names){
    col_names <- colnames(data)
    for(i in 1:length(new_names)) {
      idx_old_names <- which(col_names == old_names[i])
      col_names[idx_old_names] <- new_names[i]
    }
    colnames(data) <- col_names
    return(data)
  }
  df_month_day <- rename_columns(df_month_day, c("Rt", "Rt_plus"), c(target, target_plus))
  df_month <- rename_columns(df_month, c("Rt", "Rt_plus"), c(target, target_plus))
  df_year <- rename_columns(df_year, c("Rt", "Rt_plus"), c(target, target_plus))

  structure(
    list(
      payoff_month_day = df_month_day,
      payoff_month = df_month,
      payoff_year = df_year
    ),
    class = c("solarOptionPayoff", "list")
  )
}



#' Calibrator for solar Options
#' @inheritParams solarOption_historical
#' @param abstol The absolute convergence tolerance. Only useful for non-negative functions, as a tolerance for reaching zero.
#' @param reltol Relative convergence tolerance. The algorithm stops if it is unable to reduce the value by a factor of reltol * (abs(val) + reltol) at a step.
#' Defaults to `sqrt(.Machine$double.eps)`, typically about 1e-8.
#'
#' @examples
#' model <- Bologna
#' model_cal <- solarOption_calibrator(model, nmonths = 7, reltol=1e-3)
#' # Compare log-likelihoods
#' model$loglik
#' model_cal$loglik
#' # Compare parameters
#' model$NM_model$parameters[7,]
#' model_cal$NM_model$parameters[7,]
#' @export
solarOption_calibrator <- function(model, nmonths = 1:12, abstol = 1e-4, reltol = 1e-4, control_options = control_solarOption()){

  model_cal <- model$clone(deep = TRUE)
  # Historical monthly payoff
  payoff_call <- solarOption_historical(model_cal, put = FALSE, control_options = control_options)$payoff_month$premium
  payoff_put <- solarOption_historical(model_cal, put = TRUE, control_options = control_options)$payoff_month$premium
  # Put-call parity
  optionParity <- function(e_St, e_St_plus, strike, exercise, put = TRUE){
    put <- ifelse(put, 1, -1)
    exercise <- 1 - exercise
    price <- (e_St - e_St_plus) - exercise*strike
    return(put*price)
  }

  #min_sd <- min(c(model$NM_model$parameters$sd1, model$NM_model$parameters$sd2))*0.8
  #min_p <- min(c(model$NM_model$parameters$p1, 1-model$NM_model$parameters$p1, 0.2))
  #max_p <- max(c(model$NM_model$parameters$p1, 1-model$NM_model$parameters$p1, 0.8))
  # Loss function
  loss_function <- function(params, model_cal, nmonth){
   # if(params[3] < min_sd | params[4] < min_sd | params[5] > max_p | params[5] < min_p){
   #    return(NA)
   # }
    # Update the parameters
    model_cal$NM_model$model[[nmonth]]$update(means = params[1:2], sd = params[3:4], p = c(params[5], 1-params[5]))
    # Put pricing function
    df_put <- solarOption_model(model_cal, put = TRUE, nmonth = nmonth, control_options = control_options)
    # Put price
    price_put <- df_put$payoff_year$premium
    # Call price
    df <- df_put$payoff_month_day
    price_call <- sum(optionParity(df$GHI, df$GHI_plus, df$strike, df$exercise, put = TRUE))
    # Loss function
    loss <- abs(price_call - payoff_call[nmonth]) + abs(price_put - payoff_put[nmonth])
    message("Loss: ", loss, " Params: ", format(params, digits = 3), "\r", appendLF = FALSE)
    return(loss^2)
  }

  # Monthly calibration
  for(nmonth in nmonths){
    message("------------------------------------ Month: ", nmonth, " ------------------------------------ ")
    params <- unlist(model_cal$NM_model$parameters[nmonth, c(2:6)])
    opt <- optim(params, loss_function, model_cal = model_cal, nmonth = nmonth, control = list(abstol = abstol, reltol = reltol))
    params <- opt$par
    # Update the parameters
    model_cal$NM_model$model[[nmonth]]$update(means = params[1:2], sd = params[3:4], p = c(params[5], 1-params[5]))
  }
  # Update log-likelihood and conditional moments
  model_cal$conditional_moments()
  model_cal$unconditional_moments()
  model_cal$logLik()
  return(model_cal)
}


#' Structure payoffs
#'
#' @param payoffs object with the class `solarOptionPayoffs`. See the function \code{\link{solarOptionPayoffs}} for details.
#' @param type method used for computing the premium. If `model`, the default will be used the analytic model,
#' otherwise with `scenarios` the monte carlo scenarios stored inside the `model$scenarios$P`.
#' @param exact_daily_premium when `TRUE` the historical premium is computed as daily average among all the years.
#' Otherwise the monthly premium is computed and then divided by the number of days of the month.
#' @return The object `payoffs` with class `solarOptionPayoffs`.
#'
#' @rdname solarOption_structure
#' @name solarOption_structure
#' @export
solarOption_structure <- function(payoffs, type = "model", put = TRUE, exact_daily_premium = TRUE){

  option_type <- ifelse(put, "put", "call")
  type <- match.arg(type, choices = c("scenarios", "model"))
  payoff <- payoffs[[option_type]][c("historical", type)]

  # Yearly Premium
  df_year <- dplyr::tibble(
    side =  payoff$historical$payoff_year$side,
    ndays =  payoff$historical$payoff_year$ndays,
    # Premiums for the option
    premium = payoff$historical$payoff_year$premium,
    premium_P = payoff[[type]]$P$payoff_year$premium,
    premium_Q = payoff[[type]]$Q$payoff_year$premium,
    premium_Qdw = payoff[[type]]$Qdw$payoff_year$premium,
    premium_Qup = payoff[[type]]$Qup$payoff_year$premium,
    premium_Qr = payoff[[type]]$Qr$payoff_year$premium,
    # Probabilities of exercise the option
    exercise = payoff$historical$payoff_year$exercise,
    exercise_P = payoff[[type]]$P$payoff_year$exercise,
    exercise_Q = payoff[[type]]$Q$payoff_year$exercise,
    exercise_Qup = payoff[[type]]$Qup$payoff_year$exercise,
    exercise_Qdw = payoff[[type]]$Qdw$payoff_year$exercise,
    exercise_Qr = payoff[[type]]$Qr$payoff_year$exercise,
  )

  # Monthly Premium
  df_month <- dplyr::tibble(
    Month = payoff$historical$payoff_month$Month,
    side = payoff$historical$payoff_month$side,
    n = payoff$historical$payoff_month$ndays,
    # Premiums for the option
    premium = payoff$historical$payoff_month$premium,
    premium_P = payoff[[type]]$P$payoff_month$premium,
    premium_Q = payoff[[type]]$Q$payoff_month$premium,
    premium_Qup = payoff[[type]]$Qup$payoff_month$premium,
    premium_Qdw = payoff[[type]]$Qdw$payoff_month$premium,
    premium_Qr = payoff[[type]]$Qdw$payoff_month$premium,
    # Probabilities of exercise the option
    exercise = payoff$historical$payoff_month$exercise,
    exercise_P = payoff[[type]]$P$payoff_month$exercise,
    exercise_Q = payoff[[type]]$Q$payoff_month$exercise,
    exercise_Qup = payoff[[type]]$Qup$payoff_month$exercise,
    exercise_Qdw = payoff[[type]]$Qdw$payoff_month$exercise,
    exercise_Qr = payoff[[type]]$Qr$payoff_month$exercise,
  )

  # Monthly Daily Premium
  df_month_day_mean <- dplyr::tibble(
    Month = payoff$historical$payoff_month$Month,
    side = payoff$historical$payoff_month$side,
    n = payoff$historical$payoff_month$ndays,
    # Premiums for the option
    premium = payoff$historical$payoff_month$daily_premium,
    premium_P = payoff[[type]]$P$payoff_month$daily_premium,
    premium_Q = payoff[[type]]$Q$payoff_month$daily_premium,
    premium_Qup = payoff[[type]]$Qup$payoff_month$daily_premium,
    premium_Qdw = payoff[[type]]$Qdw$payoff_month$daily_premium,
    premium_Qr = payoff[[type]]$Qr$payoff_month$daily_premium,
  )

  # Monthly-Day Premium
  df_month_day <- dplyr::tibble(
    Month = payoff$historical$payoff_month_day$Month,
    Day = payoff$historical$payoff_month_day$Day,
    side = payoff$historical$payoff_month_day$side,
    # Premiums for the option
    premium = payoff$historical$payoff_month_day$premium,
    premium_P = payoff[[type]]$P$payoff_month_day$premium,
    premium_Q = payoff[[type]]$Q$payoff_month_day$premium,
    premium_Qup = payoff[[type]]$Qup$payoff_month_day$premium,
    premium_Qdw = payoff[[type]]$Qdw$payoff_month_day$premium,
    # Probabilities of exercise the option
    exercise = payoff$historical$payoff_month_day$exercise,
    exercise_P = payoff[[type]]$P$payoff_month_day$exercise,
    exercise_Q = payoff[[type]]$Q$payoff_month_day$exercise,
    exercise_Qup = payoff[[type]]$Qup$payoff_month_day$exercise,
    exercise_Qdw = payoff[[type]]$Qdw$payoff_month_day$exercise,
    n = payoff$hist$payoff_month_day$n
  )

  # Historical Daily Premiums
  df_payoff <- dplyr::select(payoff$historical$payoff, -exercise, -GHI, -strike)
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

  payoffs[[option_type]][[type]]$structured <- list(payoff = df_payoff,
                                                    payoff_year = df_year,
                                                    payoff_month = df_month,
                                                    payoff_month_day = df_month_day,
                                                    payoff_cum = df_cum)

  return(payoffs)
}


