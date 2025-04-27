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

#' solarOptionPayoff
#'
#' @param data df_payoff
#' @param leap_year control params
#' @keywords internal
#' @export
#' @rdname solarOptionPayoff
#' @name solarOptionPayoff
#' @return An object of the class `solarOptionPayoff`.
solarOptionPayoff <- function(data, leap_year = FALSE){
  # Reorder variables
  df_payoff <- dplyr::select(data, side, date, Year, Month, Day, Rt, strike, premium, payoff, exercise)
  # Include or not 29-th of February from computation
  if (!leap_year) {
    df_payoff <- dplyr::filter(df_payoff, !(Month == 2 & Day == 29))
  }
  # Aggregation by Year and Month
  df_year_month <- df_payoff %>%
    dplyr::group_by(Year, Month, side) %>%
    dplyr::reframe(ndays = dplyr::n(),
                   payoff = sum(payoff),
                   premium = sum(premium),
                   exercise = mean(exercise),
                   Rt = sum(Rt),
                   strike = sum(strike))

  # Aggregation by Month and Day
  df_month_day <- df_payoff %>%
    dplyr::group_by(Month, Day, side) %>%
    dplyr::reframe(premium = mean(premium),
                   payoff = mean(payoff),
                   exercise = mean(exercise),
                   Rt = mean(Rt),
                   strike = mean(strike))

  # Aggregation by Month
  df_month <- df_year_month %>%
    dplyr::group_by(Month, side)  %>%
    dplyr::reframe(ndays = dplyr::n(),
                   payoff = mean(payoff),
                   premium = mean(premium),
                   exercise = mean(exercise),
                   Rt = mean(Rt),
                   strike = mean(strike))

  # Aggregation by Year
  df_year <- df_month %>%
    dplyr::group_by(side)  %>%
    dplyr::reframe(ndays = sum(ndays),
                   payoff = sum(payoff),
                   premium = sum(premium),
                   exercise = mean(exercise),
                   Rt = sum(Rt),
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
    class = c("solarOptionPayoff", "list")
  )
}



#' Control for hedging
#'
#' @param n_panels numeric, number of meters squared of solar panels.
#' @param efficiency numeric, mean efficiency of the solar panels.
#' @param PUN numeric, mean efficiency of the solar panels.
#' @param tick numeric, conversion tick for the monetary payoff of a contract.
#' @param n_contract numeric, number of contracts
#'
#' @rdname control_hedging
#' @name control_hedging
#' @export
control_hedging <- function(n_panels = 1, efficiency = 1, PUN = 1, tick = 1, n_contracts = 1){
  structure(
    list(
      n_panels = n_panels,
      efficiency = efficiency,
      PUN = PUN,
      tick = tick,
      n_contracts = n_contracts
    ),
    class = c("control", "list")
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
#' @param control_hedging numeric, list of hedging parameters.
#'
#' @rdname solarOption_contracts
#' @name solarOption_contracts
#' @keywords 2beRevised
#' @export
solarOption_contracts  <- function(payoff, type = "model", premium = "Q", put = TRUE, nyear = 2021, control = control_hedging()){

  # Control parameters for hedging
  tick = control$tick
  efficiency = control$efficiency
  n_panels = control$n_panels
  PUN = control$PUN

  # Match option type
  option_type <- ifelse(put, "put", "call")
  # Match the type of computation
  type <- match.arg(type, choices = c("model", "scenarios"))
  # Extract historical payoff
  payoff_hist <- payoff[[option_type]]$historical$payoff
  # Extract the payoff
  payoff <- payoff[[option_type]][[type]]
  # Match the type of premium
  premium <- match.arg(premium, choices = names(payoff))
  # Extract daily premium
  payoff <- payoff[[premium]]$payoff_month_day
  # Select only relevant column
  payoff <- dplyr::select(payoff, Month, Day, premium)
  # Merge realized GHI and premium
  df_hedged <- dplyr::left_join(payoff_hist, payoff, by = c("Month", "Day"))
  # Loss function depending on the number of contracts
  loss_function <- function(n_contracts, df_hedged){
    # Compute hedged cash-flows
    df_ <- dplyr::mutate(df_hedged, hedged = PUN * n_panels * efficiency * GHI + tick * n_contracts * (payoff - premium))
    # Exclude nyear and the following from loss computation
    loss <- sd(dplyr::filter(df_, Year < nyear)$hedged, na.rm = TRUE)
    return(loss)
  }
  # Optimize the number of contracts
  opt <- optim(par = 10, fn = loss_function, method = "Brent",
               lower = 1, upper = n_panels * efficiency*10, df_hedged = df_hedged)

  list(
    nyear = nyear + 1,
    n_contracts = trunc(opt$par),
    type = type,
    premium = premium,
    sd = opt$value,
    control = control
  )
}


#' solarOptionPayoff
#'
#' @param model solarModel
#' @param control_options control list, see \code{\link{control_solarOption}} for more details.
#' @keywords old
#' @return An object of the class `solarOptionPayoffs`.
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

#' Structure payoffs
#'
#' @param payoffs object with the class `solarOptionPayoffs`. See the function \code{\link{solarOptionPayoffs}} for details.
#' @param type method used for computing the premium. If `model`, the default will be used the analytic model,
#' otherwise with `scenarios` the monte carlo scenarios stored inside the `model$scenarios$P`.
#' @param exact_daily_premium when `TRUE` the historical premium is computed as daily average among all the years.
#' Otherwise the monthly premium is computed and then divided by the number of days of the month.
#' @keywords old
#' @return The object `payoffs` with class `solarOptionPayoffs`.
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
