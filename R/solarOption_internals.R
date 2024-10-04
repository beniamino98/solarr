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
#' @keywords 2beRevised
#' @export
solarOption_contracts  <- function(payoff, type = "mod", premium = "Q", nyear = 2021, tick = 0.06, efficiency = 0.2, n_panels = 2000, pun = 0.06){

  # Extract historical payoff
  payoff_hist <- payoff$hist$payoff
  # Match the type of computation
  type <- match.arg(type, choices = c("mod", "sim"))
  # Extract the payoff
  payoffs <- payoff[[type]]
  # Match the type of premium
  premium <- match.arg(premium, names(payoffs))
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
    loss <- sd(dplyr::filter(df_, Year < nyear)$hedged, na.rm = TRUE)
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


#' Test errors solar Option model
#'
#' @examples
#' model <- Bologna
#' solarOption_model_test(model)
#' solarOption_model_test(model, nmonths = 6)
#' @export
solarOption_model_test <- function(model, nmonths = 1:12, control_options = control_solarOption()){

  payoff_model <- solarOption_model(model, nmonths = nmonths, control_options = control_options)$payoff_month
  payoff_hist <- solarOption_historical(model, nmonths = nmonths, control_options = control_options)$payoff_month

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


#' Implied expected return at maturity
#'
#' @keywords 2beRevised
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

