#' Discount factor function
#'
#' @param r level of yearly constant risk-free rate
#' @export
discount_factor <- function(r = 0.03){
  function(tau){
    (1 + r)^(-tau)
  }
}

#' Aggregate payoff
#'
#' @keywords internal
aggregate_payoff <- function(df_month_day, r = 0.03){

  #' @examples
  #' nmonth = 1:12
  #' t_now = "2012-01-01"
  #' t_hor = "2013-01-01"
  function(t_now, t_hor){
    t_now = as.Date(t_now)
    t_hor = as.Date(t_hor)
    seq_dates <- seq.Date(t_now, t_hor, by = "1 day")
    df_payoff_date <- dplyr::tibble(
      date = seq_dates,
      Month = lubridate::month(date),
      Day = lubridate::day(date),
    ) %>%
      dplyr::left_join(df_month_day)
    df_payoff_date %>%
      dplyr::reframe(
        premium = sum(premium),
        exercise = mean(exercise),
        ndays = n()
      ) %>%
      dplyr::mutate(
        disc = (1 + r)^(ndays/365),
        disc_premium = premium*disc,
        t_now = t_now, t_hor = t_hor) %>%
      dplyr::select(t_now, t_hor, dplyr::everything())
  }
}


#' Method Print for Location object
#' @keywords internal
#' @noRd
print.Location <- function(x){
  msg_1 <- paste0("Place: ", x$place, " \n (Lat: ", x$coords$lat, "; Lon: ", x$coords$lon, ") \n")
  msg_2 <- paste0("From: ", min(x$data$date), " - ", max(x$data$date), " (Nobs: ", nrow(x$data), ")")
  cat(paste0(msg_1, msg_2))
}


#'  Compute optimal number of contracts
#'
#' @export
optimal_n_contracts <- function(model, type = "model", premium = "Qr", nyear = 2021,
                                tick = 0.06, efficiency = 0.2, n_panels = 2000, pun = 0.06){

  # Extract historical payoff
  df_hist <- model$payoffs$hist$payoff
  # Extract daily premium
  df_premium <- model$payoffs[[type]][[premium]]$payoff_month_day
  # Select only relevant column
  df_premium <- dplyr::select(df_premium, Month, Day, premium)
  # Dataset with GHI and premium
  df_hedged <- dplyr::left_join(df_hist, df_premium, by = c("Month", "Day"))
  # Loss function depending on the number of contracts
  loss_function <- function(n_contracts, df_hedged){
    df_ <- dplyr::mutate(df_hedged, hedged = pun*n_panels*efficiency*GHI + tick*n_contracts*(payoff - premium))
    df_ <- dplyr::filter(df_, Year < nyear)
    sd(df_$hedged)
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
