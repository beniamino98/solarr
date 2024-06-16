cams_solar_radiation_ts <- NULL
.onLoad <- function(libname, pkgname) {
  module <- reticulate::import_from_path(module = "cams_solar_radiation_ts", path = system.file("py", package = packageName()))
  cams_solar_radiation_ts <<- module$cams_solar_radiation_ts
}

#' Perform normality tests
#'
#' @param x vector
#' @param p_value p.value
#'
#' @return a tibble
#'
#' @rdname normality_test
#' @name normality_test
#' @export

normality_test <- function(x = NULL, p_value = 0.05){
  suppressWarnings(
    dplyr::bind_rows(
      broom::tidy(nortest::ad.test(x)),
      broom::tidy(nortest::cvm.test(x)),
      broom::tidy(nortest::lillie.test(x)),
      broom::tidy(nortest::pearson.test(x)),
      broom::tidy(nortest::sf.test(x)),
      broom::tidy(stats::shapiro.test(x))
    ) %>%
      dplyr::mutate(
        H0 = ifelse(p.value <= p_value, "Rejected", "Non Rejected")
      ) %>%
      dplyr::select(method, statistic, p.value, H0)
  )
}

#' Discount factor function
#'
#' @param r level of yearly constant risk-free rate
#' @export
discount_factor <- function(r = 0.03){
  function(tau){
    (1 + r)^(-tau)
  }
}


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
#'
#' @examples
#' object <- Location("Roma")
#' object
#' @export
print.Location <- function(x){
  msg_1 <- paste0("Place: ", x$place, " \n (Lat: ", x$coords$lat, "; Lon: ", x$coords$lon, ") \n")
  msg_2 <- paste0("From: ", min(x$data$date), " - ", max(x$data$date), " (Nobs: ", nrow(x$data), ")")
  cat(paste0(msg_1, msg_2))
}
