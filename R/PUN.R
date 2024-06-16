# Function for computing the monthly mean unique national price for Italy
#' PUN
#' @name PUN
#' @rdname PUN
#' @description Function that computes the mean PUN.
#' @param month Reference month.
#' @param year Reference year.
#' @param file path
#' @return numeric, price in euros of a kwh
#' @export

PUN <- function(nyear = NULL, nmonth = NULL, file = "data/df_GME_day.RData"){

  if (!exists("df_GME_day", envir = .GlobalEnv)){
    load(file, envir = .GlobalEnv)
  }
  if (!is.null(nmonth) & is.null(nyear)) {
    df_gme <- dplyr::filter(df_GME_day, Month == nmonth)
  } else if(is.null(nmonth) & !is.null(nyear)) {
    df_gme <- dplyr::filter(df_GME_day, Year == nyear)
  } else if(!is.null(nmonth) & !is.null(nyear)) {
    df_gme <- dplyr::filter(df_GME_day, Year == nyear, Month == nmonth)
  } else {
    df_gme <- df_GME_day
  }
  return(mean(df_gme$PUN/1000))
}
