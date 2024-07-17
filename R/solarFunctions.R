#' Conversion in Radiant or Degrees
#'
#' Convert an angle in radiant into an angle in degrees.
#'
#' @param x numeric vector, angles in radiant or degrees.
#' @examples
#' # convert 0.34 radiant in degrees
#' from_radiant_to_degree(0.34)
#' @aliases from_radiant_to_degree
#' @aliases from_degree_to_radiant
#' @rdname radiant
#' @return numeric vector.
#' @export
from_radiant_to_degree <- function(x){
  assertive::assert_is_numeric(x)
  x <- as.numeric(x)
  y <- x*180/base::pi
  return(y)
}

#' @examples
#' # convert 19.48 degree in radiant
#' from_degree_to_radiant(19.48)
#'
#' @rdname radiant
#' @export
from_degree_to_radiant <- function(x){
  assertive::assert_is_numeric(x)
  x <- as.numeric(x)
  y <- x*base::pi/180
  return(y)
}

#' Detect the season
#'
#' Detect the season from a vector of dates.
#'
#' @param day_date vector of dates in the format `%YYYY-%MM-%DD`.
#' @examples
#' detect_season("2040-01-31")
#' detect_season(c("2040-01-31", "2023-04-01", "2015-09-02"))
#' @name detect_season
#' @rdname detect_season
#' @return a character vector containing the correspondent season. Can be `spring`, `summer`, `autumn`, `winter`.
#' @export
detect_season <- function(day_date = NULL){
  # first day of the year
  date_start <- as.Date("2000-01-01")
  # start dates for the seasons
  date_spring <- as.Date("2000-03-23")
  date_summer <- as.Date("2000-06-23")
  date_autumn <- as.Date("2000-09-23")
  date_winter <- as.Date("2000-12-23")
  # last day of the year
  date_end <- as.Date("2000-12-31")

  seasons <- c()
  for(i in 1:length(day_date)){
    day <- as.Date(day_date[i]) # day date
    day_ <- lubridate::day(day) # number of day of the month
    month_ <- lubridate::month(day) # number of month of the year
    day <- as.Date(paste0("2000-", month_, "-", day_)) # standardized year
    # Season for the selected day
    seasons[i] <- dplyr::case_when(day < date_spring & day >= date_start ~ "Winter",
                                   day >= date_spring & day < date_summer ~ "Spring",
                                   day >= date_summer & day < date_autumn ~ "Summer",
                                   day >= date_autumn & day < date_winter ~ "Autumn",
                                   day >= date_winter & day <= date_end ~ "Winter",
                                   TRUE ~ NA_character_)
  }
  return(seasons)
}



#' Is leap year?
#'
#' Check if an year is leap (366 days) or not (365 days).
#'
#' @param day_date dates vector in the format `%YYYY-%MM-%DD`.
#' @examples
#' is_leap_year("2024-02-01")
#' is_leap_year(c("2024-10-01", "2025-10-01"))
#' is_leap_year("2029-02-01")
#' @name is_leap_year
#' @rdname is_leap_year
#' @return Boolean. `TRUE` if it is a leap year, `FALSE` otherwise.
#' @export
is_leap_year <- function(day_date){
  day_date <- lubridate::as_date(day_date)
  out <- c()
  for(i in 1:length(day_date)){
    date_year <- lubridate::year(day_date[i])
    if (date_year %% 4 == 0) {
      out[i] <- TRUE
    } else {
      out[i] <- FALSE
    }
  }
  return(out)
}



#' Number of day
#'
#' Compute the number of day of the year given a vector of dates.
#'
#' @param day_date dates vector in the format `%YYYY-%MM-%DD`.
#'
#' @examples
#' number_of_day("2040-01-31")
#' number_of_day(c("2040-01-31", "2023-04-01", "2015-09-02"))
#' @name number_of_day
#' @rdname number_of_day
#' @return Numeric vector with the number of the day during the year.
#' @export
number_of_day <- function(day_date = NULL){

  if (is.null(day_date)) {
    warning('The argument "day_date" is NULL while it require a character of the form %YYYY-%MM-%DD')
    return(NULL)
  }
  n_of_day <- c()
  for(i in 1:length(day_date)){
    day <- lubridate::as_date(day_date[i]) # selected date
    # Check if is a valid date
    if (is.na(day)) {
      n_of_day[i] <- NA_integer_
      next
    }
    # create a sequence of dates for selected year
    selected_year <- lubridate::year(day) # selected year
    start_date <- as.Date(paste0(selected_year, "-01-01")) # start dates
    end_date <- as.Date(paste0(selected_year, "-12-31")) # end dates
    year_seq_date <- seq.Date(start_date, end_date, 1)
    # extract the number of the day
    n_of_day[i] <- which(year_seq_date == day)
  }
  names(n_of_day) <- day_date
  return(n_of_day)
}



#' Solar time adjustment
#'
#' Compute the time adjustment for a date.
#'
#' @param day_date vector of dates in the format `%YYYY-%MM-%DD`
#' @param day_end end date, if it is not NULL will be end date.
#' @examples
#' solar_time_adjustment("2040-01-31")
#' solar_time_adjustment(c("2040-01-31", "2023-04-01", "2015-09-02"))
#' @return a numeric vector containing the time adjustment in seconds.
#' @name solar_time_adjustment
#' @rdname solar_time_adjustment
#' @export
solar_time_adjustment <- function(day_date = NULL, day_end = NULL){

  if (!is.null(day_end)) {
    day_date <- seq.Date(as.Date(day_date), as.Date(day_end), 1)
  } else {
    day_date <- as.Date(day_date)
  }
  time_adj <- c()
  n <- number_of_day(day_date)  # number of days for all dates
  for(i in 1:length(day_date)){
    # B is an angle so it must be expressed in radiant.
    # Alternatively we can compute it in degree and convert it:
    # B = from_degree_to_radiant(360*(n-1)/365)
    B <- 2*base::pi*(n[i] - 1)/365
    # Time adjustment in minutes
    time_adj[i] <- 229.2*(0.000075 + 0.001868*cos(B) - 0.032077*sin(B) - 0.014615*cos(2*B) - 0.04089*sin(2*B))
    # Time adjustment in seconds
    time_adj[i] <- lubridate::seconds(time_adj[i]*60)
  }
  out <- dplyr::tibble(date = day_date, time_adjustment = time_adj)
  return(out)
}



#' Solar time constant
#'
#' Compute the solar constant for a date.
#'
#' @param day_date vector of dates in the format `YYYY-MM-DD`.
#' @param day_end end date, if it is not `NULL` will be end date.
#' @param method method used for computation, can be `cooper` or `spencer`.
#' @examples
#' solar_time_constant("2040-01-31")
#' solar_time_constant(c("2040-01-31", "2023-04-01", "2015-09-02"))
#' @name solar_time_constant
#' @rdname solar_time_constant
#' @return a numeric vector containing the solar constant.
#' @export
solar_time_constant <- function(day_date = NULL, day_end = NULL, method = "spencer"){

  # Mean solar constant
  G0 <- 1367
  # Method of computation
  method <- match.arg(method, choices = c("cooper", "spencer"))
  # Sequence of dates
  day_date <- as.Date(day_date)
  if (!is.null(day_end)) {
    day_end <- as.Date(day_end[1])
    day_date <- seq.Date(day_date[1], day_end, by = "1 day")
  }
  # Number of the day for all dates
  n <- number_of_day(day_date)
  G0n <- rep(G0, length(day_date))
  for(i in 1:length(day_date)){
    # Method cooper
    if (method == "cooper") {
      B <- (2*base::pi)*(n[i])/365
      G0n[i] <- G0n[i]*(1 + 0.033 * cos(B))
    # Method spencer
    } else if (method == "spencer") {
      B <- (n[i] - 1)*(2*base::pi/365)
      G0n[i] <- G0n[i]*(1.000110 + 0.034221*cos(B) + 0.001280*sin(B) + 0.000719*cos(2*B) + 0.000077*sin(2*B))
    }
  }
  output <- dplyr::tibble(date = day_date, method = method, G0 = G0, G0n = G0n)
  return(output)
}



#' Solar time declination
#'
#' Compute the solar declination for different dates.
#'
#' @param day_date vector of dates in the format `%YYYY-%MM-%DD`
#' @param day_end end date, if it is not NULL will be end date.
#' @param method method used for computation, can be `cooper` or `spencer`.
#' @examples
#' solar_time_declination("2040-01-01", day_end = "2040-12-31")
#' solar_time_declination(c("2040-01-31", "2023-04-01", "2015-09-02"))
#' @name solar_time_declination
#' @rdname solar_time_declination
#' @return a numeric vector containing the solar declination in minutes.
#' @export
solar_time_declination <- function(day_date = NULL, day_end = NULL, method = c("cooper", "spencer")){

  # Mean solar constant
  G0 <- 1367
  # Method of computation
  method <- match.arg(method, choices = c("cooper", "spencer"))
  # Sequence of dates
  day_date <- as.Date(day_date)
  if (!is.null(day_end)) {
    day_end <- as.Date(day_end[1])
    day_date <- seq.Date(day_date[1], day_end, by = "1 day")
  }
  # Number of the day for all dates
  n <- number_of_day(day_date)

  declination <- rep(0, length(day_date))
  for(i in 1:length(day_date)){
    # Method cooper
    if (method == "cooper") {
      declination[i] <- sin(2*base::pi*(284 + n[i])/365)*23.45
    # Method spencer
    } else if (method == "spencer") {
      B <- (n[i] - 1)*(2*base::pi/365)
      declination[i] <- (180/base::pi)*(0.006918 - 0.399912*cos(B) + 0.070257*sin(B) - 0.006758*cos(2*B))
    }
  }
  out <- dplyr::tibble(date = day_date, method = method, G0 = G0, declination = declination)
  return(out)
}



#' Solar angle minimum and maximum
#'
#' Compute the solar angle for a latitude in different dates.
#'
#' @param lat integer, latitude.
#' @param day_date vector of dates in the format `%YYYY-%MM-%DD`-
#' @param day_end end date, if it is not NULL will be end date.
#' @param method method used for computation of solar declination, can be `cooper` or `spencer`.
#' @examples
#' solar_angle_minmax(55.3, "2040-01-01", day_end = "2040-12-31")
#' solar_angle_minmax(55.3, c("2040-01-31", "2023-04-01", "2015-09-02"))
#' @name solar_angle_minmax
#' @rdname solar_angle_minmax
#' @return a tibble.
#' @export
solar_angle_minmax <- function(lat = NULL, day_date = Sys.Date(), day_end = NULL, method = "cooper"){

  # Sequence of dates
  day_date <- as.Date(day_date)
  if (!is.null(day_end)) {
    day_end <- as.Date(day_end[1])
    day_date <- seq.Date(day_date[1], day_end, by = "1 day")
  }
  # Latitude from degrees to radiant
  phi <- from_degree_to_radiant(lat)

  out <- list()
  for(i in 1:length(day_date)){
    # Solar declination
    declination <- solar_time_declination(day_date = day_date[i], method = method)
    declination <- from_degree_to_radiant(declination$declination)
    # Solar angle
    omega <- acos(-tan(declination)*tan(phi))
    # Solar angle at sunrise
    omega_min <- -from_radiant_to_degree(omega)
    # Solar angle at sunset
    omega_max <-  from_radiant_to_degree(omega)
    # Solar altitude
    alpha_max <- asin(sin(declination)*sin(phi) + cos(declination)*cos(phi))
    alpha_max <- from_radiant_to_degree(alpha_max)
    # Solar zenit angle
    zenit_angle_max <- 90 - alpha_max
    # Number of sun hours
    sun_hours <- lubridate::dhours(omega_max*2/15)
    # Time adjustment in seconds
    time_adj <- solar_time_adjustment(day_date = day_date[i])$time_adjustment
    # Output
    out[[i]] <- dplyr::tibble(date = day_date[i],
                              lat = lat,
                              declination = from_radiant_to_degree(declination),
                              omega_min = omega_min,
                              omega_max = omega_max,
                              alpha_max = alpha_max,
                              zenit_angle_max = zenit_angle_max,
                              sun_hours = sun_hours,
                              time_adj = time_adj)
  }
  out <- dplyr::bind_rows(out)
  return(out)
}

#' Solar movements
#'
#' Compute the solar angle for a latitude in different times of the day.
#'
#' @param lat latitude
#' @param lon longitude
#' @param day_date_time vector of dates in the format `%YYYY-%MM-%DD HH:MM:SS`
#' @param day_time_end end date, if it is not NULL will be end date.
#' @param method method used for computation of solar declination, can be `cooper` or `spencer`.
#' @examples
#' solar_movements(44.23, 11.20, day_date_time = "2040-01-01", day_time_end = "2040-01-03")
#' @name solar_movements
#' @rdname solar_movements
#' @return a numeric vector containing the time adjustment in minutes.
#' @export
solar_movements <- function(lat = NULL, lon = NULL, day_date_time = NULL, day_time_end = NULL, method = "spencer"){

  # Sequence of dates
  if (!is.null(day_time_end)) {
    day_date_time <- seq.POSIXt(as.POSIXct(day_date_time), as.POSIXct(day_time_end), "1 hour")
  } else {
    day_date_time <- as.POSIXct(day_date_time)
  }

  solar_constants <- list()
  for(i in 1:length(day_date_time)){

    # Convert the datetime into the local hour (lh)
    lh <- as.POSIXlt(day_date_time[i])
    # Convert in a date
    day_date_std <- as.Date(lh)
    # Detect if it is Legal Hour
    #   - from 27 Mar - 30 Oct (TRUE)
    #   - from 30 Oct - 27 Mar (FALSE)
    is_legal_hour <- (lubridate::day(lh) >= 27 & lubridate::month(lh)>= 3) & (lubridate::day(lh) <= 30 & lubridate::month(lh) <= 10)
    if (is_legal_hour) {
      solar_hour <- lh - lubridate::dhours(1)
      lh <- lh - lubridate::dhours(1)
    }  else {
      solar_hour <- lh
    }

    # Solar constants for the selected day
    solar_range <- solar_angle_minmax(lat = lat, day_date = day_date_std)
    # Solar hour for the selected day-time
    solar_hour <- solar_hour + (4*(lon-15) + as.numeric(solar_range$time_adj)/(60*60))
    # Solar angle for the selected day-time
    # Middle date
    middle_date <- as.POSIXlt(paste0(day_date_std, " 12:00:00 CET"))
    # Difference in minutes from middledate
    min_difference <- as.numeric(difftime(solar_hour, middle_date, units = "mins"))
    # Solar angle
    solar_range$omega <- (15*min_difference)/(60*60)

    # Extra information
    solar_range$local_date <- day_date_time[i]
    solar_range$solar_date <- solar_hour
    solar_range$lon <- lon

    # Solar constants
    phi <- from_degree_to_radiant(lat)
    delta <- from_degree_to_radiant(solar_range$declination)
    omega <- from_degree_to_radiant(solar_range$omega)
    # Incidence angle
    solar_range$cosI <- sin(delta)*sin(phi) + cos(delta)*cos(phi)*cos(omega)
    # Solar altitude (Alpha)
    solar_range$alpha <- from_radiant_to_degree(asin(solar_range$cosI))
    # Zenit angle
    zenit_angle <- from_degree_to_radiant(90 - solar_range$alpha)
    solar_range$zenit_angle <- from_radiant_to_degree(zenit_angle)
    # Azimut angle
    solar_azimut <- ifelse(omega > 0, 1, -1)*abs((cos(zenit_angle)*sin(phi) - sin(delta))/(sin(zenit_angle)*cos(phi)))
    solar_range$solar_azimut <- from_radiant_to_degree(solar_azimut)
    # G0 for an horizontal surface
    solar_range$G0w <- solar_range$cosI*solar_time_constant(day_date = day_date_std, method = method)$G0n
    # Select and reorder relevant variables
    solar_range <- dplyr::select(solar_range, local_date, solar_date, lat, lon,
                                 cosI, G0w, declination, omega_min, omega, omega_max, alpha_max, alpha, zenit_angle_max,
                                 zenit_angle_max, solar_azimut, time_adj, sun_hours)

    # Night corrector: solar radiation = 0
    if (solar_range$omega >= solar_range$omega_min & solar_range$omega <= solar_range$omega_max) {
      solar_range$is_night = 0
    } else {
      solar_range$is_night = 1
      solar_range$G0w = 0
      solar_range$cosI = 0
    }
    solar_constants[[i]] <- solar_range
  }
  output <- dplyr::bind_rows(solar_constants)
  return(output)
}

#' Solar extraterrestrial radiation
#'
#' Compute the solar angle for a latitude in different times of the day.
#'
#' @param lat latitude
#' @param day_date vector of dates in the format `%YYYY-%MM-%DD`
#' @param day_end end date, if it is not NULL will be end date.
#' @param method method used for computation of solar declination, can be `cooper` or `spencer`.
#' @examples
#' solar_extraterrestrial_radiation(42.23, "2022-05-01", day_end = "2022-05-31")
#' @name solar_extraterrestrial_radiation
#' @rdname solar_extraterrestrial_radiation
#' @return a numeric vector containing the time adjustment in minutes.
#' @export
solar_extraterrestrial_radiation <- function(lat = NULL, day_date = Sys.Date(), day_end = NULL, method = "spencer"){

  # Sequence of dates
  if (!is.null(day_end)) {
    day_date <- seq.Date(as.Date(day_date), as.Date(day_end), 1)
  } else {
    day_date <- as.Date(day_date)
  }

  output_ls <- list()
  for(i in 1:length(day_date)){

    # Time-adjusted solar constant
    G0n <- solar_time_constant(day_date = day_date[i], method = method)$G0n
    # Solar angles
    solar_angles <- solar_angle_minmax(lat = lat, day_date = day_date[i])
    # Solar constants
    phi <- from_degree_to_radiant(lat)
    delta <- from_degree_to_radiant(solar_angles$declination)
    omega_max <-from_degree_to_radiant(solar_angles$omega_max)
    # H0: daily irradiantce J/m^2/day
    cosZ <- cos(phi)*cos(delta)*sin(omega_max) + (base::pi*solar_angles$omega_max/180)*sin(phi)*sin(delta)
    solar_angles$G0n <- G0n
    solar_angles$cosZ <- cosZ
    solar_angles$time_factor <- 24*3600/base::pi
    solar_angles$G0 <- solar_angles$G0n*solar_angles$cosZ*solar_angles$time_factor
    output_ls[[i]] <- solar_angles
  }
  output <- dplyr::bind_rows(output_ls)
  return(output)
}

#' Solar clear sky hourly
#'
#' Compute the clear sky radiation for hourly data.
#'
#' @param cosZ cosine angle of incidence
#' @param G0 extraterrestrial radiation
#' @param altitude altitude in meters.
#' @param clime correction for different climes, can be `No Correction`, `Summer`, `Winter`, `Subartic Summer`, `Tropical`.
#'
#' @examples
#' solar_clearsky_hourly(cosZ = 0.4, G0 = 4, altitude = 2.5, clime = "No Correction")
#'
#' @name solar_clearsky_hourly
#' @rdname solar_clearsky_hourly
#' @return a numeric vector containing the time adjustment in minutes.
#' @export
solar_clearsky_hourly <- function(cosZ = NULL, G0 = NULL, altitude = 2.5, clime = "No Correction"){

  # correction for different climes
  clime <- match.arg(clime[1], choices = c("No Correction", "Summer", "Winter", "Subartic Summer", "Tropical"))

  # altitude must be converted from metre to km
  altitude <- altitude/1000
  if (altitude > 2.5) {
    a0_star <- 0.6*(1-exp(-0.214*(altitude - 1.12)))
  } else {
    a0_star <- 0.4237 - 0.00821*(6.0 - altitude)^2
  }

  a1_star <- 0.5055 - 0.00595*(6.5 - altitude)^2
  a2_star <- 0.2711 - 0.01858*(2.5 - altitude)^2

  # Correction for Climetypes
  if(clime == "No Correction"){

    a0 <- a0_star
    a1 <- a1_star
    a2 <- a2_star

  } else if (clime == "Summer"){

    a0 <- a0_star*0.97
    a1 <- a1_star*0.99
    a2 <- a2_star*1.02

  } else if (clime == "Winter"){

    a0 <- a0_star*1.03
    a1 <- a1_star*1.01
    a2 <- a2_star*1.00

  } else if (clime == "Subartic Summer"){

    a0 <- a0_star*0.99
    a1 <- a1_star*0.99
    a2 <- a2_star*1.01

  } else if (clime == "Tropical"){

    a0 <- a0_star*0.95
    a1 <- a1_star*0.98
    a2 <- a2_star*1.02
  }

  output <- dplyr::tibble(tau_beam = a0 + a1*exp(-a2/cosZ), tau_diffuse = 0.271 - 0.294*tau_beam)
  output <- dplyr::mutate(output,
                          tau_beam = ifelse(tau_beam > 1 | tau_beam  < 0, 0, tau_beam),
                          tau_diffuse = ifelse(tau_diffuse > 1 | tau_diffuse  < 0, 0, tau_diffuse))

  skymax <- G0*output$tau_beam + G0*output$tau_diffuse
  return(skymax)
}


