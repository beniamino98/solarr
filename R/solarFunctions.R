# Compute the spectral distribution for a blackbody (sun)
SpectralDistribution <- function(lambda = NULL, measure = c("nanometer", "micrometer")){

  # choose the measure: nanometer or mircrometer
  measure <- match.arg(measure, choices = c("nanometer", "micrometer"))

  Ts <- 5777          # constant, black-body temperature in (kelvin)
  Rs <- 6.955*10^(8)  # constant, radius of the sun in (meter)
  R  <- 1.495*10^(11) # constant, sun-earth distance in (meter)
  C1 <- 3.742*10^(8)  # constant, in (W micro-meter^4)/meter^2
  C2 <- 1.439*10^(4)  # constant, in (micro-meter*Kelvin)

  # if "lambda" is "nanometer", the final energy will be in (W/m2 x nano-meter)
  if (measure == "nanometer") {
    lambda <- lambda/1000  # from micrometer to nanometer
    spectra <- (Rs/R)^(2)*(C1/(lambda^(5)*(exp(C2/lambda*Ts) - 1))) # Plank's Law
    return(spectra/1000)
  # if "lambda" is "micrometer", the final energy will be in (W/m2 x micro-meter)
  } else if(measure == "micrometer") {
    spectra <- (Rs/R)^(2)*(C1/(lambda^(5)*(exp(C2/lambda*Ts) - 1))) # Plank's Law
    return(spectra)
  }
}

#' Conversion in Radiant or Degrees
#'
#' Convert an angle in radiant into an angle in degrees.
#'
#' @param x numeric vector, angles in radiant or degrees.
#'
#' @return numeric numeric vector.
#'
#' @examples
#' # convert 0.34 radiant in degrees
#' from_radiant_to_degree(0.34)
#'
#' @name from_radiant_to_degree
#' @rdname radiant
#' @export
from_radiant_to_degree <- function(x = NULL){
  assertive::assert_is_numeric(x)
  x <- as.numeric(x)
  y <- x*180/base::pi
  return(y)
}

#' @examples
#' # convert 19.48 degree in radiant
#' from_degree_to_radiant(19.48)
#'
#' @name from_degree_to_radiant
#' @rdname radiant
#' @export
from_degree_to_radiant <- function(x = NULL){
  assertive::assert_is_numeric(x)
  x <- as.numeric(x)
  y <- x*base::pi/180
  return(y)
}

#' Detect Season
#'
#' Detect the season from a vector of dates
#'
#' @param day_date vector of dates in the format `%YYYY-%MM-%DD`.
#'
#' @return a character vector containing the correspondent season.
#' Can be `spring`, `summer`, `autumn`, `winter`.
#'
#' @examples
#' detect_season("2040-01-31")
#' detect_season(c("2040-01-31", "2023-04-01", "2015-09-02"))
#'
#' @name detect_season
#' @rdname detect_season
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
#' @return Boolean. `TRUE` if it is a leap year, `FALSE` otherwise.
#'
#' @examples
#' is_leap_year("2024-02-01")
#' is_leap_year(c("2024-10-01", "2025-10-01"))
#' is_leap_year("2029-02-01")
#' @name is_leap_year
#' @rdname is_leap_year
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

#' Number of Day
#'
#' Compute the number of day of the year given a vector of dates.
#'
#' @param day_date dates vector in the format `%YYYY-%MM-%DD`.
#'
#' @return Numeric vector with the number of the day during the year. Can vary from `1` up to `365` or `366`.
#'
#' @examples
#' # detect the number of the day in 2040-01-31
#' number_of_day("2040-01-31")
#'
#' # detect the number of the day for a vector of dates
#' number_of_day(c("2040-01-31", "2023-04-01", "2015-09-02"))
#'
#' @name number_of_day
#' @rdname number_of_day
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

#' solar_time_adjustment
#'
#' Compute the time adjustment for a date.
#'
#' @param day_date vector of dates in the format `%YYYY-%MM-%DD`
#' @param day_end end date, if it is not NULL will be end date.
#' @return a numeric vector containing the time adjustment in seconds.
#' @examples
#'
#' # detect the season in 2040-01-31
#' solar_time_adjustment("2040-01-31")
#'
#' # detect the season in a vector of dates
#' solar_time_adjustment(c("2040-01-31", "2023-04-01", "2015-09-02"))
#'
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

#' solar_time_constant
#'
#' Compute the solar constant for a date.
#'
#' @param day_date vector of dates in the format `YYYY-MM-DD`.
#' @param day_end end date, if it is not `NULL` will be end date.
#' @param method method used for computation, can be `cooper` or `spencer`.
#' @return a numeric vector containing the solar constant.
#' @examples
#' solar_time_constant("2040-01-31")
#' solar_time_constant(c("2040-01-31", "2023-04-01", "2015-09-02"))
#' @name solar_time_constant
#' @rdname solar_time_constant
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

#' solar_time_declination
#' @name solar_time_declination
#' @description Compute the solar declination for different dates.
#' @param day_date vector of dates in the format `%YYYY-%MM-%DD`
#' @param day_end end date, if it is not NULL will be end date.
#' @param method method used for computation, can be `cooper` or `spencer`.
#' @return a numeric vector containing the solar declination in minutes.
#' @examples
#'
#' # detect the season in 2040-01-31
#' solar_time_declination("2040-01-01", day_end = "2040-12-31")
#'
#' # detect the season in a vector of dates
#' solar_time_declination(c("2040-01-31", "2023-04-01", "2015-09-02"))
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

#' solar_angle_minmax
#' @name solar_angle_minmax
#' @description Compute the solar angle for a latitude in different dates.
#' @param lat integer, latitude.
#' @param day_date vector of dates in the format `%YYYY-%MM-%DD`
#' @param day_end end date, if it is not NULL will be end date.
#' @param method method used for computation of solar declination, can be `cooper` or `spencer`.
#' @return a tibble.
#' @examples
#'
#' solar_angle_minmax(55.3, "2040-01-01", day_end = "2040-12-31")
#'
#' solar_angle_minmax(55.3, c("2040-01-31", "2023-04-01", "2015-09-02"))
#'
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

#' solar_movements
#' @name solar_movements
#' @description compute the solar angle for a latitude in different times of the day.
#' @param lat latitude
#' @param day_date vector of dates in the format `%YYYY-%MM-%DD`
#' @param day_end end date, if it is not NULL will be end date.
#' @param method method used for computation of solar declination, can be `cooper` or `spencer`.
#' @return a numeric vector containing the time adjustment in minutes.
#' @examples
#'
#' # detect the season in 2040-01-31
#' solar_movements(44.23, 11.20, day_date_time = "2040-01-01", day_time_end = "2040-01-03")
#'
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

#' solar_extraterrestrial_radiation
#' @name solar_extraterrestrial_radiation
#' @description compute the solar angle for a latitude in different times of the day.
#' @param lat latitude
#' @param day_date vector of dates in the format `%YYYY-%MM-%DD`
#' @param day_end end date, if it is not NULL will be end date.
#' @param method method used for computation of solar declination, can be `cooper` or `spencer`.
#' @return a numeric vector containing the time adjustment in minutes.
#' @examples
#'
#' # detect the season in 2040-01-31
#' solar_extraterrestrial_radiation(42.23, "2022-05-01", day_end = "2022-05-31")
#'
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

#' solar_clearsky_hourly
#' @name solar_clearsky_hourly
#' @description Compute the clearsky radiation for hourly data.

#' @return a numeric vector containing the time adjustment in minutes.
#' @examples
#'
#' # detect the season in 2040-01-31
#' solar_clearsky(cosZ = 1, altitude = 2.5, clime = c("No Correction", "Summer", "Winter", "Subartic Summer", "Tropical"))
#'
#' @export

solar_clearsky_hourly <- function(cosZ = NULL, G0 = NULL, altitude = 2.5, clime = c("No Correction", "Summer", "Winter", "Subartic Summer", "Tropical")){

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

#' solar_extraterrestrial_radiation_hourly
#' @name solar_extraterrestrial_radiation_hourly
#' @description Compute the extraterrestrial hourly total radiation on an horizontal surface. Note that hottel clear sky max model is included in computation.
#' @param lat latitude
#' @param day_date vector of dates in the format `%YYYY-%MM-%DD`
#' @param day_end end date, if it is not NULL will be end date.
#' @param method method used for computation of solar declination, can be `cooper` or `spencer`.
#' @return a numeric vector containing the time adjustment in minutes.
#' @examples
#'
#' # detect the season in 2040-01-31
#' solar_extraterrestrial_radiation_hourly(44.23, 11.20, day_date_time = "2040-01-01 00:00:00", day_time_end = "2040-01-03 00:00:00")
#'
#' @export

solar_extraterrestrial_radiation_hourly <- function(lat = NULL, lon = NULL,
                                                    day_date_time = NULL, day_time_end = NULL, altitude = 2.5, clime = "No Correction"){

  day_date_time <- as.POSIXct(day_date_time)
  if(!is.null(day_time_end)){
    day_time_end <- as.POSIXct(day_time_end)
    day_date_time <- seq.POSIXt(day_date_time, day_time_end, by = "1 hour")
  }

  output_ls <- list()
  for(i in 1:length(day_date_time)){

    # time dependent solar costant
    G0n <- solar_time_constant(day_date = day_date_time[i], method = method)$G0n

    # last hour
    last_hour <- day_date_time[i] - lubridate::hours(1)

    if(is.na(last_hour)){
      last_hour <- day_date_time[i] - lubridate::hours(1)
    }

    last_hour_mov <- solar_movements(lat = lat, lon = lon, day_date_time = last_hour, day_time_end = NULL, method = method)
    act_hour_mov <- solar_movements(lat = lat, lon = lon,  day_date_time = day_date_time[i], method = method)

    # Conversion in radiant
    phi <- from_degree_to_radiant(lat)
    omega1 <- from_degree_to_radiant(last_hour_mov$omega)
    omega2 <- from_degree_to_radiant(act_hour_mov$omega)
    delta <- from_degree_to_radiant(act_hour_mov$declination)

    # H0: Daily Irradiantce J/m^2/day
    cosZ <- cos(phi)*cos(delta)*(sin(omega2) - sin(omega1)) + (base::pi*(act_hour_mov$omega - last_hour_mov$omega)/180)*sin(phi)*sin(delta)

    season <- detect_season(as.Date(last_hour))
    if (season == "Winter") {
      clime <- "Winter"
    } else if(season == "Summer") {
      clime <- "Summer"
    } else {
      clime = "No Correction"
    }

    is_night <- ifelse(cosZ >= 1 | cosZ <= 0, 0, 1)

    cosZ <- cosZ*is_night

    clearsky_max <- solar_clearsky_hourly(cosZ = cosZ, altitude = altitude, clime = clime)*is_night

    output_ls[[i]] <- dplyr::tibble(date = day_date_time[i],
                                    is_night = is_night,
                                    lat = lat,
                                    lon = lon,
                                    declination = act_hour_mov$declination,
                                    omega1 = last_hour_mov$omega,
                                    omega2 = act_hour_mov$omega,
                                    G0n = G0n,
                                    time_factor = 12*3600/base::pi,
                                    cosZ = cosZ,
                                    G0 = G0n * cosZ*time_factor,
                                    tau_beam = clearsky_max$tau_beam*is_night,
                                    tau_diffuse = clearsky_max$tau_diffuse*is_night,
                                    clear_sky_beam = tau_beam*G0,
                                    clear_sky_diffuse = tau_diffuse*G0)
  }

  output <- dplyr::bind_rows(output_ls)
  return(output)
}


