#' Detect the season
#'
#' Detect the season from a vector of dates
#'
#' @param x vector of dates in the format `YYYY-MM-DD`.
#' @param invert logica, when `TRUE` the seasons will be inverted.
#'
#' @examples
#' detect_season("2040-01-31")
#' detect_season(c("2040-01-31", "2023-04-01", "2015-09-02"))
#'
#' @name detect_season
#' @rdname detect_season
#'
#' @return a character vector containing the correspondent season. Can be `spring`, `summer`, `autumn`, `winter`.
#' @export
detect_season <- function(x, invert = FALSE){
  # first day of the year
  date_start <- as.Date("2000-01-01")
  # start dates for the seasons
  date_spring <- as.Date("2000-03-23")
  date_summer <- as.Date("2000-06-23")
  date_autumn <- as.Date("2000-09-23")
  date_winter <- as.Date("2000-12-23")
  # last day of the year
  date_end <- as.Date("2000-12-31")

  seasons <- dplyr::tibble(
    date = as.Date(x),
    Month = lubridate::month(date),
    Day = lubridate::day(date)
  )
  # Standardized year
  seasons <- dplyr::mutate(seasons, date = as.Date(paste0("2000-", Month, "-", Day)))
  seasons <- dplyr::mutate(seasons,
                           season = dplyr::case_when(
                             date < date_spring  & date >= date_start ~ "Winter",
                             date >= date_spring & date < date_summer ~ "Spring",
                             date >= date_summer & date < date_autumn ~ "Summer",
                             date >= date_autumn & date < date_winter ~ "Autumn",
                             date >= date_winter & date <= date_end ~ "Winter"))
  if (invert) {
    seasons <- dplyr::mutate(seasons,
                             season = dplyr::case_when(
                               season == "Winter" & invert ~ "Summer",
                               season == "Summer" & invert ~ "Winter",
                               season == "Spring" & invert ~ "Autumn",
                               season == "Autumn" & invert ~ "Spring",
                               TRUE ~ season))
  }
  return(seasons$season)
}

#' Is leap year?
#'
#' Check if a given year is leap (366 days) or not (365 days).
#'
#' @param x numeric value or dates vector in the format `YYYY-MM-DD`.
#'
#' @examples
#' is_leap_year("2024-02-01")
#' is_leap_year(c(2023:2030))
#' is_leap_year(c("2024-10-01", "2025-10-01"))
#' is_leap_year("2029-02-01")
#' @name is_leap_year
#' @rdname is_leap_year
#' @return Boolean. `TRUE` if it is a leap year, `FALSE` otherwise.
#' @export
is_leap_year <- function(x){
  if (is.numeric(x)) {
    date_year <- x
  } else {
    date_year <- lubridate::year(lubridate::as_date(x))
  }
  return(date_year %% 4 == 0)
}

#' Number of day
#'
#' Compute the number of day of the year given a vector of dates.
#'
#' @param x dates vector in the format `YYYY-MM-DD`.
#'
#' @examples
#' number_of_day("2040-01-31")
#' number_of_day(c("2040-01-31", "2023-04-01", "2015-09-02"))
#' number_of_day(c("2029-02-28", "2029-03-01", "2020-12-31"))
#' number_of_day(c("2020-02-29", "2020-03-01", "2020-12-31"))
#' @name number_of_day
#' @rdname number_of_day
#' @return Numeric vector with the number of the day during the year.
#' @export
number_of_day <- function(x){

  if (missing(x)) {
    warning('The argument "x" is missing while it require a number or a character of the form %YYYY-%MM-%DD')
    return(NULL)
  }
  # Check if `x` is already numeric
  if (is.numeric(x)) {
    x[x > 365] <- x[x > 365] %% 365
    return(x)
  }

  # Standard sequence of dates
  seq_date_year <- seq.Date(as.Date("2010-01-01"), as.Date("2010-12-31"), by = "1 day")
  seq_date_year_leap <- seq.Date(as.Date("2012-01-01"), as.Date("2012-12-31"), by = "1 day")

  i <- 1
  n_of_day <- c()
  day_date <- lubridate::as_date(x)
  for(i in 1:length(day_date)){
    day <- day_date[i]
    # Check if is a valid date
    if (is.na(day)) {
      x[i] <- NA_integer_
      next
    }
    # Extract month and day
    m <- lubridate::month(day)
    d <- lubridate::day(day)
    # Check leap year
    if (is_leap_year(day)) {
      # Approximation for 29-02
      if (m == 2 & d == 29){
        n_of_day[i] <- 59.5
      } else {
        day <- as.Date(paste0("2012-", m,"-", d))
        # Extract the number of the day
        n_of_day[i] <- which(seq_date_year_leap == day)
        # Rescale after 29-02 for leap years
        if (m > 2) {
          n_of_day[i] <- n_of_day[i]-1
        }
      }
    } else {
      day <- as.Date(paste0("2010-", m,"-", d))
      # extract the number of the day
      n_of_day[i] <- which(seq_date_year == day)
    }
  }
  names(n_of_day) <- x
  return(n_of_day)
}

#' Solar seasonal functions
#'
#' @name seasonalSolarFunctions
#' @rdname seasonalSolarFunctions
#' @examples
#' sf <- seasonalSolarFunctions$new()
#' sf$angle_minmax("2022-01-01", 44)
#' sf$H0(1:365, 44)
#' @export
seasonalSolarFunctions <- R6::R6Class("seasonalSolarFunctions",
                                public = list(
                                  #' @description
                                  #' Initialize a `seasonalSolarFunctions` object
                                  #' @param method character, method type for computations. Can be `cooper` or `spencer`.
                                  initialize = function(method = "spencer"){
                                    # Method of computation
                                    private$method_ <- match.arg(method, choices = c("spencer", "cooper"))
                                  },
                                  #' @description
                                  #' Extract or update the method used for computations.
                                  #' @param x character, method type. Can be `cooper` or `spencer`.
                                  #' @return When `x` is missing it return a character containing the method that
                                  #' is actually used.
                                  method = function(x){
                                    if (missing(x)) {
                                      return(private$method_)
                                    } else {
                                      # Old method
                                      old_method <- private$method_
                                      # New method
                                      private$method_ <- match.arg(x, choices = c("spencer", "cooper"))
                                      return(private$method_)
                                    }
                                  },
                                  #' @description
                                  #' Seasonal adjustment parameter.
                                  #' @param n number of the day of the year
                                  #' @details The function computes
                                  #' \deqn{B(n) = \frac{2\pi}{365} n}
                                  B = function(n){
                                    (2*base::pi*n)/365
                                  },
                                  #' @description
                                  #' Convert angles in radiant into an angles in degrees.
                                  #' @param x numeric vector, angles in radiant.
                                  #' @details The function computes:
                                  #' \deqn{\frac{x 180}{\pi}}
                                  degree = function(x){
                                    x*180/base::pi
                                  },
                                  #' @description
                                  #' Convert angles in degrees into an angles in radiant
                                  #' @param x numeric vector, angles in degrees.
                                  #' @details The function computes:
                                  #' \deqn{\frac{x \pi}{180}}
                                  radiant = function(x){
                                    x*base::pi/180
                                  },
                                  #' @description
                                  #' Compute solar time adjustment in seconds
                                  #' @param n number of the day of the year
                                  #' @details The function computes
                                  #' \deqn{229.2(0.000075 + 0.001868 \cos(B) - 0.032077\sin(B) - 0.014615\cos(2B) - 0.04089\sin(2B))}
                                  time_adjustment = function(n){
                                    n <- number_of_day(n)
                                    B <- self$B(n-1)
                                    # Time adjustment in minutes
                                    time_adj <- 229.2*(0.000075 + 0.001868*cos(B) - 0.032077*sin(B) - 0.014615*cos(2*B) - 0.04089*sin(2*B))
                                    # Time adjustment in seconds
                                    time_adj <- lubridate::seconds(time_adj*60)
                                    return(time_adj)
                                  },
                                  #' @description
                                  #' Compute solar constant
                                  #' @param n number of the day of the year
                                  #' @details If the selected method is `cooper`, the function computes:
                                  #' \deqn{G_{0,n} = G_0 (1 + 0.033\cos(B))}
                                  #' otherwise when it is `spencer` it computes:
                                  #' \deqn{G_{0,n} = G_0 (1.000110 + 0.034221\cos(B) + 0.001280\sin(B) + 0.000719\cos(2B) + 0.000077\sin(2B))}
                                  G0n = function(n){
                                    n <- number_of_day(n)
                                    # Method cooper
                                    if (private$method_ == "cooper") {
                                      G0n <- private$..G0*(1 + 0.033 * cos(self$B(n)))
                                      # Method spencer
                                    } else if (private$method_ == "spencer") {
                                      B <- self$B(n-1)
                                      G0n <- private$..G0*(1.000110 + 0.034221*cos(B) + 0.001280*sin(B) + 0.000719*cos(2*B) + 0.000077*sin(2*B))
                                    }
                                    return(G0n)
                                  },
                                  #' @description
                                  #' Compute solar declination
                                  #' @param n number of the day of the year
                                  #' @details If the selected method was `cooper`, the function computes:
                                  #' \deqn{\delta(n) = 23.45 \sin \left(\frac{2 \pi (284 + n)}{365}\right)}
                                  #' otherwise when it is `spencer` it computes:
                                  #' \deqn{\delta(n) = \frac{180}{\pi}(0.006918 - 0.399912\cos(B) + 0.070257\sin(B) - 0.006758\cos(2B))}
                                  declination = function(n) {
                                    n <- number_of_day(n)
                                    if (private$method_ == "cooper") {
                                      declination <- sin(2*base::pi*(284 + n)/365)*23.45
                                      declination <- self$degree(declination)
                                    } else if (private$method_ == "spencer") {
                                      B <- self$B(n-1)
                                      declination <- (180/base::pi)*(0.006918 - 0.399912*cos(B) + 0.070257*sin(B) - 0.006758*cos(2*B))
                                    }
                                    return(declination)
                                  },
                                  #' @description
                                  #' Compute solar angle at sunset in degrees
                                  #' @param n number of the day of the year
                                  #' @param lat latitude in degrees.
                                  #' @details The function computes
                                  #' \deqn{\cos^{-1}(-\tan(\delta(n))\tan(\phi))}
                                  solar_angle = function(n, lat){
                                    # Latitude from degrees to radiant
                                    phi <- self$radiant(lat)
                                    declination <- self$radiant(self$declination(n))
                                    self$degree(acos(-tan(declination)*tan(phi)))
                                  },
                                  #' @description
                                  #' Compute solar altitude in degrees
                                  #' @param n number of the day of the year
                                  #' @param lat latitude in degrees.
                                  #' @details The function computes
                                  #' \deqn{\sin^{-1}(-\sin(\delta(n))\sin(\phi) + \cos(\delta(n))\cos(\phi))}
                                  solar_altitude = function(n, lat){
                                    # Latitude from degrees to radiant
                                    phi <- self$radiant(lat)
                                    declination <- self$radiant(self$declination(n))
                                    alpha_max <- asin(sin(declination)*sin(phi) + cos(declination)*cos(phi))
                                    self$degree(alpha_max)
                                  },
                                  #' @description
                                  #' Compute number of sun hours
                                  #' @param n number of the day of the year
                                  #' @param lat latitude in degrees.
                                  #' @details The function computes
                                  sun_hours = function(n, lat){
                                    # Solar angle at sunset
                                    omega <- self$solar_angle(n, lat)
                                    # Number of sun hours
                                    sun_hours <- lubridate::dhours(omega*2/15)
                                    return(sun_hours)
                                  },
                                  #' @description
                                  #' Compute the solar angle for a latitude in different dates.
                                  #' @param n number of the day of the year
                                  #' @param lat latitude in degrees.
                                  angle_minmax = function(n, lat) {
                                    # Latitude from degrees to radiant
                                    phi <- self$radiant(lat)
                                    # Solar declination
                                    declination <- self$declination(n)
                                    # Solar angle at sunset
                                    omega_max <-  self$solar_angle(n, lat)
                                    # Solar angle at sunrise
                                    omega_min <- -omega_max
                                    # Solar altitude
                                    alpha_max <- self$solar_altitude(n, lat)
                                    # Solar zenit angle
                                    zenit_angle <- 90 - alpha_max
                                    # Number of sun hours
                                    sun_hours <- self$sun_hours(n, lat)
                                    # Time adjustment in seconds
                                    time_adjustment <- self$time_adjustment(n)
                                    dplyr::tibble(n = n,
                                                  lat = lat,
                                                  declination = declination,
                                                  omega_min = omega_min,
                                                  omega_max = omega_max,
                                                  alpha_max = alpha_max,
                                                  zenit_angle = zenit_angle,
                                                  sun_hours = sun_hours,
                                                  time_adjustment = time_adjustment)

                                  },
                                  #' @description
                                  #' Compute the incidence angle
                                  #' @param n number of the day of the year
                                  #' @param lat latitude in degrees.
                                  cosZ = function(n, lat){
                                    # Solar constants
                                    phi <- self$radiant(lat)
                                    delta <- self$radiant(self$declination(n))
                                    omega_max <- self$radiant(self$solar_angle(n, lat))
                                    # H0: daily irradiantce J/m^2/day
                                    cosZ <- cos(phi)*cos(delta)*sin(omega_max) + (base::pi*self$degree(omega_max)/180)*sin(phi)*sin(delta)
                                    return(cosZ)
                                  },
                                  #' @description
                                  #' Compute the solar extraterrestrial radiation
                                  #' @param n number of the day of the year
                                  #' @param lat latitude in degrees.
                                  H0 = function(n, lat){
                                    dplyr::tibble(
                                      lat = lat,
                                      n = n,
                                      G0n = self$G0n(n),
                                      cosZ = self$cosZ(n, lat),
                                      G0 = G0n*cosZ*(24*3600)/base::pi,
                                      H0 = G0/(3600000)
                                    )
                                  },
                                  #' @description
                                  #' Compute the solar hour
                                  #' @param x datehour
                                  #'
                                  solar_hour = function(x){
                                    # Convert in a datetime
                                    date_h <- as.POSIXlt(x)
                                    # Convert in date
                                    date_d <- as.Date(date_h)
                                    # Detect if it is Legal Hour
                                    #   - from 27 Mar - 30 Oct (TRUE)
                                    #   - from 30 Oct - 27 Mar (FALSE)
                                    day_ <- lubridate::day(date_d)
                                    mon_ <- lubridate::month(date_d)
                                    is_legal_hour <- (day_ >= 27 & mon_ >= 3) & (day_ <= 30 & mon_ <= 10)
                                    if (is_legal_hour) {
                                      solar_hour <- date_h - lubridate::dhours(1)
                                    }  else {
                                      solar_hour <- date_h
                                    }
                                    # Solar hour for the selected day-time
                                    solar_hour <- solar_hour + (4*(lon-15) + as.numeric(self$time_adjustment(date_d))/(60*60))
                                    return(solar_hour)
                                  },
                                  #' @description
                                  #' Compute the solar angle
                                  #' @param x datehour
                                  #'
                                  omega = function(x){
                                    # Convert in date
                                    date_d <- as.Date(x)
                                    solar_hour <- self$solar_hour(x)
                                    # Solar angle for the selected day-time
                                    # Middle date
                                    middle_date <- as.POSIXlt(paste0(date_d, " 12:00:00 CET"))
                                    # Difference in minutes from middledate
                                    min_difference <- as.numeric(difftime(solar_hour, middle_date, units = "mins"))
                                    # Solar angle
                                    (15*min_difference)/(60*60)
                                  },
                                  #' @description
                                  #' Hottel clearsky
                                  #' @param cosZ solar incidence angle
                                  #' @param G0 solar constant
                                  #' @param altitude altitude in km
                                  #' @param clime clime correction
                                  clearsky = function(cosZ = NULL, G0 = NULL, altitude = 2.5, clime = "No Correction"){
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

                                    # Correction for Clime types
                                    a <- c(a0_star, a1_star, a2_star)
                                    if (clime == "Summer") {
                                      correction_factor <- c(0.97, 0.99, 1.02)
                                      a <- a*correction_factor
                                    } else if (clime == "Winter") {
                                      correction_factor <- c(1.03, 1.01, 1.00)
                                      a <- a*correction_factor
                                    } else if (clime == "Subartic Summer"){
                                      correction_factor <- c(0.99, 0.99, 1.01)
                                      a <- a*correction_factor
                                    } else if (clime == "Tropical"){
                                      correction_factor <- c(0.95, 0.98, 1.02)
                                      a <- a*correction_factor
                                    }
                                    a0 <- a[1]
                                    a1 <- a[2]
                                    a2 <- a[3]
                                    output <- dplyr::tibble(tau_beam = a0 + a1*exp(-a2/cosZ), tau_diffuse = 0.271 - 0.294*tau_beam)
                                    output <- dplyr::mutate(output,
                                                            tau_beam = ifelse(tau_beam > 1 | tau_beam  < 0, 0, tau_beam),
                                                            tau_diffuse = ifelse(tau_diffuse > 1 | tau_diffuse  < 0, 0, tau_diffuse))
                                    skymax <- G0*output$tau_beam + G0*output$tau_diffuse
                                    return(skymax)
                                  }
                                ),
                                private = list(
                                  ..G0 = 1367,
                                  method_ = ""
                                ),
                                active = list(
                                  #' @field G0 solar constant, i,e, `1367`.
                                  G0 = function(){
                                    private$..G0
                                  }
                                ))
