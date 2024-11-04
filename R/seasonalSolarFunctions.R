#' Solar seasonal functions
#'
#' @examples
#' dates <- seq.Date(as.Date("2022-01-01"), as.Date("2022-12-31"), 1)
#' # Seasonal functions object
#' sf <- seasonalSolarFunctions$new()
#'
#' # Adjustment parameter
#' sf$B(number_of_day(dates))
#'
#' # Time adjustment in minutes
#' sf$E(dates)
#'
#' # Declination
#' sf$declination(dates)
#'
#' # Solar constant
#' sf$G0
#'
#' # Solar constant adjusted
#' sf$G0n(dates)
#'
#' # Extraterrestrial radiation
#' sf$H0(dates, 43)
#'
#' # Number of hours of sun
#' sf$sun_hours(dates, 43)
#'
#' # Sunset hour angle
#' sf$sunset_hour_angle(dates, 43)
#'
#' @name seasonalSolarFunctions
#' @rdname seasonalSolarFunctions
#'
#' @references Duffie, Solar Engineering of Thermal Processes Fourth Edition.
#' @export
seasonalSolarFunctions <- R6::R6Class("seasonalSolarFunctions",
                                public = list(
                                  #' @field legal_hour Logical, when `TRUE` the clock time will be corrected for the legal hour.
                                  legal_hour = TRUE,
                                  #' @description
                                  #' Initialize a `seasonalSolarFunctions` object
                                  #' @param method character, method type for computations. Can be `cooper` or `spencer`.
                                  #' @param legal_hour Logical, when `TRUE` the clock time will be corrected for the legal hour.
                                  initialize = function(method = "spencer", legal_hour = TRUE){
                                    # Method of computation
                                    private$method_ <- match.arg(method, choices = c("spencer", "cooper"))
                                    self$legal_hour <- legal_hour
                                  },
                                  #' @description
                                  #' Extract or update the method used for computations.
                                  #' @param x character, method type. Can be `cooper` or `spencer`.
                                  #' @return When `x` is missing it return a character containing the method that
                                  #' is actually used.
                                  update_method = function(x){
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
                                  #' Compute the time adjustment in minutes.
                                  #' @param n number of the day of the year
                                  #' @details The function implement Eq. 1.5.3 from Duffie (4th edition), i.e.
                                  #' \deqn{E = 229.2(0.000075 + 0.001868 \cos(B) - 0.032077\sin(B) - 0.014615\cos(2B) - 0.04089\sin(2B))}
                                  #' @return The time adjustment in minutes.
                                  E = function(n){
                                    n <- number_of_day(n)
                                    B <- self$B(n-1)
                                    # Time adjustment in minutes
                                    E <- 229.2*(0.000075 + 0.001868*cos(B) - 0.032077*sin(B) - 0.014615*cos(2*B) - 0.04089*sin(2*B))
                                    # Time adjustment in seconds
                                    return(lubridate::dminutes(E))
                                  },
                                  #' @description
                                  #' Compute the solar time from a clock time.
                                  #' @param x datetime, clock hour.
                                  #' @param lon longitude of interest in degrees.
                                  #' @param lon_sd longitude of the Local standard meridian in degrees.
                                  #' @details The function implement Eq. 1.5.2 from Duffie (4th edition), i.e.
                                  #' \deqn{solartime = clocktime + 4 (lon_s - lon) + E(n)}
                                  #' @return A datetime object
                                  solar_time = function(x, lon, lon_sd = 15){
                                    # Convert in a datetime
                                    date_h <- as.POSIXlt(x)
                                    # Adjust the date for legal hour
                                    if (self$legal_hour) {
                                      start_legal_hour <- as.Date(paste0(lubridate::year(date_h[1]), "-03-10"))
                                      end_legal_hour <- as.Date(paste0(lubridate::year(date_h[1]), "-10-30"))
                                      #   - from 27 Mar - 30 Oct (TRUE)
                                      #   - from 30 Oct - 27 Mar (FALSE)
                                      is_legal_hour <- (as.Date(date_h) >= start_legal_hour) & (as.Date(date_h) <= end_legal_hour)
                                      date_h <- dplyr::if_else(is_legal_hour, date_h - lubridate::dhours(1), date_h)
                                    }
                                    # Solar hour for the selected day-time
                                    date_h - lubridate::dminutes(4*(lon_sd-lon)) + self$E(date_h)
                                  },
                                  #' @description
                                  #' Compute the solar angle for a specific hour of the day.
                                  #' @param x datetime, clock hour.
                                  #' @param lon longitude of interest in degrees.
                                  #' @param lon_sd longitude of the Local standard meridian in degrees.
                                  #' @return An angle in degrees
                                  hour_angle = function(x, lon, lon_sd = 15){
                                    # Solar hour
                                    LST <- self$solar_time(x, lon, lon_sd)
                                    # Local solar time
                                    solar_hour <- lubridate::hour(LST) + lubridate::minute(LST)/60
                                    # Solar angle in degrees
                                    15*(solar_hour - 12)
                                  },
                                  #' @description
                                  #' Compute the incidence angle
                                  #' @param x datetime, clock hour.
                                  #' @param lat latitude of interest in degrees.
                                  #' @param lon longitude of interest in degrees.
                                  #' @param lon_sd longitude of the Local standard meridian in degrees.
                                  #' @param beta altitude
                                  #' @param gamma orientation
                                  #' @return An angle in degrees
                                  incidence_angle = function(x, lat, lon, lon_sd = 15, beta = 0, gamma = 0){
                                    # Altitude of the surface in radiant
                                    beta = self$radiant(beta)
                                    # Orientation of the surface in radiant
                                    gamma = self$radiant(gamma)
                                    # Latitude in radiant
                                    phi = self$radiant(lat)
                                    # Solar hour angle
                                    omega = self$radiant(self$hour_angle(x, lon, lon_sd))
                                    # Declination in radiant
                                    delta = self$radiant(self$declination(x))
                                    # Components
                                    T_ = sin(delta)*(sin(phi)*cos(beta) - cos(phi)*sin(beta)*cos(gamma))
                                    U_ = cos(delta)*(cos(phi)*cos(beta) + sin(phi)*sin(beta)*cos(gamma))
                                    V_ = cos(delta)*sin(beta)*sin(gamma)
                                    # Cosine of the angle of incidence
                                    cos_theta_z = T_ + U_*cos(omega) + V_*sin(omega)
                                    # Angle of incidence
                                    self$degree(acos(cos_theta_z))
                                  },
                                  #' @description
                                  #' Compute the solar azimuth angle for a specific time of the day.
                                  #' @param x datetime, clock hour.
                                  #' @param lat latitude of interest in degrees.
                                  #' @param lon longitude of interest in degrees.
                                  #' @param lon_sd longitude of the Local standard meridian in degrees.
                                  #' @details The function implement Eq. 1.6.6 from Duffie (4th edition), i.e.
                                  #' \deqn{\gamma_s = sign(\omega) \left|\cos^{-1}\left( \frac{\cos \theta_z \sin \phi - \sin \delta}{\sin \theta_z \cos \phi} \right) \right|}
                                  #' @return The solar azimut angle in degrees
                                  azimut_angle = function(x, lat, lon, lon_sd = 15){
                                    # Declination in radiant
                                    delta <- self$radiant(self$declination(x))
                                    # Latitude in radiant
                                    phi = self$radiant(lat)
                                    # Solar hour angle
                                    omega = self$radiant(self$hour_angle(x, lon, lon_sd))
                                    # The angle of incidence is the zenith angle of the sun
                                    theta_z = self$radiant(self$incidence_angle(x, lat, lon, lon_sd, beta = 0, gamma = 0))
                                    # Azimut angle
                                    gamma_s <- sign(omega)*abs(acos((cos(theta_z)*cos(phi) - sin(delta))/(sin(theta_z)*cos(phi))))
                                    self$degree(gamma_s)
                                  },
                                  #' @description
                                  #' Compute the solar constant adjusted for the day of the year.
                                  #' @param n number of the day of the year
                                  #' @details When method is `cooper` the function implement Eq. 1.4.1a from Duffie (4th edition), i.e.
                                  #' \deqn{G_{0,n} = G_0 (1 + 0.033\cos(B))}
                                  #' otherwise when it is `spencer` it implement Eq. 1.4.1b from Duffie (4th edition):
                                  #' \deqn{G_{0,n} = G_0 (1.000110 + 0.034221\cos(B) + 0.001280\sin(B) + 0.000719\cos(2B) + 0.000077\sin(2B))}
                                  #' @return The solar constant in \eqn{W/m^2} for the day n.
                                  G0n = function(n){
                                    n <- number_of_day(n)
                                    if (private$method_ == "cooper") {
                                      # Method cooper
                                      G0n <- private$..G0*(1 + 0.033 * cos(self$B(n)))
                                    } else if (private$method_ == "spencer") {
                                      # Method spencer
                                      B <- self$B(n-1)
                                      G0n <- private$..G0*(1.000110 + 0.034221*cos(B) + 0.001280*sin(B) + 0.000719*cos(2*B) + 0.000077*sin(2*B))
                                    }
                                    return(G0n)
                                  },
                                  #' @description
                                  #' Compute solar declination in degrees.
                                  #' @param n number of the day of the year
                                  #' @details When method is `cooper` the function implement Eq. 1.6.1a from Duffie (4th edition), i.e.
                                  #' \deqn{\delta(n) = 23.45 \sin \left(\frac{2 \pi (284 + n)}{365}\right)}
                                  #' otherwise when it is `spencer` it implement Eq. 1.6.1b from Duffie (4th edition):
                                  #' \deqn{\delta(n) = \frac{180}{\pi}(0.006918 - 0.399912\cos(B) + 0.070257\sin(B) - 0.006758\cos(2B))}
                                  #' @return The solar declination in degrees.
                                  declination = function(n){
                                    n <- number_of_day(n)
                                    if (private$method_ == "cooper") {
                                      declination <- sin(2*base::pi*(284 + n)/365)*23.45
                                    } else if (private$method_ == "spencer") {
                                      B <- self$B(n-1)
                                      declination <- (180/base::pi)*(0.006918 - 0.399912*cos(B) + 0.070257*sin(B) - 0.006758*cos(2*B))
                                    }
                                    return(declination)
                                  },
                                  #' @description
                                  #' Compute the solar extraterrestrial radiation
                                  #' @param n number of the day of the year
                                  #' @param lat latitude of interest in degrees.
                                  #' @return Extraterrestrial radiation on an horizontal surface in kilowatt hour for metres squared for day.
                                  H0 = function(n, lat){
                                    # Latitude in radiant
                                    phi <- self$radiant(lat)
                                    # Declination in radiant
                                    delta <- self$radiant(self$declination(n))
                                    # Sunset hour angle in degrees
                                    omega_s <- self$sunset_hour_angle(n, lat)
                                    # Extraterrestrial radiation in daily joules/m^2 per day
                                    Gn <- (self$G0n(n)*(24*3600))/base::pi
                                    H0 <- Gn*(cos(phi)*cos(delta)*sin(self$radiant(omega_s)) + ((base::pi*omega_s)/180)*sin(phi)*sin(delta))
                                    # Convert in (kilowatt hour)/m^2 per day
                                    H0/3600000
                                  },
                                  #' @description
                                  #' Compute solar angle at sunset in degrees
                                  #' @param n number of the day of the year
                                  #' @param lat Numeric, latitude of interest in degrees.
                                  #' @details The function implement Eq. 1.6.10 from Duffie (4th edition), i.e.
                                  #' \deqn{\omega_s = \cos^{-1}(-\tan(\delta(n))\tan(\phi))}
                                  #' @return The sunset hour angle in degrees.
                                  sunset_hour_angle = function(n, lat){
                                    # Latitude from degrees to radiant
                                    phi <- self$radiant(lat)
                                    # Declination from degrees to radiant
                                    declination <- self$radiant(self$declination(n))
                                    # Sunset hour angle in degrees
                                    self$degree(acos(-tan(declination)*tan(phi)))
                                  },
                                  #' @description
                                  #' Compute number of sun hours for a day n.
                                  #' @param n number of the day of the year.
                                  #' @param lat Numeric, latitude of interest in degrees.
                                  #' @details The function implement Eq. 1.6.11 from Duffie (4th edition), i.e.
                                  #' \deqn{\frac{2}{15} \omega_s}
                                  sun_hours = function(n, lat){
                                    # Number of sun hours
                                    sun_hours <- self$sunset_hour_angle(n, lat)*(2/15)
                                    return(lubridate::dhours(sun_hours))
                                  },
                                  #' @description
                                  #' Compute solar altitude in degrees
                                  #' @param n number of the day of the year
                                  #' @param lat Numeric, latitude of interest in degrees.
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
                                  #' Compute the solar angle for a latitude in different dates.
                                  #' @param x datetime, clock hour.
                                  #' @param lat Numeric, latitude of interest in degrees.
                                  #' @param lon Numeric, longitude of interest in degrees.
                                  #' @param lon_sd Numeric, longitude of the Local standard meridian in degrees.
                                  #' @param by Character, time step. Default is `1 min`.
                                  solar_angles = function(x, lat, lon, lon_sd, by = "1 min"){
                                    day_date <- as.Date(x[1])
                                    start_date <- as.POSIXct(paste0(day_date, " 00:00:00"))
                                    end_date <- as.POSIXct(paste0(day_date+1, " 00:00:00"))
                                    day_date_seq <- seq.POSIXt(start_date, end_date, by = by)
                                    # Latitude from degrees to radiant
                                    phi <- self$radiant(lat)
                                    # Solar declination
                                    declination <- self$declination(day_date)
                                    # Solar angle at sunset
                                    omega_max <-  self$sunset_hour_angle(day_date, lat)
                                    # Solar angle at sunrise
                                    omega_min <- -omega_max
                                    # Number of sun hours
                                    sun_hours <- self$sun_hours(day_date, lat)
                                    # Time adjustment in seconds
                                    E <- self$E(day_date)
                                    # Solar angle
                                    omega <- self$hour_angle(day_date_seq, lon, lon_sd)
                                    # Solar time
                                    solartime <- self$solar_time(day_date_seq, lon, lon_sd)
                                    # Incidence angle
                                    theta_z <- self$incidence_angle(day_date_seq, lat, lon, lon_sd = lon_sd, beta = 0, gamma = 0)

                                    dplyr::tibble(date = day_date,
                                                  clocktime = day_date_seq,
                                                  solartime = solartime,
                                                  lat = lat,
                                                  lon = lon,
                                                  omega = omega,
                                                  declination = declination,
                                                  omega_min = omega_min,
                                                  omega_max = omega_max,
                                                  sun_hours = sun_hours,
                                                  theta_z = theta_z,
                                                  E = E) %>%
                                      dplyr::filter(omega > omega_min & omega < omega_max)
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
                                    } else if (clime == "Tropical") {
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
