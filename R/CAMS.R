#' CAMS Data (Location object)
#' @param place place
#'
#' @examples
#' object <- CAMS("Bologna")
#'
#' @export
CAMS <- function(place, year_max = lubridate::year(Sys.Date()), date_min = NULL, date_max = NULL){

  #' @examples
  #' place <- "Bologna"
  scale_GHI <- 1000
  # Match a location in the dataset
  place <- match.arg(place, choices = names(CAMS_data), several.ok = FALSE)
  data_place <- solarr::CAMS_data[[place]]
  # Filter for minimum date
  if (!is.null(date_min)) {
    data_place <- dplyr::filter(data_place, date >= date_min)
  }
  # Filter for maximum date
  if (!is.null(date_max)) {
    data_place <- dplyr::filter(data_place, date <= date_max)
  }
  # # Filter for maximum year
  data_place <- dplyr::filter(data_place, Year <= year_max)
  # Create dataset
  data <- dplyr::mutate(data_place, GHI = GHI/scale_GHI, clearsky = clearsky_GHI/scale_GHI)
  # Add seasonal variables
  data$Year <- lubridate::year(data$date) # year
  data$Month <- lubridate::month(data$date) # month of the year
  data$Day <- lubridate::day(data$date) # day of the month
  data$n <- solarr::number_of_day(data$date) # number of the day of the year
  data <- dplyr::select(data, date, n, Year, Month, Day, GHI, clearsky)
  # Extraterrestrial radiation for an year with 366 days
  H0_extra <- solar_extraterrestrial_radiation(lat = attr(solarr::CAMS_data[[place]], "lat"),
                                               day_date = "2020-01-01", day_end = "2020-12-31")
  # Month, Day and rescale G0 in kWh/m2
  H0_extra <- dplyr::mutate(H0_extra, H0 = G0/(3600*scale_GHI), omega = 2*base::pi/365,
                            Day = lubridate::day(date), Month = lubridate::month(date))
  H0_extra$n <- solarr::number_of_day(H0_extra$date)
  H0_extra <- dplyr::select(H0_extra, -lat, -date)
  H0_extra <- dplyr::select(H0_extra, Month, Day, n, H0, dplyr::everything())
  data <- dplyr::left_join(data, dplyr::select(H0_extra, Month, Day, H0), by = c("Month", "Day"))

  structure(
    list(
      place = data_place$place[1],
      coords = list(lat = attr(solarr::CAMS_data[[place]], "lat"),
                    lon = attr(solarr::CAMS_data[[place]], "lon"),
                    alt = attr(solarr::CAMS_data[[place]], "alt")),
      seasonal_data = H0_extra,
      data = data
    ),
    class = c("Location", "list")
  )
}

#' get CAMS data
#'
#'
#' @export
getCAMS <- function(place, lat, lon, alt, from = "2005-01-01", to = Sys.Date()){

  to <- as.character(to)
  from <- as.character(from)

  # File Name in Temp (as csv file)
  file_name <- paste0("cams-", place)
  # Create in a temporary directory a temporary file csv
  temp <- base::tempfile(pattern = file_name, fileext = ".csv")
  # Download the file
  cams_solar_radiation_ts(latitude = lat, longitude = lon, start = from, end = to, altitude = alt, filename = temp)
  # Read the file
  response <- readr::read_csv(temp, show_col_types = FALSE)
  # Metadata
  meta <- readr::read_delim(temp,  skip = 0, n_max = 40, delim = ":", show_col_types = FALSE, progress = FALSE)
  meta <- dplyr::bind_cols(variable = c("from", "to", "lat", "lon", "alt", "tz"), meta[9:14, 2]) %>%
    dplyr::mutate(variable = stringr::str_trim(variable, side = "both")) %>%
    tidyr::spread("variable", ` utf-8`) %>%
    dplyr::mutate(alt = as.numeric(alt),
                  lat = as.numeric(lat),
                  lon = as.numeric(lon),
                  to = as.POSIXct(to),
                  from = as.POSIXct(from),
                  tz = dplyr::case_when(
                    tz == " Universal time (UT)" ~ "UTC",
                    TRUE ~ tz
                  ))
  # Solar data
  data <- readr::read_delim(temp, skip = 42, delim = ";", show_col_types = FALSE, progress = FALSE)
  # Standard col names
  colnames(data) <- c("date", "TOA", "clearsky_GHI", "clearsky_BHI", "clearsky_DHI", "clearsky_BNI", "GHI", "BHI", "DHI", "BNI", "reliability")
  # Add metadata
  data$lat <- meta$lat
  data$place <- place
  data$lon <- meta$lon
  data$alt <- meta$alt
  # Add date elements
  data$date <- as.POSIXct(purrr::map_chr(strsplit(data$date, "/"), ~.x[1]), tz = meta$tz)
  data$Year <- lubridate::year(data$date)
  data$Month <- lubridate::month(data$date)
  data$Day <- lubridate::day(data$date)
  data$date <- as.Date(data$date)
  # Reorder variables
  data <- dplyr::select(data, date, place, lat, lon, alt, Year, Month, Day, TOA,
                        GHI, clearsky_GHI, BHI, clearsky_BHI, DHI, clearsky_DHI, BNI, clearsky_BNI, reliability)
  # Unlink the connection created with temp
  base::unlink(temp)
  return(data)
}


#' update CAMS data
#'
#'
#' @export
updateCAMS <- function(place, lat, lon, alt, from = "2005-01-01", to = Sys.Date()-1, CAMS_data = solarr::CAMS_data, quiet = FALSE){

  if (missing(place)) {
    places <- names(CAMS_data)
  } else {
    places <- place
  }

  bar_msg <- paste0(rep("-", 30), collapse = "")
  for(i in 1:length(places)){
    place <- places[i]
    if (place %in% names(CAMS_data)) {
      if (!quiet) message(bar_msg," ", "Updating existing place: ", place, " ", bar_msg)
      data <- CAMS_data[[place]]
      lat <- attr(CAMS_data[[place]], "lat")
      lon <- attr(CAMS_data[[place]], "lon")
      alt <- attr(CAMS_data[[place]], "alt")
      from = as.character(max(data$date))
      to = as.character(Sys.Date())
    } else {
      if (!quiet) message(bar_msg," ", "Importing new place: ", place, " ", bar_msg)
      lat <- lat[i]
      lon <- lon[i]
      alt <- alt[i]
      from = as.character(from)
      to = as.character(as.Date(to)-1)
    }
    # Get new data
    newdata <- getCAMS(place, lat, lon, alt, from = from, to = to)
    newdata <- dplyr::select(newdata, -lat, -lon, -alt)
    newdata <- newdata[!is.nan(newdata$GHI),]
    if (place %in% names(CAMS_data)){
      CAMS_data[[place]] <- dplyr::bind_rows(data, newdata)
      CAMS_data[[place]] <- CAMS_data[[place]][!duplicated(CAMS_data[[place]]$date),]
      locations[locations$place == place,]$to <<- max(CAMS_data[[place]]$date)
      locations[locations$place == place,]$nobs <<- nrow(CAMS_data[[place]])
    } else {
      CAMS_data[[place]] <- newdata
      attr(CAMS_data[[place]], "lat") <- lat
      attr(CAMS_data[[place]], "lon") <- lon
      attr(CAMS_data[[place]], "alt") <- alt
      locations <<- dplyr::bind_rows(
        locations,
        dplyr::tibble(
          place = place,
          lat = lat,
          lon = lon,
          alt = alt,
          from = min(CAMS_data[[place]]$date),
          to = max(CAMS_data[[place]]$date),
          nobs = nrow(CAMS_data[[place]]),
        )
      )
    }
    locations <<- dplyr::arrange(locations, lat)
    if (!quiet) message(bar_msg," ", "Done.", place, " ", bar_msg)
  }
  return(CAMS_data)
}

