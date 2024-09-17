#' Spatial model object
#'
#' @param locations grid of locations
#' @param solarModels list of `solarModel` objects
#'
#' @rdname spatialModel
#' @name spatialModel
#' @export
spatialModel <- function(locations, solarModels){

  structure(
    list(
      locations = locations,
      models = solarModels,
      params_models = list()
    ),
    class = c("spatialModel", "list")
  )
}

#' Find the n-closest neighborhoods of a point
#'
#' @param object a `spatialModel` object
#' @param lat numeric, latitude of the point.
#' @param lon numeric, longitude of the point.
#' @param n number of neighborhoods
#' @inheritParams IDW
#'
#' @rdname spatialModel_neighborhoods
#' @name spatialModel_neighborhoods
#' @export
spatialModel_neighborhoods <- function(object, lat, lon, n = 4, beta = 2, d0){

  # Extract grid of locations
  best_locations <- object$locations
  # Add target latitude and longitude
  best_locations$lat_E <- lat
  best_locations$lon_E <- lon
  # Compute distances from target
  best_locations <- dplyr::mutate(best_locations, dist = havDistance(lat, lon, lat_E, lon_E))
  # Arrange by ascending distances
  best_locations <- head(dplyr::arrange(best_locations, dist), n = n)
  # Compute the weights
  weight <- IDW(beta, d0)
  best_locations$wgt <- weight(best_locations$dist, normalize = TRUE)
  best_locations <- dplyr::mutate(best_locations, wgt = ifelse(dist == 0, 1, wgt))
  best_locations$wgt <- best_locations$wgt/sum(best_locations$wgt)
  # Output
  return(best_locations)
}

#' Compute all possible states
#'
#' @inheritParams spatialModel_neighborhoods
#'
#' @rdname spatialModel_combinations
#' @name spatialModel_combinations
#'
#' @export
spatialModel_combinations <- function(object, lat, lon){

  # Find the n-closest neighborhoods of a point
  nb <- spatialModel_neighborhoods(object, lat, lon, n = 4)
  # Weights in vector notation
  w <- matrix(nb$wgt, ncol = 1)
  # Extract neighborhoods models
  models <- object$models[nb$place]
  # Possible combinations in a grid with 4 reference locations
  possible_combinations <- data.frame(
    State = 1:16,
    Month = rep(1, 16),
    Location_A = c(1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0),
    Location_B = c(1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0),
    Location_C = c(1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0),
    Location_D = c(1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0)
  )

  # Main loop
  m <- 1
  combinations <- list()
  for(m in 1:12){
    # Initialize the dataset
    combinations[[m]] <- dplyr::as_tibble(possible_combinations)
    combinations[[m]]$probs <- 0
    combinations[[m]]$mean <- 0
    combinations[[m]]$sd <- 0
    combinations[[m]]$n <- 0
    combinations[[m]]$cv <- as.list(1:16)
    # Extract data
    df_A <- dplyr::filter(models[[1]]$data, Month == m & isTrain)
    df_B <- dplyr::filter(models[[2]]$data, Month == m & isTrain)
    df_C <- dplyr::filter(models[[3]]$data, Month == m & isTrain)
    df_D <- dplyr::filter(models[[4]]$data, Month == m & isTrain)
    # Select only relevant variables
    df_A <- dplyr::select(df_A, date, B_A = "B", u_A = "u")
    df_B <- dplyr::select(df_B, date, B_B = "B", u_B = "u")
    df_C <- dplyr::select(df_C, date, B_C = "B", u_C = "u")
    df_D <- dplyr::select(df_D, date, B_D = "B", u_D = "u")
    # Create a unique dataset
    df <- dplyr::left_join(df_A, df_B, by = "date") %>%
      dplyr::left_join(df_C, by = "date") %>%
      dplyr::left_join(df_D, by = "date")
    i <- 1
    for(i in 1:16){
      # Extract the states
      states <- unlist(possible_combinations[i,][-c(1:2)])
      # Compute joint probabilities
      combinations[[m]]$probs[i] <- mean(df$B_A == states[1] & df$B_B == states[2] & df$B_C == states[3] & df$B_D == states[4], na.rm = TRUE)
      # Extract the realized joint series
      df_u <- dplyr::filter(df, B_A == states[1] & B_B == states[2] & B_C == states[3] & B_D == states[4])
      df_u <- dplyr::select(df_u, dplyr::contains("u"))
      # Number of observations for the estimates
      combinations[[m]]$n[i] <- nrow(df_u)
      # If the number of cases is too low set probability to zero
      if (nrow(df_u) <= 10) {
        combinations[[m]]$probs[i] <- 0
        next
      }
      # Variance-covariance matrix
      combinations[[m]]$cv[[i]] <- cov(df_u)
      # Expected value for the i-th state
      combinations[[m]]$mean[i] <- sum(dplyr::summarise_all(df_u, mean)*w)
      combinations[[m]]$sd[i] <- sqrt(t(w) %*% combinations[[m]]$cv[[i]] %*% w)[1]
    }
    # Normalize the probabilities
    combinations[[m]]$probs <- combinations[[m]]$probs/sum(combinations[[m]]$probs)
    combinations[[m]]$Month <- m
    colnames(combinations[[m]]) <- c("State", "Month", names(models), "probs", "mean", "sd","n", "cv", "mu2")
  }
  dplyr::bind_rows(combinations)
}


#' Interpolate the solar radiation for a location
#'
#' @inheritParams spatialModel_neighborhoods
#' @param day_date day date for interpolation. If missing all the available dates will be used.
#' @param quiet logical
#'
#' @rdname spatialModel_interpolate_GHI
#' @name spatialModel_interpolate_GHI
#'
#' @export
spatialModel_interpolate_GHI <- function(object, lat, lon, n = 4, day_date, quiet = FALSE, ...){

  if (!quiet) message("Lat: ", lat, " Lon: ", lon, "\r", appendLF = FALSE)
  # Known points in the grid
  known_lat <- round(object$locations$lat, digits = 3)
  knonw_lon <- round(object$locations$lon, digits = 3)
  if (lon < min(knonw_lon) || lon > max(knonw_lon) || lat < min(known_lat) || lat > max(known_lat)){
    msg <- paste0("Point (", lat, ", ", lon, ") outside the grid ", "(",
                  min(known_lat), " - ", max(known_lat), " | ", min(knonw_lon), " - ", max(knonw_lon), ")")
    if(!quiet) message(msg)
    return(dplyr::tibble(date = lubridate::as_date(NA), lat = lat, lon = lon, GHI = NA, interp = FALSE))
  }
  idx_known_location <- which(rowSums(c(lat, lon) == cbind(known_lat, knonw_lon)) == 2)
  # If the point is known it is returned the known model
  if (!purrr::is_empty(idx_known_location)){
    id_location <- object$locations[idx_known_location,]$place
    model <- object$models[[id_location]]
    interp_GHI <- model$data
    if (!missing(day_date)) {
      day_date <- as.Date(day_date)
      interp_GHI <-dplyr::filter(interp_GHI, date == day_date)
      interp_GHI$interp <- FALSE
    }
  } else {
    # Detect n-neighborhoods in the grid
    nb <- spatialModel_neighborhoods(object, lat, lon, n = n, ...)
    # Extract neighborhoods models
    models <- object$models[nb$place]
    if (!missing(day_date)) {
      day_date <- as.Date(day_date)
      l_data <- purrr::map(models, ~dplyr::filter(.x$data, date == day_date))
      interp_GHI <- purrr::map2_df(l_data, nb$wgt, ~dplyr::bind_cols(dplyr::select(.x, date, GHI), wgt = .y))
      interp_GHI <- dplyr::summarise(interp_GHI, date = day_date, GHI = sum(GHI*wgt))
    } else {
      # Interpolate the realized GHI
      interp_GHI <- purrr::map2(models, nb$wgt, ~dplyr::bind_cols(dplyr::select(.x$data, date, GHI), wgt = .y))
      interp_GHI <- dplyr::bind_cols(purrr::map(interp_GHI, ~.x$GHI*.x$wgt))
      interp_GHI$GHI <- rowSums(interp_GHI)
      interp_GHI$date <- models[[1]]$data$date
    }
    interp_GHI$interp <- TRUE
  }
  # Add interpolated coordinates
  interp_GHI$lat <- lat
  interp_GHI$lon <- lon
  # Select only the time series of GHI
  dplyr::select(interp_GHI, date, lat, lon, GHI, interp)
}


#' Compute a solar model for a location
#'
#' @inheritParams spatialModel_neighborhoods
#' @param quiet logical
#'
#' @rdname spatialModel_interpolate
#' @name spatialModel_interpolate
#'
#' @export
spatialModel_interpolate <- function(object, lat, lon, n = 4, quiet = FALSE, ...){

  # Known points in the grid
  known_lat <- round(object$locations$lat, digits = 3)
  knonw_lon <- round(object$locations$lon, digits = 3)
  if (lon < min(knonw_lon) || lon > max(knonw_lon) || lat < min(known_lat) || lat > max(known_lat)){
    msg <- paste0("Point (", lat, ", ", lon, ") outside the grid ", "(",
                  min(known_lat), " - ", max(known_lat), " | ", min(knonw_lon), " - ", max(knonw_lon), ")")
    if(!quiet) message(msg)
    return(dplyr::tibble(date = lubridate::as_date(NA), lat = lat, lon = lon, GHI = NA, interp = FALSE))
  }
  # If the point is known it is returned the known model
  idx_lon_lat <- which(lat == known_lat & lon == knonw_lon)
  if (!purrr::is_empty(idx_lon_lat)){
    id_location <- object$locations$place[idx_lon_lat]
    return(object$models[[id_location]])
  }

  # Detect n-neighborhoods in the grid
  nb <- spatialModel_neighborhoods(object, lat, lon, n = n, ...)
  # Interpolate the realized GHI
  interp_GHI <- spatialModel_interpolate_GHI(object, lat = lat, lon = lon, n = n, quiet = quiet, ...)
  # Initialize a model that will be updated
  model <- object$models[nb$place][[1]]
  # Update place label
  model$place <- paste0(nb$place, collapse = "-")
  # Update coordinates
  model$coords <- list(lat = lat, lon = lon, alt = NA)
  # Add the interpolated GHI
  data <- dplyr::select(model$data, -GHI)
  data <- dplyr::left_join(data, dplyr::select(interp_GHI, date, GHI), by = "date")
  data <- dplyr::select(data, date:Day, GHI, dplyr::everything())
  data <- dplyr::filter(data, !is.na(GHI) & !is.na(date))
  # Update the dataset
  model$data <- data

  # Compute extraterrestrial radiation at the given location
  seasonal_H0 <- solar_extraterrestrial_radiation(lat = lat, day_date = "2020-01-01", day_end = "2020-12-31")
  seasonal_H0 <- dplyr::bind_cols(n = 1:nrow(model$seasonal_data), dplyr::select(seasonal_H0, H0 = "G0"))
  seasonal_H0$H0 <- seasonal_H0$H0/(3600*1000)
  model$seasonal_data$H0 <-  seasonal_H0$H0
  # Update H0 insiede clearsky model
  model$seasonal_model_Ct$.__enclos_env__$private$seasonal_data <- seasonal_H0
  # Predict the best parameters for the target location
  params <- spatialParameters_predict(object$params_models, lat = lat, lon = lon, as_tibble = FALSE)[[1]]
  # Update the model parameters
  model_upd <- solarModel_update(model, params)
  # Update model time series
  model_upd <- solarModel_filter(model_upd)
  return(model_upd)
}
