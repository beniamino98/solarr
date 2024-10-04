#' Spatial model object
#'
#' @rdname spatialModel
#' @name spatialModel
#' @export
spatialModel <- R6::R6Class("spatialModel",
                            public = list(
                              #' @description
                              #' Initialize the spatial model
                              #' @param locations A tibble with columns `place`, `lat`, `lon`, `from`, `to`, `nobs`.
                              #' @param models A list of `solarModel` objects
                              #' @param paramsModels A `spatialParameters` object.
                              #' @param beta numeric, used in exponential and power functions.
                              #' @param d0 numeric, used only in exponential function.
                              #' @param quiet logical
                              initialize = function(locations, models, paramsModels, beta = 2, d0, quiet = FALSE){
                                private$..locations <- locations
                                private$..models <- models
                                private$..parameters <- paramsModels
                                private$weight <- IDW(beta, d0)
                                private$quiet <- quiet
                              },
                              #' @description
                              #' Find the n-closest neighborhoods of a point
                              #' @param lat numeric, latitude of a point in the grid.
                              #' @param lon numeric, longitude of a point in the grid.
                              #' @param n number of neighborhoods
                              neighborhoods = function(lat, lon, n = 4){
                                # Extract the neighborhood given a pair of coordinates
                                nb <- spatialModel_neighborhoods(private$..locations, lat, lon, n = n, weight = private$weight)
                                # Output
                                return(nb)
                              },
                              #' @description
                              #' Check if a point is already in the spatial grid
                              #' @param lat numeric, latitude of a location.
                              #' @param lon numeric, longitude of a location.
                              #' @return `TRUE` when the point is a known point and `FALSE` otherwise.
                              is_known_location = function(lat, lon){
                                locations <- private$..locations
                                idx_known_location <- c()
                                if (self$is_inside_limits(lat, lon)) {
                                  # Known points in the grid
                                  known_lat <- round(locations$lat, digits = 3)
                                  knonw_lon <- round(locations$lon, digits = 3)
                                  # Check if the point is a known point in the grid
                                  idx_known_location <- which(rowSums(c(lat, lon) == cbind(known_lat, knonw_lon)) == 2)
                                }
                                if (purrr::is_empty(idx_known_location)) {
                                  return(FALSE)
                                } else {
                                  return(TRUE)
                                }
                              },
                              #' @description
                              #' Get a known model in the grid from place or coordinates.
                              #' @param place character, id of the location.
                              #' @param lat numeric, latitude of a location.
                              #' @param lon numeric, longitude of a location.
                              gridModel = function(place, lat, lon){
                                locations <- private$..locations
                                if (missing(place)) {
                                  if (self$is_known_location(lat, lon)) {
                                    # Known points in the grid
                                    known_lat <- round(locations$lat, digits = 3)
                                    knonw_lon <- round(locations$lon, digits = 3)
                                    # Check if the point is a known point in the grid
                                    idx_known_location <- which(rowSums(c(lat, lon) == cbind(known_lat, knonw_lon)) == 2)
                                    # Get the model place
                                    place <- locations[idx_known_location,]$place
                                  }
                                }
                                # Get the model at the location
                                return(private$..models[[place]])
                              },
                              #' @description
                              #' Check if a point is inside the limits of the spatial grid.
                              #' @param lat numeric, latitude of a location.
                              #' @param lon numeric, longitude of a location.
                              #' @return `TRUE` when the point is inside the limits and `FALSE` otherwise.
                              is_inside_limits = function(lat, lon){
                                locations <- private$..locations
                                # Known points in the grid
                                known_lat <- round(locations$lat, digits = 3)
                                knonw_lon <- round(locations$lon, digits = 3)
                                # Check if the point is outside the limit of the grid
                                if (lon < min(knonw_lon) || lon > max(knonw_lon) || lat < min(known_lat) || lat > max(known_lat)) {
                                  msg <- paste0("Point (", lat, ", ", lon, ") outside the limits of the grid ",
                                                "(", min(known_lat), " - ", max(known_lat), " | ", min(knonw_lon), " - ", max(knonw_lon), ")")
                                  if(!private$quiet) message(msg)
                                  return(FALSE)
                                }
                                return(TRUE)
                              },
                              #' @description
                              #' Perform the bilinear interpolation for a target variable.
                              #' @param lat numeric, latitude of the location to be interpolated.
                              #' @param lon numeric, longitude of the location to be interpolated.
                              #' @param target character, name of the target variable to interpolate.
                              #' @param day_date date for interpolation, if missing all the available dates will be used.
                              #' @param n number of neighborhoods to use for interpolation.
                              interpolator = function(lat, lon, target = "GHI", n = 4, day_date){

                                if (!private$quiet) message("Interpolating ", target, " (Lat: ", lat, " Lon: ", lon, ")\r", appendLF = FALSE)
                                # Initialize the output dataset
                                interp_data <- dplyr::tibble(x1 = lat, x2 = lon, x3 = lubridate::as_date(NA), x4 = NA, x5 = FALSE)
                                colnames(interp_data) <- c("lat", "lon", "date", target, "interpolated")
                                # Check if the point is outside the limit of the grid
                                if (!self$is_inside_limits(lat, lon)) {
                                  return(interp_data)
                                }
                                # Check if the point is a known point in the grid
                                if (self$is_known_location(lat, lon)){
                                  # Get the data at the location
                                  data <- self$gridModel(lat = lat, lon = lon)$data
                                  # Filter for specific dates
                                  if (!missing(day_date)) {
                                    data <- dplyr::filter(data, date %in% day_date)
                                  }
                                  # Returned data are NOT interpolated
                                  interp_data <- as.list(interp_data)
                                  interp_data$date <- data$date
                                  interp_data[[target]] <- data[[target]]
                                  interp_data$interpolated <- FALSE
                                  interp_data <- dplyr::bind_cols(interp_data)
                                  return(interp_data)
                                }
                                # Detect n-neighborhoods in the grid
                                nb <- self$neighborhoods(lat, lon, n = n)
                                # Extract neighborhoods models
                                nb_models <- private$..models[nb$place]
                                # Bilinear interpolation
                                spatial_interp <- spatialModel_interpolator(nb_models, target = target, n = n, weights = nb$wgt, day_date)
                                # Add latitudes and longitudes
                                interp_data <- dplyr::bind_cols(interp_data[,c(1,2)], spatial_interp)
                                return(interp_data)
                              },
                              #' @description
                              #' Interpolator function for a `solarModel` object
                              #' @param lat numeric, latitude of a point in the grid.
                              #' @param lon numeric, longitude of a point in the grid.
                              #' @param n number of neighborhoods
                              solarModel = function(lat, lon, n = 4){
                                # Check if the point is outside the limit of the grid
                                if (!self$is_inside_limits(lat, lon)) {
                                  return(invisible(NULL))
                                }
                                # Check if the point is a known point in the grid
                                if (self$is_known_location(lat, lon)) {
                                  # Get the model at the location
                                  model <- self$gridModel(lat = lat, lon = lon)
                                  return(model)
                                }
                                # Detect n-neighborhoods in the grid
                                nb <- self$neighborhoods(lat, lon, n = n)
                                # Initialize a model to be updated
                                model_inter <- private$..models[nb$place][[1]]$clone(deep = TRUE)
                                # Update coordinates
                                model_inter$coords <- list(lat = lat, lon = lon, alt = NA)
                                # Update place label
                                model_inter$place <- paste0(nb$place, collapse = "-")
                                # Interpolate the realized GHI
                                model_inter$.__enclos_env__$private$..data[["GHI"]] <- self$interpolator(lat, lon, target = "GHI", n = n)$GHI
                                # Interpolate the realized Clearsky
                                model_inter$.__enclos_env__$private$..data[["clearsky"]] <- self$interpolator(lat, lon, target = "clearsky", n = n)$clearsky
                                # Update H0 inside clear sky model
                                H0 <- seasonalSolarFunctions$new("spencer")$H0(model_inter$seasonal_data$n, lat)$H0
                                model_inter$.__enclos_env__$private$..seasonal_model_Ct$seasonal_data$H0 <- H0
                                # Predict the best parameters for the target location
                                params <- private$..parameters$predict(lat = lat, lon = lon, as_tibble = FALSE)[[1]]
                                # Update the model parameters
                                model_inter$update(params)
                                model_inter$filter()
                                model_inter
                              },
                              #' @description
                              #' Compute monthly moments for mixture with 16 components
                              #' @param lat numeric, latitude of a point in the grid.
                              #' @param lon numeric, longitude of a point in the grid.
                              #' @param nmonths numeric, months to consider
                              #' @param nobs.min numeric, minimum number of joint states under which the state is considered with 0 probability.
                              combinations = function(lat, lon, nmonths = 1:12, nobs.min = 3){
                                # Find the n-closest neighborhoods of a point
                                nb <- self$neighborhoods(lat, lon, n = 4)
                                # Extract neighborhoods models
                                nb_models <- private$..models[nb$place]
                                # Compute mixture combinations
                                combinations <- spatialModel_combinations(nb_models, nmonths = nmonths, weights = nb$wgt, nobs.min = nobs.min)
                                # output
                                return(combinations)
                              }
                            ),
                            private = list(
                              # Grid of points
                              ..locations = NA,
                              # Models for each point in the locations grid
                              ..models = NA,
                              # Spatial parameters models
                              ..parameters = NA,
                              # Slot for weighting function
                              weight = NA,
                              quiet = FALSE
                            ),
                            active = list(
                              #' @field models list of `solarModel` objects
                              models = function(){
                                private$..models
                              },
                              #' @field locations dataset with all the locations.
                              locations = function(){
                                private$..locations
                              },
                              #' @field parameters `spatialParameters` object
                              parameters = function(){
                                private$..parameters
                              }
                            )
)
