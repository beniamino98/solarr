#' Spatial Grid
#'
#' Create a grid from a range of latitudes and longitudes.
#'
#' @return a tibble with two columns `lat` and `lon`.
#'
#' @examples
#' grid <- spatialGrid$new()
#' grid$make_grid()
#' grid$grid
#' grid$weights <- IDW(beta = 2)
#' grid$is_known_point(49.95, 12.15)
#' grid$known_point(c(44.85, 44.9), c(12.15, 11.2))
#' grid$is_inside_bounds(44.8, 10.9)
#' grid$neighborhoods(44.9, 12.1)
#' grid$neighborhoods(44.95, 12.15)
#' filter(grid$grid, lat == 44.95 & lon == 12.15)
#'
#' @rdname spatialGrid
#' @name spatialGrid
#' @export
spatialGrid <- R6::R6Class("spatialGrid",
                           public = list(
                             weights = NA,
                             initialize = function(lat_min = 43.7, lat_max = 45.1, lon_min = 9.2, lon_max = 12.7, lat_by = 0.1, lon_by = 0.1, weights = IDW(2)){
                               # Latitudes bounds
                               private$lat_min <- lat_min
                               private$lat_max <- lat_max
                               private$lat_by <- lat_by
                               # Sequence of latitudes
                               private$lats = round(seq(lat_min, lat_max, by = lat_by), digits = 5)
                               # Longitudes bounds
                               private$lon_min <- lon_min
                               private$lon_max <- lon_max
                               private$lon_by <- lon_by
                               # Sequence of longitudes
                               private$lons = round(seq(lon_min, lon_max, by = lon_by), digits = 5)
                               # Weighting function
                               self$weights <- weights
                             },
                             make_grid = function(labels){
                               # Matrix of latitudes and longitudes
                               grid <- dplyr::tibble()
                               for(lat in private$lats) {
                                 grid_lat <- dplyr::tibble(lat = lat, lon = private$lons)
                                 grid <- dplyr::bind_rows(grid, grid_lat)
                               }
                               # Assign a unique ID for each place
                               if (!missing(labels)) {
                                 if (length(labels) == nrow(grid)) {
                                   grid <- dplyr::bind_cols(ID = labels, grid)
                                 }
                               } else {
                                 labels = paste0("ID_", 1:nrow(grid))
                                 grid <- dplyr::bind_cols(id = labels, grid)
                               }
                               private$..grid <- grid
                             },
                             #' @description
                             #' Check if a point is inside the bounds of the spatial grid.
                             #' @param lat numeric, latitude of a location.
                             #' @param lon numeric, longitude of a location.
                             #' @return `TRUE` when the point is inside the limits and `FALSE` otherwise.
                             is_inside_bounds = function(lat, lon){
                               coords <- dplyr::tibble(lat = lat, lon = lon)
                               condition <- rep(TRUE, nrow(coords))
                               for(i in 1:nrow(coords)) {
                                 condition[i] <- coords$lon[i] >= private$lon_min & coords$lon[i] <= private$lon_max
                                 condition[i] <- condition[i] & coords$lat[i] >= private$lat_min & coords$lat[i] <= private$lat_max
                               }
                              return(condition)
                             },
                             #' @description
                             #' Check if a point is already in the spatial grid
                             #' @param lat numeric, latitude of a location.
                             #' @param lon numeric, longitude of a location.
                             #' @return `TRUE` when the point is a known point and `FALSE` otherwise.
                             is_known_point = function(lat, lon){
                               coords <- dplyr::tibble(lat = lat, lon = lon)
                               condition <- rep(NA, nrow(coords))
                               for(i in 1:nrow(coords)) {
                                 condition_ <- dplyr::filter(self$grid, lat == coords$lat[i] & lon == coords$lon[i])
                                 condition[i] <- nrow(condition_) != 0
                               }
                               return(condition)
                             },
                             #' @description
                             #' Return the ID and coordinates of a point that is already in the spatial grid
                             #' @param lat numeric, latitude of a location.
                             #' @param lon numeric, longitude of a location.
                             known_point = function(lat, lon){
                               coords <- dplyr::tibble(lat = lat, lon = lon)
                               output <- dplyr::tibble()
                               for(i in 1:nrow(coords)) {
                                 if (self$is_known_point(coords$lat[i], coords$lon[i])) {
                                   df_known_point <- dplyr::filter(self$grid, lat == coords$lat[i] & lon == coords$lon[i])
                                   output <- dplyr::bind_rows(output, df_known_point)
                                 }
                               }
                               return(output)
                             },
                             #' @description
                             #' Find the n-closest neighborhoods of a point
                             #' @param lat numeric, latitude of a point in the grid.
                             #' @param lon numeric, longitude of a point in the grid.
                             #' @param n number of neighborhoods
                             neighborhoods = function(lat, lon, n = 4){
                               if (self$is_known_point(lat, lon)) {
                                 nb <- self$known_point(lat, lon)
                                 # Add target latitude and longitude
                                 nb$lat_E <- lat
                                 nb$lon_E <- lon
                                 # No interpolation if it is a known point
                                 nb$wgt <- 1
                               } else {
                                 # Grid of locations
                                 nb <- self$grid
                                 # Add target latitude and longitude
                                 nb$lat_E <- lat
                                 nb$lon_E <- lon
                                 # Compute distances from target
                                 nb <- dplyr::mutate(nb, dist = havDistance(lat, lon, lat_E, lon_E))
                                 # Arrange by ascending distances
                                 nb <- head(dplyr::arrange(nb, dist), n = n)
                                 # Compute the weights
                                 nb$wgt <- self$weights(nb$dist)
                                 # Normalize the weights
                                 nb$wgt <- nb$wgt/sum(nb$wgt)
                               }
                               # Output
                               return(nb)
                             }
                           ),
                           private = list(
                             ..grid = NA,
                             lats = NA,
                             lons = NA,
                             lat_min = NA,
                             lat_max = NA,
                             lon_min = NA,
                             lon_max = NA,
                             lon_by = NA,
                             lat_by = NA
                           ),
                           active = list(
                             grid = function(){
                               private$..grid
                             }
                           )
                          )

#' Haversine distance
#'
#' Compute the Haversine distance between two points.
#'
#' @param lat_1 numeric, latitude of first point.
#' @param lon_1 numeric, longitude of first point.
#' @param lat_2 numeric, latitude of second point.
#' @param lon_2 numeric, longitude of second point.
#'
#' @examples
#' havDistance(43.3, 12.1, 43.4, 12.2)
#' havDistance(43.35, 12.15, 43.4, 12.2)
#' @name havDistance
#' @rdname havDistance
#' @return Numeric vector the distance in kilometers.
#' @export
havDistance <- function(lat_1, lon_1, lat_2, lon_2){
  r <- 6371 # Radius of the earth in kilometres
  dist <- c()
  # Solar functions
  sf <- seasonalSolarFunctions$new()
  # Compute distances
  for(i in 1:length(lat_1)){
    # Coordinates
    coord_A <- c(lat_1[i], lon_1[i])
    coord_B <- c(lat_2[i], lon_2[i])
    # Convert in radiant
    point_A <- sf$radiant(coord_A)
    point_B <- sf$radiant(coord_B)
    # Latitudes
    phi <- c(point_A[1], point_B[1])
    delta_phi <- phi[2] - phi[1]
    # Longitudes
    lam <- c(point_A[2], point_B[2])
    delta_lam <- lam[2] - lam[1]
    # Compute distance
    dist[i] <- 2*r*asin(sqrt(0.5*(1 - cos(delta_phi) + cos(phi[1])*cos(phi[2])*(1-cos(delta_lam)))))
  }
  return(dist)
}


#' Inverse Distance Weighting Functions
#'
#' Return a distance weighting function
#'
#' @param beta parameter used in exponential and power functions.
#' @param d0 parameter used only in exponential function.
#' @details When the parameter `d0` is not specified the function returned will be of power type otherwise of exponential type.
#'
#' @examples
#' # Power weighting
#' IDW_pow <- IDW(2)
#' IDW_pow(c(2, 3,10))
#' IDW_pow(c(2, 3,10), normalize = TRUE)
#' # Exponential weighting
#' IDW_exp <- IDW(2, d0 = 5)
#' IDW_exp(c(2, 3,10))
#' IDW_exp(c(2, 3,10), normalize = TRUE)
#' @export
IDW <- function(beta, d0){
  # Power function
  IDW.pow <- function(beta) {
    function(d, normalize = FALSE){
      w <- 1/(d^beta)
      if (normalize) {
        w <- w/sum(w)
      }
      return(w)
    }
  }
  # Exponential function
  IDW.exp <- function(beta, d0) {
    function(d, normalize = FALSE){

      w <- exp(-(d/d0)^beta)
      if (normalize) {
        w <- w/sum(w)
      }
      return(w)
    }
  }
  # Output function
  if (missing(d0)) {
    IDW.pow(beta)
  } else {
    IDW.exp(beta, d0)
  }
}




