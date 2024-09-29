#' Spatial Grid
#'
#' Create a grid from a range of latitudes and longitudes.
#'
#' @param range_lat vector with latitudes. Only the minimum and maximum values are considered.
#' @param range_lon vector with longitudes. Only the minimum and maximum values are considered.
#' @param by step for longitudes and latitudes. If two values are specified the first will be used for
#' latitudes and the second for longitudes
#'
#' @return a tibble with two columns `lat` and `lon`.
#' @examples
#' spatialGrid(lat = c(43.7, 43.8), lon = c(12.5, 12.7), by = 0.1)
#' spatialGrid(lat = c(43.7, 43.75, 43.8), lon = c(12.6, 12.6, 12.7), by = c(0.05, 0.01))
#'
#' @rdname spatialGrid
#' @name spatialGrid
#' @export
spatialGrid <- function(lat = c(43.7, 45.1), lon = c(9.2, 12.7), by = c(0.1, 0.1)) {
  # Compute grid bounds
  range_lat <- range(lat)
  range_lon <- range(lon)
  # Grid of latitudes
  by_lat <- ifelse(length(by) == 1, by, by[1])
  grid_lat <- seq(range_lat[1], range_lat[2], by = by_lat)
  # Grid of longitudes
  by_lon <- ifelse(length(by) == 1, by, by[2])
  grid_lon <- seq(range_lon[1], range_lon[2], by = by_lon)
  # Matrix of latitudes and longitudes
  grid_lat_lon <- dplyr::tibble()
  for(lat in grid_lat) {
    df_grid_lat <- dplyr::tibble(lat = lat, lon = grid_lon)
    grid_lat_lon <- dplyr::bind_rows(grid_lat_lon, df_grid_lat)
  }
  return(grid_lat_lon)
}


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
  # Compute distances
  for(i in 1:length(lat_1)){
    # Coordinates
    coord_A <- c(lat_1[i], lon_1[i])
    coord_B <- c(lat_2[i], lon_2[i])
    # Convert in radiant
    point_A <- seasonalSolarFunctions$new()$radiant(coord_A)
    point_B <- seasonalSolarFunctions$new()$radiant(coord_B)
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


#' Inverse Distance Weighting Function
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
      if (normalize) {
        w <- 1/(d^beta)
        w/sum(w)
      } else {
        1/(d^beta)
      }
    }
  }
  # Exponential function
  IDW.exp <- function(beta, d0) {
    function(d, normalize = FALSE){
      if (normalize) {
        w <- exp(-(d/d0)^beta)
        w/sum(w)
      } else {
        exp(-(d/d0)^beta)
      }
    }
  }
  # Output function
  if (missing(d0)) {
    IDW.pow(beta)
  } else {
    IDW.exp(beta, d0)
  }
}
