#' Print method for the class `solarModel`
#' @param object an object of the class `solarModel`.
#' @keywords internal
#' @export
#' @noRd
print.solarModel <- function(object){
  msg_1 <- paste0("Place: ", object$place, " \n (Lat: ", object$coords$lat, "; Lon: ", object$coords$lon, ") \n")
  msg_2 <- paste0("From: ", min(object$data$date), " - ", max(object$data$date), " (Nobs: ", nrow(object$data), ")")
  cat(paste0(msg_1, msg_2))
}

#' Print method for the class `solarOptionPayoff`
#' @param object an object of the class `solarOptionPayoff`.
#' @keywords internal
#' @export
#' @noRd
print.solarOption <- function(object){

  msg_1 <- paste0("------------------------ Solar Option Payoff ------------------------ \n")
  msg_2 <- paste0("Yearly payoff: ", format(object$payoff_year$premium, digits = 5), "\n")
  msg_3 <- paste0("Monthly payoffs: \n ",
                  "  Jan: ", format(object$payoff_month$premium[1], digits = 3),
                  "  Feb: ", format(object$payoff_month$premium[2], digits = 3),
                  "  Mar: ", format(object$payoff_month$premium[3], digits = 3), "\n",
                  "  Apr: ", format(object$payoff_month$premium[4], digits = 3),
                  "  May: ", format(object$payoff_month$premium[5], digits = 3),
                  "  Jun: ", format(object$payoff_month$premium[6], digits = 3), "\n",
                  "  Jul: ", format(object$payoff_month$premium[7], digits = 3),
                  "  Ago: ", format(object$payoff_month$premium[8], digits = 3),
                  "  Sep: ", format(object$payoff_month$premium[9], digits = 3), "\n",
                  "  Oct: ", format(object$payoff_month$premium[10], digits = 3),
                  "  Nov: ", format(object$payoff_month$premium[11], digits = 3),
                  "  Dec: ", format(object$payoff_month$premium[12], digits = 3), "\n")
  cat(paste0(msg_1, msg_2, msg_3))
}

#' Print method for the class `spatialModel`
#' @param object an object of the class `spatialModel`.
#' @keywords internal
#' @export
#' @noRd
print.spatialModel <- function(object){

  range_lat <- range(object$locations$lat)
  range_lon <- range(object$locations$lon)
  n_points <- nrow(object$locations)
  n_models <- length(object$models)
  n_params <- length(object$params_models$models)

  msg_1 <- paste0("Spatial model: ", "Points", " (", n_points, ")",
                  " - Models ", "(", n_models, ")",
                  " - Parameters", " (", n_params, ") \n")
  msg_2 <- paste0("    Latitude: ", range_lat[1], " - ", range_lat[2], " - ",
                  "Longitude: ", range_lon[1], " - ", range_lon[2], " \n")
  cat(paste0(msg_1, msg_2))
}

#' Print method for the class `spatialTransform`
#' @param object an object of the class `spatialTransform`.
#' @keywords internal
#' @export
#' @noRd
print.solarTransform <- function(object){
  alpha_ <- format(object$params$alpha, digits = 3, scientific = FALSE)
  beta_ <- format(object$params$beta, digits = 3, scientific = FALSE)
  msg_1 <- paste0("Solar Transform in Range for Xt: ", object$params$Xt_min, " - ", object$params$Xt_max, "\n")
  msg_2 <- paste0("Xt: ",  alpha_, " + ", beta_, " exp(-exp(Yt)) \n")
  msg_3 <- paste0("Yt: log(log(",  beta_, ") - log(Xt-", alpha_, "))")
  cat(msg_1, msg_2, msg_3)
}
