#' Print method for the class `solarModel`
#'
#' @param object an object of the class \code{\link{solarModel_spec}} or \code{\link{solarModel}}.
#'
#' @keywords internal
#' @noRd
#' @export
print.solarModelSpec <- function(object){
  model_type <- class(object)[1]
  # Complete data specifications
  data <- object$dates$data
  # Train data specifications
  train <- object$dates$train
  train$perc <- format(train$perc*100, digits = 4)
  # Test data specifications
  test <- object$dates$test
  test$perc <- format(test$perc*100, digits = 4)
  msg_0 <- paste0("--------------------- ", model_type, " (", object$place, ") ", "--------------------- \n")
  msg_1 <- paste0("Target: ", object$target, " \n Lat: ", object$coords$lat, "\n Lon: ", object$coords$lon, "\n Alt: ", object$coords$alt, " \n")
  msg_1 <- paste0("Target: ", object$target, " \n Coordinates: (Lat: ", object$coords$lat, ", Lon: ", object$coords$lon, ", Alt: ", object$coords$alt, ") \n")
  msg_2 <- paste0(" Dates: ", data$from, " - ", data$to, "\n Observations: ", data$nobs, "\n")
  msg_3 <- paste0("---------------------------------------------------------------\n")
  msg_4 <- paste0("Train dates: ", train$from, " - ", train$to, " (", train$nobs, " points ~ ", train$perc, "%)", "\n")
  msg_5 <- paste0("Test  dates: ", test$from, " - ", test$to, " (", test$nobs, " points ~ ", test$perc, "%)", "\n")
  msg_6 <- paste0("---------------------------------------------------------------\n")
  msg_7 <- paste0("Likelihood: ", format(object$loglik, digits = 8), "\n")
  cat(paste0(msg_0, msg_1, msg_2, msg_3, msg_4, msg_5, msg_6, msg_7))
}

#' Print method for the class `solarModel`
#'
#' @param object an object of the class \code{\link{solarModel}}.
#'
#' @keywords internal
#' @noRd
#' @export
print.solarModel <- function(object){
  NextMethod(object)
}

#' Print method for the class `seasonalModel`
#'
#' @param object an object of the class  \code{\link{seasonalModel}}.
#'
#' @keywords internal
#' @noRd
#' @export
print.seasonalModel <- function(object){
  cat(paste0("----------------------- ", class(object)[1], " ----------------------- \n"))
  msg_1 <- paste0(" - Order: ", object$order, "\n - Period: ", object$period, "\n")
  n_external_reg <- ncol(object$seasonal_data)-1
  if (n_external_reg == 0) {
    msg_2 <- paste0("- External regressors: ", n_external_reg, "\n")
  } else {
    msg_2 <- paste0("- External regressors: ", n_external_reg, " (", names(object$seasonal_data)[-1], ")\n")
  }
  cat(c(msg_1, msg_2))
  cat(paste0("--------------------------------------------------------------\n"))
  print(object$model)
}


#' Print method for the class `solarOption`
#'
#' @param object an object of the class  \code{\link{solarOption}}.
#'
#' @keywords internal
#' @noRd
#' @export
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
#'
#' @param object an object of the class \code{\link{spatialModel}}.
#'
#' @keywords internal
#' @noRd
#' @export
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

#' Print method for the class `solarTransform`
#'
#' @param object an object of the class \code{\link{solarTransform}}.
#'
#' @keywords internal
#' @noRd
#' @export
print.solarTransform <- function(object){
  alpha_ <- format(object$params$alpha, digits = 3, scientific = FALSE)
  beta_ <- format(object$params$beta, digits = 3, scientific = FALSE)
  msg_1 <- paste0("Solar Transform in Range for Xt: ", object$params$Xt_min, " - ", object$params$Xt_max, "\n")
  msg_2 <- paste0("Xt: ",  alpha_, " + ", beta_, " exp(-exp(Yt)) \n")
  msg_3 <- paste0("Yt: log(log(",  beta_, ") - log(Xt-", alpha_, "))")
  cat(msg_1, msg_2, msg_3)
}

#' Print method for the class `spatialScenarioSpec`
#'
#' @param object an object of the class \code{\link{spatialScenario_spec}}.
#'
#' @keywords internal
#' @noRd
#' @export
print.spatialScenarioSpec <- function(object){
  cat("------------------------- spatialScenarioSpec -------------------------\n")
  msg_1 <- paste0(" Number of models: ", length(object$spec), "\n")
  msg_2 <- paste0("Dates: ", object$from, " - ", object$to, "\n")
  msg_3 <- paste0("Number of simulations: ", object$nsim, "\n",
                  " - Residuals: (", object$residuals, ")", "\n",
                  " - Filter: (", object$filter, ")")
  cat(msg_1, msg_2, msg_3)
}

#' Print method for the class `solarScenarioSpec`
#'
#' @param object an object of the class \code{\link{solarScenario_spec}}.
#'
#' @keywords internal
#' @noRd
#' @export
print.solarScenarioSpec <- function(object){
  msg_0 <- paste0("--------------------- ", "solarScenarioSpec", " (", object$place, ") ", "--------------------- \n")
  msg_1 <- paste0("Target: ", object$target, " \n Lat: ", object$coords$lat, "\n Lon: ", object$coords$lon, "\n Alt: ", object$coords$alt, " \n")
  msg_2 <- paste0("---------------------------------------------------------------\n")
  msg_3 <- paste0(" Dates: ", object$from, " - ", object$to, " Observations: ", nrow(object$sim), "\n")
  msg_4 <- paste0("Number of simulations: ", object$nsim, "\n",
                  " - Residuals: (", !purrr::is_empty(object$residuals), ")", "\n",
                  " - Filter: (", !purrr::is_empty(object$simulations), ")")
  cat(paste0(msg_0, msg_1, msg_2, msg_3, msg_4))
}


