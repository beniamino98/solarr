#' Print method for the class `solarModel`
#'
#' @param object an object of the class \code{\link{solarModel_spec}}.
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


#' Print method for the class `solarOptionPayoff`
#'
#' @param object an object of the class `solarOptionPayoff`.
#'
#' @keywords internal
#' @noRd
#' @export
print.solarOptionPayoff <- function(object){
  msg_1 <- paste0("------------------------ \033[1;35mSolar Option Payoffs\033[0m (", object$payoff_year$side, ") ------------------------ \n")
  msg_2 <- paste0("Yearly payoff: \033[1;31m", format(object$payoff_year$premium, digits = 5), "\033[0m\n")

  # Monthly premiums
  premiums <- object$payoff_month$premium
  # Count the number of integers on the left
  n_integers <- purrr::map_dbl(as.integer(premiums), ~length(stringr::str_split(.x, "")[[1]]))
  n_integers <- ifelse(n_integers == 1, 2, 1)
  # Format the number accordingly
  premiums <- purrr::map2_chr(premiums, n_integers, ~format(.x, digits = 3, nsmall = .y))

  msg_3 <- paste0("Monthly payoffs: \n ",
                  " Jan: \033[1;32m", premiums[1], "\033[0m",
                  "   Feb: \033[1;32m", premiums[2], "\033[0m",
                  "   Mar: \033[1;32m", premiums[3], "\033[0m\n",
                  "  Apr: \033[1;32m", premiums[4], "\033[0m",
                  "   May: \033[1;32m", premiums[5], "\033[0m",
                  "   Jun: \033[1;32m", premiums[6], "\033[0m\n",
                  "  Jul: \033[1;32m", premiums[7], "\033[0m",
                  "   Ago: \033[1;32m", premiums[8], "\033[0m",
                  "   Sep: \033[1;32m", premiums[9], "\033[0m\n",
                  "  Oct: \033[1;32m", premiums[10], "\033[0m",
                  "   Nov: \033[1;32m", premiums[11], "\033[0m",
                  "   Dec: \033[1;32m", premiums[12], "\033[0m\n")
  cat(paste0(msg_1, msg_2, msg_3))
}


#' Print method for the class `solarOptionPayoffs`
#'
#' @param object an object of the class \code{\link{solarOptionPayoffs}}.
#'
#' @keywords internal
#' @noRd
#' @export
print.solarOptionPayoffs <- function(object){

  msg_title <- paste0("-------------- Solar Option Payoffs -------------- \n")
  idx_not_NA <- which(purrr::map_lgl(object$call$scenarios, ~!is.na(.x[1])))
  msg_1 <- paste0("Call scenarios payoffs: (", paste0(names(object$call$scenarios)[idx_not_NA], collapse = ", "), ") \n")
  idx_not_NA <- which(purrr::map_lgl(object$put$scenarios, ~!is.na(.x[1])))
  msg_2 <- paste0("Put scenarios payoffs: (", paste0(names(object$put$scenarios)[idx_not_NA], collapse = ", "), ") \n")
  msg_line <- paste0("------------------------------------------------- \n")
  idx_not_NA <- which(purrr::map_lgl(object$call$model, ~!is.na(.x[1])))
  msg_3 <- paste0("Call model payoffs: (", paste0(names(object$call$model)[idx_not_NA], collapse = ", "), ") \n")
  idx_not_NA <- which(purrr::map_lgl(object$put$model, ~!is.na(.x[1])))
  msg_4 <- paste0("Put model payoffs: (", paste0(names(object$put$model)[idx_not_NA], collapse = ", "), ") \n")

  cat(c(msg_title, msg_1, msg_2, msg_line, msg_3, msg_4))
}


#' Print method for the class `solarOptionBoot`
#'
#' @param object an object of the class  \code{\link{solarOption_historical_bootstrap}}.
#'
#' @keywords internal
#' @noRd
#' @export
print.solarOptionBoot <- function(object){
  msg_1 <- paste0("------------------------ \033[1;35mSolar Option Bootstrap\033[0m (", object$payoff_year$side, ") ------------------------ \n")
  msg_2 <- paste0("Yearly payoff: \033[1;31m", format(object$payoff_year$premium, digits = 5), "\033[0m\n")
  msg_3 <- paste0("Confidence interval ", "(", names(object$payoff_year$premium_dw), ")",
                  "\033[1;31m ", format(object$payoff_year$premium_dw, digits = 5), "\033[0m", "\n",
                  "Confidence interval ", "(", names(object$payoff_year$premium_up), ")",
                  "\033[1;31m ", format(object$payoff_year$premium_up, digits = 5), "\033[0m \n")

  format_premiums <- function(premiums){
    # Count the number of integers on the left
    n_integers <- purrr::map_dbl(as.integer(premiums), ~length(stringr::str_split(.x, "")[[1]]))
    n_integers <- ifelse(n_integers == 1, 2, 1)
    # Format the number accordingly
    premiums <- purrr::map2_chr(premiums, n_integers, ~format(.x, digits = 3, nsmall = .y))
    premiums
  }

  # Monthly premiums
  premiums <- format_premiums(object$payoff_month$premium)
  # Monthly premiums (up)
  premiums_up <- format_premiums(object$payoff_month$premium_up)
  # Monthly premiums (dw)
  premiums_dw <- format_premiums(object$payoff_month$premium_dw)
  msg_line <- paste0(paste0(rep("-", 78), collapse = ""), "\n")
  msg_4 <- paste0("Monthly payoffs: \n ",
                  " Jan: \033[1;32m", premiums[1], "\033[0m", paste0(" (\033[1;31m", premiums_dw[1], "\033[0m", " - ", "\033[1;31m", premiums_up[1], "\033[0m)"),
                  "   Feb: \033[1;32m", premiums[2], "\033[0m", paste0(" (\033[1;31m", premiums_dw[2], "\033[0m", " - ", "\033[1;31m", premiums_up[2], "\033[0m)"),
                  "   Mar: \033[1;32m", premiums[3], "\033[0m", paste0(" (\033[1;31m", premiums_dw[3], "\033[0m", " - ", "\033[1;31m", premiums_up[3], "\033[0m) \n"),
                  "  Apr: \033[1;32m", premiums[4], "\033[0m", paste0(" (\033[1;31m", premiums_dw[4], "\033[0m", " - ", "\033[1;31m", premiums_up[4], "\033[0m)"),
                  "   May: \033[1;32m", premiums[5], "\033[0m", paste0(" (\033[1;31m", premiums_dw[5], "\033[0m", " - ", "\033[1;31m", premiums_up[5], "\033[0m)"),
                  "   Jun: \033[1;32m", premiums[6], "\033[0m", paste0(" (\033[1;31m", premiums_dw[6], "\033[0m", " - ", "\033[1;31m", premiums_up[6], "\033[0m) \n"),
                  "  Jul: \033[1;32m", premiums[7], "\033[0m", paste0(" (\033[1;31m", premiums_dw[7], "\033[0m", " - ", "\033[1;31m", premiums_up[7], "\033[0m)"),
                  "   Ago: \033[1;32m", premiums[8], "\033[0m", paste0(" (\033[1;31m", premiums_dw[8], "\033[0m", " - ", "\033[1;31m", premiums_up[8], "\033[0m)"),
                  "   Sep: \033[1;32m", premiums[9], "\033[0m", paste0(" (\033[1;31m", premiums_dw[9], "\033[0m", " - ", "\033[1;31m", premiums_up[9], "\033[0m) \n"),
                  "  Oct: \033[1;32m", premiums[10], "\033[0m", paste0(" (\033[1;31m", premiums_dw[10], "\033[0m", " - ", "\033[1;31m", premiums_up[10], "\033[0m)"),
                  "   Nov: \033[1;32m", premiums[11], "\033[0m", paste0(" (\033[1;31m", premiums_dw[11], "\033[0m", " - ", "\033[1;31m", premiums_up[11], "\033[0m)"),
                  "   Dec: \033[1;32m", premiums[12], "\033[0m", paste0(" (\033[1;31m", premiums_dw[12], "\033[0m", " - ", "\033[1;31m", premiums_up[12], "\033[0m)"))
  cat(paste0(msg_1, msg_2, msg_3, msg_line, msg_4))
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

