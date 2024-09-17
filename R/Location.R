#' Generate a location
#'
#' @rdname Location
#' @name Location
#' @keywords internal
Location <- function(place, nsim = 50, by = "1 month",
                     exact_daily_premium = FALSE,
                     measures = c("Q", "Qdw", "Qup"),
                     control_model = control_solarModel(),
                     control_options = control_solarOption(),
                     control_esscher = control_solarEsscher(),
                     seed = 1){
  # New environment
  object <- new.env(parent = .GlobalEnv)
  # Add place
  object$place <- place
  # Number of simulations
  object$nsim <- nsim
  # Structure cams data
  spec <- solarModel_spec(object$place, from = as.Date("2005-01-01"), to = as.Date("2023-11-08"),
                          CAMS_data = solarr::CAMS_data, control_model = control_model)

  # 1) fit a solar model
  object$model <- solarModel(spec)
  # Calibrate Esscher bounds and optimal parameters
  # Esscher parametric function, bounds and optimal theta for change of measure
  object$model <- solarEsscher_bounds(model = object$model, control_esscher = control_esscher, control_options = control_options)

  # 2) Simulate scenarios under different measures
  # Simulate P-scenarios
  object$model$scenarios$P <- solarModel_scenario(object$model, by = by, nsim = nsim, lambda = 0, seed = seed,
                                                  from = control_options$from, to = control_options$to)
  # Aggregate P-Payoff
  object$model$payoffs$sim$P <- solarOption_scenario(sim = object$model$scenarios$P$sim, control_options = control_options)

  for(measure in measures){
    # From P to Q measure (simulations)
    object$model$scenarios[[measure]] <- solarModel_scenario(object$model, by = by, nsim = nsim,
                                                    lambda = object$model$esscher$params[[measure]], seed = seed,
                                                    from = control_options$from, to = control_options$to)
    # Aggregate payoff (simulations)
    object$model$payoffs$sim[[measure]] <- solarOption_scenario(sim = object$model$scenarios[[measure]]$sim, control_options = control_options)
  }

  # Add options control
  object$model$payoffs$control <- control_options
  # Structure the payoff
  if (length(measures) == 3) {
    object$model <- solarOption_structure(object$model, type = "sim", exact_daily_premium = exact_daily_premium)
  }
  object$model <- solarOption_structure(object$model, type = "model", exact_daily_premium = exact_daily_premium)

  return(object)
}
