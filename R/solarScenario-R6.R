#' Scenario specification
#'
#' R6 class used to store the data and model components needed to simulate
#' scenarios from a fitted `solarModel` object.
#'
#' @rdname solarScenario_spec
#' @name solarScenario_spec
#' @keywords solarScenario
#' @note Version 1.0.3
#' @export
solarScenario_spec <- R6::R6Class("solarScenario_spec",
                             public = list(
                               #' @field location Tibble containing the name of the reference location and its coordinates.
                               location = dplyr::tibble(),
                               #' @field target Character, target variable.
                               target = "",
                               #' @field quiet Logical. If `TRUE`, suppress selected messages and warnings.
                               quiet = FALSE,
                               #' @field control Named list with date range and simulation setup values.
                               control = list(from = "", to = "", i_start = 1, exclude_known = FALSE),
                               #' @field simulations List of simulated scenario data sets.
                               simulations = list(),
                               #' @description
                               #' Initialize a `solarScenario_spec` object.
                               #' @param model A fitted `solarModel` object.
                               #' @param from Date or character value coercible to `Date`. First date of the scenario period.
                               #' @param to Date or character value coercible to `Date`. Last date of the scenario period.
                               #' @param theta Numeric value stored in the simulation data as `theta`.
                               #' @param exclude_known Logical. If `TRUE`, empirical data stored in the specification are restricted to dates from `from` to `to`.
                               #' @param seed Integer seed stored for scenario simulation.
                               #' @param quiet Logical. If `TRUE`, suppress selected messages and warnings.
                               initialize = function(model, from = "2010-01-01", to = "2010-12-31", theta = 0, exclude_known = FALSE, seed = 1, quiet = FALSE){
                                 # Model's specification
                                 spec <- model$spec
                                 # Reference location
                                 place <- spec$place
                                 # Model's data
                                 data <- model$data
                                 # Maximum number of lags to start
                                 i_start <- max(c(spec$mean.model$order, spec$variance.model$order)) + 1
                                 # Initial date
                                 from <- as.Date(from)
                                 # End date
                                 to <- as.Date(to)
                                 # Selected columns
                                 cols_emp <- c("date", "n", "Year", "Month", "Day", "GHI", "clearsky",
                                               "Xt", "Yt", "Yt_tilde", "Yt_tilde_hat", "eps", "eps_tilde", "sigma", "u_tilde", "B", "z1", "z2")
                                 cols_sim <- c(cols_emp, "Ct", "Yt_bar", "GHI_bar", "sigma_bar", "Yt_tilde_uncond", "sigma_uncond", "mu1", "mu2", "sd1", "sd2", "p1")

                                 # Initialize a dataset
                                 max_date_from <- max(data$date)
                                 max_date_to <- max_date_from - i_start
                                 if (max_date_to >= to) {
                                   df_emp <- dplyr::filter(data, date >= (from - lubridate::days(i_start-1)) & date <= to)
                                   df_emp <- dplyr::bind_cols(place = place, df_emp)
                                 } else if (max_date_to >= from & max_date_from >= from) {
                                   df_emp <- dplyr::filter(data, date >= (from - lubridate::days(i_start-1)))
                                   df_new_emp <- dplyr::tibble(date = seq.Date(max(df_emp$date) + 1, to, by = "1 day"))
                                   df_emp <- dplyr::bind_rows(df_emp, df_new_emp)
                                   df_emp <- dplyr::mutate(df_emp,
                                                           Year = lubridate::year(date),
                                                           Month = lubridate::month(date),
                                                           Day = lubridate::day(date))
                                   df_emp$n <- solarr::number_of_day(df_emp$date)
                                   # Add seasonal variables
                                   df_emp <- df_emp[, cols_emp]
                                   df_emp <- dplyr::left_join(df_emp, model$seasonal_data, by = c("Month", "Day", "n"))
                                 } else {
                                   msg <- paste0("The maximum date for starting a simulation is: ", max_date_from)
                                   if (!quiet) warning(msg)
                                   return(NULL)
                                 }
                                 # Initialize simulation dataset
                                 df_sim <- df_emp[, cols_sim]
                                 # Initialize lambda
                                 df_sim$theta <- theta
                                 # Filter df_emp to be in [from - to] dates
                                 df_emp <- df_emp[, cols_emp]
                                 if (exclude_known) {
                                   df_emp <- dplyr::filter(df_emp, date >= from & date <= to)
                                 }
                                 # Reference location
                                 self$location <- spec$location
                                 # Reference target variable
                                 self$target <- spec$target
                                 # Verbose messages
                                 self$quiet <- quiet
                                 # Control
                                 self$control[["from"]] <- from
                                 self$control[["to"]] <- to
                                 self$control[["exclude_known"]] <- exclude_known
                                 self$control[["i_start"]] <- i_start
                                 # Store data
                                 private[["..data_sim"]] <- df_sim
                                 private[["..data_emp"]] <- df_emp
                                 # Store only required parameters
                                 # ARMA model
                                 ARMA <- spec$mean.model
                                 private[["..ARMA"]] <- list(arOrder = max(c(ARMA$arOrder, 1)),
                                                             maOrder = max(c(ARMA$maOrder, 1)),
                                                             intercept = ARMA$intercept,
                                                             phi = ARMA$phi,
                                                             theta = ARMA$theta)
                                 # GARCH model
                                 GARCH <- spec$variance.model
                                 private[["..GARCH"]] <- list(archOrder = max(c(1, GARCH$archOrder)),
                                                              garchOrder = max(c(1, GARCH$garchOrder)),
                                                              omega = GARCH$omega,
                                                              alpha = GARCH$alpha,
                                                              beta = GARCH$beta)
                                 # Transform function
                                 private[["..transform"]] <- spec$transform
                                 # Marginal probabilities
                                 private[["..prob"]] <- spec$mixture.model$prob$clone(TRUE)
                               },
                               #' @description
                               #' Marginal mixture probabilities.
                               #' @param nmonth Integer vector of month numbers.
                               #' @return A numeric vector or object returned by the stored probability model's `predict()` method.
                               prob = function(nmonth){
                                 private[["..prob"]]$predict(nmonth)
                               },
                               #' @description
                               #' Print a summary of the scenario specification.
                               #' @return Invisibly returns `NULL`.
                               print = function(){
                                 cat("---------------- solarScenario Specification ----------------  \n")
                                 cat(paste0("Place: ", self$place, "\n"))
                                 cat(paste0("From - To: ", self$control[["from"]], " - ", self$control[["to"]], "\n"))
                                 cat(paste0("Number of simulations: ", self$nsim, "\n"))
                                 cat(paste0("Version: ", private$version, "\n"))
                               }
                             ),
                             private = list(
                               version = "1.0.2",
                               ..data_sim = NA,
                               ..data_emp = NA,
                               ..ARMA = NA,
                               ..GARCH = NA,
                               ..transform = NA,
                               ..prob = NA
                             ),
                             active = list(
                               #' @field place Character value with the reference place.
                               place = function(){
                                 self$location$place
                               },
                               #' @field ARMA List of stored ARMA model parameters used for simulation.
                               ARMA = function(){
                                 private$..ARMA
                               },
                               #' @field GARCH List of stored GARCH model parameters used for simulation.
                               GARCH = function(){
                                 private$..GARCH
                               },
                               #' @field transform Transformation object copied from the model specification.
                               transform = function(){
                                 private$..transform
                               },
                               #' @field data_sim Tibble containing the data used as the simulation template.
                               data_sim = function(){
                                 private$..data_sim
                               },
                               #' @field emp Tibble containing empirical data stored in the specification.
                               emp = function(){
                                 private$..data_emp
                               },
                               #' @field nsim Integer number of stored simulations.
                               nsim = function(){
                                 length(self$simulations)
                               },
                               #' @field scenarios Nested tibble built from stored simulations.
                               scenarios = function() {
                                 simulations <- self$simulations
                                 if (purrr::is_empty(simulations)) {
                                   return(simulations)
                                   stop("The slot `self$simulations` is empty! Consider running `simSpec <- solarScenario_filter(scenario)` before!")
                                 }
                                 dplyr::bind_rows(simulations) %>%
                                   dplyr::group_by(date, Year, Month, Day) %>%
                                   tidyr::nest() %>%
                                   dplyr::ungroup()
                               }
                             ))

#' Solar scenario simulation
#'
#' R6 class used to simulate residuals and store filtered scenario paths from a
#' `solarScenario_spec` object.
#'
#' @rdname solarScenario
#' @name solarScenario
#' @keywords solarScenario
#' @note Version 1.0.3
#' @export
solarScenario <- R6::R6Class("solarScenario",
                             public = list(
                               #' @field seed Integer seed used when simulating residuals.
                               seed = 1,
                               #' @field residuals Tibble or list-like object containing simulated residuals waiting to be filtered.
                               residuals = list(),
                               #' @description
                               #' Initialize a `solarScenario` object.
                               #' @param simSpec A `solarScenario_spec` object.
                               #' @param seed Integer seed used for residual simulation.
                               initialize = function(simSpec, seed = 1){
                                 # Random seed
                                 self$seed <- seed
                                 # Model specification
                                 private[["..spec"]] <- simSpec$clone(TRUE)
                               },
                               #' @description
                               #' Simulate residuals for scenario generation.
                               #' @param nsim Integer number of residual paths to simulate.
                               #' @return Updates the `residuals` field.
                               simulate_residuals = function(nsim = 1){
                                 # First time simulation
                                 if (purrr::is_empty(self$residuals)) {
                                   self$residuals <- solarScenario_residuals(self$spec, nsim = nsim, seed = self$seed)
                                   self$residuals$filter <- TRUE
                                 } else {
                                   self$seed <- self$seed + 1
                                   new_residuals <- solarScenario_residuals(self$spec, nsim = nsim, seed = self$seed)
                                   new_residuals$filter <- TRUE
                                   self$residuals <- dplyr::bind_rows(self$residuals, new_residuals)
                                 }
                               },
                               #' @description
                               #' Filter simulated residuals into scenario paths.
                               #' @param all Logical. If `TRUE`, filter all stored residuals; otherwise, filter residuals marked for filtering.
                               #' @return Appends filtered scenario paths to `spec$simulations` and updates residual flags.
                               filter = function(all = FALSE){
                                 if (all) {
                                   index_filter <- 1:nrow(self$residuals)
                                 } else {
                                   index_filter <- which(self$residuals$filter)
                                 }
                                 if (!purrr::is_empty(index_filter)) {
                                   residuals <- self$residuals[index_filter,]
                                   simulations <- solarScenario_filter(self$spec, residuals)
                                   private$..spec$simulations <- append(private$..spec$simulations, simulations)
                                   self$residuals$filter <- FALSE
                                 } else {
                                   cli::cli_alert_warning("All the residuals have been already filtered!
                                                           Use `all = TRUE` to filter them again or simulate new residuals with $simulate_residuals()!")
                                 }
                               }
                             ),
                             private = list(
                               version = "1.0.2",
                               ..spec = NA
                             ),
                             active = list(
                               #' @field spec Stored `solarScenario_spec` object.
                               spec = function(){
                                 private$..spec
                               }
                             ))






