#' R6 Class implementation for `solarScenario_spec`
#'
#' Specification of scenarios of a `solarModel`
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
                               #' @field quiet Logical
                               quiet = FALSE,
                               #' @field control Named list
                               control = list(from = "", to = "", i_start = 1, exclude_known = FALSE),
                               #' @field simulations List with simulations
                               simulations = list(),
                               #' @description
                               #' Initialize a `solarScenario_spec`
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
                               #' Marginal mixture probabilities
                               #' @param nmonth Integer, number of month.
                               prob = function(nmonth){
                                 private[["..prob"]]$predict(nmonth)
                               },
                               #' @description
                               #' Print method for the class `solarScenario_spec`
                               print = function(){
                                 cat("---------------- solarScenario Specification ----------------  \n")
                                 cat(paste0("Place: ", self$place, "\n"))
                                 cat(paste0("From - To: ", from, " - ", to, "\n"))
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
                               place = function(){
                                 self$location$place
                               },
                               ARMA = function(){
                                 private$..ARMA
                               },
                               GARCH = function(){
                                 private$..GARCH
                               },
                               transform = function(){
                                 private$..transform
                               },
                               data_sim = function(){
                                 private$..data_sim
                               },
                               emp = function(){
                                 private$..data_emp
                               },
                               nsim = function(){
                                 length(self$simulations)
                               },
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

#' R6 Class implementation for `solarScenario`
#'
#' Specification of scenarios of a `solarModel`
#'
#' @rdname solarScenario
#' @name solarScenario
#' @keywords solarScenario
#' @note Version 1.0.3
#' @export
solarScenario <- R6::R6Class("solarScenario",
                             public = list(
                               seed = 1,
                               residuals = list(),
                               initialize = function(simSpec, seed = 1){
                                 # Random seed
                                 self$seed <- seed
                                 # Model specification
                                 private[["..spec"]] <- simSpec$clone(TRUE)
                               },
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
                               spec = function(){
                                 private$..spec
                               }
                             ))







