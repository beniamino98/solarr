#' Control function for a `solarModel` object
#'
#' @description
#' Control function for a `solarModel` object that contains all the setups used for the estimation.
#'
#' @examples
#' control <- solarModel_spec$new()
#' @rdname solarModel_spec
#' @name solarModel_spec
#' @keywords solarModel
#' @note Version 1.0.3
#' @export
solarModel_spec <- R6::R6Class("solarModel_spec",
                               public = list(
                                 #' @field outliers list with outliers.
                                 outliers = NULL,
                                 #' @field place Character, optional name of the location considered.
                                 place = "",
                                 #' @field coords A named list with the name and the coordinates of the location considered. Contains:
                                 #' \describe{
                                 #'  \item{lat}{Numeric, reference latitude in degrees.}
                                 #'  \item{lon}{Numeric, reference longitude in degrees.}
                                 #'  \item{alt}{Numeric, reference altitude in metres.}
                                 #'}
                                 coords = dplyr::tibble(lat = NA, lon = NA, alt = NA),
                                 #' @description
                                 #' Initialize a `solarModel_spec` object.
                                 initialize = function(){
                                   self$set_params()
                                   self$set_transform()
                                   self$set_clearsky()
                                   self$set_seasonal.mean()
                                   self$set_mean.model()
                                   self$set_seasonal.variance()
                                   self$set_variance.model()
                                   self$set_mixture.model()
                                 },
                                 #' @description
                                 #' Generic controls
                                 #' @param stochastic_clearsky Logical, when `TRUE` the clear sky will be considered stochastic.
                                 #' @param clearsky_threshold Numeric, parameter > 1, used to scale up CAMS clearsky to avoid that clear sky radiaion and global horizontal radiation are equal.
                                 #' @param quiet Logical. When `TRUE` the function will not display any message. The dafault if `TRUE`.
                                 set_params = function(stochastic_clearsky = FALSE, clearsky_threshold = 1.01, quiet = FALSE){
                                   private$..stochastic_clearsky <- stochastic_clearsky
                                   private$..clearsky_threshold <- clearsky_threshold
                                   private$..quiet <- quiet
                                 },
                                 #' @description
                                 #' Control parameters for the `solarTransform`. See \code{\link{solarTransform}} for more details.
                                 #' @param min_pos Integer, position of the minimum. For example when `2` the minimum is the second lowest value.
                                 #' @param max_pos Integer, position of the maximum. For example when `3` the maximum is the third greatest value.
                                 #' @param delta transform params
                                 #' @param link Character, link function.
                                 #' @param threshold Numeric. Threshold used to estimate the transformation parameters \deqn{\alpha} and \deqn{\beta}.
                                 #' The default is `0.01`. See \code{\link{solarTransform}} for more details.
                                 set_transform = function(min_pos = 1, max_pos = 1, link = "invgumbel", delta = 0.05, threshold = 0.01){
                                   # Initialize transform object
                                   private$..transform <- solarTransform$new(alpha = 0, beta = 1, link = link)
                                   # Add extra control parameters
                                   private$..transform$control[["min_pos"]] <- min_pos
                                   private$..transform$control[["max_pos"]] <- max_pos
                                   private$..transform$control[["delta"]] <- delta
                                   private$..transform$control[["threshold"]] <- threshold
                                 },
                                 #' @description
                                 #' List with specification's parameters of the clear sky model.
                                 #' @param spec Named list, model's specification. See the function \code{\link{seasonalClearsky_spec}} for more details.
                                 #' @param control Named list, control parameters. See the function \code{\link{control_seasonalClearsky}} for more details.
                                 set_clearsky = function(spec = seasonalClearsky_spec(), control = control_seasonalClearsky()){
                                   private$..seasonal_model_Ct <- seasonalClearsky$new(spec = spec, control = control)
                                 },
                                 #' @description
                                 #' List with specification's parameters of the seasonal mean \eqn{\bar{Y}_t} for \eqn{Y_t}.
                                 #' @param order Integer. Specify the order of the seasonal mean \deqn{\bar{Y}_t}. The default is `1`.
                                 #' @param period Integer, seasonal periodicity, the default is `365`.
                                 #' @param include.trend Logical. When `TRUE` an yearly trend \deqn{t} will be included in the seasonal model, otherwise will be excluded. The default is `FALSE`.
                                 #' @param include.intercept Logical. When `TRUE` the intercept \deqn{a_0} will be included in the seasonal model, otherwise will be excluded. The default is `TRUE`.
                                 #' @param monthly.mean Logical. When `TRUE` a vector of 12 monthly means will be computed on the deseasonalized series \deqn{\tilde{Y}_t = Y_t - \bar{Y}_t}
                                 #'  and it is subtracted to ensure that the time series is centered around zero for all the months. The dafault if `TRUE`.
                                 set_seasonal.mean = function(order = 1, period = 365, include.trend = FALSE, include.intercept = TRUE, monthly.mean = FALSE){
                                   # 1) Initialization of the seasonal model
                                   base_formula <- ifelse(include.intercept, "Yt ~ 1", "Yt ~ -1")
                                   formula <- ifelse(include.trend, paste0(base_formula, " + n"), base_formula)
                                   private$..seasonal.mean <- seasonalModel$new(formula = formula, order = order, period = period)
                                   # **************************************************** #
                                   # 2) Standardize the parameters names
                                   new_coefs_names <- private$..seasonal.mean$coefs_names
                                   if(include.intercept) {
                                     new_coefs_names[1] <- "a_0"
                                   }
                                   idx_sincos <- stringr::str_detect(new_coefs_names, "sin|cos")
                                   if(!purrr::is_empty(idx_sincos)) {
                                     new_coefs_names[idx_sincos] <- paste0("a_", new_coefs_names[idx_sincos])
                                   }
                                   # Update parameters name inside the R6 object
                                   private$..seasonal.mean$update_coefs_names(new_coefs_names)
                                   # **************************************************** #
                                   # 3) Store extra control parameters
                                   private$..seasonal.mean$control[["include.intercept"]] <- include.intercept
                                   private$..seasonal.mean$control[["include.trend"]] <- include.trend
                                   private$..seasonal.mean$control[["monthly.mean"]] <- monthly.mean
                                 },
                                 #' @description
                                 #' List with specification's parameters of the ARMA model for deseasonalized series \eqn{\tilde{Y}_t = Y_t - \bar{Y}_t}.
                                 #' @param arOrder Integer. An integer specifying the order of the AR component. The default is `1`.
                                 #' @param maOrder Integer. An integer specifying the order of the MA component. The default is `0`.
                                 #' @param include.intercept Logical. When `TRUE` the intercept \deqn{\phi_0} will be included in the seasonal model, otherwise will be excluded. The default is `FALSE`.
                                 set_mean.model = function(arOrder = 1, maOrder = 0, include.intercept = FALSE){
                                   private$..mean.model = ARMA_modelR6$new(arOrder = arOrder, maOrder = maOrder, include.intercept = include.intercept)
                                 },
                                 #' @description
                                 #' List with specification's parameters of the seasonal variance \eqn{\bar{\sigma}_t} for ARMA's residuals \eqn{e_t}
                                 #' @param order Integer. Specify the order of the seasonality of the seasonal variance. The default is `1`.
                                 #' @param period Integer, seasonal periodicity, the default is `365`.
                                 #' @param include.trend Logical. When `TRUE` an yearly trend \deqn{t} will be included in the seasonal model, otherwise will be excluded. The default is `FALSE`.
                                 #' @param correction Logical. When `TRUE` the parameters of seasonal variance are corrected to ensure
                                 #'  that the standardize the residuals have exactly a unitary variance. The dafault if `TRUE`.
                                 #' @param monthly.mean Logical. When `TRUE` a vector of 12 monthly std. deviations will be computed
                                 #'  on the standardized residuals  \deqn{\tilde{\varepsilon}_t} and used to standardize the time series
                                 #'  such that it has unitary variance for all the months. The default if `TRUE`.
                                 set_seasonal.variance = function(order = 1, period = 365, include.trend = FALSE, correction = FALSE, monthly.mean = FALSE){
                                   # 1) Initialization of the seasonal model
                                   formula <- ifelse(include.trend, "eps2 ~ 1 + n", "eps2 ~ 1")
                                   # Initialize the seasonal model
                                   private$..seasonal.variance <- seasonalModel$new(formula = formula, order = order, period = period)
                                   # **************************************************** #
                                   # 2) Standardize the parameters names
                                   new_coefs_names <- private$..seasonal.variance$coefs_names
                                   new_coefs_names["intercept"] <- "c_0"
                                   # Add delta on standard coefficients
                                   idx_sincos <- stringr::str_detect(new_coefs_names, "sin|cos")
                                   new_coefs_names[idx_sincos] <- paste0("c_", new_coefs_names[idx_sincos])
                                   # Update parameters name inside the R6 object
                                   private$..seasonal.variance$update_coefs_names(new_coefs_names)
                                   # **************************************************** #
                                   # 3) Store extra control parameters
                                   private$..seasonal.variance$control[["include.intercept"]] <- TRUE
                                   private$..seasonal.variance$control[["include.trend"]] <- include.trend
                                   private$..seasonal.variance$control[["correction"]] <- correction
                                   private$..seasonal.variance$control[["sigma20"]] <- 1
                                   private$..seasonal.variance$control[["monthly.mean"]] <- monthly.mean
                                 },
                                 #' @description
                                 #' List with specification's parameters of the GARCH variance \eqn{\sigma_t} for deseasonalized residuals \eqn{\tilde{e}_t = e_t/\bar{\sigma}_t}.
                                 #' @param archOrder Integer. An integer specifying the order of the ARCH component. The default is `1`.
                                 #' @param garchOrder Integer. An integer specifying the order of the GARCH component. The default is `1`.
                                 set_variance.model = function(archOrder = 1, garchOrder = 1){
                                   # 2) Store extra control parameters
                                   if (archOrder == 0 & garchOrder == 0){
                                     # 1) Initialize a GARCH model
                                     private$..variance.model <- sGARCH$new(archOrder = 1, garchOrder = 1, mode = "unitOmega")
                                     private$..variance.model$.__enclos_env__$private$..archOrder <- 0
                                     private$..variance.model$.__enclos_env__$private$..garchOrder <- 0
                                     # Set garch_variance = FALSE
                                     private$..variance.model$control[["garch_variance"]] <- FALSE
                                   } else {
                                     # 1) Initialize a GARCH model
                                     private$..variance.model <- sGARCH$new(archOrder = archOrder, garchOrder = garchOrder, mode = "unitOmega")
                                     # Set garch_variance = TRUE
                                     private$..variance.model$control[["garch_variance"]] <- TRUE
                                   }
                                 },
                                 #' @description
                                 #' List with specification's parameters of the Gaussian mixture model for GARCH residuals \eqn{u_t = \tilde{e}_t/\sigma_t}.
                                 #' @param abstol Numeric. Absolute level for convergence of the EM-algorithm. The default is `1e-20`.
                                 #' @param match.expectation Logical, when `TRUE` the mixture parameters ensures that the expected value is matched.
                                 #' @param match.variance Logical, when `TRUE` the mixture parameters ensures that the variance is matched.
                                 #' @param match.empiric Logical, when `TRUE` and `match.expectation = TRUE` or  `match.variance = TRUE` the mixture parameters
                                 #' will be estimated ensuring that mean and variance matches the empirical parameters. Otherwise if `FALSE` and
                                 #'  `match.expectation = TRUE` or `match.variance = TRUE` the target expectation will be zero and the target variance 1.
                                 #' @param method Character, package used to fit the parameters. Can be `mclust` or `mixtools`.
                                 #' @param maxit Integer. Maximum number of iterations for EM-algorithm. The default is `5000`.
                                 #' @param maxrestarts Integer. Maximum number of restarts when EM-algorithm does not converge. The default is `500`.
                                 set_mixture.model = function(abstol = 1e-20, match.expectation = FALSE, match.variance = FALSE,
                                                              match.empiric = FALSE, method = "mclust", maxit = 5000, maxrestarts = 500){
                                   # 1) Initialize a solarMixture
                                   private$..mixture.model <- solarMixture$new(components = 2,
                                                                               method = method,
                                                                               maxit = maxit,
                                                                               maxrestarts = maxrestarts,
                                                                               abstol = abstol)
                                   # 2) Store extra control parameters
                                   private$..mixture.model$control[["match.expectation"]] <- match.expectation
                                   private$..mixture.model$control[["match.variance"]] <- match.variance
                                   private$..mixture.model$control[["match.empiric"]] <- match.empiric
                                 },
                                 #' @description
                                 #' Specification function for a `solarModel`
                                 #' @param place Character, name of an element in the `CAMS_data` list.
                                 #' @param target Character, target variable to model. Can be `GHI` or `clearsky`.
                                 #' @param min_date Character. Date in the format `YYYY-MM-DD`. Minimum date for the complete data. If `missing` will be used the minimum data available.
                                 #' @param max_date Character. Date in the format `YYYY-MM-DD`. Maximum date for the complete data. If `missing` will be used the maximum data available.
                                 #' @param from Character. Date in the format `YYYY-MM-DD`. Starting date to use for training data.
                                 #' If `missing` will be used the minimum data available after filtering for `min_date`.
                                 #' @param to character. Date in the format `YYYY-MM-DD`. Ending date to use for training data.
                                 #' If `missing` will be used the maximum data available after filtering for `max_date`.
                                 #' @param data data for the selected location.
                                 specification = function(place, target = "GHI", min_date, max_date, from, to, data){
                                   # Match the target variable to model
                                   target <- match.arg(target, choices = c("GHI", "clearsky"))
                                   # Extract CAMS data for the selected location
                                   if (missing(data)){
                                     # Match a location in the dataset
                                     place <- match.arg(place, choices = names(CAMS_data), several.ok = FALSE)
                                     data <- solarr::CAMS_data[[place]]
                                   }
                                   # Minimum date for the complete data
                                   if (missing(min_date) || is.null(min_date) || is.na(min_date)) {
                                     min_date <- min(data$date, na.rm = TRUE)
                                   } else {
                                     min_date <- as.Date(min_date)
                                   }
                                   # Maximum date for the complete data
                                   if (missing(max_date) || is.null(max_date) || is.na(max_date)) {
                                     max_date <- max(data$date, na.rm = TRUE)
                                   } else {
                                     max_date <- as.Date(max_date)
                                   }
                                   # Minimum date for train data
                                   if (missing(from) || is.null(from) || is.na(from)) {
                                     from <- min(data$date, na.rm = TRUE)
                                   } else {
                                     from <- as.Date(from)
                                   }
                                   # Maximum date for train data
                                   if (missing(to) || is.null(to) || is.na(to)) {
                                     to <- max(data$date, na.rm = TRUE)
                                   } else {
                                     to <- as.Date(to)
                                   }
                                   # Filter for min and max dates the complete dataset
                                   data <- dplyr::filter(data, date >= min_date & date <= max_date)
                                   # Increase clearsky to avoid NaN
                                   data$clearsky <- data$clearsky * self$clearsky_threshold
                                   # It may happen in CAMS data that clear sky value is just a little bit greater than GHI (~10-3)
                                   # Therefore, before using such time series it is convenient to impute GHI value such that they
                                   # became equal to the given CAMS clear sky
                                   # Detect and impute outliers
                                   outliers <- clearsky_outliers(data$GHI, data$clearsky, date = data$date, threshold = 0, quiet = self$quiet)
                                   # Update the time series of GHI with adjusted values
                                   data$GHI <- outliers$x
                                   # Label for data used for estimation
                                   data <- dplyr::mutate(data,
                                                         isTrain = ifelse(date >= from & date <= to, TRUE, FALSE),
                                                         weights = ifelse(isTrain, 1, 0))
                                   # Add the normalized weights
                                   data$weights <-  data$weights / sum(data$weights)
                                   # Train observations and percentage
                                   nobs_train <- length(data$isTrain[data$isTrain])
                                   # Compute percentage of train obs. on total obs.
                                   perc_train <- nobs_train / nrow(data)
                                   # Train observations and percentage
                                   nobs_test <- length(data$isTrain[!data$isTrain])
                                   # Compute percentage of test obs. on total obs.
                                   perc_test <- nobs_test / nrow(data)
                                   # Model dates
                                   model_dates = list(data = list(from = min_date, to = max_date, nobs = nrow(data), perc = 1),
                                                      train = list(from = from, to = to, nobs = nobs_train, perc = perc_train),
                                                      test = list(from = to, to = max_date, nobs = nobs_test, perc = perc_test))
                                   # Store the data
                                   self$place = attr(data, "place")
                                   self$coords = attr(data, "coords")
                                   private$..data = data
                                   private$..dates = model_dates
                                   private$..target = target
                                   # Update latitude for clear sky
                                   private$..seasonal_model_Ct$lat <- self$coords$lat
                                 },
                                 #' @description
                                 #' Update the parameters inside object
                                 #' @param params List of parameters. See the slot `$coefficients` for a template.
                                 update = function(params){
                                   if (!missing(params)) {
                                     if (!is.null(dim(params))) {
                                       params <- solarModel_match_params(params, self$coefficients)
                                     }
                                     # Update transform parameters
                                     private$..transform$update(params$params$alpha, params$params$beta)
                                     # Update clear sky model
                                     private$..seasonal_model_Ct$update(unlist(params$seasonal_model_Ct))
                                     # Update seasonal mean model
                                     private$..seasonal.mean$update(unlist(params$seasonal_model_Yt))
                                     # Update AR mean model
                                     private$..mean.model$update(unlist(params$ARMA))
                                     # Update seasonal variance model
                                     private$..seasonal.variance$update(unlist(params$seasonal_variance))
                                     # Update GARCH variance model
                                     private$..variance.model$update(unlist(params$GARCH))
                                     # ***************** Update Gaussian Mixture model *****************
                                     means <- cbind(mu1 = unlist(params$NM_mu_up), mu2 = unlist(params$NM_mu_dw))
                                     sd <- cbind(sd1 = unlist(params$NM_sd_up), sd2 = unlist(params$NM_sd_dw))
                                     p <- cbind(p1 = unlist(params$NM_p_up), p2 = 1 - unlist(params$NM_p_up))
                                     private$..mixture.model$update(means = means, sd = sd, p = p)
                                   } else {
                                     message("`params` is missing nothing to update!")
                                   }
                                 },
                                 #' @description
                                 #' Print method for `solarModel_spec` class.
                                 print = function(){
                                   # Seasonal mean order
                                   sm_order <- self$seasonal.mean$order
                                   # ARMA order
                                   arOrder <- self$mean.model$arOrder
                                   maOrder <- self$mean.model$maOrder
                                   # Seasonal variance order
                                   sv_order <- self$seasonal.variance$order
                                   # GARCH order
                                   archOrder <- self$variance.model$archOrder
                                   garchOrder <- self$variance.model$garchOrder
                                   green <- function(x) paste0("\033[1;32m", x, "\033[0m")
                                   red <- function(x) paste0("\033[1;31m", x, "\033[0m")
                                   msg_col <- function(x) ifelse(x, green(x), red(x))

                                   if (!is.na(self$place)) {
                                     # Complete data specifications
                                     data <- self$dates$data
                                     # Train data specifications
                                     train <- self$dates$train
                                     train$perc <- format(train$perc*100, digits = 4)
                                     # Test data specifications
                                     test <- self$dates$test
                                     test$perc <- format(test$perc*100, digits = 4)
                                     msg_0 <- paste0("--------------------- ", "solarModel", " (\033[4;35m", self$place, "\033[0m) ", "--------------------- \n")
                                     msg_1 <- paste0("\033[1;35m Target \033[0m: ", self$target, " \n\033[1;35m Coordinates\033[0m: ",
                                                     "(\033[1;35mLat\033[0m: ", self$coords$lat,
                                                     ", \033[1;35mLon\033[0m: ", self$coords$lon,
                                                     ", \033[1;35mAlt\033[0m: ", self$coords$alt, ") \n")
                                     msg_2 <- paste0("\033[1;35m Observations\033[0m: ", data$nobs, "\n")
                                     msg_3 <- paste0("---------------------------------------------------------------\n")
                                     msg_4 <- paste0("\033[1;34m Dates\033[0m: ", data$from, " - ", data$to, "\n")
                                     msg_5 <- paste0("  - \033[1;34mTrain\033[0m: ", train$from, " - ", train$to, " (", train$nobs, " points ~ ", train$perc, "%)", "\n")
                                     msg_6 <- paste0("  - \033[1;34mTest\033[0m: ", test$from, " - ", test$to, " (", test$nobs, " points ~ ", test$perc, "%)", "\n")
                                     cat(paste0(msg_0, msg_1, msg_2, msg_3, msg_4, msg_5, msg_6))
                                   }
                                   msg_0 <- "------------------------\033[4;35m Specification \033[0m------------------------ \n"
                                   msg_1 <- paste0(" - Stochastic clearsky: ", msg_col(self$stochastic_clearsky), "\n")
                                   cat(c(msg_0, paste0(self$model_name, "\n"), msg_1))
                                   msg_0 <- "-------------------------\033[4;35m Mean Models \033[0m------------------------- \n"
                                   msg_1 <- paste0(" - Trend: ", msg_col(self$seasonal.mean$control$include.trend), "\n")
                                   msg_2 <- paste0(" - Intercept: ", msg_col(self$seasonal.mean$control$include.intercept), "\n")
                                   msg_3 <- paste0(" - Monthly correction: ", msg_col(self$seasonal.mean$control$monthly.mean), "\n")
                                   cat(c(msg_0, msg_1, msg_2, msg_3))
                                   msg_0 <- "-----------------------\033[4;35m Variance Models \033[0m----------------------- \n"
                                   msg_1 <- paste0(" - Trend: ", msg_col(self$seasonal.variance$control$include.trend), "\n")
                                   msg_2 <- paste0(" - Correction: ", msg_col(self$seasonal.variance$control$correction), "\n")
                                   msg_3 <- paste0(" - Monthly correction: ", msg_col(self$seasonal.variance$control$monthly.mean), "\n")
                                   cat(c(msg_0, msg_1, msg_2, msg_3))
                                   msg_0 <- "------------------------\033[4;35m Mixture Model \033[0m------------------------ \n"
                                   msg_1 <- paste0(" - Method: ", self$mixture.model$method, "\n")
                                   msg_2 <- paste0(" - Match Expectation: ", msg_col(self$mixture.model$control$match.expectation), "\n")
                                   msg_3 <- paste0(" - Match Variance: ", msg_col(self$mixture.model$control$match.variance), "\n")
                                   msg_4 <- paste0(" - Match Empiric: ", msg_col(self$mixture.model$control$match.empiric), "\n")
                                   cat(c(msg_0, msg_1, msg_2, msg_3, msg_4))
                                 }
                               ),
                               private = list(
                                 version = "1.0.3",
                                 ..data = NA,
                                 ..dates = NA,
                                 ..target = NA,
                                 # ****************************
                                 ..transform = list(),
                                 ..seasonal_model_Ct = list(),
                                 ..seasonal.mean = list(),
                                 ..mean.model = list(),
                                 ..seasonal.variance = list(),
                                 ..variance.model = list(),
                                 ..mixture.model = list(),
                                 ..hessian = NA,
                                 ..jacobian = NA,
                                 # ****************************
                                 ..stochastic_clearsky = FALSE,
                                 ..clearsky_threshold = 1.01,
                                 ..quiet = FALSE
                               ),
                               active = list(
                                 #' @field target Character, name of the target variable to model. Can be `"GHI"` or `"clearsky"`.
                                 target = function(){
                                   private$..target
                                 },
                                 #' @field location A Tibble with the name and the coordinates of the location considered. Contains:
                                 #' \describe{
                                 #'  \item{place}{Character, selected location}
                                 #'  \item{target}{Target variable}
                                 #'  \item{lat}{Numeric, reference latitude in degrees.}
                                 #'  \item{lon}{Numeric, reference longitude in degrees.}
                                 #'  \item{alt}{Numeric, reference altitude in metres.}
                                 #'}
                                 location = function(){
                                   dplyr::bind_cols(place = self$place, target = self$target, dplyr::bind_rows(self$coords))
                                 },
                                 #' @field model_name Character, standard model's name.
                                 model_name = function(){
                                   paste0(self$transform$link, "-ARMA(", self$mean.model$arOrder, ", ", self$mean.model$maOrder, ")-",
                                          "GARCH(", self$variance.model$archOrder, ", ", self$variance.model$garchOrder, ")")
                                 },
                                 #' @field dates A named list, with three sub-lists: `data` containing the information on the complete dataset,
                                 #' `train` containing the information on the train dataset and `test` containing the information on the test dataset.
                                 #' Each sub-list is structured as follows:
                                 #' \describe{
                                 #'  \item{from}{Character date, minmum date in the dataset.}
                                 #'  \item{to}{Character date, maximum date in the dataset.}
                                 #'  \item{nobs}{Integer scalar, number of observations contained in the dataset between `from` and `to`.}
                                 #'  \item{perc}{Numeric scalar, percentage of data in the dataset with respect to the complete data.}
                                 #'}
                                 dates = function(){
                                   private$..dates
                                 },
                                 #' @field data Tibble, dataset with CAMS solar radiation data.
                                 data = function(){
                                   private$..data
                                 },
                                 #' @field transform A \code{\link{solarTransform}} object with the transformation applied to the data.
                                 transform = function(){
                                   private$..transform
                                 },
                                 #' @field seasonal_model_Ct A \code{\link{clearskySeasona}} object with the clear-sky model applied to the data.
                                 seasonal_model_Ct = function(){
                                   private$..seasonal_model_Ct
                                 },
                                 #' @field seasonal.mean A \code{\link{seasonalModel}} object with the seasonal model applied to the data.
                                 seasonal.mean = function(){
                                   private$..seasonal.mean
                                 },
                                 #' @field mean.model An \code{\link{ARMA_modelR6}} object with the ARMA model applied to the data.
                                 mean.model = function(){
                                   private$..mean.model
                                 },
                                 #' @field seasonal.variance A \code{\link{seasonalModel}} object with the seasonal variance applied to the data.
                                 seasonal.variance = function(){
                                   private$..seasonal.variance
                                 },
                                 #' @field variance.model An \code{\link{sGARCH}} object with the GARCH model applied to the data.
                                 variance.model = function(){
                                   private$..variance.model
                                 },
                                 #' @field mixture.model A \code{\link{solarMixture}} object with the monthly Gaussian Mixture models applied to the data.
                                 mixture.model = function(){
                                   private$..mixture.model
                                 },
                                 #' @field coefficients Get the model parameters as a named list.
                                 coefficients = function(){
                                   # 1. Clear sky seasonal model
                                   coefs_names <- c()
                                   seasonal_model_Ct <- as.list(self$seasonal_model_Ct$coefficients)
                                   # 2. Seasonal model Yt
                                   seasonal_model_Yt <- as.list(self$seasonal.mean$coefficients)
                                   # 3. AR model
                                   ARMA <- as.list(self$mean.model$coefficients)
                                   # 4. Seasonal variance model
                                   seasonal_variance <- as.list(self$seasonal.variance$coefficients)
                                   # 5. GARCH variance model
                                   GARCH <- as.list(self$variance.model$coefficients)
                                   # 6. Gaussian mixture model
                                   NM_model <- self$mixture.model$coefficients
                                   NM_mu_up <- dplyr::bind_rows(setNames(NM_model$mu1, paste0("mu1_", 1:12)))
                                   NM_mu_dw <- dplyr::bind_rows(setNames(NM_model$mu2, paste0("mu2_", 1:12)))
                                   NM_sd_up <- dplyr::bind_rows(setNames(NM_model$sd1, paste0("sd1_", 1:12)))
                                   NM_sd_dw <- dplyr::bind_rows(setNames(NM_model$sd2, paste0("sd2_", 1:12)))
                                   NM_p_up <- dplyr::bind_rows(setNames(NM_model$p1, paste0("p1_", 1:12)))
                                   # Output list of parameter for each model
                                   params <- list(
                                     location = self$location,
                                     params = dplyr::bind_cols(alpha = self$transform$alpha, beta = self$transform$beta),
                                     seasonal_model_Ct = dplyr::bind_cols(seasonal_model_Ct),
                                     seasonal_model_Yt = dplyr::bind_cols(seasonal_model_Yt),
                                     ARMA = dplyr::bind_cols(ARMA),
                                     seasonal_variance = dplyr::bind_cols(seasonal_variance),
                                     GARCH = dplyr::bind_cols(GARCH),
                                     NM_mu_up = NM_mu_up,
                                     NM_mu_dw = NM_mu_dw,
                                     NM_sd_up = NM_sd_up,
                                     NM_sd_dw = NM_sd_dw,
                                     NM_p_up = NM_p_up)

                                   return(params)
                                 },
                                 #' @field std.errors Get the model parameters as a named list.
                                 std.errors = function(){
                                   # 1. Clear sky seasonal model
                                   coefs_names <- c()
                                   seasonal_model_Ct <- as.list(self$seasonal_model_Ct$std.errors)
                                   # 2. Seasonal model Yt
                                   seasonal_model_Yt <- as.list(self$seasonal.mean$std.errors)
                                   # 3. AR model
                                   ARMA <- as.list(self$mean.model$std.errors)
                                   # 4. Seasonal variance model
                                   seasonal_variance <- as.list(self$seasonal.variance$std.errors)
                                   # 5. GARCH variance model
                                   GARCH <- as.list(self$variance.model$std.errors)
                                   # 6. Gaussian mixture model
                                   NM_model <- self$mixture.model$std.errors
                                   NM_mu_up <- dplyr::bind_rows(setNames(NM_model$mu1, paste0("mu1_", 1:12)))
                                   NM_mu_dw <- dplyr::bind_rows(setNames(NM_model$mu2, paste0("mu2_", 1:12)))
                                   NM_sd_up <- dplyr::bind_rows(setNames(NM_model$sd1, paste0("sd1_", 1:12)))
                                   NM_sd_dw <- dplyr::bind_rows(setNames(NM_model$sd2, paste0("sd2_", 1:12)))
                                   NM_p_up <- dplyr::bind_rows(setNames(NM_model$p1, paste0("p1_", 1:12)))
                                   # Output list of parameter for each model
                                   params <- list(
                                     location = self$location,
                                     params = dplyr::bind_cols(alpha = self$transform$alpha, beta = self$transform$beta),
                                     seasonal_model_Ct = dplyr::bind_cols(seasonal_model_Ct),
                                     seasonal_model_Yt = dplyr::bind_cols(seasonal_model_Yt),
                                     ARMA = dplyr::bind_cols(ARMA),
                                     seasonal_variance = dplyr::bind_cols(seasonal_variance),
                                     GARCH = dplyr::bind_cols(GARCH),
                                     NM_mu_up = NM_mu_up,
                                     NM_mu_dw = NM_mu_dw,
                                     NM_sd_up = NM_sd_up,
                                     NM_sd_dw = NM_sd_dw,
                                     NM_p_up = NM_p_up)

                                   return(params)
                                 },
                                 #' @field hessian Get the model hessian
                                 hessian = function(){
                                   private$..hessian
                                 },
                                 #' @field jacobian Get the model jacobian
                                 jacobian = function(){
                                   private$..jacobian
                                 },
                                 #' @field clearsky_threshold Numeric, parameter > 1, used to scale up CAMS clearsky.
                                 clearsky_threshold = function(){
                                   private$..clearsky_threshold
                                 },
                                 #' @field stochastic_clearsky Logical, when `TRUE` the clear sky is considered stochastic.
                                 stochastic_clearsky = function(){
                                   private$..stochastic_clearsky
                                 },
                                 #' @field quiet Logical. When `TRUE` the function will not display any message. The dafault if `TRUE`.
                                 quiet = function(){
                                   private$..quiet
                                 }
                               ))

