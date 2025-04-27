#' Control function for a `solarModel` object
#'
#' Control function for a `solarModel` object that contains all the setups used for the estimation.
#'
#' @examples
#' control <- control_solarModel$new()
#'
#' @rdname control_solarModel
#' @name control_solarModel
#' @keywords internal solarModel control
#' @export
solarModel_spec <- R6::R6Class("solarModel_spec",
                               public = list(
                                 initialize = function(){
                                   self$set_params()
                                   self$set_clearsky()
                                   self$set_seasonal.mean()
                                   self$set_mean.model()
                                   self$set_seasonal.variance()
                                   self$set_variance.model()
                                   self$set_mixture.model()
                                 },
                                 #' @description Specification function for a `solarModel`
                                 #' @param place Character, name of an element in the `CAMS_data` list.
                                 #' @param target Character, target variable to model. Can be `GHI` or `clearsky`.
                                 #' @param min_date Character. Date in the format `YYYY-MM-DD`. Minimum date for the complete data. If `missing` will be used the minimum data available.
                                 #' @param max_date Character. Date in the format `YYYY-MM-DD`. Maximum date for the complete data. If `missing` will be used the maximum data available.
                                 #' @param from Character. Date in the format `YYYY-MM-DD`. Starting date to use for training data.
                                 #' If `missing` will be used the minimum data available after filtering for `min_date`.
                                 #' @param to character. Date in the format `YYYY-MM-DD`. Ending date to use for training data.
                                 #' If `missing` will be used the maximum data available after filtering for `max_date`.
                                 #' @param CAMS_data named list with radiation data for different locations.
                                 specification = function(place, target = "GHI", min_date, max_date, from, to, CAMS_data = solarr::CAMS_data){
                                   # Match the target variable to model
                                   target <- match.arg(target, choices = c("GHI", "clearsky"))
                                   # Match a location in the dataset
                                   place <- match.arg(place, choices = names(CAMS_data), several.ok = FALSE)
                                   # Extract CAMS data for the selected location
                                   data <- CAMS_data[[place]]

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
                                   private$..place = attr(data, "place")
                                   private$..coords = attr(data, "coords")
                                   private$..data = data
                                   private$..dates = model_dates
                                   private$..target = target
                                 },
                                 #' Generic controls
                                 #' @param threshold Numeric. Threshold used to estimate the transformation parameters \deqn{\alpha} and \deqn{\beta}.
                                 #' The default is `0.01`. See \code{\link{solarTransform}} for more details.
                                 #' @param quiet Logical. When `TRUE` the function will not display any message. The dafault if `TRUE`.
                                 set_params = function(stochastic_clearsky = FALSE, min_pos = 1, max_pos = 1, transform_delta = 0.05,
                                                       clearsky_threshold = 1.01, threshold = 0.01, quiet = FALSE){
                                   private$..stochastic_clearsky = stochastic_clearsky
                                   private$..min_pos = min_pos
                                   private$..max_pos = max_pos
                                   private$..transform_delta = transform_delta
                                   private$..clearsky_threshold = clearsky_threshold
                                   private$..threshold = threshold
                                   private$..quiet = quiet
                                 },
                                 #' @description List of control parameters for the clear sky model.
                                 #' @param control See \code{\link{control_seasonalClearsky}} for more details.
                                 set_clearsky = function(control = control_seasonalClearsky()){
                                   private$..clearsky <- control
                                 },
                                 #' @description List of parameters to control the seasonal mean.
                                 #' @param order Integer. Specify the order of the seasonal mean \deqn{\bar{Y}_t}. The default is `1`.
                                 #' @param period Integer, seasonal periodicity, the default is `365`.
                                 #' @param include.trend Logical. When `TRUE` an yearly trend \deqn{t} will be included in the seasonal model, otherwise will be excluded. The default is `FALSE`.
                                 #' @param include.intercept Logical. When `TRUE` the intercept \deqn{a_0} will be included in the seasonal model, otherwise will be excluded. The default is `TRUE`.
                                 #' @param monthly.mean Logical. When `TRUE` a vector of 12 monthly means will be computed on the deseasonalized series \deqn{\tilde{Y}_t = Y_t - \bar{Y}_t}
                                 #'  and it is subtracted to ensure that the time series is centered around zero for all the months. The dafault if `TRUE`.
                                 set_seasonal.mean = function(order = 1, period = 365, include.trend = FALSE, include.intercept = TRUE, monthly.mean = TRUE){
                                   private$..seasonal.mean = list(order = order, period = period,
                                                                  include.trend = include.trend, include.intercept = include.intercept,
                                                                  monthly.mean = monthly.mean)
                                 },
                                 #' @description List of parameters to control the AR model.
                                 #' @param arOrder Integer. An integer specifying the order of the AR component. The default is `1`.
                                 #' @param maOrder Integer. An integer specifying the order of the MA component. The default is `0`.
                                 #' @param include.intercept Logical. When `TRUE` the intercept \deqn{\phi_0} will be included in the seasonal model, otherwise will be excluded. The default is `FALSE`.
                                 set_mean.model = function(arOrder = 1, maOrder = 0, include.intercept = FALSE){
                                   private$..mean.model = list(arOrder = arOrder, maOrder = maOrder, include.intercept = include.intercept)
                                 },
                                 #' @description List of parameters to control the seasonal variance.
                                 #' @param order Integer. Specify the order of the seasonality of the seasonal variance. The default is `1`.
                                 #' @param period Integer, seasonal periodicity, the default is `365`.
                                 #' @param include.trend Logical. When `TRUE` an yearly trend \deqn{t} will be included in the seasonal model, otherwise will be excluded. The default is `FALSE`.
                                 #' @param correction Logical. When `TRUE` the parameters of seasonal variance are corrected to ensure
                                 #'  that the standardize the residuals have exactly a unitary variance. The dafault if `TRUE`.
                                 #' @param monthly.mean Logical. When `TRUE` a vector of 12 monthly std. deviations will be computed
                                 #'  on the standardized residuals  \deqn{\tilde{\varepsilon}_t} and used to standardize the time series
                                 #'  such that it has unitary variance for all the months. The default if `TRUE`.
                                 set_seasonal.variance = function(order = 1, period = 365, include.trend = FALSE, correction = TRUE, monthly.mean = TRUE){
                                   private$..seasonal.variance = list(order = order, period = period, correction = correction,
                                                                      include.trend = include.trend, monthly.mean = monthly.mean)
                                 },
                                 #' @description List of parameters for the GARCH model.
                                 #' @param arOrder Integer. An integer specifying the order of the ARCH component. The default is `1`.
                                 #' @param maOrder Integer. An integer specifying the order of the GARCH component. The default is `1`.
                                 #' @param garch_variance Logical. When `TRUE` the GARCH model will be used to standardize the residuals otherwise will be excluded. The dafault if `TRUE`.
                                 set_variance.model = function(archOrder = 1, garchOrder = 1, garch_variance = TRUE){
                                   private$..garch_variance = garch_variance
                                   # Modify only if variance is not false
                                   if (private$..garch_variance != FALSE) {
                                     private$..variance.model = rugarch::ugarchspec(variance.model = list(garchOrder = c(archOrder,garchOrder), variance.targeting = 1),
                                                                                    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE))
                                     if (archOrder == 0 & garchOrder == 0){
                                       private$..garch_variance <- FALSE
                                     } else {
                                       private$..garch_variance <- TRUE
                                     }
                                   } else {
                                     private$..variance.model = rugarch::ugarchspec(variance.model = list(garchOrder = c(0,0), variance.targeting = 1),
                                                                                    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE))
                                   }
                                 },
                                 #' @description List of parameters for the Gaussian mixture model.
                                 #' @param abstol Numeric. Absolute level for convergence of the EM-algorithm. The default is `1e-20`.
                                 #' @param match.moments description
                                 #' @param maxit Integer. Maximum number of iterations for EM-algorithm. The default is `100`.
                                 set_mixture.model = function(abstol = 1e-20, match.moments = TRUE, B = 0, method = "mclust", maxit = 5000){
                                   method <- match.arg(method, choices = c("mclust", "mixtools"))
                                   private$..mixture.model = list(abstol = abstol,
                                                                  B = B,
                                                                  method = method,
                                                                  match.moments = match.moments,
                                                                  maxit = maxit)
                                 },
                                 #' @description
                                 #' Print method for `control_solarModel` class.
                                 print = function(){
                                   # Seasonal mean order
                                   sm_order <- self$seasonal.mean$order
                                   # ARMA order
                                   arOrder <- self$mean.model$arOrder
                                   maOrder <- self$mean.model$maOrder
                                   # Seasonal variance order
                                   sv_order <- self$seasonal.variance$order
                                   # GARCH order
                                   par <- self$variance.model@model$modelinc
                                   par_names <- names(par[par == 1])
                                   archOrder <- sum(stringr::str_detect(par_names, "alpha"))
                                   garchOrder <- sum(stringr::str_detect(par_names, "beta"))

                                   if (!is.na(self$place)) {
                                     # Complete data specifications
                                     data <- self$dates$data
                                     # Train data specifications
                                     train <- self$dates$train
                                     train$perc <- format(train$perc*100, digits = 4)
                                     # Test data specifications
                                     test <- self$dates$test
                                     test$perc <- format(test$perc*100, digits = 4)
                                     msg_0 <- paste0("--------------------- ", "solarModel", " (", self$place, ") ", "--------------------- \n")
                                     msg_1 <- paste0("Target: ", self$target, " \n Coordinates: (Lat: ", self$coords$lat, ", Lon: ", self$coords$lon, ", Alt: ", self$coords$alt, ") \n")
                                     msg_2 <- paste0(" Dates: ", data$from, " - ", data$to, "\n Observations: ", data$nobs, "\n")
                                     msg_3 <- paste0("---------------------------------------------------------------\n")
                                     msg_4 <- paste0("Train dates: ", train$from, " - ", train$to, " (", train$nobs, " points ~ ", train$perc, "%)", "\n")
                                     msg_5 <- paste0("Test  dates: ", test$from, " - ", test$to, " (", test$nobs, " points ~ ", test$perc, "%)", "\n")
                                     cat(paste0(msg_0, msg_1, msg_2, msg_3, msg_4, msg_5))
                                   }
                                   cat("-------------------- control Parameters  -------------------- \n")
                                   model_name <- paste0("S(", sm_order, ", ", sv_order, ")-ARMA(", arOrder, ", ", maOrder, ")")
                                   if (self$garch_variance) {
                                     model_name <- paste0(model_name,"-GARCH(", archOrder, ", ", archOrder, ")")
                                   }
                                   cat(paste0(model_name, "\n"))
                                   cat(paste0(" - Stochastic clearsky: ", self$stochastic_clearsky, "\n"))
                                   cat(paste0(" - Monthly means correction: ", self$seasonal.mean$monthly.mean, "\n"))
                                   cat(paste0(" - Monthly variance correction: ", self$seasonal.variance$monthly.mean, "\n"))
                                 }
                               ),
                               private = list(
                                 ..place = NA,
                                 ..coords = NA,
                                 ..data = NA,
                                 ..dates = NA,
                                 ..target = NA,
                                 ..clearsky = control_seasonalClearsky(),
                                 ..seasonal.mean = list(),
                                 ..mean.model = list(),
                                 ..seasonal.variance = list(),
                                 ..variance.model = list(),
                                 ..mixture.model = list(),
                                 ..garch_variance = FALSE,
                                 ..clearsky_threshold = 1.01,
                                 ..threshold = 0.01,
                                 ..quiet = FALSE,
                                 ..stochastic_clearsky = FALSE,
                                 ..min_pos = 1,
                                 ..max_pos = 1,
                                 ..transform_delta = 0.05
                               ),
                               active = list(
                                 place = function(){
                                   private$..place
                                 },
                                 coords = function(){
                                   private$..coords
                                 },
                                 data = function(){
                                   private$..data
                                 },
                                 dates = function(){
                                   private$..dates
                                 },
                                 target = function(){
                                   private$..target
                                 },
                                 clearsky = function(){
                                   private$..clearsky
                                 },
                                 seasonal.mean = function(){
                                   private$..seasonal.mean
                                 },
                                 mean.model = function(){
                                   private$..mean.model
                                 },
                                 seasonal.variance = function(){
                                   private$..seasonal.variance
                                 },
                                 variance.model = function(){
                                   private$..variance.model
                                 },
                                 mixture.model = function(){
                                   private$..mixture.model
                                 },
                                 garch_variance = function(){
                                   private$..garch_variance
                                 },
                                 clearsky_threshold = function(){
                                   private$..clearsky_threshold
                                 },
                                 threshold = function(){
                                   private$..threshold
                                 },
                                 stochastic_clearsky = function(){
                                   private$..stochastic_clearsky
                                 },
                                 quiet = function(){
                                   private$..quiet
                                 },
                                 min_pos = function(){
                                   private$..min_pos
                                 },
                                 max_pos = function(){
                                   private$..max_pos
                                 },
                                 transform_delta = function(){
                                   private$..transform_delta
                                 }
                               ))


#' Compute the conditional moments
#' @keywords internal solarModel
#' @export
solarModel_conditional_moments = function(data, control){
  if (!control$garch_variance){
    data$sigma <- 1
  }
  # Compute conditional moments
  data <- dplyr::mutate(data,
                        # Conditional expectation Yt
                        e_Yt = Yt_bar + Yt_tilde_hat + Yt_tilde_uncond,
                        # Conditional std. deviation Yt
                        sd_Yt = sigma * sigma_bar * sigma_uncond,
                        # Conditional moments Yt (state up)
                        M_Y1 = e_Yt + sd_Yt * mu1,
                        S_Y1 = sd_Yt * sd1,
                        # Conditional moments Yt (state dw)
                        M_Y0 = e_Yt + sd_Yt * mu2,
                        S_Y0 = sd_Yt * sd2)
  # Store only relevant variables
  data <- dplyr::select(data, date, Year, Month, Day, e_Yt, sd_Yt, M_Y1, S_Y1, M_Y0, S_Y0)
  return(data)
}

#' Compute the unconditional moments
#' @keywords internal solarModel
#' @export
solarModel_unconditional_moments = function(seasonal_data, ARMA, GARCH){
  # Arguments
  data <- seasonal_data
  # ARMA variance
  arma_variance <- sqrt(ARMA$variance)
  # Compute unconditional moments
  data <- dplyr::mutate(data,
                        # Unconditional expectation Yt
                        e_Yt = ARMA$intercept + Yt_bar + Yt_tilde_uncond,
                        # Unconditional std. deviation Yt
                        sd_Yt = sigma_bar * sigma_uncond * GARCH$vol,
                        # Unconditional moments Yt (state up)
                        M_Y1 = e_Yt + sd_Yt * mu1 * arma_variance,
                        S_Y1 = sd_Yt * sd1 * arma_variance,
                        # Unconditional moments Yt (state dw)
                        M_Y0 = e_Yt + sd_Yt * mu2 * arma_variance,
                        S_Y0 = sd_Yt * sd2 * arma_variance)
  # Store only relevant variables
  data <- dplyr::select(data, Month, Day, e_Yt, sd_Yt, M_Y1, S_Y1, M_Y0, S_Y0)
  return(data)
}

#' Moments from a solarModel object
#'
#' @examples
#' model = solarModel$new(spec)
#' model$fit()
#' t_now = "2022-01-30"
#' t_hor = "2022-05-01"
#' data = model$data
#' ARMA = model$ARMA
#' GARCH  = model$GARCH
#' NM_model = model$NM_model
#' transform = model$transform
#' solarModel_moments(t_now, t_hor, data, ARMA, GARCH, NM_model, transform, quiet = FALSE)
#' @keywords internal solarModel
#' @export
solarModel_moments = function(t_now, t_hor, data, ARMA, GARCH, NM_model, transform, quiet = FALSE){
  if(!quiet) message("Forecast: ", t_hor, " given ", t_now, "\r", appendLF = FALSE)
  # Conditioning date
  t_now <- as.Date(t_now)
  # Horizon date
  t_hor <- as.Date(t_hor)
  # Maximum order of AR / GARCH
  lag_max <- max(c(ARMA$order, GARCH$order))
  # Filter data between (t_now - lag_max + 1) and t_hor
  data <- dplyr::filter(data, date >= (t_now - lag_max + 1) & date <= t_hor)
  # Extract mixture moments
  mom <- NM_model$moments
  nm <- NM_model$coefficients
  # Unknown data till maturity
  df_tT <- dplyr::left_join(data[-c(1:(lag_max)),], mom, by = "Month")
  # Known data at time t_now (used as vector only in state-space forecast)
  df_t <- data[c(1:lag_max),]
  # Must be reversed ordered, i.e. most recent (t, t-1, t-2, ...)
  df_t <- df_t[order(df_t$date, decreasing = TRUE),]
  # Forecasting horizon
  h <- nrow(df_tT)
  # *******************************************************************************
  #  1) Forecast GARCH moments
  # *******************************************************************************
  # Second moment GARCH variance (exact)
  df_tT$e_sigma2 <- e_sigma2_h_mix(h - 1, GARCH$omega, GARCH$alpha, GARCH$beta, e_x2 = df_tT$m2[-h], df_tT$sigma[1]^2)
  # Second moment of GARCH std. dev (exact)
  df_tT$e_sigma4 <- e_sigma4_h_mix(h - 1, GARCH$omega, GARCH$alpha, GARCH$beta, e_x2 = df_tT$m2[-h], e_x4 = df_tT$m4[-h], df_tT$sigma[1]^4)
  # Variance of GARCH (exact)
  df_tT$v_sigma2 <- df_tT$e_sigma4 - df_tT$e_sigma2^2
  # First moment of GARCH std.dev (approximated)
  df_tT$e_sigma1 <-  df_tT$e_sigma2^(1/2) - (1/8) * df_tT$v_sigma2 / sqrt(df_tT$e_sigma2)^3
  # Third moment of GARCH std.dev (approximated)
  df_tT$e_sigma3 <-  df_tT$e_sigma2^(3/2) + (3/8) * df_tT$v_sigma2 / sqrt(df_tT$e_sigma2)
  # Variance of GARCH (exact)
  df_tT$v_sigma <- df_tT$e_sigma2 - df_tT$e_sigma1^2
  # Update the mixture variances
  # df_tT$sigma2_1 <- df_tT$e_sigma2 * df_tT$sd1^2 + df_tT$v_sigma * df_tT$mu1^2
  # df_tT$sigma2_2 <- df_tT$e_sigma2 * df_tT$sd2^2 + df_tT$v_sigma * df_tT$mu2^2
  # Central moments
  # df_tT$delta_1 <- df_tT$mu1 - df_tT$mean
  # df_tT$delta_2 <- df_tT$mu2 - df_tT$mean
  # df_tT$variance <- (df_tT$sigma2_1 + df_tT$delta_1^2) * df_tT$p1 + (df_tT$sigma2_2 + df_tT$delta_2^2) * df_tT$p2
  # Mixture third moment
  # df_tT$m3 <- (3 * df_tT$delta_1 * df_tT$sigma2_1 + df_tT$delta_1^3) * df_tT$p1 + (3 * df_tT$delta_2 * df_tT$sigma2_2 + df_tT$delta_2^3) * df_tT$p2
  # *******************************************************************************
  #  2)  Forecast seasonal variance
  # *******************************************************************************
  df_tT$sigma_bar <- df_tT$sigma_bar * df_tT$sigma_uncond
  # *******************************************************************************
  #  3)  Forecast mean and variance of Yt_tilde
  # *******************************************************************************
  # Companion matrix
  A <- ARMA$Phi
  # Residuals vector for mean
  b <- matrix(ARMA$b, ncol = 1)
  # Residuals matrix for variance
  B <- b %*% t(b)
  # Extract first component
  e1 <- matrix(rep(1, length(b)), ncol = 1)
  e1[-1] <- 0
  # Conditioning values
  Y0 <- df_t$Yt_tilde[1:ARMA$order[1]]
  eps0 <- c()
  if (ARMA$order[2] > 0){
    eps0 <- df_t$eps[1:ARMA$order[2]]
  }
  # State vector
  Xt <- c(Y0, eps0)
  # Initialize the variables
  df_tT$psi_j <- df_tT$psi2_j <- df_tT$psi3_j <- NA
  df_tT$intercept <- 0
  for(j in 1:h){
    # Pre compute the powers
    A_pow_j <- pow_matrix(A, j)
    A_pow_hj <- pow_matrix(A, h - j)
    # Compute weights for expectations
    df_tT$psi_j[j] <- t(e1) %*% (A_pow_hj %*% (b * df_tT$sigma_bar[j] * df_tT$e_sigma1[j] * df_tT$mean[j]))
    # Intercept contribution
    df_tT$intercept[j] <- t(e1) %*% (A_pow_j %*% (e1 * (ARMA$intercept)))
    # Variance
    df_tT$psi2_j[j] <- t(e1) %*% (A_pow_hj %*% (df_tT$e_sigma2[j] * df_tT$sigma_bar[j]^2 * B * df_tT$variance[j]) %*% t(A_pow_hj)) %*% e1
    # Skewness (third moment)
    df_tT$psi3_j[j] <- (t(e1) %*% (pow_matrix(A_pow_hj, 3) %*% (b * df_tT$e_sigma3[j] * df_tT$sigma_bar[j]^3 * df_tT$m3[j])))
    # Forecasted value
    df_tT$Yt_tilde_hat[j] <- ARMA$intercept + sum(df_tT$intercept[1:j]) + t(e1) %*% A_pow_j %*% Xt
  }
  df_tT$psi_j[h] <- 0
  # Forecasted mean value
  df_tT$e_Yt_tilde <- df_tT$Yt_tilde_hat + cumsum(df_tT$psi_j)
  df_tT$e_Yt <- df_tT$e_Yt_tilde + df_tT$Yt_tilde_uncond + df_tT$Yt_bar
  # Forecasted multinomial variance
  df_tT$Sigma2 <- cumsum(df_tT$psi2_j)
  # Forecasted multinomial third moment
  df_tT$Omega <- cumsum(df_tT$psi3_j)
  # Last value
  df_T <- tail(df_tT, 1)
  # *******************************************************************************
  #  4)  Approximate the multinomial mixture with a 2-component GM
  # *******************************************************************************
  # Approximated mixture means
  M_Y <- c(df_T$mu1, df_T$mu2) * df_T$sigma_bar * df_T$e_sigma1 + df_T$e_Yt
  # Approximated mixture variances
  S2_Y <- c(df_T$sd1^2, df_T$sd2^2) * df_T$sigma_bar^2 * df_T$e_sigma2 + cumsum(df_tT$psi2_j)[h-1]
  # Exact mixture conditional probabilities
  p <- c(df_T$p1, 1-df_T$p1)
  # Target second moment constraint
  m2_target <- df_T$Sigma2
  # Target third moment constraint
  m3_target <- df_T$Omega
  # Central moments
  delta_k <- M_Y - sum(M_Y*p)
  # Create the system in matrix form
  # A_ <- matrix(c(p[1], p[2], 3 * p[1] * delta_k[1], 3 * p[2] * delta_k[2]),nrow = 2, ncol = 2, byrow = TRUE)
  # b_ <- matrix(c(m2_target - p[1] * delta_k[1]^2 - p[2] * delta_k[2]^2, m3_target - p[1] * delta_k[1]^3 - p[2] * delta_k[2]^3), ncol = 1)
  # S2_Y_star <- c(solve(A_) %*% b_)
  # Solve the variances to match the second and third moment
  D <- 3 * p[1] * p[2] * (delta_k[2] - delta_k[1])
  D1 <- (m2_target - p[1] * delta_k[1]^2 - p[2] * delta_k[2]^2) * 3 * delta_k[2] * p[2] - p[2] * (m3_target - p[1] * delta_k[1]^3 - p[2] * delta_k[2]^3)
  D0 <- p[1] * (m3_target - p[1] * delta_k[1]^3 - p[2] * delta_k[2]^3)  - (m2_target - p[1] * delta_k[1]^2 - p[2] * delta_k[2]^2) * 3 * delta_k[1] * p[1]
  S2_Y_star <- c(D1/D, D0/D)
  # Ensure variances are positive
  if (any(S2_Y_star < 0)) {
    message("S_Y gives NaN for ", t_hor, " given ", t_now)
    S_Y_star <- sqrt(S2_Y)
  } else {
    S_Y_star <- sqrt(S2_Y_star)
  }
  dplyr::tibble(
    date = t_hor,
    Year = lubridate::year(t_hor),
    Month = lubridate::month(t_hor),
    Day = lubridate::day(t_hor),
    e_Yt = df_T$e_Yt,
    sd_Yt = sqrt(df_T$Sigma2),
    M_Y1 = M_Y[1],
    S_Y1 = S_Y_star[1],
    M_Y0 = M_Y[2],
    S_Y0 = S_Y_star[2],
    p1 = df_T$p1,
    Ct = df_T$Ct,
    alpha = transform$alpha,
    beta = transform$beta
  )
}

#' Match solarModel parameters in vector form
#'
#' @examples
#' model = solarModel$new(spec)
#' model$fit()
#' vec_params <- c(theta = 1, alpha1 = 10)
#' solarModel_match_params(vec_params, model$coefficients)
#' @rdname solarModel_match_params
#' @name solarModel_match_params
#' @keywords internal solarModel
#' @export
solarModel_match_params <- function(vec_params, params){
  # List of names of the original parameters
  names_params <- purrr::map(params, ~names(.x))
  # Names of the parameters in vector form
  names_vec_params <- names(vec_params)
  for(i in 1:length(vec_params)){
    # Parameter to use
    target_param <- vec_params[i]
    # Parameter name
    names_target_param <- names(target_param)
    # Check the parameter inside the list of parameters
    idx_list <- purrr::map_lgl(names_params, ~sum(.x %in% names_target_param) == 1)
    if (purrr::is_empty(which(idx_list))) {
      cli::cli_alert_warning(paste0('Parameter: "', names_target_param, '" not found in the model specification. Ignored!'))
      next
    }
    idx_list <- which(idx_list)
    # Position of the parameter inside the list
    idx_position <- names(params[[idx_list]]) == names_target_param
    # Update the parameter
    params[[idx_list]][idx_position] <- target_param[[1]]
  }
  return(params)
}


#' Produce a forecast from a solarModel object
#' @param theta Esscher parameter
#' @examples
#' model = solarModel$new(spec)
#' model$fit()
#' moments <- model$moments$conditional[14,]
#' object <- solarModel_predict(model, moments, ci = 0.01)
#' object
#' @rdname solarModel_predict
#' @name solarModel_predict
#' @export
solarModel_predict <- function(model, moments, theta = 0, ci = 0.01){
  # Moments
  df_n <- moments
  # Create datasets for the moments
  comb <- dplyr::tibble(mean = c(df_n$M_Y1, df_n$M_Y0), sd = c(df_n$S_Y1, df_n$S_Y0), probs = c(df_n$p1, 1-df_n$p1))
  # Distribution and Density of Yt
  if (theta == 0) {
    # Mixture pdf
    pdf_Yt <- function(x) dmixnorm(x, comb$mean, comb$sd, comb$probs)
    pdf_Yt_up <- function(x) dmixnorm(x, comb$mean[1], comb$sd[1], 1)
    pdf_Yt_dw <- function(x) dmixnorm(x, comb$mean[2], comb$sd[2], 1)
    # Mixture cdf
    cdf_Yt <- function(x) pmixnorm(x, comb$mean, comb$sd, comb$probs)
    cdf_Yt_up <- function(x) pmixnorm(x, comb$mean[1], comb$sd[1], 1)
    cdf_Yt_dw <- function(x) pmixnorm(x, comb$mean[2], comb$sd[2], 1)
  } else {
    # Esscher Mixture pdf
    pdf_Yt <- desscherMixture(comb$mean, comb$sd, comb$probs, theta)
    pdf_Yt_up <- desscherMixture(x, comb$mean[1], comb$sd[1], 1, theta)
    pdf_Yt_dw <- desscherMixture(x, comb$mean[2], comb$sd[2], 1, theta)
    # Esscher Mixture cdf
    cdf_Yt <- pesscherMixture(comb$mean, comb$sd, comb$probs, theta)
    cdf_Yt_up <- pesscherMixture(x, comb$mean[1], comb$sd[1], 1, theta)
    cdf_Yt_dw <- pesscherMixture(x, comb$mean[2], comb$sd[2], 1, theta)
  }
  # Expected value of Rt^q
  e_Rt_q <- function(q = 1, pdf_Yt) integrate(function(x) model$transform$GHI_y(x, df_n$Ct)^q * pdf_Yt(x), lower = -Inf, upper = Inf)$value

  # Expected values
  # Mixture
  df_n$e_Rt <- e_Rt_q(1, pdf_Yt)
  # Sunny state
  df_n$e_Rt_up <- e_Rt_q(1, pdf_Yt_up)
  # Cloudy state
  df_n$e_Rt_dw <- e_Rt_q(1, pdf_Yt_dw)

  # Variances
  # Mixture
  df_n$v_Rt <- e_Rt_q(2, pdf_Yt) - df_n$e_Rt^2
  # Sunny state
  df_n$v_Rt_up <- e_Rt_q(2, pdf_Yt_up) - df_n$e_Rt_up^2
  # Cloudy state
  df_n$v_Rt_dw <- e_Rt_q(2, pdf_Yt_dw)- df_n$e_Rt_dw^2

  # Confidence interval
  # Mixture
  df_n$ci_mix_lo <- qsolarGHI(ci, df_n$Ct, model$transform$alpha, model$transform$beta, cdf_Yt)
  df_n$ci_mix_hi <- qsolarGHI(1-ci, df_n$Ct, model$transform$alpha, model$transform$beta, cdf_Yt)
  # Sunny state
  df_n$ci_up_lo <- qsolarGHI(ci, df_n$Ct, model$transform$alpha, model$transform$beta, cdf_Yt_up)
  df_n$ci_up_hi <- qsolarGHI(1-ci, df_n$Ct, model$transform$alpha, model$transform$beta, cdf_Yt_up)
  # Cloudy state
  df_n$ci_dw_lo <- qsolarGHI(ci, df_n$Ct, model$transform$alpha, model$transform$beta, cdf_Yt_dw)
  df_n$ci_dw_hi <- qsolarGHI(1-ci, df_n$Ct, model$transform$alpha, model$transform$beta, cdf_Yt_dw)

  # Number of points for the grid
  n_points <- 100
  # Compute bounds for GHI
  lower_Rt = df_n$Ct*model$transform$bounds("Kt")[1]
  upper_Rt = df_n$Ct*model$transform$bounds("Kt")[2]
  # Grid for PDF plot
  grid_x <- seq(lower_Rt, upper_Rt, length.out = n_points+2)[-c(1,n_points+2)]
  grid <- dplyr::tibble(x = grid_x)
  # Density GHI
  pdf_Rt <- function(x, pdf_Yt) dsolarGHI(x, df_n$Ct, model$transform$alpha, model$transform$beta, pdf_Yt)
  # Density points (Mixture)
  grid$pdf_Rt_mix <- pdf_Rt(grid$x, pdf_Yt)
  # Density points (Mixture, up)
  grid$pdf_Rt_mix_up <- pdf_Rt(grid$x, pdf_Yt_up) * df_n$p1
  # Density points (Mixture, dw)
  grid$pdf_Rt_mix_dw <- pdf_Rt(grid$x, pdf_Yt_dw) * (1 - df_n$p1)

  # Add the value of the points for plotting
  # Expected value and variance
  # Mixture
  df_n$pdf_e_Rt <- pdf_Rt(df_n$e_Rt, pdf_Yt)
  # Confidence intervals
  # Mixture
  df_n$pdf_ci_mix_lo <- pdf_Rt(df_n$ci_mix_lo, pdf_Yt)
  df_n$pdf_ci_mix_hi <- pdf_Rt(df_n$ci_mix_hi, pdf_Yt)
  # Sunny state
  df_n$pdf_ci_up_lo <- pdf_Rt(df_n$ci_up_lo, pdf_Yt)
  df_n$pdf_ci_up_hi <- pdf_Rt(df_n$ci_up_hi, pdf_Yt)
  # Cloudy state
  df_n$pdf_ci_dw_lo <- pdf_Rt(df_n$ci_dw_lo, pdf_Yt)
  df_n$pdf_ci_dw_hi <- pdf_Rt(df_n$ci_dw_hi, pdf_Yt)

  # Realized GHI
  Rt <- dplyr::filter(model$data, date == df_n$date)$GHI
  df_n$Rt <- ifelse(purrr::is_empty(Rt), NA_integer_, Rt)

  structure(
    list(
      grid = grid,
      df_n = df_n,
      ci = ci
    ),
    class = c("solarModelForecast", "list")
  )
}

#' Iterate the forecast on multiple dates
#' @rdname solarModel_forecast
#' @name solarModel_forecast
#' @export
solarModel_forecast <- function(model, moments, ci = 0.1, theta = 0){
  fun <- function(df_n){
    safe_forecaster <- purrr::safely(solarModel_predict)
    smf <- safe_forecaster(model, moments = df_n, ci = ci, theta = theta)$result
    if (is.null(smf)){
      return(NULL)
    }
    return(smf[-1][[1]])
  }
  out <- purrr::map_df(1:nrow(moments), ~fun(moments[.x,]))
  return(out)
}

#' Calibrate the parameters
#' @rdname solarModel_loglik_calibrator
#' @name solarModel_loglik_calibrator
#' @export
solarModel_loglik_calibrator <- function(model, abstol = 1e-4, reltol = 1e-4, quiet = FALSE) {
  # Extract all the parameters
  coefs <- model$coefficients
  # Extract ARMA and seasonal variance parameters
  #params <- c(unlist(coefs$ARMA), unlist(coefs$seasonal_variance))
  # Remove ARMA intercept
  #if (params[1] == 0) {
  #  params <- params[-c(1)]
  #}
  params <- c()
  # Extract GARCH parameters if used
  if (model$spec$garch_variance) {
    params <- c(params, unlist(coefs$GARCH))
    # Remove GARCH intercept
    if (!is.na(model$GARCH$.__enclos_env__$private$..sigma20)) {
      params <- params[-c(which(names(params) == "omega"))]
    }
  }
  params <- params[params != 0]
  # Loss function
  loss_function <- function(params, model){
    # Clone the model
    model_cal <- model$clone(TRUE)
    # Check positivity arch parameters
    alpha <- params[stringr::str_detect(names(params), "alpha")]
    if (!purrr::is_empty(alpha)){if (any(alpha < 0)){return(1e6)}}
    # Check positivity garch parameters
    beta <- params[stringr::str_detect(names(params), "beta")]
    if (!purrr::is_empty(beta)){if (any(beta < 0)){return(1e6)}}
    # Check stationary arma
    phi <- params[stringr::str_detect(names(params), "phi")]
    if (!purrr::is_empty(phi)){if (sum(phi) > 1) {return(1e6)}}
    # Check positivity seasonal variance
    c_ <- params[stringr::str_detect(names(params), "c_")]
    if (!purrr::is_empty(c_)){if (c_[1] < sqrt(c_[2]^2 + c_[3]^2)) {return(1e6)}}
    if(!quiet) {print(params)}
    # Update the parameters
    model_cal$update(params)
    # Filter the time series
    model_cal$filter()
    # Re-fit the Gaussian Mixture
    model_cal$fit_NM_model()
    # Update the moments
    model_cal$update_moments()
    # Update log-likelihood
    model_cal$update_logLik()
    # Negative Log-likelihood
    loglik <- -model_cal$loglik
    # Loss
    loss <- loglik
    if (!quiet) {
      print(paste0("Log-likelihood: ", loglik))
    }
    return(loss)
  }
  # Optimization
  opt <- optim(params, loss_function, model = model,
               control = list(maxit = 30, abstol = abstol, reltol = reltol), hessian = TRUE)
  # Optimal parameters
  vec_params <- opt$par
  names(vec_params) <- names(params)
  # Calibrated model
  model_cal <- model$clone(TRUE)
  # Update the parameters
  model_cal$update(vec_params)
  # Filter the time series
  model_cal$filter()
  # Re-fit the Gaussian Mixture
  model_cal$fit_NM_model()
  model_cal$update_NM_classification()
  # Update the moments
  model_cal$update_moments()
  # Update log-likelihood
  model_cal$update_logLik()
  return(model_cal)
}

#' Distribution test
#'
#' Evaluate a Kolmogorov-smirnov test on the residuals of a `solarModel` model
#' object against the estimated Gaussian mixture distribution.
#'
#' @param model An object of the class `solarModel`
#' @inheritParams ks_test
#' @examples
#' model = solarModel$new(spec)
#' model$fit()
#' solarModel_test_distribution(model)
#' @rdname solarModel_test_distribution
#' @name solarModel_test_distribution
#' @export
solarModel_test_distribution <- function(model, ci = 0.05, min_quantile = 0.025, max_quantile = 0.985){
  test <- list()
  for(nmonth in 1:12){
    # Residuals
    x <- dplyr::filter(model$data, Month == nmonth)$u_tilde
    # Gaussian mixture model
    gm <- model$NM_model$model[[nmonth]]
    # Mixture CDF
    cdf_Yt <- function(x) pmixnorm(x, gm$means, gm$sd, gm$p)
    # Distribution test
    test[[nmonth]] <- ks_test(x, cdf_Yt, ci = ci, min_quantile = min_quantile, max_quantile = max_quantile)
  }
  test <- dplyr::bind_cols(Month = 1:12, dplyr::bind_rows(test))
  return(test)
}

#' Autocorrelation test
#'
#' Evaluate the autocorrelation in the components of a `solarModel` object.
#' @param model An object of the class `solarModel`
#' @param lag.max Numeric, scalar. Maximum lag to consider for the test.
#' @param ci Numeric, scalar. Minimum p-value to consider the test `"passed"`.
#' @param method Character. Type of test. Can be `"bp"` for breush-pagan or `"lb"` for Box-pierce.
#' @examples
#' model = solarModel$new(spec)
#' model$fit()
#' solarModel_test_autocorr(model, method = "lb")
#' @rdname solarModel_test_autocorr
#' @name solarModel_test_autocorr
#' @export
solarModel_test_autocorr <- function(model, lag.max = 3, ci = 0.05, method = "bp"){

  # Test Dataset
  data <- model$data
  data <- dplyr::left_join(data, model$NM_model$moments, by = c("Month"))
  # Standardized residuals
  data$z <- (data$u_tilde - data$mean) / sqrt(data$variance)
  # Compute squared residuals
  data$Yt_tilde2 <- data$Yt_tilde^2
  data$eps2 <- data$eps^2
  data$eps_tilde2 <- data$eps_tilde^2
  data$u_tilde2 <- data$u_tilde^2
  data$z2 <- data$z^2

  # Autocorrelation test Breusch-Godfrey
  bpagan_test <- function(target, data, lag.max, ci, expected = "Not-rejected"){
    formula <- as.formula(paste0(target, "~ 1"))
    test <- lmtest::bgtest(formula, order = lag.max, data = data)
    test <- broom::tidy(test)
    test$target <- target
    test <- dplyr::select(test, target, statistic, p.value, H0 = "method", lags = "parameter")
    test$H0 <- ifelse(test$p.value > ci, "Not-rejected", "Rejected")
    test$p.value <- round(test$p.value, digits = 5)
    test$result <- ifelse(test$H0 == expected, "passed", "Not-passed")
    return(test)
  }
  # Autocorrelation test Box.pierce
  ljungbox_test <- function(target, data, lag.max, ci, expected = "Not-rejected"){
    formula <- as.formula(paste0(target, "~ 1"))
    test <- Box.test(data[[target]], lag = lag.max)
    test <- broom::tidy(test)
    test$target <- target
    test <- dplyr::select(test, target, statistic, p.value, H0 = "method", lags = "parameter")
    test$H0 <- ifelse(test$p.value > ci, "Not-rejected", "Rejected")
    test$p.value <- round(test$p.value, digits = 5)
    test$result <- ifelse(test$H0 == expected, "passed", "Not-passed")
    return(test)
  }
  # Residuals
  targets <- c("Yt_tilde", "eps", "eps_tilde", "u_tilde", "z")
  expected <- c("Rejected", "Not-rejected", "Not-rejected", "Not-rejected", "Not-rejected")
  # Squared residuals
  targets <- c(targets, "Yt_tilde2", "eps2", "eps_tilde2", "u_tilde2", "z2")
  expected <- c(expected, "Rejected", "Rejected", "Rejected", "Not-rejected", "Not-rejected")
  names(expected) <- targets

  tests <- list()
  for(target in targets){
    if (method == "bp"){
      tests[[target]] <- bpagan_test(target, data, lag.max, ci, expected = expected[target])
    } else if (method == "lb") {
      tests[[target]] <- ljungbox_test(target, data, lag.max, ci, expected = expected[target])
    }
  }
  tests = bind_rows(tests)
  return(tests)
}


#' Autocorrelation and Distribution tests
#'
#' Evaluate a Kolmogorov-Smirnov test on the residuals of a `solarModel` model
#' object against the estimated Gaussian mixture distribution and a Breush-pagan or Box-pierce
#' test on the residuals.
#'
#' @inheritParams solarModel_test_distribution
#' @inheritParams solarModel_test_autocorr
#' @examples
#' model = solarModel$new(spec)
#' model$fit()
#' solarModel_tests(model)
#' @rdname solarModel_tests
#' @name solarModel_tests
#' @export
solarModel_tests <- function(model, lag.max = 3, ci = 0.05, min_quantile = 0.025, max_quantile = 0.985, method = "bp"){
  # Distribution test
  df_test_distrib <- solarModel_test_distribution(model, ci, min_quantile, max_quantile)
  df_test_distrib$result <- ifelse(df_test_distrib$H0 == "Non-Rejected", "Passed", "Not-passed")
  df_test_distrib <- dplyr::bind_cols(target = "u_tilde", df_test_distrib)
  # Autocorrelation test
  df_test_autocorr <- solarModel_test_autocorr(model, lag.max, ci, method)
  structure(
    list(
      autocorr = df_test_autocorr,
      distrib = df_test_distrib
    )
  )
}

#' Best models fit
#'
#' @examples
#' solarModels_fit_best_model("Roma", 10, 0.05)
#' solarModels_fit_best_model("Bologna", 10, 0.05)
#' solarModels_fit_best_model("Palermo", 10, 0.05)
#'
#' @export
solarModels_fit_best_model <- function(place, lag.max = 5, ci = 0.01, control_model = control_solarModel(), ...){
  control_model$garch_variance <- FALSE
  control_model$mean.model$arOrder <- 1
  control_model$mean.model$maOrder <- 0
  narch <- ngarch <- 0
  condition <- TRUE
  while(condition) {
    # Model specification
    spec <- solarModel_spec(place, ..., control_model = control_model)
    # Initialize the model
    model <- solarModel$new(spec)
    # Model fit
    model$fit()
    # Model names
    model_name <- paste0("ARMA(", model$ARMA$order[1], ", ", model$ARMA$order[2],")-GARCH(",model$GARCH$order[1], ", ", model$GARCH$order[2],")")
    # Evaluate the tests
    tests <- solarModel_tests(model, lag.max = lag.max, ci = ci)
    condition_autocorr <- sum(tests$autocorr$result[c(2,3,4,5)] == "passed") != 4
    if (condition_autocorr) {
      cli::cli_alert_danger(paste0("Autocorrelation test not passed for model: ", model_name))
      if (control_model$mean.model$maOrder == control_model$mean.model$arOrder) {
        control_model$mean.model$arOrder <- control_model$mean.model$arOrder + 1
      } else if (control_model$mean.model$maOrder < control_model$mean.model$arOrder){
        #message("Increase MA component by 1")
        control_model$mean.model$maOrder <- control_model$mean.model$maOrder + 1
      }
    } else {
      cli::cli_alert_success(paste0("All autocorrelation test for model ", model_name, " passed!"))
    }
    condition_autocorr_var <- sum(tests$autocorr$result[c(9, 10)] == "passed") != 2 & !condition_autocorr
    if (condition_autocorr_var) {
      control_model$garch_variance <- TRUE
      cli::cli_alert_danger(paste0("Autocorrelation test in variance not passed for model: ", model_name))
      if (narch == ngarch) {
        #message("Increase ARCH component by 1")
        ngarch <- ifelse(ngarch == 0, 1, ngarch)
        control_model$variance.model <- rugarch::ugarchspec(variance.model = list(garchOrder = c(narch + 1, ngarch), variance.targeting = 1),
                                                      mean.model = list(armaOrder = c(0, 0), include.mean = FALSE))
        narch <- narch + 1
      } else if (ngarch < narch){
        #message("Increase GARCH component by 1")
        control_model$variance.model <- rugarch::ugarchspec(variance.model = list(garchOrder = c(narch, ngarch +1), variance.targeting = 1),
                                                      mean.model = list(armaOrder = c(0, 0), include.mean = FALSE))
        ngarch <- ngarch + 1
      }
    } else {
      cli::cli_alert_success(paste0("All autocorrelation test in variance for model ", model_name, " passed!"))
    }
    condition <- condition_autocorr | condition_autocorr_var
  }
  return(model)
}
