#' Solar Model in R6 Class
#'
#' @description
#' The `solarModel` class is an implementation of a comprehensive solar model that includes fitting seasonal models,
#' detecting outliers, performing transformations, and applying time-series models such as AR and GARCH. This model
#' is specifically designed to predict solar radiation data, and it uses seasonal and Gaussian Mixture models to
#' capture the underlying data behavior.
#'
#' @details
#' The `solarModel` class allows for the step-by-step fitting and transformation of solar radiation data, from clear sky
#' models to GARCH models for residual analysis. It utilizes various private and public methods to fit the seasonal
#' clearsky model, compute risk drivers, detect outliers, and apply time-series models.
#'
#' @examples
#'
#' # Model specification
#' spec <- solarModel_spec$new()
#' spec$set_mean.model(arOrder = 1, maOrder = 1)
#' spec$specification("Bologna")
#' spec
#' # Model fit
#' Bologna <- solarModel$new(spec)
#' Bologna$fit()
#'
#' Bologna <- solarModel_QMLE(Bologna)
#' solarOption_model(Bologna, Bologna$moments$conditional[1:365,])
#'
#' # save(spec, file = "data/Bologna.RData")
#'
#' # Extract and update the parameters
#' model <- Bologna$clone(TRUE)
#' params <- model$coefficients
#' model$update(params)
#' model$filter()
#'
#' # Fit a model with the realized clear sky
#' spec$control$stochastic_clearsky <- TRUE
#' # Initialize a new model
#' model <- solarModel$new(spec)
#' # Model fit
#' model$fit()
#'
#' # Fit a model for the clearsky
#' spec_Ct <- spec
#' spec_Ct$control$stochastic_clearsky <- FALSE
#' spec_Ct$target <- "clearsky"
#' # Initialize a new model
#' model <- solarModel$new(spec)
#' # Model fit
#' model$fit()
#'
#' @rdname solarModel
#' @name solarModel
#' @keywords solarModel
#' @note Version 1.0.4
#' @export
solarModel <- R6::R6Class("solarModel",
                          # ====================================================================================================== #
                          #                                             Public slots
                          # ====================================================================================================== #
                          public = list(
                            #' @field interpolated Logical values, when TRUE the time-series and the parameters are interpolated.
                            interpolated = FALSE,
                            #' @description
                            #' Initialize a `solarModel`
                            #' @param spec an object with class `solarModelSpec`. See the function \code{\link{solarModel_spec}} for details.
                            initialize = function(spec){
                              # Seasonal data by month and day for an year with 366 days
                              dates <- seq.Date(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "1 day")
                              seasonal_data <- dplyr::tibble(date = dates)
                              seasonal_data <- dplyr::mutate(seasonal_data,
                                                             Month = lubridate::month(date),
                                                             Day = lubridate::day(date),
                                                             n = number_of_day(date))
                              seasonal_data <- dplyr::select(seasonal_data, -date)
                              seasonal_data <- dplyr::arrange(seasonal_data, Month, Day)
                              # **************************************************** #
                              # Private components
                              private$..data <- dplyr::select(spec$data, -H0)
                              private$..data$loglik <- 0
                              private$..seasonal_data <- seasonal_data
                              private$..monthly_data <- dplyr::tibble(Month = 1:12, Yt_tilde_uncond = 0, sigma_uncond = 1)
                              private$..loglik <- NA
                              private$..spec <- spec$clone(TRUE)
                            },
                            # ***************************************************************************** #
                            #' @description
                            #' Initialize and fit a \code{\link{solarModel}} object given the specification contained in `$control`.
                            fit = function(){
                              # 1) Clearsky
                              self$fit_seasonal_model_Ct()
                              # 2) Risk-driver
                              self$compute_risk_drivers()
                              # 3) Solar transform
                              self$fit_transform()
                              # 4) Seasonal mean
                              self$fit_seasonal_mean()
                              #    Center the mean to be exactly equal to zero
                              self$fit_monthly_mean()
                              # 5) ARMA model
                              self$fit_mean_model()
                              # 6) Seasonal variance
                              self$fit_seasonal_variance()
                              # 7) GARCH variance
                              self$fit_variance_model()
                              # 8) Mixture model
                              self$fit_mixture_model()
                              self$update_classification()
                              # Update the moments
                              self$update_moments()
                              # Update the log-likelihoods
                              self$update_logLik()
                            },
                            #' @description
                            #' Initialize and fit a \code{\link{seasonalClearsky}} model given the specification contained in `$control`.
                            fit_seasonal_model_Ct = function(){
                              # Arguments
                              data <- dplyr::filter(private$..data, isTrain & weights != 0)
                              # Reference variables
                              x = data[[self$spec$target]]
                              dates = data[["date"]]
                              lat = self$spec$coords$lat
                              clearsky = data[["clearsky"]]
                              # **************************************************** #
                              # Initialize a seasonal model for clear sky radiation
                              private$..spec$seasonal_model_Ct$fit(x = x,
                                                                   dates = dates,
                                                                   lat = lat,
                                                                   clearsky = clearsky)
                              # Optimize the model when method is not OLS
                              private$..spec$.__enclos_env__$private$..seasonal_model_Ct <- clearsky_optimizer(x = x,
                                                                                                        dates = dates,
                                                                                                        clearsky = clearsky,
                                                                                                        seasonal_model_Ct = private$..spec$seasonal_model_Ct)

                              # **************************************************** #
                              # Add seasonal clear sky max to seasonal data
                              private$..data[["Ct"]] <- self$spec$seasonal_model_Ct$predict(newdata = self$data)
                            },
                            #' @description
                            #' Compute the risk drivers and impute the observation that are greater or equal to the clear sky level.
                            compute_risk_drivers = function(){
                              # Arguments
                              target <- self$spec$target
                              control <- self$spec
                              data <- self$data
                              transform <- self$spec$transform
                              # **************************************************** #
                              if (control$stochastic_clearsky) {
                                # Risk driver
                                data$Xt <- transform$X(data[[target]], data$clearsky)
                                # Detect and impute outliers
                                outliers <- clearsky_outliers(data[["GHI"]], data$clearsky, date = data$date, quiet = control$quiet)
                                # Update GHI
                                data[["GHI"]] <- outliers$x
                                # Risk driver
                                data[["Xt"]] <- transform$X(data[["GHI"]], data$clearsky)
                              } else {
                                # Detect and impute outliers
                                outliers <- clearsky_outliers(data[["GHI"]], data$Ct, date = data$date, quiet = control$quiet)
                                # Update GHI
                                data[["GHI"]] <- outliers$x
                                # Add computed risk driver in data
                                data[["Xt"]] <- transform$X(data[["GHI"]], data$Ct)
                              }
                              # **************************************************** #
                              # Add computed risk driver in data
                              private$..data[["Xt"]] <- data$Xt
                              # Store outliers data
                              private$..spec$outliers <- outliers
                            },
                            #' @description
                            #' Fit the parameters of the \code{\link{solarTransform}} object.
                            #' @param detect_outliers Logical, when true outliers are detected
                            fit_transform = function(detect_outliers = TRUE){
                              # Arguments
                              control <- self$spec$transform$control
                              # Full sample
                              data <- private$..data
                              # Train data
                              data_train <- dplyr::filter(data, isTrain & weights != 0)
                              outliers <- private$..spec$outliers
                              # **************************************************** #
                              # Transformation parameters
                              params <- private$..spec$transform$fit(data_train$Xt, control$threshold, control$min_pos, control$max_pos)
                              # Update transform parameters
                              private$..spec$transform$update(params$alpha, params$beta)
                              # Rescale the minimum to avoid that extreme values
                              # breaks the ARMA-GARCH routines (especially the minimum)
                              if (detect_outliers || is.null(outliers$index_type$transform)) {
                                idx_Xt_min <- which(data$Xt <= params$Xt_min)
                                idx_Xt_max <- which(data$Xt >= params$Xt_max)
                                # Store the index of the imputed values
                                outliers$index_type$transform <- c(idx_Xt_min = idx_Xt_min, idx_Xt_max = idx_Xt_max)
                                # Update outliers index and dates
                                outliers$index <- unique(c(outliers$index, outliers$index_type$transform))
                                outliers$date <- data$date[outliers$index]
                                # Store the updated outliers data
                                private$..spec$outliers <- outliers
                              } else {
                                idx_Xt_min <- outliers$index_type$transform[1]
                                idx_Xt_max <- outliers$index_type$transform[2]
                              }
                              # Rescale the minimum and maximum
                              data$Xt[idx_Xt_min] <- params$Xt_min * (1 + control$delta)
                              data$Xt[idx_Xt_max] <- params$Xt_max * (1 - control$delta)
                              # **************************************************** #
                              # Compute the transformed variable
                              data$Yt <- self$spec$transform$Y(self$spec$transform$X_prime(data[["Xt"]]))
                              # Add Yt to private data
                              private$..data[["Yt"]] <- data$Yt
                              # Store the updated outliers
                              private$outliers <- outliers
                            },
                            #' @description
                            #' Fit a \code{\link{seasonalModel}} the transformed variable (`Yt`) and compute deseasonalized series (`Yt_tilde`).
                            fit_seasonal_mean = function(){
                              # Arguments
                              target <- self$spec$target
                              data <- private$..data
                              # **************************************************** #
                              # Train data
                              data_train <- dplyr::filter(data, isTrain & weights != 0)
                              # Fit the seasonal model
                              private$..spec$seasonal.mean$fit(data_train)
                              # **************************************************** #
                              # Fitted seasonal mean
                              data$Yt_bar <- self$spec$seasonal.mean$predict(newdata = data)
                              # Fitted deseasonalized transformed variable
                              data$Yt_tilde <- data$Yt - data$Yt_bar
                              # Fitted seasonal mean of target variable
                              target_bar <- paste0(target, "_bar")
                              data[[target_bar]] <- self$spec$transform$iRY(data$Yt_bar, data$Ct)
                              # **************************************************** #
                              # Update private components
                              private$..data[["Yt_bar"]]   <- data$Yt_bar
                              private$..data[["Yt_tilde"]] <- data$Yt_tilde
                              private$..data[[target_bar]] <- data[[target_bar]]
                            },
                            #' @description
                            #' Correct the deseasonalized series (`Yt_tilde`) by subtracting its monthly mean (`Yt_tilde_uncond`).
                            fit_monthly_mean = function(){
                              # Unconditional mean for Yt_tilde
                              if (self$spec$seasonal.mean$control$monthly.mean) {
                                data <- private$..data
                                # Train data
                                train_data <- dplyr::filter(data, isTrain & weights != 0)
                                # Compute monthly unconditional mean
                                monthly_mean <- train_data %>%
                                  dplyr::group_by(Month) %>%
                                  dplyr::summarise(Yt_tilde_uncond = mean(Yt_tilde)) %>%
                                  dplyr::ungroup()
                                # Add unconditional mean to the dataset
                                data <- dplyr::left_join(data, monthly_mean, by = c("Month"))
                                # Update Yt_tilde
                                data$Yt_tilde <- data$Yt_tilde - data$Yt_tilde_uncond
                                # **************************************************** #
                                # Update private components
                                # Add unconditional mean to monthly data
                                private$..monthly_data[["Yt_tilde_uncond"]] <- monthly_mean$Yt_tilde_uncond
                                # Add Yt_tilde to data
                                private$..data[["Yt_tilde"]] <- data$Yt_tilde
                              }
                            },
                            #' @description
                            #' Fit an ARMA model (`Yt_tilde`) and compute the residuals (`eps`).
                            fit_mean_model = function(){
                              # Arguments
                              data <- private$..data
                              # **************************************************** #
                              # Train data
                              data_train <- dplyr::filter(data, isTrain)# & weights != 0)
                              # Fit ARMA model
                              private$..spec$mean.model$fit(data_train$Yt_tilde)
                              # **************************************************** #
                              # Fitted Yt_tilde
                              data$Yt_tilde_hat <- self$spec$mean.model$filter(data$Yt_tilde)
                              # Fitted residuals
                              data$eps <- data$Yt_tilde - data$Yt_tilde_hat
                              # **************************************************** #
                              # Update private components
                              private$..data[["Yt_tilde_hat"]] <- data$Yt_tilde_hat
                              private$..data[["eps"]] <- data$eps
                            },
                            #' @description
                            #' Fit a \code{\link{seasonalModel}} on AR squared residuals (`eps`) and compute deseasonalized residuals `eps_tilde`.
                            fit_seasonal_variance = function(){
                              # Dataset
                              data <- private$..data
                              # **************************************************** #
                              # Train data
                              data_train <- dplyr::filter(data, isTrain & weights != 0)
                              # Squared residuals
                              data_train$eps2 <- data_train$eps^2
                              # Fit seasonal model
                              private$..spec$seasonal.variance$fit(data = data_train)
                              # **************************************************** #
                              # Fitted seasonal standard deviation
                              data$sigma_bar <- sqrt(self$spec$seasonal.variance$predict(newdata = data))
                              # Fitted seasonal-standardized residuals
                              data$eps_tilde <- data$eps / data$sigma_bar
                              # **************************************************** #
                              # Fitted seasonal standard deviation
                              private$..data[["sigma_bar"]] <- data$sigma_bar
                              # Compute standardized residuals
                              private$..data[["eps_tilde"]] <- data$eps_tilde
                              # **************************************************** #
                              # Compute monthly corrective variance to ensure unitary monthly variance
                              self$fit_monthly_variance()
                              # Correction of the parameters to ensure unitary variance
                              self$correct_seasonal_variance()
                            },
                            #' @description
                            #' Correct the standardized series (`eps_tilde`) by subtracting its monthly mean (`sigma_uncond`).
                            fit_monthly_variance = function(){
                              # Condition
                              condition <- self$spec$seasonal.variance$control$monthly.mean
                              # Unconditional variance for eps_tilde
                              if (condition) {
                                # Full dataset
                                data <- private$..data
                                # Train data
                                train_data <- dplyr::filter(data, isTrain & weights != 0)
                                # Compute monthly unconditional mean
                                monthly_data <- train_data %>%
                                  dplyr::group_by(Month) %>%
                                  dplyr::summarise(sigma_uncond = sd(eps_tilde, na.rm = TRUE)) %>%
                                  dplyr::ungroup()
                                # **************************************************** #
                                # Add unconditional variance to monthly data
                                private$..monthly_data[["sigma_uncond"]] <- monthly_data$sigma_uncond
                                # Updated standardized residuals
                                private$..data[["eps_tilde"]] <- data$eps / (data$sigma_bar * self$data$sigma_uncond)
                              }
                            },
                            #' @description
                            #' Correct the parameters of the seasonal variance to ensure a unitary variance
                            correct_seasonal_variance = function(){
                              # Condition
                              condition <- self$spec$seasonal.variance$control$monthly.mean
                              # Correction to ensure unitary variance
                              if (condition) {
                                # Full dataset
                                data <- self$data
                                # Train data
                                train_data <- dplyr::filter(data, isTrain & weights != 0)
                                # Store correction factor
                                correction <- var(data_train$eps / (data_train$sigma_bar * data_train$sigma_uncond))
                                # Correct parameters to ensure unitary variance
                                std.errors <- self$spec$seasonal.variance$std.errors
                                # Update correction
                                private$..spec$.__enclos_env__$private$..seasonal.variance$control[["sigma20"]] <- correction
                                private$..spec$seasonal.variance$update(self$spec$seasonal.variance$coefficients * correction)
                                private$..spec$seasonal.variance$update_std.errors(std.errors * correction)
                                # **************************************************** #
                                # Update private components
                                private$..data[["sigma_bar"]] <- sqrt(self$spec$seasonal_variance$predict(newdata = data))
                                private$..data[["eps_tilde"]] <- data$eps / (private$..data[["sigma_bar"]] * data$sigma_uncond)
                              }
                            },
                            #' @description
                            #' Fit a `GARCH` model on the deseasonalized residuals (`eps_tilde`).
                            #' Compute the standardized (`u`) and monthly deseasonalized residuals (`u_tilde`).
                            fit_variance_model = function(){
                              # Arguments
                              control <- self$spec
                              data <- private$..data
                              # **************************************************** #
                              # Train data
                              data_train <- dplyr::filter(data, isTrain & weights != 0)
                              # Control for garch variance
                              if (sum(self$spec$variance.model$order) > 0) {
                                # Fit the model
                                private$..spec$variance.model$fit(data_train$eps_tilde, data_train$weights)
                                # Fitted GARCH std. deviation
                                data$sigma <- sqrt(self$spec$variance.model$filter(data$eps_tilde))
                                # Fitted GARCH standardized residuals
                                data$u_tilde <- data$eps_tilde / data$sigma
                              } else {
                                # Fitted std. deviation
                                data$sigma <- 1
                                # Not compute standardized residuals
                                data$u_tilde <- data$eps_tilde
                              }
                              # **************************************************** #
                              # Update private components
                              private$..data[["sigma"]] <- data$sigma
                              private$..data[["u_tilde"]] <- data$u_tilde
                            },
                            #' @description
                            #' Initialize and fit a `solarMixture` object.
                            fit_mixture_model = function(){
                              # Arguments
                              control <- self$spec$mixture.model$control
                              outliers <- private$outliers
                              # Train data
                              data_train <- dplyr::filter(private$..data, isTrain & weights != 0)
                              # **************************************************** #
                              # Match moments
                              if (control$match.expectation | control$match.variance) {
                                # Target empirical moments
                                target <- data_train %>%
                                  group_by(Month) %>%
                                  summarise(mu_target = mean(u_tilde),
                                            var_target = var(u_tilde))
                                # Match or not empirical moments
                                if (!control$match.empiric) {
                                  target$mu_target <- 0
                                  target$var_target <- 1
                                }
                                # Match or not expectation
                                if (!control$match.expectation) {
                                  target$mu_target <- NA
                                }
                                # Match or not variance
                                if (!control$match.variance) {
                                  target$var_target <- NA
                                }
                                # Store Gaussian Mixture parameters
                                private$..spec$mixture.model$fit(x = data_train$u_tilde, date = data_train$date, weights = data_train$weights,
                                                                 mu_target = target$mu_target, var_target = target$var_target)
                              } else {
                                # Store Gaussian Mixture parameters
                                private$..spec$mixture.model$fit(x = data_train$u_tilde, date = data_train$date, weights = data_train$weights)
                              }
                            },
                            #' @description
                            #' Update the parameters inside object
                            #' @param params List of parameters. See the slot `$coefficients` for a template.
                            update = function(params){
                              # Update the parameters
                              private$..spec$update(params)
                              # Set Log-likelihood to NA
                              private$..loglik <- NA
                            },
                            #' @description
                            #' Update the moments inside object
                            update_moments = function(){
                              # Update conditional moments
                              private$..moments$conditional <- solarMoments_conditional(self$data, theta = 0)
                            },
                            #' @description
                            #' Update the log-likelihood inside object
                            update_logLik = function(){
                              # Compute the log-likelihoods
                              log.likelihoods <- self$logLik()
                              # **************************************************** #
                              # Update the log-likelihoods
                              private$..data[["loglik"]] <- log.likelihoods
                              # Update total log-likelihood
                              private$..loglik <- sum(log.likelihoods)
                            },
                            #' @description
                            #' Update the clear sky and risk drivers
                            update_risk_drivers = function(){
                              # Update clear sky
                              private$..data[["Ct"]] <- self$spec$seasonal_model_Ct$predict(newdata = self$data)
                              # Update risk driver
                              self$compute_risk_drivers()
                              # Update solar transform and compute Yt
                              self$fit_transform(detect_outliers = FALSE)
                            },
                            #' @description
                            #' Update the classification of the Bernoulli random variable.
                            #' @param filter Logical, when `TRUE` before the classification will be runned the command
                            #' `filter` to update the mixture classification.
                            update_classification = function(filter = FALSE){
                              # Update mixture classification
                              if (filter) {
                                private$..spec$mixture.model$filter()
                              }
                              # Full data
                              data <- private$..data
                              # **************************************************** #
                              # Ensure that no columns with B name are included
                              data <- data[, !(colnames(data) %in% c("B", "z1","z2"))]
                              # Classify the series
                              df_m <- list()
                              for(nmonth in 1:12){
                                df_m[[nmonth]] <- dplyr::filter(data, Month == nmonth)
                                df_nm <- private$..spec$mixture.model$model[[nmonth]]$classify(df_m[[nmonth]]$u_tilde)
                                df_m[[nmonth]]$B <- df_nm$B1
                                df_m[[nmonth]]$z1 <- df_nm$z1
                                df_m[[nmonth]]$z2 <- df_nm$z2
                              }
                              # Add the data
                              df_m <- dplyr::select(dplyr::bind_rows(df_m), date, B, z1, z2)
                              data <- dplyr::left_join(data, df_m, by = "date")
                              # **************************************************** #
                              # Private components
                              # Update fitted series of Bernoulli
                              private$..data[["B"]] <- data$B
                              # Update fitted series of Gaussians components
                              private$..data[["z1"]] <- data$z1
                              private$..data[["z2"]] <- data$z2
                            },
                            # ***************************************************************************** #
                            #' @description
                            #' Filter the time series when new parameters are supplied in the method `$update(params)`.
                            #' @param fit Logical, when `FALSE`, if in the model's specification, the monthly mean and variances will be re estimated and the seasonal variance corrected
                            #' such that the total variance of the deseasonalized residuals is zero.
                            #' @return Update the slots `$data`, `$seasonal_data`, `$monthly_data`
                            filter = function(fit = FALSE){
                              # Arguments
                              target <- self$spec$target
                              # Full dataset (old)
                              data <- self$data
                              # **************************************************** #
                              # Update seasonal clear sky to seasonal data
                              private$..data[["Ct"]] <- self$spec$seasonal_model_Ct$predict(newdata = data)
                              # Add computed risk driver in data
                              self$compute_risk_drivers()
                              # Compute the transformed variable
                              self$fit_transform(detect_outliers = FALSE)
                              # Update seasonal mean of Yt
                              private$..data[["Yt_bar"]] <- self$spec$seasonal.mean$predict(newdata = data)
                              # Update seasonal mean of target variable
                              private$..data[[paste0(target, "_bar")]] <- self$spec$transform$iRY(private$..data[["Yt_bar"]], private$..data[["Ct"]])
                              # Update Yt_tilde
                              private$..data[["Yt_tilde"]] <- private$..data[["Yt"]] - private$..data[["Yt_bar"]]
                              # Fit the corrective mean
                              if (fit) {
                                self$fit_monthly_mean()
                              }
                              # Update Yt_tilde_hat
                              private$..data[["Yt_tilde_hat"]] <- self$spec$mean.model$filter(private$..data[["Yt_tilde"]])
                              # Update ARMA residuals
                              private$..data[["eps"]] <- private$..data[["Yt_tilde"]] - private$..data[["Yt_tilde_hat"]]
                              # **************************************************** #
                              # Update seasonal std. deviation
                              private$..data[["sigma_bar"]] <- sqrt(self$spec$seasonal.variance$predict(newdata = data))
                              # Update standardized residuals
                              private$..data[["eps_tilde"]] <- private$..data[["eps"]] / private$..data[["sigma_bar"]]
                              if (fit) {
                                # Compute corrective monthly variance
                                self$fit_monthly_variance()
                                # Correct the seasonal variance
                                self$correct_seasonal_variance()
                              }
                              # **************************************************** #
                              # Update Garch standard deviation
                              if (sum(self$spec$variance.model$order) != 0) {
                                private$..data[["sigma"]] <- sqrt(self$spec$variance.model$filter(private$..data[["eps_tilde"]]))
                              } else {
                                private$..data[["sigma"]] <- 1
                              }
                              # Fitted standardized residuals
                              private$..data[["u_tilde"]] <- private$..data[["eps_tilde"]] / private$..data[["sigma"]]
                            },
                            #' @description
                            #' Compute the log-likelihood of the model and update the slot `$loglik`.
                            #' @param moments Dataset containing the moments to use for computation.
                            #' @param target Character. Target variable to use "Yt" or "GHI".
                            #' @param quasi Logical, when `TRUE` is computed the pseudo-likelihood with Gaussian link.
                            logLik = function(moments, target = "Yt", quasi = FALSE){
                              # Default argument
                              if (missing(moments)) {
                                moments <- self$moments$conditional
                              }
                              # Add Yt and weights
                              moments <- dplyr::left_join(moments, private$..data[, c("date", target, "weights", "isTrain")], by = "date")
                              # Compute the quasi and true log-likelihood
                              moments$loglik <- 0
                              # Density: solar radiation
                              if (target == "GHI") {
                                for(i in 1:nrow(moments)){
                                  df_n <- moments[i,]

                                  if (df_n$weights == 0 & df_n$isTrain) {
                                    moments$loglik[i] <- 0
                                    next
                                  }
                                  if (quasi) {
                                    # Normal density
                                    pdf_quasi <- function(x) dnorm(x, df_n$e_Yt, df_n$sd_Yt)
                                    # Quasi Log-likelihood
                                    moments$loglik[i] <- log(dsolarGHI(df_n$GHI, df_n$Ct, df_n$alpha, df_n$beta, pdf_quasi, link = self$spec$transform$link))
                                  } else {
                                    # True mixture density
                                    pdf_Yt <- function(x) dmixnorm(x, mean = c(df_n$M_Y1, df_n$M_Y0), sd = c(df_n$S_Y1, df_n$S_Y0), alpha = c(df_n$p1, 1-df_n$p1))
                                    # True Log-likelihood
                                    moments$loglik[i] <- log(dsolarGHI(df_n$GHI, df_n$Ct, df_n$alpha, df_n$beta, pdf_Yt, link = self$spec$transform$link))
                                  }
                                }
                              } else {
                                if (quasi){
                                  # Standardize the time series
                                  moments$z <- (moments$Yt - moments$e_Yt) / moments$sd_Yt
                                  # Pseudo log-likelihood
                                  moments$loglik <- log(dnorm(moments$z) / moments$sd_Yt)
                                } else {
                                  # Standardize the time series into its components
                                  moments$z1 <- (moments$Yt - moments$M_Y1) / moments$S_Y1
                                  moments$z0 <- (moments$Yt - moments$M_Y0) / moments$S_Y0
                                  # True mixture log-likelihood
                                  moments$loglik <- log(dnorm(moments$z1) / moments$S_Y1 * moments$p1 + dnorm(moments$z0) / moments$S_Y0 * (1 - moments$p1))
                                  sum(moments$loglik)
                                  nrow(moments)
                                }
                              }
                              # **************************************************** #
                              return(moments$loglik)
                            },
                            #' @description
                            #' Print method for `solarModel` class.
                            print = function(){
                              # Complete data specifications
                              data <- self$spec$dates$data
                              # Train data specifications
                              train <- self$spec$dates$train
                              train$perc <- format(train$perc*100, digits = 3)
                              # Test data specifications
                              test <- self$spec$dates$test
                              test$perc <- format(test$perc*100, digits = 3)
                              coords <- self$spec$location
                              # **********************************************************************
                              msg0 <- paste0("--------------------- ", "solarModel", " (", self$spec$place, ") ", "---------------------")
                              cat(msg0, "\n")
                              cat(" Model:", self$spec$model_name, "\n",
                                  "Target:", self$spec$target, "\n",
                                  "Lat:", coords$lat, "Lon:", coords$lon, "Alt:", coords$alt, "\n",
                                  "Dates:", as.character(data$from), "-", as.character(data$to), "\n",
                                  "Observations:", data$nobs, "\n",
                                  "Interpolated:", self$interpolated, "\n",
                                  "Log-Likelihood:", format(self$loglik, digits = 8), "\n",
                                  paste0(rep("*", length(strsplit(msg0, "")[[1]])-1), collapse = ""), "\n",
                                  "Train (~", train$perc, "%):", as.character(train$from), "-", as.character(train$to),
                                  paste0("(", train$nobs, " points)"), "\n",
                                  "Test  (~", test$perc, "%):", as.character(test$from), "-", as.character(test$to),
                                  paste0("(", test$nobs, " points)"), "\n",
                                  paste0(rep("*", length(strsplit(msg0, "")[[1]])-1), collapse = ""), "\n",
                                  "Version:", private$version, "\n",
                                  paste0(rep("-", length(strsplit(msg0, "")[[1]])-1), collapse = ""), "\n")
                            }
                          ),
                          # ====================================================================================================== #
                          #                                             Private slots
                          # ====================================================================================================== #
                          private = list(
                            version = "1.0.4",
                            ..spec = NA,
                            ..data = NA,
                            ..seasonal_data = NA,
                            ..monthly_data = NA,
                            outliers = NA,
                            ..loglik = NA,
                            ..moments = list(conditional = NA)
                          ),
                          # ====================================================================================================== #
                          #                                             Active slots
                          # ====================================================================================================== #
                          active = list(
                            #' @field data A data frame with the fitted data, and the seasonal and monthly parameters.
                            data = function(){
                              # Seasonal data
                              seasonal_data <- dplyr::select(self$seasonal_data, -n)
                              dplyr::left_join(private$..data, seasonal_data, by = c("Month", "Day"))
                            },
                            #' @field seasonal_data A data frame containing seasonal and monthly parameters.
                            seasonal_data = function(){
                              dplyr::left_join(private$..seasonal_data, self$monthly_data, by = "Month")
                            },
                            #' @field monthly_data A data frame that contains monthly parameters.
                            monthly_data = function(){
                              if (!purrr::is_empty(self$spec$mixture.model)) {
                                dplyr::left_join(private$..monthly_data, self$spec$mixture.model$coefficients, by = "Month")
                              } else {
                                private$..monthly_data
                              }
                            },
                            #' @field loglik The log-likelihood computed on train data.
                            loglik = function(){
                              private$..loglik
                            },
                            #' @field spec A list with the specification that govern the behavior of the model's fitting process.
                            spec = function(){
                              private$..spec
                            },
                            #' @field moments Get a list containing the conditional moments.
                            moments = function(){
                              moments <- private$..moments
                              # Extra data to add
                              data_extra <- dplyr::select(self$data, date, p1, GHI_bar, Ct)
                              data_extra <- dplyr::mutate(data_extra, alpha = self$spec$transform$alpha, beta = self$spec$transform$beta)
                              # Conditional moments
                              if (length(moments$conditional) != 1) {
                                moments$conditional <- dplyr::left_join(moments$conditional, data_extra, by = c("date"))
                              }
                              return(moments)
                            },
                            #' @field coefficients Get the model parameters as a named list.
                            coefficients = function(){
                              private$..spec$coefficients
                            }
                          )
                        )

