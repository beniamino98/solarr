#' Solar Model in R6 Class
#'
#' The `solarModel` class is an implementation of a comprehensive solar model that includes fitting seasonal models,
#' detecting outliers, performing transformations, and applying time-series models such as AR and GARCH. This model
#' is specifically designed to predict solar radiation data, and it uses seasonal and Gaussian Mixture models to
#' capture the underlying data behavior.
#'
#' @description
#' The `solarModel` class allows for the step-by-step fitting and transformation of solar radiation data, from clear sky
#' models to GARCH models for residual analysis. It utilizes various private and public methods to fit the seasonal
#' clearsky model, compute risk drivers, detect outliers, and apply time-series models.
#'
#' @examples
#' # Control list
#' control <- control_solarModel(outliers_quantile = 0.005)
#' # Model specification
#' spec <- solarModel_spec("Bologna", from="2005-01-01", to="2022-01-01", control_model = control)
#' Bologna <- solarModel$new(spec)
#' # Model fit
#' Bologna$fit()
#'
#' # Extract and update the parameters
#' params <- sm$parameters
#' sm$update(params)
#' sm$filter()
#'
#' # Fit a model with the realized clear sky
#' spec$control$stochastic_clearsky <- TRUE
#' # Initialize a new model
#' model <- solarModel$new(spec)
#' #' # Model fit
#' model$fit()
#'
#' # Fit a model for the clearsky
#' spec_Ct <- spec
#' spec_Ct$control$stochastic_clearsky <- FALSE
#' spec_Ct$target <- "clearsky"
#' # Initialize a new model
#' model <- solarModel$new(spec)
#' #' # Model fit
#' model$fit()
#' @rdname solarModel
#' @name solarModel
#' @export
solarModel <- R6::R6Class("solarModel",
                            public = list(
                              #' @field place Character, optional name of the location considered.
                              place = NA,
                              #' @field target Character, name of the target variable to model. Can be `GHI` or `clearsky`.
                              target = NA,
                              #' @field dates List, with the range of dates used in the model.
                              dates = NA,
                              #' @field coords A data frame with coordinates of the location considered.
                              coords = NA,
                              #' @description
                              #' Initialize a `solarModel`
                              #' @param spec an object with class `solarModelSpec`. See the function \code{\link{solarModel_spec}} for details.
                              initialize = function(spec){
                                # **************************************************** #
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
                                # Public components
                                self$place <- spec$place
                                self$target <- spec$target
                                self$dates <- spec$dates
                                self$coords <- spec$coord
                                # **************************************************** #
                                # Private components
                                private$..data <- dplyr::select(spec$data, -H0)
                                private$..seasonal_data <- seasonal_data
                                private$..monthly_data <- dplyr::tibble(Month = 1:12, Yt_tilde_uncond = 0, sigma_m = 1)
                                private$..loglik <- dplyr::tibble(Month = 1:12, loglik = 0)
                                # private$dates <- spec$dates
                                private$..control <- spec$control
                                private$..transform <- solarTransform$new(0, 1)
                              },
                              #' @method fit solarModel
                              #' @description
                              #' Fit the model given the specification contained in `control`.
                              fit = function(){
                                self$fit_clearsky_model()
                                self$compute_risk_drivers()
                                self$fit_solar_transform()
                                self$detect_outliers_Yt()
                                self$fit_seasonal_mean()
                                self$corrective_monthly_mean()
                                self$fit_AR_model()
                                self$fit_seasonal_variance()
                                self$fit_GARCH_model()
                                self$corrective_monthly_variance()
                                self$fit_mixture_model()
                                # Compute moments
                                self$conditional_moments()
                                self$unconditional_moments()
                                # Compute log-likelihood
                                self$logLik()
                              },
                              #' @description
                              #' Fit a `seasonalClearsky` given a certain specification
                              fit_clearsky_model = function(){
                                # Private arguments
                                data <- private$..data
                                # **************************************************** #
                                # Initialize a seasonal model for clear sky radiation
                                seasonal_model_Ct <- seasonalClearsky$new(control = self$control$clearsky)
                                # Fit the seasonal model
                                seasonal_model_Ct$fit(data[[self$target]], data[["date"]], self$coords$lat, data[["clearsky"]])
                                # **************************************************** #
                                # Add seasonal clearsky max to seasonal data
                                private$..seasonal_data$Ct <- seasonal_model_Ct$predict(private$..seasonal_data$n)
                                # Store the model
                                private$..seasonal_model_Ct <- seasonal_model_Ct$clone(deep = TRUE)
                              },
                              #' @description
                              #' Compute the risk drivers and detect outliers with respect to clearsky
                              compute_risk_drivers = function(){
                                # Self Arguments
                                target <- self$target
                                control <- self$control
                                # Private arguments
                                data <- private$..data
                                # **************************************************** #
                                # Compute risk drivers
                                if (control$stochastic_clearsky) {
                                  data$clearsky <- data$clearsky*1.01
                                  # Detect and impute outliers
                                  outliers <- clearsky_outliers(data[[target]], data$clearsky, date = data$date, quiet = control$quiet)
                                  # Update GHI
                                  data[[target]] <- outliers$x
                                  data$Xt <- self$transform$iGHI(data[[target]], data$clearsky)
                                } else {
                                  # Predict seasonal clearsky
                                  data$Ct <- self$seasonal_model_Ct$predict(data$n)
                                  # Detect and impute outliers
                                  outliers <- clearsky_outliers(data[[target]], data$Ct, date = data$date, quiet = control$quiet)
                                  # Update GHI
                                  data[[target]] <- outliers$x
                                  data$Xt <- self$transform$iGHI(data[[target]], data$Ct)
                                }
                                # **************************************************** #
                                # Update target in data
                                private$..data[[target]] <- data[[target]]
                                # Update clearsky in data
                                private$..data[["clearsky"]] <- data$clearsky
                                # Add computed risk driver in data
                                private$..data[["Xt"]] <- data$Xt
                                # Store outliers data
                                private$outliers <- outliers
                              },
                              #' @description
                              #' Fit the parameters of the solar tranform
                              fit_solar_transform = function(){
                                # Self Arguments
                                control <- self$control
                                st <- self$transform
                                # Private arguments
                                data <- private$..data
                                outliers <- private$outliers
                                # **************************************************** #
                                # Transformation parameters
                                params <- st$parameters(data$Xt, control$threshold)
                                st$update(params$alpha, params$beta)
                                # Remove minimum and maximum of Xt from computations
                                outliers$index_type$transform <- c(which.min(data$Xt), which.max(data$Xt))
                                # Update outliers index and dates
                                outliers$index <- unique(c(outliers$index, outliers$index_type$transform))
                                outliers$date <- data$date[outliers$index]
                                # **************************************************** #
                                # Update solar transform
                                private$..transform <- st$clone(deep = TRUE)
                                # Update outliers
                                private$outliers <- outliers
                              },
                              #' @description
                              #' Detect the outliers that will be excluded from computations
                              detect_outliers_Yt = function(){
                                # Private arguments
                                data <- private$..data
                                outliers <- private$outliers
                                # **************************************************** #
                                # Compute Yt
                                data$Yt <- self$transform$Y(data$Xt)
                                # Exclude outliers from computations
                                idx_outliers_l <- which(data$Yt <= quantile(data$Yt, probs = self$control$outliers_quantile, na.rm = TRUE))
                                idx_outliers_r <- which(data$Yt >= quantile(data$Yt, probs = 1-self$control$outliers_quantile, na.rm = TRUE))
                                # Update outliers index and dates
                                outliers$index_type$threshold <- unique(c(idx_outliers_l, idx_outliers_r))
                                outliers$index <- unique(c(outliers$index, outliers$index_type$threshold))
                                outliers$date <- data$date[outliers$index]
                                # Rescale weights
                                if (!purrr::is_empty(outliers$index)) {
                                  data$weights[outliers$index] <- 0
                                  data$weights <- data$weights/sum(data$weights)
                                }
                                # **************************************************** #
                                # Update outliers
                                private$outliers <- outliers
                                # Update weights
                                private$..data$weights <- data$weights
                              },
                              #' @description
                              #' Fit a `seasonalModel` on `Yt` and compute deseasonalized series `Yt_tilde`.
                              fit_seasonal_mean = function(){
                                # Self Arguments
                                target <- self$target
                                control <- self$control$seasonal.mean
                                # Private arguments
                                data <- private$..data
                                # **************************************************** #
                                # Compute Yt
                                data$Yt <- self$transform$Y(data$Xt)
                                # Train data
                                data_train <- dplyr::filter(data, isTrain & weights != 0)
                                # Custom formula
                                formula_Yt <- ifelse(control$include.intercept, "Yt ~ 1", "Yt ~ -1")
                                formula_Yt <- paste0(formula_Yt, ifelse(control$include.H0, "+ H0", ""))
                                # Initialize the seasonal model for Yt
                                seasonal_model_Yt <- seasonalModel$new(order = control$seasonalOrder, period = 365)
                                # Seasonal model for Yt
                                seasonal_model_Yt$fit(formula = formula_Yt, data = data_train)
                                # Compute Yt_bar
                                data$Yt_bar <- seasonal_model_Yt$predict(data$n)
                                # Compute Yt_tilde
                                data$Yt_tilde <- data$Yt - data$Yt_bar
                                # Standardize parameters names
                                coefs_names <- c()
                                params <- seasonal_model_Yt$coefficients
                                if (control$include.intercept) {
                                  coefs_names[1] <- "a_0"
                                }
                                if (control$include.H0) {
                                  coefs_names <- c(coefs_names, "a_extra")
                                }
                                if (control$seasonalOrder > 0) {
                                  base_names <- paste0("a_", rep(c("cos_", "sin_")))
                                  coefs_names <- c(coefs_names, unlist(purrr::map(1:control$seasonalOrder, ~paste0(base_names, .x))))
                                }
                                names(params) <- coefs_names
                                seasonal_model_Yt$update(params)
                                # **************************************************** #
                                # Store seasonal model for Yt
                                private$..seasonal_model_Yt <- seasonal_model_Yt$clone(deep = TRUE)
                                # Add Yt_bar to seasonal data
                                private$..seasonal_data[["Yt_bar"]] <- private$..seasonal_model_Yt$predict(private$..seasonal_data$n)
                                # Add Yt to data
                                private$..data[["Yt"]] <- data$Yt
                                # Add Yt_tilde to data
                                private$..data[["Yt_tilde"]] <- data$Yt_tilde
                                # Add impled seasonal mean of the target to seasonal data
                                private$..seasonal_data[[paste0(target, "_bar")]] <- self$transform$GHI_y(private$..seasonal_data$Yt_bar, private$..seasonal_data$Ct)
                              },
                              #' @description
                              #' Correct the deseasonalized series `Yt_tilde` by subtracting its monthly mean.
                              corrective_monthly_mean = function(){
                                # Self Arguments
                                control <- self$control$seasonal.mean
                                # Private arguments
                                data <- private$..data
                                # **************************************************** #
                                # Unconditional mean for Yt_tilde
                                if (control$monthly.mean) {
                                  # Train data
                                  train_data <- dplyr::filter(data, isTrain)
                                  # Compute monthly unconditional mean
                                  monthly_data <- train_data %>%
                                    dplyr::group_by(Month) %>%
                                    dplyr::summarise(Yt_tilde_uncond = sum(Yt_tilde*weights, na.rm = TRUE)/sum(weights, na.rm = TRUE)) %>%
                                    dplyr::ungroup()
                                  # Add unconditional mean to the dataset
                                  data <- dplyr::left_join(data, monthly_data, by = c("Month"))
                                  # Update Yt_tilde
                                  data$Yt_tilde <- data$Yt_tilde - data$Yt_tilde_uncond
                                  # **************************************************** #
                                  # Add unconditional mean to monthly data
                                  private$..monthly_data[["Yt_tilde_uncond"]] <- monthly_data$Yt_tilde_uncond
                                  # Add Yt_tilde to data
                                  private$..data[["Yt_tilde"]] <- data$Yt_tilde
                                }
                              },
                              #' @description
                              #' Fit an AR model on `Yt_tilde` and compute residuals
                              fit_AR_model = function(){
                                # Self arguments
                                control <- self$control$mean.model
                                # Private arguments
                                data <- private$..data
                                # **************************************************** #
                                # Train data
                                data_train <- dplyr::filter(data, isTrain)
                                # AR with flexible formula for multiple orders
                                AR_formula_Yt <- "I(Yt_tilde) ~ "
                                if (control$arOrder > 0){
                                  for(i in 1:control$arOrder){
                                    AR_formula_Yt <- paste0(AR_formula_Yt, " + I(dplyr::lag(Yt_tilde,", i, "))")
                                  }
                                }
                                # Custom formula for intercept
                                AR_formula_Yt <- paste0(AR_formula_Yt, ifelse(control$include.intercept, "", "-1"))
                                # Fitted AR model
                                AR_model_Yt <- lm(formula = as.formula(AR_formula_Yt), data = data_train, weights = data_train$weights)
                                # Fitted Yt_tilde
                                data$Yt_tilde_hat <- predict.lm(AR_model_Yt, newdata = data)
                                # Set the initial values as the real ones
                                if (control$arOrder > 0) {
                                  data$Yt_tilde_hat[1:control$arOrder] <- data$Yt_tilde[1:control$arOrder]
                                }
                                # Fitted AR residuals
                                data$eps <- data$Yt_tilde - data$Yt_tilde_hat
                                # Standardize parameters names
                                params <- AR_model_Yt$coefficients
                                coefs_names <- c()
                                if (control$include.intercept) {
                                  coefs_names[1] <- "phi_0"
                                }
                                if (control$arOrder > 0){
                                  coefs_names <- c(coefs_names, paste0("phi_", 1:control$arOrder))
                                }
                                names(AR_model_Yt$coefficients) <- coefs_names
                                # **************************************************** #
                                # Store AR model
                                private$..AR_model_Yt <- AR_model_Yt
                                # Add fitted Yt_tilde
                                private$..data[["Yt_tilde_hat"]] <- data$Yt_tilde_hat
                                # Add fitted residuals
                                private$..data[["eps"]] <- data$eps
                              },
                              #' @description
                              #' Fit a `seasonalModel` on `eps^2` and compute deseasonalized residuals `eps_tilde`.
                              fit_seasonal_variance = function(){
                                # Self arguments
                                control <- self$control$seasonal.variance
                                # Private arguments
                                data <- private$..data
                                # **************************************************** #
                                # Train data
                                data_train <- dplyr::filter(data, isTrain & weights != 0)

                                # Initialize the seasonal model for eps^2
                                seasonal_variance <- seasonalModel$new(order = control$seasonalOrder, period = 365)
                                # Seasonal model for Yt
                                seasonal_variance$fit(formula = "I(eps^2) ~ 1", data = data_train)
                                # Fitted seasonal standard deviation
                                data$sigma_bar <- sqrt(seasonal_variance$predict(n = data$n))
                                # Fitted standardized residuals
                                data$eps_tilde <- data$eps/data$sigma_bar
                                # Correction to ensure unitary variance
                                if (control$correction) {
                                  # Update train data
                                  data_train <- dplyr::filter(data, isTrain & weights != 0)
                                  # Correct parameters to ensure unitary variance
                                  seasonal_variance$update(seasonal_variance$coefficients*mean(data_train$eps_tilde^2, na.rm = TRUE))
                                  # Fitted seasonal standard deviation
                                  data$sigma_bar <- sqrt(seasonal_variance$predict(n = data$n))
                                  # Updated standardized residuals
                                  data$eps_tilde <- data$eps/data$sigma_bar
                                }
                                coefs_names <- c("c_0")
                                if (control$seasonalOrder > 0) {
                                  base_names <- paste0("c_", rep(c("cos_", "sin_")))
                                  coefs_names <- c(coefs_names, unlist(purrr::map(1:control$seasonalOrder, ~paste0(base_names, .x))))
                                }
                                params <- seasonal_variance$coefficients
                                names(params) <- coefs_names
                                seasonal_variance$update(params)
                                # **************************************************** #
                                # Store seasonal variance model
                                private$..seasonal_variance <- seasonal_variance$clone(deep = TRUE)
                                # Add fitted deseasonalized residuals
                                private$..data[["eps_tilde"]] <- data$eps_tilde
                                # Add Yt_bar to seasonal data
                                private$..seasonal_data[["sigma_bar"]] <- sqrt(seasonal_variance$predict(private$..seasonal_data$n))
                              },
                              #' @description
                              #' Fit a `GARCH` model on `eps_tilde` and compute standardized `u` and monthly deseasonalized residuals `u_tilde`.
                              fit_GARCH_model = function(){
                                # Self arguments
                                control <- self$control
                                # Private arguments
                                data <- private$..data
                                # **************************************************** #
                                # Train data
                                data_train <- dplyr::filter(data, isTrain)
                                # GARCH specification
                                GARCH_spec <- control$variance.model
                                # Fitted model
                                GARCH_spec <- ugarchfit_fp(GARCH_spec, data_train$eps_tilde, data_train$weights, sigma0 = 1)
                                control$variance.model <- GARCH_spec
                                # Update fitted std. deviation
                                data$sigma <- rugarch::ugarchfilter(GARCH_spec, data$eps_tilde, out.sample = 0)@filter$sigma
                                # Compute standardized residuals
                                data$u <- data$eps_tilde/data$sigma
                                # Fitted standardized residuals
                                data$u_tilde <- data$u
                                # Initialize a GARCH object
                                GARCH <- list(coefficients = GARCH_spec@model$fixed.pars, vol = 1, omega = 0, alpha = c(), beta = c(), p = 0, q = 0)
                                # Intercept
                                GARCH$omega <- GARCH$coefficients[names(GARCH$coefficients) == "omega"]
                                # ARCH parameters
                                GARCH$alpha <- GARCH$coefficients[stringr::str_detect(names(GARCH$coefficients), "alpha")]
                                # ARCH order
                                GARCH$p <- max(c(length(GARCH$alpha), 0))
                                # GARCH parameters
                                GARCH$beta <- GARCH$coefficients[stringr::str_detect(names(GARCH$coefficients), "beta")]
                                # GARCH order
                                GARCH$q <- max(c(length(GARCH$beta), 0))
                                # **************************************************** #
                                # Store GARCH specification
                                private$..control$variance.model <- GARCH_spec
                                # Store GARCH parameters
                                private$..GARCH <- GARCH
                                # Add fitted GARCH std. deviation
                                private$..data[["sigma"]] <- data$sigma
                                # Add fitted standardized residuals
                                private$..data[["u"]] <- data$u
                                # Add fitted standardized and monthly deseasonalized residuals
                                private$..data[["u_tilde"]] <- data$u_tilde
                              },
                              #' @description
                              #' Correct the standardized series `u` by dividing by its monthly std. deviation.
                              corrective_monthly_variance = function(){
                                # Self arguments
                                control <- self$control
                                # Private arguments
                                data <- private$..data
                                # **************************************************** #
                                # Monthly correction
                                if (control$seasonal.variance$monthly.mean) {
                                  monthly_data <- dplyr::filter(data, isTrain & weights != 0) %>%
                                    dplyr::group_by(Month) %>%
                                    dplyr::summarise(sigma_m = sd(u)) %>%
                                    dplyr::ungroup()
                                  # Add unconditional variance to the dataset
                                  data <- dplyr::left_join(data, monthly_data, by = c("Month"))
                                  # Fitted standardized residuals
                                  data$u_tilde <- data$eps_tilde/(data$sigma*data$sigma_m)
                                  # Remove unconditional mean to the dataset
                                  data <- dplyr::select(data, -sigma_m)
                                  # **************************************************** #
                                  # Add fitted standardized and monthly deseasonalized residuals
                                  private$..data[["u_tilde"]] <- data$u_tilde
                                  # Add corrective std. deviation to monthly data
                                  private$..monthly_data[["sigma_m"]] <- monthly_data$sigma_m
                                }
                              },
                              #' @description
                              #' Fit a `gaussianMixture` monthly model on `u_tilde` and return a series of bernoulli `B` and standardized components `z1` and `z2`.
                              fit_mixture_model = function(){
                                # Self arguments
                                control <- self$control$mixture.model
                                # Private arguments
                                data <- private$..data
                                # **************************************************** #
                                # 4) Gaussian Mixture
                                GM_model <- solarModel_mixture(data$u_tilde, date = data$date, weights = data$weights,
                                                               match_moments = control$match_moments,
                                                               abstol = control$abstol, maxit = control$maxit)
                                # Add the fitted series of Bernoulli
                                data <- dplyr::left_join(data, dplyr::select(GM_model$data, date, B), by = c("date"))
                                # Compute standardized components
                                data <- dplyr::left_join(data, GM_model$model, by = "Month")
                                data$z1 <- data$B*(data$u_tilde - data$mu_up)/data$sd_up
                                data$z2 <- (1-data$B)*(data$u_tilde - data$mu_dw)/data$sd_dw
                                # **************************************************** #
                                # Store Gaussian Mixture parameters
                                private$..NM_model <- GM_model$model
                                # Add fitted series of bernoulli and of standardized components
                                private$..data[["B"]] <- data$B
                                private$..data[["z1"]] <- data$z1
                                private$..data[["z2"]] <- data$z2
                              },
                              #' @description
                              #' Update the parameters inside object
                              #' @param params updated parameters
                              update = function(params){
                                # Update clear sky model
                                private$..seasonal_model_Ct$update(unlist(params$seasonal_model_Ct))
                                # Update seasonal model
                                private$..seasonal_model_Yt$update(unlist(params$seasonal_model_Yt))
                                # Update AR model
                                private$..AR_model_Yt$coefficients <- unlist(params$AR_model_Yt)
                                # Update seasonal variance
                                private$..seasonal_variance$update(unlist(params$seasonal_variance))
                                # Update GARCH variance
                                private$..GARCH$coefficients <- unlist(params$GARCH)
                                private$..GARCH$vol <- sqrt(private$..GARCH$coefficients[1]/(1-sum(private$..GARCH$coefficients[-1])))
                                private$..control$variance.model@model$fixed.pars <- unlist(params$GARCH)
                                # Update Gaussian Mixture
                                private$..NM_model$mu_up <- unlist(params$NM_mu_up)
                                private$..NM_model$mu_dw <- unlist(params$NM_mu_dw)
                                private$..NM_model$sd_up <- unlist(params$NM_sd_up)
                                private$..NM_model$sd_dw <- unlist(params$NM_sd_dw)
                                private$..NM_model$p_up <- unlist(params$NM_p_up)
                                private$..NM_model$p_dw <- 1 - private$..NM_model$p_up
                                # Update theoric moments
                                NM <- private$..NM_model
                                private$..NM_model$e_x <- NM$mu_up*NM$p_up + NM$mu_dw*NM$p_dw
                                private$..NM_model$v_x <- NM$p_up*NM$p_dw*(NM$mu_up-NM$mu_dw)^2 + NM$p_up*NM$sd_up^2 + NM$p_dw*NM$sd_dw^2
                                # Set Log-likelihood to NA
                                private$..loglik$loglik <- rep(NA, 12)
                              },
                              #' @description
                              #' Update the time series inside object given certain parameters
                              filter = function(){
                                # Self arguments
                                control <- self$control
                                target <- self$target
                                # **************************************************** #
                                # Update risk driver
                                self$compute_risk_drivers()
                                # Update solar transform
                                self$fit_solar_transform()
                                # Update Yt
                                private$..data[["Yt"]] <- self$transform$Y(private$..data[["Xt"]])
                                # Update Yt_tilde
                                private$..data[["Yt_tilde"]] <- private$..data[["Yt"]] - self$seasonal_model_Yt$predict(private$..data$n)
                                # Fit the corrective mean
                                self$corrective_monthly_mean()
                                # Update Yt_tilde_hat
                                private$..data[["Yt_tilde_hat"]] <-  predict.lm(self$AR_model_Yt, newdata = private$..data)
                                # Initial values as the real ones
                                if (control$mean.model$arOrder > 0) {
                                  private$..data[["Yt_tilde_hat"]][1:control$mean.model$arOrder] <- private$..data[["Yt_tilde"]][1:control$mean.model$arOrder]
                                }
                                # Update AR residuals
                                private$..data[["eps"]] <- private$..data[["Yt_tilde"]] - private$..data[["Yt_tilde_hat"]]
                                # Fitted standardized residuals
                                private$..data[["eps_tilde"]] <- private$..data[["eps"]]/sqrt(self$seasonal_variance$predict(private$..data$n))
                                # Update Garch standard deviation
                                private$..data[["sigma"]] <- rugarch::ugarchfilter(control$variance.model, private$..data[["eps_tilde"]], out.sample = 0)@filter$sigma
                                # Update fitted standard residuals
                                private$..data[["u"]] <- private$..data[["eps_tilde"]]/private$..data[["sigma"]]
                                private$..data[["u_tilde"]] <- private$..data[["u"]]
                                # Fit the corrective variance
                                self$corrective_monthly_variance()
                                # Update seasonal data
                                private$..seasonal_data[["Ct"]] <- self$seasonal_model_Ct$predict(self$seasonal_data$n)
                                private$..seasonal_data[["Yt_bar"]] <- self$seasonal_model_Yt$predict(self$seasonal_data$n)
                                private$..seasonal_data[[paste0(target, "_bar")]] <- self$transform$iY(private$..seasonal_data[["Yt_bar"]])
                                private$..seasonal_data[["sigma_bar"]] <- sqrt(self$seasonal_variance$predict(self$seasonal_data$n))
                                # Update moments
                                self$conditional_moments()
                                self$unconditional_moments()
                                self$logLik()
                              },
                              #' @description
                              #' Compute the conditional moments.
                              conditional_moments = function(){
                                # Self arguments
                                data <- self$data
                                # **************************************************** #
                                # Compute conditional moments
                                data <- dplyr::mutate(data,
                                                      # Conditional expectation Yt
                                                      e_Yt = Yt_bar + Yt_tilde_hat + Yt_tilde_uncond,
                                                      # Conditional std. deviation Yt
                                                      sd_Yt = sigma*sigma_bar*sigma_m,
                                                      # Conditional moments Yt (state up)
                                                      e_Yt_up = e_Yt + sd_Yt*mu_up,
                                                      sd_Yt_up = sd_Yt*sd_up,
                                                      # Conditional moments Yt (state dw)
                                                      e_Yt_dw = e_Yt + sd_Yt*mu_dw,
                                                      sd_Yt_dw = sd_Yt*sd_dw)
                                # Store only relevant variables
                                data <- dplyr::select(data, date, Year, Month, Day, e_Yt, sd_Yt, e_Yt_up, sd_Yt_up, e_Yt_dw, sd_Yt_dw)
                                # **************************************************** #
                                private$..moments$conditional <- data
                              },
                              #' @description
                              #' Compute the unconditional seasonal moments.
                              unconditional_moments = function(){
                                # Self arguments
                                data <- self$seasonal_data
                                # **************************************************** #
                                # Compute unconditional moments
                                data <- dplyr::mutate(data,
                                                      # Unconditional expectation Yt
                                                      e_Yt = Yt_bar + Yt_tilde_uncond,
                                                      # Unconditional std. deviation Yt
                                                      sd_Yt = sigma_bar*sigma_m,
                                                      # Unconditional moments Yt (state up)
                                                      e_Yt_up = e_Yt + sd_Yt*mu_up,
                                                      sd_Yt_up = sd_Yt*sd_up,
                                                      # Unconditional moments Yt (state dw)
                                                      e_Yt_dw = e_Yt + sd_Yt*mu_dw,
                                                      sd_Yt_dw = sd_Yt*sd_dw)
                                # Store only relevant variables
                                data <- dplyr::select(data, Month, Day, e_Yt, sd_Yt, e_Yt_up, sd_Yt_up, e_Yt_dw, sd_Yt_dw)
                                # **************************************************** #
                                private$..moments$unconditional <- data
                              },
                              #' @description
                              #' Compute the log-likelihood of the model given the parameters.
                              logLik = function(){
                                # Self arguments
                                moments <- dplyr::left_join(private$..moments$conditional, private$..data[, c("date", "Yt")], by = "date")
                                moments <- dplyr::left_join(moments, self$NM_model[, c("Month", "p_up")], by = "Month")
                                # Add weights
                                moments$weights <- self$data$weights
                                # **************************************************** #
                                # Compute log-likelihood
                                moments$loglik <- 0
                                for(i in 1:nrow(moments)){
                                  df_n <- moments[i,]
                                  if (df_n$weights == 0){
                                    moments$loglik[i] <- 0
                                    next
                                  }
                                  # Conditional mixture Pdf
                                  pdf_Yt <- function(x) dmixnorm(x, means = c(df_n$e_Yt_up, df_n$e_Yt_dw), sd = c(df_n$sd_Yt_up, df_n$sd_Yt_dw), p = c(df_n$p_up, 1-df_n$p_up))
                                  moments$loglik[i] <- log(pdf_Yt(df_n$Yt))
                                }
                                # Aggregate log-likelihood by month
                                moments <- moments %>%
                                  dplyr::group_by(Month) %>%
                                  dplyr::summarise(loglik = sum(loglik, na.rm = TRUE))
                                # **************************************************** #
                                # Update log-likelihood
                                private$..loglik$loglik <- moments$loglik
                              }
                            ),
                            private = list(
                              ..data = NA,
                              ..seasonal_data = NA,
                              ..monthly_data = NA,
                              ..coords = NA,
                              ..control = NA,
                              ..transform = NA,
                              ..combinations = NA,
                              interpolated = FALSE,
                              outliers = NA,
                              ..seasonal_model_Ct = NA,
                              ..seasonal_model_Yt = NA,
                              ..AR_model_Yt = NA,
                              ..seasonal_variance = NA,
                              ..GARCH = NA,
                              ..NM_model = NA,
                              ..loglik = NA,
                              ..moments = list(conditional = NA, unconditional = NA)
                            ),

                            active = list(
                              #' @field data Get a data frame containing the complete data with seasonal and monthly parameters.
                              data = function(){
                                dplyr::left_join(private$..data, dplyr::select(self$seasonal_data, -n), by = c("Month", "Day"))
                              },
                              #' @field seasonal_data Get a data frame containing seasonal and monthly parameters.
                              seasonal_data = function(){
                                dplyr::left_join(private$..seasonal_data, self$monthly_data, by = "Month")
                              },
                              #' @field monthly_data Get a data frame that contains monthly parameters.
                              monthly_data = function(){
                                dplyr::left_join(private$..monthly_data, self$NM_model, by = "Month")
                              },
                              #' @field loglik Get the log-likelihood of the train data.
                              loglik = function(){
                                sum(private$..loglik)
                              },
                              #' @field control A list of control parameters that govern the behavior of the model's fitting process and other configurations.
                              control = function(){
                                private$..control
                              },
                              #' @field location A data frame with coordinates of the location considered.
                              location = function(){
                                dplyr::bind_cols(place = self$place, dplyr::bind_rows(self$coords))
                              },
                              #' @field transform An object representing the transformation functions applied to the data.
                              transform = function(){
                                private$..transform
                              },
                              #' @field seasonal_model_Ct The fitted model for clear sky radiation, used for predict the maximum radiation available.
                              seasonal_model_Ct = function(){
                                private$..seasonal_model_Ct
                              },
                              #' @field seasonal_model_Yt The fitted seasonal model for the target variable.
                              seasonal_model_Yt = function(){
                                private$..seasonal_model_Yt
                              },
                              #' @field AR_model_Yt The fitted Autoregressive (AR) model for the target variable.
                              AR_model_Yt = function(){
                                private$..AR_model_Yt
                              },
                              #' @field seasonal_variance The fitted model for seasonal variance.
                              seasonal_variance = function(){
                                private$..seasonal_variance
                              },
                              #' @field GARCH A model object representing the GARCH model fitted to the residuals.
                              GARCH = function(){
                                private$..GARCH
                              },
                              #' @field NM_model A model object representing the Gaussian Mixture model fitted to the standardized residuals.
                              NM_model = function(){
                                private$..NM_model[, c(1:6)]
                              },
                              #' @field moments Get a list containing the conditional and unconditional moments.
                              moments = function(){
                                private$..moments
                              },
                              #' @field parameters Get the model parameters as a named list.
                              parameters = function(){
                                # 1. Clears sky seasonal model
                                coefs_names <- c()
                                seasonal_model_Ct <- as.list(self$seasonal_model_Ct$coefficients)
                                # 2. Seasonal model Yt
                                seasonal_model_Yt <- as.list(self$seasonal_model_Yt$coefficients)
                                # 3. AR model
                                AR_model_Yt <- as.list(self$AR_model_Yt$coefficients)
                                # 4. Seasonal variance model
                                seasonal_variance <- as.list(self$seasonal_variance$coefficients)
                                # 5. GARCH(1,1) variance model
                                GARCH <- as.list(self$GARCH$coefficients)
                                # 6. Gaussian mixture model
                                NM_mu_up <- self$NM_model$mu_up
                                names(NM_mu_up) <- paste0("mu_up_", 1:12)
                                NM_mu_dw <- self$NM_model$mu_dw
                                names(NM_mu_dw) <- paste0("mu_dw_", 1:12)
                                NM_sd_up <- self$NM_model$sd_up
                                names(NM_sd_up) <- paste0("sd_up_", 1:12)
                                NM_sd_dw <- self$NM_model$sd_dw
                                names(NM_sd_dw) <- paste0("sd_dw_", 1:12)
                                NM_p_up <- self$NM_model$p_up
                                names(NM_p_up) <- paste0("p_", 1:12)

                                # Output list of parameter for each model
                                list(
                                  location = dplyr::tibble(place = self$place,
                                                           target = self$target,
                                                           lat = self$coords$lat,
                                                           lon = self$coords$lon,
                                                           alt = self$coords$alt),
                                  params = dplyr::bind_cols(alpha = self$transform$alpha, beta = self$transform$beta),
                                  seasonal_model_Ct = dplyr::bind_cols(seasonal_model_Ct),
                                  seasonal_model_Yt = dplyr::bind_cols(seasonal_model_Yt),
                                  AR_model_Yt = dplyr::bind_cols(AR_model_Yt),
                                  seasonal_variance = dplyr::bind_cols(seasonal_variance),
                                  GARCH = dplyr::bind_cols(GARCH),
                                  NM_mu_up = dplyr::bind_rows(NM_mu_up),
                                  NM_mu_dw = dplyr::bind_rows(NM_mu_dw),
                                  NM_sd_up = dplyr::bind_rows(NM_sd_up),
                                  NM_sd_dw = dplyr::bind_rows(NM_sd_dw),
                                  NM_p_up = dplyr::bind_rows(NM_p_up)
                                )
                              }
                            )
                          )

