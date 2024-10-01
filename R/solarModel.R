#' Fit a model for solar radiation
#'
#' @param spec an object with class `solarModelSpec`. See the function \code{\link{solarModel_spec}} for details.
#'
#' @rdname solarModel
#' @name solarModel
#' @examples
#' # Control list
#' control <- control_solarModel(outliers_quantile = 0.005)
#' # Model specification
#' spec_1 <- solarModel_spec("Bologna", from="2005-01-01", to="2022-01-01", control_model = control)
#' # Model fit
#' Bologna <- solarModel(spec_1)
#' # save(Bologna, file = "data/Bologna.RData")
#'
#' # Example for Oslo
#' # Model specification
#' spec_2 <- solarModel_spec("Oslo", from="2005-01-01", to="2022-01-01", control_model = control)
#' # Model fit
#' Oslo <- solarModel(spec_2)
#' # save(Oslo, file = "data/Oslo.RData")
#'
#' # Use the realized clearsky
#' spec_1$control$stochastic_clearsky <- TRUE
#' model_B <- solarModel(spec_1)
#'
#' # Fit a model for clearsky
#' spec_Ct <- spec_1
#' spec_Ct$control$stochastic_clearsky <- FALSE
#' spec_Ct$target <- "clearsky"
#' model_Ct <- solarModel(spec_Ct)
#' @export
solarModel <- function(spec){

  # Extract control parameters
  control <- spec$control
  # Target variable
  target <- spec$target
  # Complete data
  data <- spec$data

  # Seasonal model for clear sky radiation
  seasonal_model_Ct <- seasonalClearsky(spec$data[[target]], spec$data$date, spec$coords$lat, spec$data$clearsky, control$clearsky)
  data$Ct <- seasonal_model_Ct$predict(data$n)

  # 0) Risk drivers
  st <- solarTransform$new(0, 1)
  # Compute risk drivers
  if (control$stochastic_clearsky) {
    data$clearsky <- data$clearsky*1.01
    # Detect and impute outliers
    outliers <- clearsky_outliers(data[[target]], data$clearsky, date = data$date, quiet = control$quiet)
    # Update GHI
    data[[target]] <- outliers$x
    data$Xt <- st$iGHI(data[[target]], data$clearsky)
  } else {
    # Detect and impute outliers
    outliers <- clearsky_outliers(data[[target]], data$Ct, date = data$date, quiet = control$quiet)
    # Update GHI
    data[[target]] <- outliers$x
    data$Xt <- st$iGHI(data[[target]], data$Ct)
  }

  # 1) Solar Transform
  # Remove minimum and maximum of Xt from computations
  outliers$index <- c(outliers$index, which.min(data$Xt), which.max(data$Xt))
  # Transformation parameters
  params <- st$parameters(data$Xt, control$threshold)
  #params <- st$parameters()
  st$update(params$alpha, params$beta)
  # Compute Yt
  data$Yt <- st$Y(data$Xt)
  # Exclude outliers from computations
  idx_outliers_l <- which(data$Yt <= quantile(data$Yt, probs = control$outliers_quantile, na.rm = TRUE))
  idx_outliers_r <- which(data$Yt >= quantile(data$Yt, probs = 1-control$outliers_quantile, na.rm = TRUE))
  outliers$index <- unique(c(outliers$index, idx_outliers_l, idx_outliers_r))
  # Rescale weights
  if (!purrr::is_empty(outliers$index)) {
    data$weights[outliers$index] <- 0
    data$weights <- data$weights/sum(data$weights)
  }

  # 2) Seasonal Mean
  data_train <- dplyr::filter(data, isTrain & weights != 0)
  # Custom formula
  formula_Yt <- ifelse(control$seasonal.mean$include.intercept, "Yt ~ 1", "Yt ~ -1")
  formula_Yt <- paste0(formula_Yt, ifelse(control$seasonal.mean$include.H0, "+ H0", ""))
  # Seasonal model for Yt
  seasonal_model <- seasonalModel$new(formula = formula_Yt, order = control$seasonal.mean$seasonalOrder, data = data_train)
  # Fitted seasonal mean for Yt
  data$Yt_bar <- seasonal_model$predict(data$n)
  # Fitted deseasonalized Yt
  data$Yt_tilde <- data$Yt - data$Yt_bar

  # Unconditional mean for Yt_tilde
  if (control$seasonal.mean$monthly.mean) {
    # Store monthly unconditional mean only
    monthly_data <- dplyr::filter(data, isTrain) %>%
      dplyr::group_by(Month) %>%
      dplyr::summarise(Yt_tilde_uncond = sum(Yt_tilde*weights, na.rm = TRUE)/sum(weights, na.rm = TRUE)) %>%
      dplyr::ungroup()
    # Add unconditional mean to the dataset
    data <- dplyr::left_join(data, monthly_data, by = c("Month"))
    # Deseasonalized series centered in zero
    data$Yt_tilde <- data$Yt_tilde - data$Yt_tilde_uncond
    # Remove unconditional mean to the dataset
    data <- dplyr::select(data, -Yt_tilde_uncond)
  } else {
    monthly_data <- dplyr::tibble(Month = 1:12, Yt_tilde_uncond = 0)
  }

  # 3) AR model
  # AR with flexible formula for multiple orders
  AR_formula_Yt <- "I(Yt_tilde) ~ "
  if (control$mean.model$arOrder > 0){
    for(i in 1:control$mean.model$arOrder){
      AR_formula_Yt <- paste0(AR_formula_Yt, " + I(dplyr::lag(Yt_tilde,", i, "))")
    }
  }
  # Control for intercept
  AR_formula_Yt <- paste0(AR_formula_Yt, ifelse(control$mean.model$include.intercept, "", "-1"))
  # Fitted AR model
  data_train <- dplyr::filter(data, isTrain)
  AR_model_Yt <- lm(formula = as.formula(AR_formula_Yt), data = data_train, weights = data_train$weights)
  # Fitted Yt_tilde
  data$Yt_tilde_hat <- predict.lm(AR_model_Yt, newdata = data)
  # Initial values as the real ones
  if (control$mean.model$arOrder > 0) {
    data$Yt_tilde_hat[1:control$mean.model$arOrder] <- data$Yt_tilde[1:control$mean.model$arOrder]
  }
  # Fitted AR(2) residuals
  data$eps <- data$Yt_tilde - data$Yt_tilde_hat

  # 4) Seasonal variance
  data_train <- dplyr::filter(data, isTrain & weights != 0)
  # Custom formula
  seasonal_variance <- seasonalModel$new(formula = "I(eps^2) ~ 1", order = control$seasonal.variance$seasonalOrder, data = data_train)
  # Fitted seasonal standard deviation
  data$sigma_bar <- sqrt(seasonal_variance$predict(n = data$n))
  # Fitted standardized residuals
  data$eps_tilde <- data$eps/data$sigma_bar
  # Correction to ensure unitary variance
  if (control$seasonal.variance$correction) {
    # Update train data
    data_train <- dplyr::filter(data, isTrain & weights != 0)
    # Correct parameters to ensure unitary variance
    seasonal_variance$update(seasonal_variance$coefficients*mean(data_train$eps_tilde^2, na.rm = TRUE))
    # Fitted seasonal standard deviation
    data$sigma_bar <- sqrt(seasonal_variance$predict(n = data$n))
    # Updated standardized residuals
    data$eps_tilde <- data$eps/data$sigma_bar
  }

  # 6) GARCH model
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

  # Monthly correction
  if (control$seasonal.variance$monthly.mean) {
    sd_u_m <- dplyr::filter(data, isTrain & weights != 0) %>%
      dplyr::group_by(Month) %>%
      dplyr::summarise(sigma_m = sd(u)) %>%
      dplyr::ungroup()
    monthly_data$sigma_m <- sd_u_m$sigma_m
    # Add unconditional variance to the dataset
    data <- dplyr::left_join(data, sd_u_m, by = c("Month"))
    # Fitted standardized residuals
    data$u_tilde <- data$eps_tilde/(data$sigma*data$sigma_m)
    # Remove unconditional mean to the dataset
    data <- dplyr::select(data, -sigma_m)
    # Store the corrective parameter
    monthly_data$sigma_m <- sd_u_m$sigma_m
  } else {
    monthly_data$sigma_m <- 1
    # Fitted standardized residuals
    data$u_tilde <- data$u
  }

  # 4) Gaussian Mixture
  GM_model <- solarModel_mixture(data$u_tilde, date = data$date, weights = data$weights,
                                 match_moments = control$mixture.model$match_moments,
                                 abstol = control$mixture.model$abstol, maxit = control$mixture.model$maxit)
  # Return also the fitted series of Bernoulli
  data <- dplyr::left_join(data, dplyr::select(GM_model$data, date, B), by = c("date"))
  # Store the parameters
  GM_model <- GM_model$model

  # Compute conditional parameters and likelihood
  df <- dplyr::left_join(data, monthly_data, by = c("Month"))
  df <- dplyr::left_join(df, GM_model[, c(1:7)], by = c("Month"))
  # Fitted values
  df <- dplyr::mutate(df,
                      Yt_hat = Yt_bar + Yt_tilde_uncond + Yt_tilde_hat,
                      mu_up = Yt_hat + sigma_bar*sigma_m*sigma*mu_up,
                      mu_dw = Yt_hat + sigma_bar*sigma_m*sigma*mu_dw,
                      sd_up = sigma_bar*sigma_m*sigma*sd_up,
                      sd_dw = sigma_bar*sigma_m*sigma*sd_dw)
  df <- dplyr::select(df, date, Yt, mu_up, mu_dw, sd_up, sd_dw, p_up)

  df$log.lik <- 0
  for(i in 1:nrow(df)){
    means <- c(df$mu_up[i], df$mu_dw[i])
    sd <- c(df$sd_up[i], df$sd_dw[i])
    probs <- c(df$p_up[i], 1-df$p_up[i])
    df$log.lik[i] <- dmixnorm(df$Yt[i], means, sd, probs, log = TRUE)
  }

  # Seasonal data by month and day for an year with 366 days
  reference_year <- unique(data$Year)[which(is_leap_year(paste0(unique(data$Year), "-01-01")))[1]]
  from_date <- as.Date(paste0(reference_year, "-01-01"))
  to_date <- as.Date(paste0(reference_year, "-12-31"))
  seasonal_data <- dplyr::filter(data, date >= from_date & date <= to_date)
  seasonal_data <- dplyr::select(seasonal_data, Month, Day, H0, Ct, Yt_bar, sigma_bar)
  seasonal_data[[paste0(target, "_bar")]] <- st$GHI_y(seasonal_data$Yt_bar, seasonal_data$Ct)

  # Remove seasonal variables
  data <- dplyr::select(data, -Ct, -H0, -Yt_bar, -sigma_bar)

  # Update specification data
  spec$data <- dplyr::select(data, -weights)
  # Structure model data
  model <-  list(
    seasonal_model_Ct = seasonal_model_Ct,
    seasonal_model_Yt = seasonal_model,
    AR_model_Yt = AR_model_Yt,
    seasonal_variance = seasonal_variance,
    GARCH = list(coef = GARCH_spec@model$fixed.pars, vol = 1),
    NM_model = GM_model,
    seasonal_data = seasonal_data,
    monthly_data = monthly_data,
    transform = st,
    params = params,
    log_lik = sum(df$log.lik, na.rm = TRUE),
    outliers = append(outliers[-c(1)], list(weights = data$weights)),
    control = control,
    # Extra slot: scenarios
    scenarios = list(P = NA, Q = NA, Qdw = NA, Qup = NA),
    # Extra slot: payoffs
    payoffs = list(hist = NA,
                   boot = NA,
                   model = list(P = NA, Q = NA, Qdw = NA, Qup = NA, structured = NA),
                   sim = list(P = NA, Q = NA, Qdw = NA, Qup = NA, structured = NA),
                   control = NA),
    # Extra slot: esscher parameters
    esscher = list(grid = NA, params = NA, model = NA,
                   coefficients = NA, theta = NA, theta_m = NA, implied_r = NA,
                   control = NA)
  )

  structure(
    append(spec, model),
    class = c("solarModel", "solarModelSpec", "list")
  )
}


#' Update the time series inside a `solarModel` object
#'
#' @param model object with the class `solarModel`. See the function \code{\link{solarModel}} for details.
#' @param params named list of parameters. See the function \code{\link{solarModel_parameters}} to structure the list of new parameters.
#' model.
#' @param as_tibble logical, when `TRUE` the output will be converted in a `tibble`.
#' @examples
#' model <- Bologna
#' # Extract the parameters
#' params <- solarModel_parameters(model)
#' # Make some changes on the parameters
#' params$AR_model_Yt$phi_1 <- params$AR_model_Yt$phi_1*0.95
#' # Update the parameters inside the model
#' model <- solarModel_update(model, params)
#' # Filter the time series with new parameters
#' model <- solarModel_filter(model)
#' @rdname solarModel_filter
#' @name solarModel_filter
#' @aliases solarModel_update
#' @aliases solarModel_parameters
#' @export
solarModel_filter <- function(model){

  # Extract control parameters
  control <- model$control
  target <- model$target

  # Complete data
  data <- dplyr::bind_cols(model$data, weights = model$outliers$weights)
  data <- dplyr::mutate(data, isTrain = ifelse(is.na(isTrain), FALSE, isTrain))
  # 0) Clear sky and risk drivers
  # Fitted seasonal clear sky
  data$Ct <- model$seasonal_model_Ct$predict(data$n)
  # Risk Driver (Xt is computed here)
  outliers <- clearsky_outliers(data[[target]], data$Ct, data$date, quiet = control$clearsky$quiet)
  # Update the dataset with imputed values
  data[[target]] <- outliers$x
  # 1) Solar Transform
  st <- model$transform
  # Compute risk drivers
  data$Xt <- st$iGHI(data[[target]], data$Ct)

  # Remove minimum and maximum of Xt from computations
  outliers$index <- unique(c(model$outliers$index, which.min(data$Xt), which.max(data$Xt)))
  data$weights[outliers$index] <- 0
  # Transformation parameters
  params <- st$parameters(data$Xt, control$threshold)
  # Update solar transform
  st$update(params$alpha, params$beta)
  # Compute Yt
  data$Yt <- st$Y(data$Xt)
  # Fitted seasonal mean for Yt
  data$Yt_bar <- model$seasonal_model_Yt$predict(data$n)
  # Fitted deseasonalized Yt
  data$Yt_tilde <- data$Yt - data$Yt_bar
  # Unconditional mean for Yt_tilde
  if (control$seasonal.mean$monthly.mean) {
    # Store monthly unconditional mean only
    monthly_data <- dplyr::filter(data, isTrain) %>%
      dplyr::group_by(Month) %>%
      dplyr::summarise(Yt_tilde_uncond = sum(Yt_tilde*weights, na.rm = TRUE)/sum(weights, na.rm = TRUE)) %>%
      dplyr::ungroup()
    # Add unconditional mean to the dataset
    data <- dplyr::left_join(data, monthly_data, by = c("Month"))
    # Deseasonalized series centered in zero
    data$Yt_tilde <- data$Yt_tilde - data$Yt_tilde_uncond
    # Remove unconditional mean to the dataset
    data <- dplyr::select(data, -Yt_tilde_uncond)
  } else {
    monthly_data <- dplyr::tibble(Month = 1:12, Yt_tilde_uncond = 0)
  }

  # 2) AR model
  # Fitted Yt_tilde
  data$Yt_tilde_hat <- predict.lm(model$AR_model_Yt, newdata = data)
  # Initial values as the real ones
  if (control$mean.model$arOrder > 0) {
    data$Yt_tilde_hat[1:control$mean.model$arOrder] <- data$Yt_tilde[1:control$mean.model$arOrder]
  }
  # Fitted AR(2) residuals
  data$eps <- data$Yt_tilde - data$Yt_tilde_hat

  # 3) Seasonal variance
  # Fitted seasonal standard deviation
  data$sigma_bar <- sqrt(model$seasonal_variance$predict(data$n))
  # Fitted standardized residuals
  data$eps_tilde <- data$eps/data$sigma_bar
  # Correction to ensure unitary variance
  # Update train data
  data_train <- dplyr::filter(data, isTrain & weights != 0)
  # Correction to ensure unitary variance
  model$seasonal_variance$update(model$seasonal_variance$coefficients*mean(data_train$eps_tilde^2, na.rm = TRUE))
  # Fitted seasonal standard deviation
  data$sigma_bar <- sqrt(model$seasonal_variance$predict(data$n))
  # Fitted standardized residuals
  data$eps_tilde <- data$eps/data$sigma_bar

  # 4) GARCH variance
  # Add fitted standard deviation
  data$sigma <- rugarch::ugarchfilter(control$variance.model, data$eps_tilde, out.sample = 0)@filter$sigma
  # Add fitted final residuals
  data$u <- data$eps_tilde/data$sigma
  # Monthly correction
  if (control$seasonal.variance$correction) {
    sd_u_m <- dplyr::filter(data, isTrain & weights != 0) %>%
      dplyr::group_by(Month) %>%
      dplyr::summarise(sigma_m = sd(u)) %>%
      dplyr::ungroup()
    monthly_data$sigma_m <- sd_u_m$sigma_m
    # Add unconditional variance to the dataset
    data <- dplyr::left_join(data, sd_u_m, by = c("Month"))
    # Fitted standardized residuals
    data$u_tilde <- data$eps_tilde/(data$sigma*data$sigma_m)
    # Remove unconditional mean to the dataset
    data <- dplyr::select(data, -sigma_m)
    # Store the corrective parameter
    monthly_data$sigma_m <- sd_u_m$sigma_m
  } else {
    monthly_data$sigma_m <- 1
    # Fitted standardized residuals
    data$u_tilde <- data$u
  }

  # Update seasonal data
  seasonal_data <- model$seasonal_data
  seasonal_data$Ct <- model$seasonal_model_Ct$predict(1:nrow(seasonal_data))
  seasonal_data$Yt_bar <- model$seasonal_model_Yt$predict(1:nrow(seasonal_data))
  seasonal_data$Xt_bar <- st$iY(seasonal_data$Yt_bar)
  seasonal_data[[paste0(target, "_bar")]] <- st$GHI(seasonal_data$Xt_bar, seasonal_data$Ct)
  seasonal_data$sigma_bar <- sqrt(model$seasonal_variance$predict(1:nrow(seasonal_data)))

  # Update object
  model$params <- params
  model$transform <- st
  model$monthly_data <- monthly_data
  model$seasonal_data <- seasonal_data
  model$outliers <- outliers
  model$outliers$weights <- data$weights
  model$data <- dplyr::select(data, -Ct, -Yt_bar, -sigma_bar, -weights)
  return(model)
}

#' Extract the parameters of a `solarModel`
#'
#' @rdname solarModel_filter
#' @export
solarModel_parameters <- function(model, as_tibble = FALSE){

  # Control parameters
  control <- model$control

  # 1. Clears sky seasonal model
  coefs_names <- c()
  seasonal_model_Ct <- as.list(model$seasonal_model_Ct$coefficients)
  if (control$clearsky$include.intercept) {
    coefs_names[1] <- "delta_0"
  }
  if (control$clearsky$order > 0){
    coefs_names <- c(coefs_names, paste0("delta_", 1:(2*control$clearsky$order+1)))
  }
  names(seasonal_model_Ct) <- coefs_names

  # 2. Seasonal model Yt
  coefs_names <- c()
  seasonal_model_Yt <- as.list(model$seasonal_model_Yt$coefficients)
  if (control$seasonal.mean$include.intercept) {
    coefs_names[1] <- "a_0"
  }
  if (control$seasonal.mean$seasonalOrder > 0){
    coefs_names <- c(coefs_names, paste0("a_", 1:(length(seasonal_model_Yt)-1)))
  }
  names(seasonal_model_Yt) <- coefs_names

  # 3. AR model
  AR_model_Yt <- as.list(model$AR_model_Yt$coefficients)
  coefs_names <- c()
  if (control$mean.model$include.intercept) {
    coefs_names[1] <- "phi_0"
  }
  if (control$seasonal.mean$seasonalOrder > 0){
    coefs_names <- c(coefs_names, paste0("phi_", 1:control$mean.model$arOrder))
  }
  names(AR_model_Yt) <- coefs_names

  # 4. Seasonal variance model
  seasonal_variance <- as.list(model$seasonal_variance$coefficients)
  coefs_names <- c("c_0")
  if (control$seasonal.variance$seasonalOrder > 0){
    coefs_names <- c(coefs_names, paste0("c_", 1:(2*control$seasonal.variance$seasonalOrder)))
  }
  names(seasonal_variance) <- coefs_names

  # 5. GARCH(1,1) variance model
  GARCH <- as.list(model$GARCH$coef)

  # 6. Gaussian mixture model
  NM_mu_up <- model$NM_model$mu_up
  names(NM_mu_up) <- paste0("mu_up_", 1:12)
  NM_mu_dw <- model$NM_model$mu_dw
  names(NM_mu_dw) <- paste0("mu_dw_", 1:12)
  NM_sd_up <- model$NM_model$sd_up
  names(NM_sd_up) <- paste0("sd_up_", 1:12)
  NM_sd_dw <- model$NM_model$sd_dw
  names(NM_sd_dw) <- paste0("sd_dw_", 1:12)
  NM_p_up <- model$NM_model$p_up
  names(NM_p_up) <- paste0("p_", 1:12)

  # Output list of parameter for each model
  params <- structure(
    list(
      location = dplyr::tibble(place = model$place,
                               target = model$target,
                               lat = model$coords$lat,
                               lon = model$coords$lon,
                               alt = model$coords$alt),
      params = dplyr::bind_cols(model$params[1:2]),
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
    ),
    class = c("solarParameters", "list")
  )

  if (as_tibble) {
    params <- dplyr::bind_cols(params)
  }
  return(params)
}


#' Update the parameters of a `solarModel` object
#' @rdname solarModel_filter
#' @export
solarModel_update <- function(model, params){

  # Update clear sky model
  model$seasonal_model_Ct$update(unlist(params$seasonal_model_Ct))
  # Update seasonal model
  model$seasonal_model_Yt$update(unlist(params$seasonal_model_Yt))
  # Update AR model
  model$AR_model_Yt$coefficients <- unlist(params$AR_model_Yt)
  # Update seasonal variance
  model$seasonal_variance$update(unlist(params$seasonal_variance))
  # Update GARCH variance
  model$GARCH$coef <- unlist(params$GARCH)
  model$GARCH$vol <- sqrt(model$GARCH$coef[1]/(1-sum(model$GARCH$coef[-1])))
  model$control$variance.model@model$fixed.pars <- unlist(params$GARCH)
  # Update Gaussian Mixture
  NM <- model$NM_model
  NM$mu_up <- unlist(params$NM_mu_up)
  NM$mu_dw <- unlist(params$NM_mu_dw)
  NM$sd_up <- unlist(params$NM_sd_up)
  NM$sd_dw <- unlist(params$NM_sd_dw)
  NM$p_up <- unlist(params$NM_p_up)
  NM$p_dw <- 1-NM$p_up
  # Update theoric moments
  NM$e_x <- NM$mu_up*NM$p_up + NM$mu_dw*NM$p_dw
  NM$v_x <- NM$p_up*NM$p_dw*(NM$mu_up-NM$mu_dw)^2 + NM$p_up*NM$sd_up^2 + NM$p_dw*NM$sd_dw^2
  # Update Gaussian Mixture model
  model$NM_model <- NM
  return(model)
}


#' Update Gaussian Mixture parameters for a given month
#'
#' @export
solarModel_update_GM <- function(model, params, nmonth){

  if (missing(params)) {
    warning("Params are missing")
    return(model)
  }
  if (missing(nmonth)) {
    warning("nmonth is missing")
    return(model)
  }

  model$NM_model[nmonth, "mu_up"] <- params[1]
  model$NM_model[nmonth, "mu_dw"] <- params[2]
  model$NM_model[nmonth, "sd_up"] <- params[3]
  model$NM_model[nmonth, "sd_dw"] <- params[4]
  model$NM_model[nmonth, "p_up"] <- params[5]
  model$NM_model[nmonth, "p_dw"] <- 1 - model$NM_model[nmonth, "p_up"]
  return(model)
}


#' Empiric Gaussian Mixture parameters
#'
#' @export
solarModel_empiric_GM <- function(model, match_moments = FALSE){

  # Select only train data
  data <- model$data
  data$weights <- model$outliers$weights
  data <- dplyr::filter(data, isTrain & weights != 0)
  data <- dplyr::select(data, Month, u_tilde, B)
  for(m in 1:12){
    # Monthly data
    df_m <- dplyr::filter(data, Month == m)
    u_up <- dplyr::filter(df_m, B == 1)$u_tilde
    u_dw <- dplyr::filter(df_m, B == 0)$u_tilde
    nm <- model$NM_model[m,]
    # Up moments
    n1 <- length(u_up)
    nm$mu_up <- mean(u_up, na.rm = TRUE)
    nm$sd_up <- sd(u_up, na.rm = TRUE)*sqrt(n1/(n1-1))
    nm$p_up <- mean(df_m$B, na.rm = TRUE)
    # Down moments
    n2 <- length(u_dw)
    nm$mu_dw <- mean(u_dw, na.rm = TRUE)
    nm$sd_dw <- sd(u_dw, na.rm = TRUE)*sqrt(n2/(n2-1))
    nm$p_dw <- 1 - nm$p_up
    df_m
    if (match_moments) {
      if (n1 >= n2) {
        nm$mu_dw <- (nm$e_x_hat - nm$mu_up*nm$p_up)/(1-nm$p_up)
        nm$sd_dw <- sqrt((nm$v_x_hat - nm$sd_up^2*nm$p_up)/(1-nm$p_up) - nm$p_up*(nm$mu_up - nm$mu_dw)^2)
      } else {
        nm$mu_up <- (nm$e_x_hat - nm$mu_dw*nm$p_dw)/(1-nm$p_dw)
        nm$sd_up <- sqrt((nm$v_x_hat - nm$sd_dw^2*nm$p_dw)/(1-nm$p_dw) - nm$p_dw*(nm$mu_up - nm$mu_dw)^2)
      }
    }
    # Extract parameters
    means = c(nm$mu_up, nm$mu_dw)
    sd = c(nm$sd_up, nm$sd_dw)
    p = c(nm$p_up, nm$p_dw)
    # Update log-likelihood
    pdf <- function(x) dmixnorm(x, means = means, sd = sd, p = p)
    nm$loss <- sum(log(pdf(df_m$u_tilde)), na.rm = TRUE)
    nm$nobs <- nrow(df_m)
    # Update moments
    nm$e_x <- sum(means*p)
    nm$v_x <- sum((means^2 + sd^2)*p) - nm$e_x^2
    # Update data
    params <- unlist(nm[1,c(2:7)])
    model <- solarModel_update_GM(model, params, m)
  }
  model$log_lik <- solarModel_loglik(model)
  model
}
