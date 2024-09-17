#' Control parameters for a `solarModel` object
#'
#' Control function for a solarModel
#'
#' @param clearsky.model list with control parameters, see \code{\link{control_clearskyModel}} for details.
#' @param mean.model a list of parameters. Available choices are:
#'\describe{
#'  \item{`seasonalOrder`}{An integer specifying the order of the seasonal component in the model. The default is `1`.}
#'  \item{`arOrder`}{An integer specifying the order of the autoregressive component in the model. The default is `2`.}
#'  \item{`include.intercept`}{When `TRUE` the intercept will be included in the AR model. The dafault if `FALSE`.}
#'  \item{`monthly.mean`}{When `TRUE` a set of 12 monthly means parameters will be computed from the deseasonalized time series to center it perfectly around zero.}
#'}
#' @param variance.model a list of parameters.
#' @param outliers_quantile quantile for outliers detection. If different from 0, the observations that are below the quantile at confidence levels `outliers_quantile` and
#' the observation above the quantile at confidence level 1-`outliers_quantile` will have a weight equal to zero and will be excluded from estimations.
#' @param threshold numeric, threshold for the estimation of alpha and beta.
#' @param quiet logical, when `TRUE` the function will not display any message.
#'
#' @rdname control_solarModel
#' @export
control_solarModel <- function(clearsky.model = control_clearskyModel(),
                               mean.model = list(seasonalOrder = 1, arOrder = 2, include.intercept = FALSE, monthly.mean = TRUE),
                               variance.model = list(seasonalOrder = 1, unconditional_variance = NA, match_moments = FALSE, monthly.mean = TRUE, abstol = 1e-20, maxit = 100),
                               threshold = 0.01, outliers_quantile = 0, quiet = FALSE){

  # Mean model default parameters
  mean_model = list(seasonalOrder = 1, arOrder = 2, include.intercept = FALSE, monthly.mean = FALSE)
  names_mean_model <- names(mean_model)
  for(name in names_mean_model){
    arg <- mean.model[[name]]
    if (!is.null(arg)) {
      mean_model[[name]] <- mean.model[[name]]
    }
  }
  # Variance model default parameters
  variance_model = list(seasonalOrder = 1, unconditional_variance = NA, match_moments = FALSE, monthly.mean = TRUE, abstol = 1e-20, maxit = 100)
  names_variance_model <- names(variance_model)
  for(name in names_variance_model){
    arg <- variance.model[[name]]
    if (!is.null(arg)) {
      variance_model[[name]] <- variance.model[[name]]
    }
  }

  structure(
    list(
      clearsky.model = clearsky.model,
      mean.model = mean_model,
      variance.model = variance_model,
      threshold = threshold,
      outliers_quantile = outliers_quantile,
      quiet = quiet
    ),
    class = c("control", "list")
  )
}

#' Specification function for a `solarModel`
#'
#' @param place character, name of an element in the `CAMS_data` list.
#' @param min_date character. Date in the format `YYYY-MM-DD`. Minimum date for the complete data. If `missing` will be used the minimum data available.
#' @param max_date character. Date in the format `YYYY-MM-DD`. Maximum date for the complete data. If `missing` will be used the maximum data available.
#' @param from character. Date in the format `YYYY-MM-DD`. Starting date to use for training data.
#' If `missing` will be used the minimum data available after filtering for `min_date`.
#' @param to character. Date in the format `YYYY-MM-DD`. Ending date to use for training data.
#' If `missing` will be used the maximum data available after filtering for `max_date`.
#' @param CAMS_data named list with radiation data for different locations.
#' @param control_model list with control parameters, see \code{\link{control_solarModel}} for details.
#'
#' @rdname solarModel_spec
#' @name solarModel_spec
#'
#' @export
solarModel_spec <- function(place, min_date, max_date, from, to, CAMS_data = solarr::CAMS_data, control_model = control_solarModel()){

  # Match a location in the dataset
  place <- match.arg(place, choices = names(CAMS_data), several.ok = FALSE)
  data <- CAMS_data[[place]]

  # Minimum date for the complete dataset
  if (missing(min_date) || is.null(min_date) || is.na(min_date)) {
    min_date <- min(data$date, na.rm = TRUE)
  } else {
    min_date <- as.Date(min_date)
  }
  # Maximum date for the complete dataset
  if (missing(max_date) || is.null(max_date) || is.na(max_date)) {
    max_date <- max(data$date, na.rm = TRUE)
  } else {
    max_date <- as.Date(max_date)
  }

  # Filter for min and max dates the complete dataset
  data <- dplyr::filter(data, date >= min_date & date <= max_date)

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

  # Extract seasonal data
  seasonal_data <- attr(data, "seasonal")
  # Add number of the day
  seasonal_data$n <- number_of_day(seasonal_data$date)

  structure(
    list(
      place = attr(data, "place"),
      coords = attr(data, "coords"),
      seasonal_data = seasonal_data,
      data = data,
      from = from,
      to = to,
      min_date = min_date,
      max_date = max_date,
      control = control_model
    ),
    class = c("solarModelSpec", "list")
  )
}

#' Fit a model for solar radiation
#'
#' @param spec an object with class `solarModelSpec`. See the function \code{\link{solarModel_spec}} for details.
#'
#' @rdname solarModel
#' @name solarModel
#' @examples
#' control <- control_solarModel(outliers_quantile = 0.0005)
#' spec <- solarModel_spec("Palermo", from="2005-01-01", to="2024-01-01", control_model = control)
#' model <- solarModel(spec)
#' @export
solarModel <- function(spec){

  # Extract control parameters
  control <- spec$control
  control_cs <- control$clearsky.model
  # Complete data
  data <- spec$data
  # Label for data used for estimation
  data <- dplyr::mutate(data, isTrain = ifelse(date >= spec$from & date <= spec$to, TRUE, FALSE))
  # - Compute weights for each observation
  n_train <- nrow(dplyr::filter(data, isTrain))
  data$weights <- c(rep(1/n_train, n_train), rep(0, nrow(data)-n_train))

  # Extract the seasonal model for Ct
  data_train <- dplyr::filter(data, isTrain)
  seasonal_data <- dplyr::select(spec$seasonal_data, n, H0)
  seasonal_model_Ct <- clearskySeasonal(data_train, seasonal_data, control_cs$method, control_cs)
  data$Ct <- seasonal_model_Ct$predict(data$n)
  # Detect and impute outliers
  outliers <- clearsky_outliers(data$GHI, data$Ct, date = data$date, quiet = control$clearsky.model$quiet)
  # Update GHI
  data$GHI <- outliers$x

  # 1) Solar Transform
  st <- solarTransform$new()
  # Compute risk drivers
  data$Xt <- st$GHI(data$GHI, data$Ct, inverse = TRUE)
  # Remove minimum and maximum of Xt from computations
  outliers$index <- c(outliers$index, which.min(data$Xt), which.max(data$Xt))
  # Transformation parameters
  params <- st$parameters(data$Xt, control$threshold)
  st$update(params$alpha, params$beta)
  # Compute Yt
  data$Yt <- st$Yt(data$Xt)
  # - Exclude outliers from computations
  idx_outliers_l <- which(data$Yt <= quantile(data$Yt, probs = control$outliers_quantile))
  idx_outliers_r <- which(data$Yt >= quantile(data$Yt, probs = 1-control$outliers_quantile))
  outliers$index <- unique(c(outliers$index, idx_outliers_l, idx_outliers_r))

  # - Rescale weights
  if (!purrr::is_empty(outliers$index)) {
    data$weights[outliers$index] <- 0
    data$weights <- data$weights/sum(data$weights)
  }

  # 2) Seasonal Mean
  data_train <- dplyr::filter(data, isTrain & weights != 0)
  # Seasonal model for Yt
  seasonal_model_Yt <- seasonalModel$new(formula = "Yt ~ 1", order = control$mean.model$seasonalOrder, data = data_train)
  # Fitted seasonal mean for Yt
  data$Yt_bar <- seasonal_model_Yt$predict(data$n)
  # Fitted deseasonalized Yt
  data$Yt_tilde <- data$Yt - data$Yt_bar
  # Unconditional mean for Yt_tilde
  if (control$mean.model$monthly.mean) {
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

  # 3) AR-seasonal-GARCH(1,1) model
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

  # 3-b) Seasonal variance
  data_train <- dplyr::filter(data, isTrain & weights != 0)
  seasonal_variance <- seasonalModel$new(formula = "I(eps^2) ~ 1", order = control$variance.model$seasonalOrder, data = data_train)
  # Fitted seasonal standard deviation
  data$sigma_bar <- sqrt(seasonal_variance$predict(n = data$n))
  # Fitted standardized residuals
  data$eps_tilde <- data$eps/data$sigma_bar
  # Update train data
  data_train <- dplyr::filter(data, isTrain & weights != 0)
  # Correction to ensure unitary variance
  seasonal_variance$update(seasonal_variance$coefficients*mean(data_train$eps_tilde^2, na.rm = TRUE))
  # Fitted seasonal standard deviation
  data$sigma_bar <- sqrt(seasonal_variance$predict(n = data$n))
  # Updated standardized residuals
  data$eps_tilde <- data$eps/data$sigma_bar

  # 3-c) GARCH model
  data_train <- dplyr::filter(data, isTrain)
  # Variance specification
  GARCH_spec <- rugarch::ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1,1), external.regressors = NULL),
    mean.model = list(armaOrder = c(0,0), include.mean = FALSE), distribution.model = "norm")
  # Fitted model
  GARCH_model <- rugarch::ugarchfit(data = data_train$eps_tilde, spec = GARCH_spec, out.sample = 0)
  # Specify unconditional variance
  sigma2_0 <- 1
  if (!is.na(control$variance.model$unconditional_variance)) {
    sigma2_0 <- control$variance.model$unconditional_variance
  }

  init_params <- c(0.4, 0.05, 0.3)
  for(i in 1:50) {
    # Alternative fit
    GARCH_alternative <- GARCH_optimize(params = init_params, x = data_train$eps_tilde,
                                        weights = data_train$weights, sigma2 = sigma2_0)
    check_convergence <- any(is.na(GARCH_alternative$fitted$sigma))
    if (check_convergence) {
      message("GARCH routine do not converge... iteration ", i, "/50")
    } else {
      if(i > 1) message("GARCH routine converged in iteration ", i, "/50")
      break
    }
    init_params <- c(1, runif(1, 0.01, 0.1), runif(1, 0.1, 0.5))
    init_params[1] <- 1-init_params[2]-init_params[3]
  }

  # Check Log-likelihood improvement
  if (GARCH_alternative$logLik > sum(-GARCH_model@fit$log.likelihoods)) {
    GARCH_model@fit$coef <- GARCH_alternative$params
    GARCH_model@fit$sigma <- GARCH_alternative$fitted$sigma
  }
  # Constraint to omega0 to match the unconditional variance
  GARCH_model@fit$coef[1] <- sigma2_0*(1-GARCH_model@fit$coef[2]-GARCH_model@fit$coef[3])
  # Update fitted variance and residuals
  GARCH_model_filter <- GARCH_filter(data$eps_tilde, params = GARCH_model@fit$coef, x0 = NA, sigma0 = GARCH_model@fit$sigma[1])
  # Update model fitted variance
  GARCH_model@fit$var <- GARCH_model_filter$sigma^2
  # Update model fitted std. deviation
  GARCH_model@fit$sigma <- GARCH_model_filter$sigma
  # Update model fitted residuals
  GARCH_model@fit$residuals <- GARCH_model_filter$u
  # Add fitted standard deviation
  data$sigma <- GARCH_model@fit$sigma
  # Add fitted final residuals
  data$u <- data$eps_tilde/data$sigma
  # Monthly correction
  if (control$variance.model$monthly.mean) {
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
  GM_model <- solar_monthly_mixture(data$u_tilde, date = data$date, weights = data$weights,
                                    match_moments = control$variance.model$match_moments,
                                    abstol = control$variance.model$abstol, maxit = control$variance.model$maxit)
  # Return also the fitted series of Bernoulli
  data <- dplyr::left_join(data, dplyr::select(GM_model$data, date, B), by = c("date"))
  # Store the parameters
  GM_model <- GM_model$model

  # Seasonal data by month and day for an year with 366 days
  reference_year <- unique(data$Year)[which(is_leap_year(paste0(unique(data$Year), "-01-01")))[1]]
  from_date <- as.Date(paste0(reference_year, "-01-01"))
  to_date <- as.Date(paste0(reference_year, "-12-31"))
  seasonal_data <- dplyr::filter(data, date >= from_date & date <= to_date)
  seasonal_data$Xt_bar <- st$Yt(seasonal_data$Yt_bar, inverse = TRUE)
  seasonal_data$GHI_bar <- st$GHI(seasonal_data$Xt_bar, seasonal_data$Ct)
  seasonal_data <- dplyr::select(seasonal_data, Month, Day, Ct, GHI_bar, Xt_bar, Yt_bar, sigma_bar)
  # Add extraterrestrial radiation (H0) and sun hours
  seasonal_data <- dplyr::left_join(seasonal_data, dplyr::select(spec$seasonal_data, Month, Day, H0, sun_hours), by = c("Month", "Day"))
  # Remove seasonal variables
  data <- dplyr::select(data, -Ct, -H0, -Yt_bar, -sigma_bar)

  # Structure model data
  structure(
    list(
      data = dplyr::select(data, -weights),
      place = spec$place,
      coords = spec$coords,
      seasonal_model_Ct = seasonal_model_Ct,
      seasonal_model_Yt = seasonal_model_Yt,
      AR_model_Yt = AR_model_Yt,
      seasonal_variance = seasonal_variance,
      GARCH = list(coef = GARCH_model@fit$coef, vol = sqrt(sigma2_0)),
      seasonal_data = seasonal_data,
      monthly_data = monthly_data,
      NM_model = GM_model,
      transform = st,
      params = params,
      log_lik = sum(GM_model$loss),
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
    ),
    class = c("solarModel", "list")
  )
}


#' Monthly Gaussian mixture with two components
#'
#' @param x arg
#' @param date arg
#' @param weights arg
#' @param match_moments arg
#' @param ... arg
#'
#' @rdname solar_monthly_mixture
#' @name solar_monthly_mixture
#' @export
solar_monthly_mixture <- function(x, date, weights, match_moments = FALSE, prior_p, ...){

  data <- dplyr::tibble(date = date, Month = lubridate::month(date), eps = x, w = weights)
  # Gaussian Mixture parameters
  GM_model <- list()
  # Monthly data
  data_months <- list()
  for(m in 1:12){
    data_months[[m]] <- dplyr::filter(data, Month == m)
    w <- data_months[[m]]$w
    # Monthly data
    eps <- data_months[[m]]$eps
    # Initial parameters
    mu_0 <- quantile(eps, c(0.2, 0.8), na.rm = TRUE)
    sd_0 <- sd(eps, na.rm = TRUE)
    init_means <- c(mu1 = mu_0[1], mu2 = mu_0[2])
    init_sd <- c(sd1 = sd_0, sd2 = sd_0)
    init_p <- c(0.5, 0.5)
    if (!missing(prior_p)){
      pri_p <- unlist(prior_p[m,])
    } else {
      pri_p <- rep(NA, 2)
    }
    # Fitted model
    gm <- gaussianMixture(eps, means = init_means, components = 2, sd = init_sd, p = init_p, weights = w,
                          prior_p = pri_p, ...)

    # Compute the theoretical expected value
    e_u <- sum(gm$par$mean*gm$par$p)
    # Compute the theoric variance
    v_u <- sum((gm$par$mean^2 + gm$par$sd^2)*gm$par$p) - e_u^2
    # Compute the sample mean
    e_u_hat <- mean(eps[w!=0], na.rm = TRUE)
    # Compute the sample variance
    v_u_hat <- var(eps[w!=0], na.rm = TRUE)
    # Fitted parameters
    df_par <- dplyr::bind_cols(purrr::map(gm$par, ~dplyr::bind_rows(.x)))
    # Add fitted Bernoulli series
    if (df_par$mu1 > df_par$mu2) {
      data_months[[m]]$B <- gm$fitted$B1
    } else {
      data_months[[m]]$B <- gm$fitted$B2
    }
    # Monthly data
    GM_model[[m]] <- dplyr::tibble(Month = m,
                                   df_par,
                                   loss = gm$log_lik,
                                   nobs = length(eps[w!=0]),
                                   e_x = e_u,
                                   v_x = v_u,
                                   e_x_hat = e_u_hat,
                                   v_x_hat = v_u_hat)
  }

  # Reorder the variables
  GM_model <- dplyr::bind_rows(GM_model)
  GM_model <- dplyr::mutate(GM_model,
                            mu_up = dplyr::case_when(
                              mu1 > mu2 ~ mu1,
                              TRUE ~ mu2),
                            mu_dw = dplyr::case_when(
                              mu1 > mu2 ~ mu2,
                              TRUE ~ mu1),
                            sd_up = dplyr::case_when(
                              mu1 > mu2 ~ sd1,
                              TRUE ~ sd2),
                            sd_dw = dplyr::case_when(
                              mu1 > mu2 ~ sd2,
                              TRUE ~ sd1),
                            p_up = dplyr::case_when(
                              mu1 > mu2 ~ p1,
                              TRUE ~ p2),
                            p_dw = 1 - p_up)

  GM_model <- dplyr::select(GM_model, Month, mu_up:p_dw, loss, nobs, e_x, v_x, e_x_hat, v_x_hat)
  # Fitted series
  data_months <- dplyr::bind_rows(data_months)

  structure(
    list(
      model = GM_model,
      data = data_months
    )
  )
}


#' Simulate trajectories
#'
#' Simulate trajectories of solar radiation with a `solarModel` object.
#'
#' @inheritParams solarModel_scenario
#'
#' @examples
#' spec <- solarModel_spec("Ferrara", from="2005-01-01", to="2020-01-01")
#' model <- solarModel(spec)
#' solarModel_simulate(model, from = "2010-01-01", to = "2010-12-31", nsim = 1)
#'
#' @rdname solarModel_simulate
#' @name solarModel_simulate
#' @export
solarModel_simulate <- function(model, from = "2010-01-01", to = "2010-12-31", nsim = 1, lambda = 0, seed = 1, exclude_known = FALSE, quiet = FALSE){

  # Number of lags to consider
  i_start <- model$control$mean.model$arOrder+1
  data <- model$data
  place <- model$place

  # AR(2) model (GHI)
  AR_model_Yt <- model$AR_model_Yt
  NM_model <- model$NM_model
  GARCH <- model$GARCH
  params <- model$params

  esscher_probability <- function(params = c(0,0,1,1,0.5), df_n, h = 0){
    params <- list(
      mu_up = df_n$Yt_bar + df_n$Yt_tilde_uncond + df_n$Yt_tilde_hat + df_n$sigma*df_n$sigma_bar*params[1],
      mu_dw = df_n$Yt_bar + df_n$Yt_tilde_uncond + df_n$Yt_tilde_hat + df_n$sigma*df_n$sigma_bar*params[2],
      sd_up = params[3]*df_n$sigma_bar*df_n$sigma,
      sd_dw = params[4]*df_n$sigma_bar*df_n$sigma,
      p_up = params[5]
    )
    params <- unlist(params)
    num <- params[5]*exp(h*params[1] + (h^2*params[3]^2)/2)
    den <- (1-params[5])*exp(h*params[2] + (h^2*params[4]^2)/2)
    num/(num + den)
  }

  # Initial date
  from <- as.Date(from)
  # End date
  to <- as.Date(to)
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
    df_emp <- dplyr::select(df_emp, -dplyr::any_of(colnames(seasonal_data)[-c(1:2)]))
    df_emp <- dplyr::mutate(df_emp,
                            Year = lubridate::year(date),
                            Month = lubridate::month(date),
                            Day = lubridate::day(date))
    df_emp$n <- solarr::number_of_day(df_emp$date)
  } else {
    msg <- paste0("The maximum date for starting a simulation is: ", max_date_from)
    if (!quiet) warning(msg)
    return(model)
  }
  # Add seasonal data
  seasonal_data <- dplyr::select(model$seasonal_data, Month, Day, Ct, Yt_bar, GHI_bar, sigma_bar)
  df_emp <- dplyr::left_join(df_emp, seasonal_data, by = c("Month", "Day"))
  df_emp <- dplyr::left_join(df_emp, model$monthly_data, by = c("Month"))

  # Garch parameters
  omega0 <- GARCH$coef[1]
  omega1 <- GARCH$coef[2]
  omega2 <- GARCH$coef[3]

  # Initialize the template dataset
  df_sim_init <- dplyr::left_join(df_emp, NM_model[, c("Month","mu_up", "mu_dw", "sd_up", "sd_dw", "p_up")], by = "Month")
  # Filter df_emp to be in [from - to] dates
  if (exclude_known) {
    df_emp <- dplyr::filter(df_emp, date >= from & date <= to)
  }
  # Initialize Market risk premium
  df_sim_init$lambda <- lambda
  # Remove redundant variables
  df_sim_init <- dplyr::select(df_sim_init, -clearsky, -isTrain)

  j <- 1
  simulations <- list()
  for(j in 1:nsim){
    # Initialize dataset for storing the simulation
    df_sim <- df_sim_init
    # Simulate Normal mixture (ut)
    df_sim$u <- rnorm(nrow(df_sim), mean = 0, sd = 1)
    if (!quiet) message("Simulation: ", j, "/", nsim, " (", round(j/nsim*100, 4), " %) \r", appendLF = FALSE)
    set.seed(seed)
    i <- i_start
    for(i in i_start:nrow(df_sim)){

      # Simulated GARCH standard deviation
      df_sim$sigma[i] <- sqrt(omega0 + omega1*df_sim$eps_tilde[i-1]^2 + omega2*df_sim$sigma[i-1]^2)
      # Simulated Yt_tilde
      df_sim$Yt_tilde_hat[i] <- predict(AR_model_Yt, newdata = df_sim[(i-i_start):i,])[i_start]

      if (lambda != 0) {
        params <- c(df_sim$mu_up[i], df_sim$mu_dw[i], df_sim$sd_up[i], df_sim$sd_dw[i], df_sim$p_up[i])
        df_sim$p_up[i] <- esscher_probability(params, df_n = df_sim[i,], h = lambda)
      }

      # Bernoulli jump
      df_sim$B[i] <- rbinom(1, 1, df_sim$p_up[i])
      # Simulated Esscher parameter
      df_sim$lambda[i] <- df_sim$B[i]*(df_sim$sd_up[i]^2)*df_sim$lambda[i] + (1-df_sim$B[i])*(df_sim$sd_dw[i]^2)*df_sim$lambda[i]
      df_sim$lambda[i] <- df_sim$lambda[i]*(df_sim$sigma[i]*df_sim$sigma_bar[i])^2

      # Simulated normal mixture
      df_sim$u_tilde[i] <- (df_sim$mu_up[i] + df_sim$sd_up[i]*df_sim$u[i])*df_sim$B[i] + (1-df_sim$B[i])*(df_sim$mu_dw[i] + df_sim$sd_dw[i]*df_sim$u[i])
      # Simulated standardized monthly residuals
      df_sim$u[i] <- df_sim$sigma_m[i]*df_sim$u_tilde[i]
      # Simulated standardized residuals
      df_sim$eps_tilde[i] <- df_sim$sigma[i]*df_sim$u[i]
      # Simulated AR residuals
      df_sim$eps[i] <- df_sim$eps_tilde[i]*df_sim$sigma_bar[i]
      # Simulated Yt_tilde
      df_sim$Yt_tilde[i] <- df_sim$Yt_tilde_hat[i] + df_sim$eps[i]
      # Simulated Yt
      df_sim$Yt[i] <- df_sim$Yt_bar[i] + df_sim$Yt_tilde[i] + df_sim$Yt_tilde_uncond[i] + df_sim$lambda[i]
      # Simulated Xt
      df_sim$Xt[i] <- model$transform$Yt(df_sim$Yt[i], inverse = TRUE)
      # Simulated GHI
      df_sim$GHI[i] <- model$transform$GHI(df_sim$Xt[i], df_sim$Ct[i])
    }
    # Remove redundant variables
    df_sim <- dplyr::select(df_sim, -mu_up, -mu_dw, -sd_up, -sd_dw, -p_up, -Ct, -Yt_bar, -sigma_bar, -sigma_m, -Yt_tilde_hat, -Yt_tilde_uncond, -lambda)
    # Remove initial values
    if (exclude_known) {
      df_sim <- dplyr::filter(df_sim, date >= from & date <= to)
    }
    # Store simulations
    simulations[[j]] <- dplyr::bind_cols(seed = seed, df_sim)
    # Update seed
    seed <- seed + j
  }

  # Remove redundant variables
  df_emp <- dplyr::select(df_emp, -clearsky, -Ct, -Yt_bar, -sigma_bar, -sigma_m, -Yt_tilde_hat, -Yt_tilde_uncond, -isTrain)

  structure(
    list(
      sim = simulations,
      emp = df_emp,
      params = list(seed = seed, from = from, to = to, nsim = nsim, lambda = lambda)
    ),
    class = c("solarModelSimulation", "list")
  )
}

#' Simulate multiple scenarios
#'
#' Simulate multiple scenarios of solar radiation with a `solarModel` object.
#'
#' @param model object with the class `solarModel`. See the function \code{\link{solarModel}} for details.
#' @param from character, start Date for simulations in the format `YYYY-MM-DD`.
#' @param to character, end Date for simulations in the format `YYYY-MM-DD`.
#' @param by character, steps for multiple scenarios, e.g. `1 day` (day-ahead simulations), `15 days`, `1 month`, `3 months`, ecc.
#' For each step are simulated `nsim` scenarios.
#' @param nsim integer, number of simulations.
#' @param lambda numeric, Esscher parameter.
#' @param seed scalar integer, starting random seed.
#' @param quiet logical
#'
#' @examples
#' spec <- solarModel_spec("Ferrara", from="2005-01-01", to="2020-01-01")
#' model <- solarModel(spec)
#' solarModel_scenario(model, from = "2010-01-01", to = "2010-12-31", nsim = 2, by = "1 month")
#'
#' @rdname solarModel_scenario
#' @name solarModel_scenario
#' @export
solarModel_scenario <- function(model, from = "2010-01-01", to = "2010-12-31", by = "1 month", nsim = 1, lambda = 0, seed = 1, quiet = FALSE){

  idx_date <- seq.Date(as.Date(from), as.Date(to), by = by)
  scenarios <- list()
  df_emp <- dplyr::tibble()
  n_scenario <- length(idx_date)
  j <- 2
  for(j in 2:n_scenario){
    if (!quiet) {
      # To report progress
      pb <- txtProgressBar(min = 1,            # Minimum value of the progress bar
                           max = n_scenario,   # Maximum value of the progress bar
                           style = 3,          # Progress bar style (also available style = 1 and style = 2)
                           width = 50,         # Progress bar width. Defaults to getOption("width")
                           char = "#")
      setTxtProgressBar(pb, j)
    }
    sim <- solarModel_simulate(model, from = idx_date[j-1], to = idx_date[j]-1, nsim = nsim,
                               lambda = lambda, seed = seed, exclude_known = TRUE, quiet = TRUE)
    df_emp <- dplyr::bind_rows(df_emp, sim$emp)
    scenarios[[j]] <- dplyr::bind_rows(sim$sim)
    seed <- sim$params$seed + 1
  }
  if (!quiet) close(pb)
  df_emp <- df_emp[!duplicated(df_emp),]
  df_sim <- dplyr::bind_rows(scenarios) %>%
    dplyr::group_by(date, Year, Month, Day, n)%>%
    tidyr::nest() %>%
    dplyr::ungroup()

  structure(
    list(
      sim = df_sim,
      emp = df_emp,
      seed = seed
    ),
    class = c("scenarioSolarModel", "list")
  )
}

#' Extract the parameters of a `solarModel`
#'
#' @param model object with the class `solarModel`. See the function \code{\link{solarModel}} for details.
#' @param as_tibble logical, when `TRUE` the output will be converted in a tibble.
#' @examples
#' solarModel_parameters(Bologna)
#'
#' @return a named list with all the parameters
#' @export
solarModel_parameters <- function(model, as_tibble = FALSE){

  control <- model$control

  # 1. Clears sky model
  seasonal_model_Ct <- as.list(model$seasonal_model_Ct$coefficients)
  if (control$clearsky.model$include.intercept) {
    names(seasonal_model_Ct) <- c("delta_0", "delta_1", paste0("delta_", 2:(2*control$clearsky.model$order + 1)))
  } else {
    names(seasonal_model_Ct) <- c(paste0("delta_1", "delta_", 2:(2*control$clearsky.model$order + 1)))
  }
  # 2. Seasonal mean model
  seasonal_model_Yt <- as.list(model$seasonal_model_Yt$coefficients)
  names(seasonal_model_Yt) <- c("a_0", paste0("a_", 1:(2*control$mean.model$seasonalOrder)))
  # 3. AR model
  AR_model_Yt <- as.list(model$AR_model_Yt$coefficients)
  if (control$mean.model$include.intercept) {
    names(AR_model_Yt) <- c("phi_0", paste0("phi_", 1:control$mean.model$arOrder))
  } else {
    names(AR_model_Yt) <- c(paste0("phi_", 1:control$mean.model$arOrder))
  }
  # 4. Seasonal variance model
  seasonal_variance <- as.list(model$seasonal_variance$coefficients)
  names(seasonal_variance) <- c("c_0", paste0("c_", 1:(2*control$variance.model$seasonalOrder)))
  # 5. GARCH(1,1) variance model
  GARCH <- as.list(model$GARCH$coef)
  names(GARCH) <- paste0("omega_", 0:2)
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
                               lat = model$coords$lat,
                               lon = model$coords$lon),
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
#'
#' @param model `solarModel` object
#' @param params named list of parameters. See the function \code{\link{solarModel_parameters}} to structure the list of new parameters.
#' @export
solarModel_update <- function(model, params){

  # Extract control parameters
  control <- model$control
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
  model$GARCH$vol <- sqrt(model$GARCH$coef[1]/(1-sum(model$GARCH$coef[2:3])))
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
  # Update empirical moments
  moments_u <- model$data %>%
    dplyr::group_by(Month) %>%
    dplyr::summarise(e_x_hat = mean(u_tilde), v_x_hat = var(u_tilde))
  NM$e_x_hat <- moments_u$e_x_hat
  NM$v_x_hat <- moments_u$v_x_hat
  # Update Gaussian Mixture model
  model$NM_model <- NM

  return(model)
}


#' Update the time series of a `solarModel` object
#'
#' @param model `solarModel` object
#'
#' @export
solarModel_filter <- function(model){

  # Extract control parameters
  control <- model$control
  # Complete data
  data <- dplyr::bind_cols(model$data, weights = model$outliers$weights)
  data <- dplyr::mutate(data, isTrain = ifelse(is.na(isTrain), FALSE, isTrain))
  # 0) Clear sky and risk drivers
  # Predict seasonal clear sky
  data$Ct <- model$seasonal_model_Ct$predict(data$n)
  # Risk Driver (Xt is computed here)
  outliers <- clearsky_outliers(data$GHI, data$Ct, data$date, quiet = control$clearsky.model$quiet)
  # Update the dataset with imputed values
  data$GHI <- outliers$x
  # 1) Solar Transform
  st <- solarTransform$new()
  # Compute risk drivers
  data$Xt <- st$GHI(data$GHI, data$Ct, inverse = TRUE)
  # Remove minimum and maximum of Xt from computations
  outliers$index <- c(model$outliers$index, which.min(data$Xt), which.max(data$Xt))
  # Transformation parameters
  params <- st$parameters(data$Xt, control$threshold)
  # Update solar transform
  st$update(params$alpha, params$beta)
  # Compute Yt
  data$Yt <- st$Yt(data$Xt)
  # Fitted seasonal mean for Yt
  data$Yt_bar <- model$seasonal_model_Yt$predict(data$n)
  # Fitted deseasonalized Yt
  data$Yt_tilde <- data$Yt - data$Yt_bar
  # Unconditional mean for Yt_tilde
  if (control$mean.model$monthly.mean) {
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
  # Update fitted variance and residuals
  GARCH_model_filter <- GARCH_filter(data$eps_tilde, params = model$GARCH$coef)
  # Add fitted standard deviation
  data$sigma <- GARCH_model_filter$sigma
  # Add fitted final residuals
  data$u <- data$eps_tilde/data$sigma
  # Monthly correction
  if (control$variance.model$monthly.mean) {
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
  seasonal_data$Xt_bar <- st$Yt(seasonal_data$Yt_bar, inverse = TRUE)
  seasonal_data$GHI_bar <- st$GHI(seasonal_data$Xt_bar, seasonal_data$Ct)
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

#' Compute the log-likelihood of a `solarModel` object
#'
#' @param model `solarModel` object
#' @param nmonths months to consider
#' @export
solarModel_loglik <- function(model, target = c("Yt", "GHI"), nmonths = 1:12){

  target <- match.arg(target, choices = c("Yt", "GHI"))

  data <- dplyr::filter(model$data, Month %in% nmonths)
  data <- dplyr::left_join(data, model$seasonal_data, by = c("Month", "Day"))
  data <- dplyr::left_join(data, model$monthly_data, by = c("Month"))
  data <- dplyr::left_join(data, model$NM_model[,c(1:7)], by = c("Month"))
  data <- dplyr::mutate(data,
                        mu_up = Yt_bar + Yt_tilde_hat + Yt_tilde_uncond + sigma*sigma_bar*sigma_m*mu_up,
                        mu_dw = Yt_bar + Yt_tilde_hat + Yt_tilde_uncond + sigma*sigma_bar*sigma_m*mu_dw,
                        sd_up = sigma*sd_up*sigma_bar*sigma_m,
                        sd_dw = sigma*sd_dw*sigma_bar*sigma_m)

  data$log_lik <- 0
  for(i in 1:nrow(data)){
    df_n <- data[i,]
    # Conditional mixture Pdf
    pdf_Yt <- desscherMixture(means = c(df_n$mu_up, df_n$mu_dw), sd = c(df_n$sd_up, df_n$sd_dw), p = c(df_n$p_up, 1-df_n$p_up))
    if (target == "GHI"){
      data$log_lik[i] <- log(dsolarGHI(df_n$GHI, df_n$Ct, model$params$alpha, model$params$beta, pdf_Yt))
    } else if (target == "Yt") {
      data$log_lik[i] <- pdf_Yt(df_n$Yt, log = TRUE)
    }
  }
  model$log_lik <- sum(data$log_lik, na.rm = TRUE)
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
  m <- 11
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
    pdf <- dmixnorm(means = means, sd = sd, p = p)
    nm$loss <- sum(pdf(df_m$u_tilde, log = TRUE), na.rm = TRUE)
    nm$nobs <- nrow(df_m)
    # Update moments
    nm$e_x <- sum(means*p)
    nm$v_x <- sum((means^2 + sd^2)*p) - nm$e_x^2
    # Update data
    params <- unlist(nm[1,c(2:7)]) # + unlist(model$NM_model[1,c(2:7)]))/2
    model <- solarModel_update_GM(model, params, m)
  }
  model$log_lik <- sum(model$NM_model$loss)
  model
}

#' Calibrator for solar Options
#'
#' @export
solarModel_calibrator <- function(model, nmonths = 1:12, control_options = control_solarOption()){

  # Historical monthly payoff
  control_call <- control_options
  control_call$put <- FALSE
  payoff_call <- solarOption_historical(model, control_options = control_call)$payoff_month$premium
  control_put <- control_options
  control_put$put <- TRUE
  payoff_put <- solarOption_historical(model, control_options = control_put)$payoff_month$premium

  min_sd <- 0 # min(model$NM_model$sd_up)*0.9
  # Loss function
  loss_function <- function(params, model, nmonth){

    if(params[3] < min_sd | params[4] < min_sd | params[5] > 0.8 | params[5] < 0.2){
      return(NA)
    }

    # Updated model
    upd_mod <- solarModel_update_GM(model, params, nmonth)
    # Model price
    mod_call <- solarOption_model(upd_mod, nmonth = nmonth, control_options = control_call)$payoff_year$premium
    mod_put <- solarOption_model(upd_mod, nmonth = nmonth, control_options = control_put)$payoff_year$premium
    loss <- abs(mod_call - payoff_call[nmonth]) + abs(mod_put - payoff_put[nmonth])
    message("Loss: ", loss, " Params: ", format(params, digits = 3), "\r", appendLF = FALSE)
    loss <- -loss^2
    return(loss^2)
  }

  # Monthly calibration
  model_cal <- model
  for(nmonth in nmonths){
    message("------------------------------------ Month: ", nmonth, " ------------------------------------ ")
    params <- unlist(model$NM_model[nmonth, c(2:6)])
    opt <- optim(params, loss_function, model = model, nmonth = nmonth)
    model_cal <- solarModel_update_GM(model_cal, opt$par, nmonth)
  }
  return(model_cal)
}

#' Stationarity and distribution test (Gaussian mixture) for a `solarModel`
#'
#' @export
solarModel_test_residuals <- function(model, nrep = 500, ci = 0.05, min_quantile = 0.015, max_quantile = 0.985, seed = 1){

  distribution_test <- list()
  stationary_test <- list()
  for(nmonth in 1:12){
    # Residuals
    x <- dplyr::filter(model$data, Month == nmonth)$u_tilde
    # Gaussian mixture parameters
    nm <- model$NM_model[nmonth,]
    means = c(nm$mu_up,  nm$mu_dw)
    sd = c(nm$sd_up,  nm$sd_dw)
    p = c(nm$p_up,  nm$p_dw)

    # Mixture CDF
    cdf_GM <- function(x) pmixnorm(x, means = means, sd = sd, p = p)
    # Distribution test
    distribution_test[[nmonth]] <- dplyr::bind_cols(Month = nmonth,
                                                    ks_test(x, cdf_GM, ci = ci, min_quantile = min_quantile, max_quantile = max_quantile))
    # Stationary test
    iid_test <- function(seed) ks_ts_test(x, ci = ci, min_quantile = min_quantile, max_quantile = max_quantile, seed = seed)
    iid_tests <- purrr::map_df(1:nrep, ~iid_test(.x))
    Rejected <- mean(iid_tests$H0 == "Non-Rejected")
    stationary_test[[nmonth]] <- dplyr::bind_cols(Month = nmonth, Rejected = 1-Rejected, Non_Rej = Rejected)
    stationary_test[[nmonth]] <- dplyr::bind_cols(stationary_test[[nmonth]], iid_test(seed))
  }

  structure(
    list(
      distribution = dplyr::bind_rows(distribution_test),
      stationary = dplyr::bind_rows(stationary_test)
    )
  )
}
