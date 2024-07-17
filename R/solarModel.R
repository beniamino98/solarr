#' Control parameters for a `solarModel` object
#'
#' @param clearsky.model list with control parameters, see `clearskyModel_control()`.
#' @param mean.model a list of parameters.
#' @param variance.model a list of parameters.
#' @param threshold Threshold for the estimation of alpha and beta.
#' @param quiet logical, when `TRUE` the function will not display any message.
#'
#' @rdname solarModel_control
#' @export
solarModel_control <- function(clearsky.model = clearskyModel_control(),
                               mean.model = list(seasonalOrder = 1, arOrder = 2, include.intercept = FALSE),
                               variance.model = list(seasonalOrder = 1, match_moments = FALSE, algo = "em"),
                               threshold = 0.001, quiet = FALSE){

  # Mean model default parameters
  mean_model = list(seasonalOrder = 1, arOrder = 2, include.intercept = FALSE)
  names_mean_model <- names(mean_model)
  for(name in names_mean_model){
    arg <- mean.model[[name]]
    if (!is.null(arg)) {
      mean_model[[name]] <- mean.model[[name]]
    }
  }
  # Variance model default parameters
  variance_model = list(seasonalOrder = 1, match_moments = FALSE, algo = "em")
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
      quiet = quiet
    ),
    class = c("control", "list")
  )
}

#' Specification for a `solarModel` object
#'
#' @param place character, name for the selected location in `CAMS_data` list.
#' @param year_max integer, maximum year in the dataset
#' @param from character. Date in the format `YYYY-MM-DD`. Minimum date in the data in `CAMS_data`. If `NULL` will be used the maximum available.
#' @param to character. Date in the format `YYYY-MM-DD`. Maximum date in the data in `CAMS_data`. If `NULL` will be used the maximum available.
#' @param CAMS_data list with radiation data for different locations.
#' @param control list with control parameters, see `control_solarModel()`.
#'
#' @rdname solarModel_spec
#' @name solarModel_spec
#'
#' @export
solarModel_spec <- function(place, year_max = NULL, from = NULL, to = NULL, CAMS_data = solarr::CAMS_data, control = solarModel_control()){

  # Match a location in the dataset
  place <- match.arg(place, choices = names(CAMS_data), several.ok = FALSE)
  data <- CAMS_data[[place]]

  # Filter for minimum date
  if (!is.null(from)) {
    data <- dplyr::filter(data, date >= from)
  }
  # Filter for maximum date
  if (!is.null(to)) {
    data <- dplyr::filter(data, date <= to)
  }
  # Filter for maximum year
  if (is.null(year_max)) {
    year_max <- lubridate::year(Sys.Date())
  }
  data <- dplyr::filter(data, Year <= year_max)
  seasonal_data <- attr(data, "seasonal")
  seasonal_data$n <- number_of_day(seasonal_data$date)

  structure(
    list(
      place = attr(data, "place"),
      coords = attr(data, "coords"),
      seasonal_data = seasonal_data,
      data = data,
      control = control
    ),
    class = c("Location", "list")
  )
}

#' Fit a model for solar radiation
#'
#' @param object object with class `solarModel`. See the function `solarModel_spec()`. For example  `solarModel_spec("Bologna")`.
#' @param ... additional parameters for the function `fit_dnorm_mix_em()`.
#' @rdname solarModel_fit
#' @name solarModel_fit
#'
#' @export
solarModel_fit <- function(object, ...){

  # Extract data from object
  data <- object$data  # dataset
  control <- object$control

  # Risk Driver: Xt is computed here
  obj <- clearskyModel_fit(data, seasonal_data = dplyr::select(object$seasonal_data, n, H0),
                           control = control$clearsky.model)
  # Update the dataset with imputed values
  data <- obj$data
  # Extract the seasonal model for Ct
  seasonal_model_Ct <- obj$seasonal_model_Ct

  # ---- 1) Solar Transform ----
  # Upper and lower bounds: alpha_ and beta_
  epsilon <- min(data$Xt)*control$threshold
  min_Xt <- min(data$Xt)
  max_Xt <- max(data$Xt)
  alpha_ <- min_Xt - epsilon
  beta_ <- max_Xt - min_Xt + 2*epsilon
  # Vector of parameters
  params <- list(alpha = alpha_, beta = beta_, min_Xt = min_Xt, max_Xt = max_Xt, epsilon = epsilon)

  # Transformation functions
  Xt = function(x) alpha_ + beta_*exp(-exp(x))
  Yt = function(x) log(log(beta_) - log(x - alpha_))
  GHI = function(Ct, x) Ct*(1 - x)
  nday <- function(n) purrr::map_int(n, ~ifelse(is.character(.x) | lubridate::is.Date(.x), number_of_day(.x), .x))
  Ct = function(n) {clearskyModel_predict(seasonal_model_Ct, n = nday(n))}
  Yt_bar = function(n) {seasonalModel_predict(seasonal_model_Yt, n = nday(n))}
  GHI_bar = function(n) {Ct(n)*(1 - Xt(Yt_bar(n)))}

  # ---- 2) Seasonal Mean ----
  # Compute Yt
  data$Yt <- Yt(data$Xt)
  # Impute outliers
  data$Yt[data$Yt < quantile(data$Yt, probs = 0.001)] <- quantile(data$Yt, probs = 0.001)
  data$Xt <- Xt(data$Yt)
  data$GHI <- GHI(data$Ct, data$Xt)
  # Seasonal model for Yt
  seasonal_model_Yt <- seasonalModel_fit(formula = "Yt ~ 1", order = control$mean.model$seasonalOrder, period = 365, data = data)
  # Fitted seasonal mean for Yt
  data$Yt_bar <- seasonalModel_predict(seasonal_model_Yt, n = data$n)
  # Fitted deseasonalized Yt
  data$Yt_tilde <- data$Yt - data$Yt_bar
  # Fitted seasonal mean (GHI)
  data$GHI_bar <- GHI(data$Ct, Xt(data$Yt_bar))

  # ---- 3) AR model ----
  # AR with flexible formula for multiple orders
  AR_formula_Yt <- "Yt_tilde ~ "
  if (control$mean.model$arOrder > 0){
    for(i in 1:control$mean.model$arOrder){
      AR_formula_Yt <- paste0(AR_formula_Yt, " + I(dplyr::lag(Yt_tilde,", i, "))")
    }
  }
  # Control for intercept
  AR_formula_Yt <- paste0(AR_formula_Yt, ifelse(control$mean.model$include.intercept, "", "-1"))
  # Fitted AR model
  AR_model_Yt <- lm(formula = as.formula(AR_formula_Yt), data = data)
  # Fitted Yt_tilde
  data$Yt_tilde_hat <- predict.lm(AR_model_Yt, newdata = data)
  # Initial values as the real ones
  if (control$mean.model$arOrder > 0) {
    data$Yt_tilde_hat[1:control$mean.model$arOrder] <- data$Yt_tilde[1:control$mean.model$arOrder]
  }
  # Fitted Yt
  data$Yt_hat <- data$Yt_tilde_hat + data$Yt_bar
  # Fitted risk driver
  data$Xt_hat <- Xt(data$Yt_hat)
  # Fitted GHI
  data$GHI_hat <- GHI(data$Ct, data$Xt_hat)
  # Fitted residuals
  data$eps <- data$Yt_tilde - data$Yt_tilde_hat

  # ---- 4) Seasonal variance ----
  # Fitted model
  seasonal_variance <- seasonalModel_fit(formula = "I(eps^2) ~ 1", order = control$variance.model$seasonalOrder, period = 365, data = data)
  # Fitted seasonal standard deviation
  data$sigma_bar <- sqrt(seasonalModel_predict(seasonal_variance, n = data$n))
  # Fitted standardized residuals
  data$eps_tilde <- data$eps/data$sigma_bar

  # ---- 5) GARCH variance ----
  # Variance specification
  GARCH_spec <- rugarch::ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1,1), external.regressors = NULL),
    mean.model = list(armaOrder = c(0,0), include.mean = FALSE), distribution.model = "norm")
  # Fitted model
  GARCH_model <- rugarch::ugarchfit(data = data$eps_tilde, spec = GARCH_spec, out.sample = 0)
  # GARCH_model@fit$coef[1] <- 1 - GARCH_model@fit$coef[2] - GARCH_model@fit$coef[3]
  # Fitted variance
  data$sigma2 <- GARCH_model@fit$var
  # Fitted standard deviation
  data$sigma <- GARCH_model@fit$sigma
  # Fitted final residuals
  data$ut <- data$eps_tilde/data$sigma

  # ---- 6) Gaussian Mixture ----
  NM_model <- fit_dnorm_mix_monthly(x = data$ut, date = data$date, match_moments = control$variance.model$match_moments,
                                    algo = control$variance.model$algo, ...)

  # Seasonal data by month and day for an year with 366 days
  seasonal_data <- dplyr::filter(data, date >= as.Date("2016-01-01") & date <= as.Date("2016-12-31"))
  seasonal_data$Xt_bar <- Xt(seasonal_data$Yt_bar)
  seasonal_data <- dplyr::select(seasonal_data, Month, Day, Ct, GHI_bar, Xt_bar, Yt_bar, sigma_bar)
  seasonal_data <- dplyr::left_join(seasonal_data, object$seasonal_data, by = c("Month", "Day"))

  # Creat new object data
  structure(
    list(
      data = data,
      place = object$place,
      coords = object$coords,
      seasonal_model_Ct = seasonal_model_Ct,
      seasonal_model_Yt = seasonal_model_Yt,
      AR_model_Yt = AR_model_Yt,
      seasonal_variance = seasonal_variance,
      GARCH = list(model = GARCH_model,
                   spec = GARCH_spec,
                   vol = sqrt(GARCH_model@fit$coef[1]/(1 - GARCH_model@fit$coef[2] - GARCH_model@fit$coef[3]))),
      seasonal_data = seasonal_data,
      NM_model = NM_model,
      Xt = Xt,
      Yt = Yt,
      GHI = GHI,
      nday = nday,
      Ct = Ct,
      Yt_bar = Yt_bar,
      GHI_bar = GHI_bar,
      params = params,
      log_lik = sum(NM_model$loss),
      control = control,
      # Extra slot: scenarios
      scenarios = list(P = NA, Q = NA, P_Q = NA, P_Qr = NA, P_Qdw = NA, P_Qup = NA),
      # Extra slot: payoffs
      payoffs = list(hist = NA,
                     boot = NA,
                     model = list(P = NA, Q = NA, Qr = NA, Qdw = NA, Qup = NA, structured = NA),
                     sim = list(P = NA, Q = NA, Qr = NA, Qdw = NA, Qup = NA, structured = NA),
                     control = NA),
      # Extra slot: esscher parameters
      esscher = list(params = NA, model = NA, control = NA)
    ),
    class = c("solarModel", "Location", "list")
  )
}

#' Simulate trajectories
#'
#' Simulate trajectories of solar radiation with a `solarModel` object.
#'
#' @param from date, starting date for simulations.
#' @param to date, end date for simulations.
#' @param nsim integer, number of simulations.
#' @param lambda numeric, Esscher parameter. When `rf = FALSE`, the input parameter `lambda` will be transformed in negative.
#' @param vol numeric, unconditional mean of GARCH(1,1) standard deviation. If `NA` will be used the estimated one.
#' @param rf logical. When `TRUE` the AR(2) component will be set to zero.
#' @param seed scalar integer, starting random seed.
#' @param quiet logical
#'
#' @rdname solarModel_simulate
#' @name solarModel_simulate
#' @export
solarModel_simulate <- function(object, from = "2010-01-01", to = "2010-12-31", nsim = 1, lambda = 0, vol = NA, rf = FALSE, seed = 1, quiet = FALSE){

  # Number of lags to consider
  i_start <- object$control$mean.model$arOrder+1
  data <- object$data
  place <- object$place
  seasonal_data <- object$seasonal_data
  # AR(2) model (GHI)
  AR_model_Yt <- object$AR_model_Yt
  NM_model <- object$NM_model
  GARCH <- object$GARCH
  params <- object$params

  # Initial date
  from <- as.Date(from)
  # End date
  to <- as.Date(to)
  # Initialize a dataset
  max_date_from <- max(data$date)
  max_date_to <- max_date_from - i_start
  if (max_date_to >= to) {
    df_emp <- dplyr::filter(data, date >= (from - lubridate::days(i_start)) & date <= to)
    df_emp <- dplyr::bind_cols(place = place, df_emp)
  } else if (max_date_to >= from & max_date_from >= from) {
    df_emp <- dplyr::filter(data, date >= (from - lubridate::days(i_start)))
    df_new_emp <- dplyr::tibble(date = seq.Date(max(df_emp$date) + 1, to, by = "1 day"))
    df_emp <- dplyr::bind_rows(df_emp, df_new_emp)
    df_emp <- dplyr::select(df_emp, -dplyr::any_of(colnames(seasonal_data)[-c(1:2)]))
    df_emp <- dplyr::mutate(df_emp,
                            Year = lubridate::year(date),
                            Month = lubridate::month(date),
                            Day = lubridate::day(date))
    df_emp <- dplyr::left_join(dplyr::bind_cols(place = place, df_emp),
                               dplyr::select(seasonal_data, Month, Day, Ct, GHI_bar, Yt_bar, sigma_bar), by = c("Month", "Day"))
    df_emp$n <- solarr::number_of_day(df_emp$date)
  } else {
    msg <- paste0("The maximum date for starting a simulation is: ", max_date_from)
    if (!quiet) warning(msg)
    return(object)
  }
  # Garch parameters
  omega1 <- GARCH$model@fit$coef[2]
  omega2 <- GARCH$model@fit$coef[3]
  # "vol" set the level for the unconditional mean
  if (is.na(vol)) {
    vol <- GARCH$vol
    omega0 <- (vol^2)*(1 - omega1 - omega2)
  } else {
    omega0 <- GARCH$model@fit$coef[1]
  }
  # Initialize the template dataset
  df_sim_init <- dplyr::left_join(df_emp, NM_model[,c(1:6)], by = "Month")
  # Filter df_emp to be in [from - to] dates
  df_emp <- dplyr::filter(df_emp, date >= from & date <= to)
  # Initialize simulation variables
  df_sim_init$z <- rep(0, nrow(df_sim_init))
  # Initialize Market risk premium
  df_sim_init$lambda <- lambda
  df_sim_init$vol <- vol

  j <- 1
  simulations <- list()
  for(j in 1:nsim){
    # Initialize dataset for storing the simulation
    df_sim <- df_sim_init
    # Simulate Normal mixture (ut)
    df_sim$ut <- rnorm(nrow(df_sim), mean = 0, sd = 1)
    if (!quiet) message("Simulation: ", j, "/", nsim, " (", round(j/nsim*100, 4), " %) \r", appendLF = FALSE)
    set.seed(seed)
    for(i in i_start:nrow(df_sim)){
      # Bernoulli jump
      df_sim$z[i] <- rbinom(1, 1, df_sim$p_up[i])
      if (df_sim$z[i] == 1) {
        # Simulated Esscher parameter (up)
        df_sim$lambda[i] <- (df_sim$sd_up[i]^2)*df_sim$lambda[i]
        # Simulated normal mixture (up)
        df_sim$ut[i] <- df_sim$mu_up[i] + df_sim$sd_up[i]*df_sim$ut[i] # (up)
      } else {
        # Simulated Esscher parameter (down)
        df_sim$lambda[i] <- (df_sim$sd_dw[i]^2)*df_sim$lambda[i]
        # Simulated normal mixture (down)
        df_sim$ut[i] <- df_sim$mu_dw[i] + df_sim$sd_dw[i]*df_sim$ut[i] # (down)
      }
      # Simulated GARCH variance
      df_sim$sigma2[i] <- omega0 + omega1*df_sim$eps_tilde[i-1]^2  + omega2*df_sim$sigma2[i-1]
      # Simulated GARCH standard deviation
      df_sim$sigma[i] <- sqrt(df_sim$sigma2[i])
      # Simulated standardized residuals
      df_sim$eps_tilde[i] <- df_sim$sigma[i]*df_sim$ut[i]
      # Simulated AR residuals
      df_sim$eps[i] <- df_sim$eps_tilde[i]*df_sim$sigma_bar[i]
      # Simulated Esscher parameter
      df_sim$lambda[i] <- df_sim$lambda[i]*(df_sim$sigma[i]*df_sim$sigma_bar[i])^2
      # Simulated Yt_tilde
      df_sim$Yt_tilde_hat[i] <- predict(AR_model_Yt, newdata = df_sim[(i-i_start):i,])[i_start]
      # Risk-free dynamic has Yt_tilde_hat = 0
      df_sim$Yt_tilde_hat[i] <- ifelse(rf, 0, df_sim$Yt_tilde_hat[i])
      # Simulated Yt_tilde
      df_sim$Yt_tilde[i] <- df_sim$Yt_tilde_hat[i] + df_sim$eps[i]
      # Risk-free dynamic has plus lambda
      df_sim$Yt_tilde[i] <- ifelse(rf, df_sim$Yt_tilde[i] + df_sim$lambda[i], df_sim$Yt_tilde[i] - df_sim$lambda[i])
      # Simulated Yt
      df_sim$Yt[i] <- df_sim$Yt_bar[i] + df_sim$Yt_tilde[i]
      # Simulated Xt
      df_sim$Xt[i] <- params$alpha + params$beta*exp(-exp(df_sim$Yt[i]))
      # Simulated GHI
      df_sim$GHI[i] <- df_sim$Ct[i]*(1 - df_sim$Xt[i])
    }
    df_sim <- dplyr::select(df_sim, -mu_up, -mu_dw, -sd_up, -sd_dw, -p_up)
    df_sim <- dplyr::select(df_sim, -clearsky, -Yt_bar, -sigma_bar, -sigma2, -Yt_hat, -GHI_hat, -Xt_hat, -Yt_tilde_hat)
    df_sim <- dplyr::filter(df_sim, date >= from & date <= to)
    simulations[[j]] <- dplyr::bind_cols(seed = seed, df_sim)
    seed <- seed + j
  }

  df_emp <- dplyr::select(df_emp, -clearsky, -Yt_bar, -sigma_bar, -sigma2, -Yt_hat, -GHI_hat, -Xt_hat, -Yt_tilde_hat)
  structure(
    list(
      sim = simulations,
      emp = df_emp,
      params = list(seed = seed, from = from, to = to, nsim = nsim, lambda = lambda, vol = vol, rf = rf)
    ),
    class = c("solarModelSimulation", "list")
  )
}


#' Simulate multiple scenarios
#'
#' Simulate multiple scenarios of solar radiation with a `solarModel` object.
#'
#' @param from character, start Date for simulations in the format `YYYY-MM-DD`.
#' @param to character, end Date for simulations in the format `YYYY-MM-DD`.
#' @param by character, steps for multiple scenarios, e.g. `1 day` (day-ahead simulations), `15 days`, `1 month`, `3 months`, ecc.
#' For each step are simulated `nsim` scenarios.
#' @param nsim integer, number of simulations.
#' @param lambda numeric, Esscher parameter. When `rf = FALSE`, the input parameter `lambda` will be transformed in negative.
#' @param vol numeric, unconditional mean of GARCH(1,1) standard deviation. If `NA` will be used the estimated one.
#' @param rf logical. When `TRUE` the AR(2) component will be set to zero.
#' @param seed scalar integer, starting random seed.
#' @param quiet logical
#'
#' @rdname solarModel_scenario
#' @name solarModel_scenario
#' @export
solarModel_scenario <- function(object, from = "2010-01-01", to = "2010-12-31", by = "1 month", nsim = 1, lambda = 0, vol = NA, rf = FALSE, seed = 1, quiet = FALSE){

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
    sim <- solarModel_simulate(object, from = idx_date[j-1], to = idx_date[j]-1, nsim = nsim,
                               lambda = lambda, vol = vol, rf = rf, seed = seed, quiet = TRUE)
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
