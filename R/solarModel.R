#' Solar Model
#'
#' @param object Location object, `CAMS("Bologna")`
#' @param control control settings, `control.solarModel()`.
#'
#' @export


solarModel <- function(object, control = control.solarModel()){

  #' @examples
  #' object <- object
  #' control <- control_model

  # Clearsky Model
  # Risk Driver: Xt is computed here
  object <- clearsky.seasonalModel(object, control = control$clearsky.model)
  # Initialize dataset
  data <- object$data
  # Solar Transform
  # Upper and lower bounds: alpha_ and beta_
  epsilon <- min(data$Xt)*control$threshold
  min_Xt <- min(data$Xt)
  max_Xt <- max(data$Xt)
  alpha_ <- min_Xt - epsilon
  beta_ <- max_Xt - min_Xt + 2*epsilon
  params = list(alpha = alpha_, beta = beta_, min_Xt = min_Xt, max_Xt = max_Xt, epsilon = epsilon)

  # Transform Functions
  Xt = function(x) alpha_ + beta_*exp(-exp(x))
  Yt = function(x) log(log(beta_) - log(x - alpha_))
  GHI = function(Ct, x) Ct*(1 - x)
  # Compute Yt
  data$Yt <- Yt(data$Xt)
  # Impute outliers
  data$Yt[data$Yt < quantile(data$Yt, probs = 0.001)] <- quantile(data$Yt, probs = 0.001)
  data$Xt <- Xt(data$Yt)
  data$GHI <- GHI(data$Ct, data$Xt)
  # Seasonal model for Yt
  seasonal_model_Yt <- seasonalModel(formula = "Yt ~ 1", order = control$mean.model$seasonalOrder, period = 365, data = data)
  # Fitted seasonal mean for Yt
  data$Yt_bar <- predict.seasonalModel(seasonal_model_Yt, n = data$n)
  # Fitted deseasonalized Yt
  data$Yt_tilde <- data$Yt - data$Yt_bar
  # Fitted seasonal mean (GHI)
  data$GHI_bar <- GHI(data$Ct, Xt(data$Yt_bar))

  # AR model
  # AR flexible Formula
  AR_formula_Yt <- "Yt_tilde ~ "
  if (control$mean.model$arOrder > 0){
    for(i in 1:control$mean.model$arOrder){
      AR_formula_Yt <- paste0(AR_formula_Yt, " + I(dplyr::lag(Yt_tilde,", i, "))")
    }
  }
  AR_formula_Yt <- paste0(AR_formula_Yt, ifelse(control$mean.model$include.intercept, "", "-1"))

  # Fitted model
  AR_model_Yt <- lm(formula = as.formula(AR_formula_Yt), data = data)
  # Fitted Yt_tilde
  data$Yt_tilde_hat <- predict(AR_model_Yt, newdata = data)
  if (control$mean.model$arOrder > 0){
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

  # Seasonal variance model
  # Fitted model
  seasonal_variance <- seasonalModel(formula = "I(eps^2) ~ 1", order = control$variance.model$seasonalOrder, period = 365, data = data)
  # Fitted seasonal standard deviation
  data$sigma_bar <- sqrt(predict.seasonalModel(seasonal_variance, n = data$n))
  # Fitted standardized residuals
  data$eps_tilde <- data$eps/data$sigma_bar

  # GARCH Model
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
  # Fitted residuals
  data$ut <- data$eps_tilde/data$sigma
  # Normal Mixture Model
  NM_model <- solarModel.monthly_mixture(data, loss = control$loss, match_moments = control$variance.model$match_moments)

  # Seasonal data by month and day for an year with 366 days
  seasonal_data <- dplyr::filter(data, date >= as.Date("2016-01-01") & date <= as.Date("2016-12-31"))
  seasonal_data$Xt_bar <- Xt(seasonal_data$Yt_bar)
  seasonal_data <- dplyr::select(seasonal_data, Month, Day, Ct, GHI_bar, Xt_bar, Yt_bar, sigma_bar)
  seasonal_data <- dplyr::left_join(seasonal_data, object$seasonal_data, by = c("Month", "Day"))

  # Update object data
  structure(
    list(
      data = dplyr::select(data, -H0),
      place = object$place,
      coords = object$coords,
      seasonal_model_Ct = object$seasonal_model_Ct,
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
      params = params,
      log_lik = ifelse(control$loss == "ml", sum(NM_model$loss), NA),
      fitted = TRUE,
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



#' Solar Normal Mixture Model
#'
#' @param data dataset with at least a column with `Month` and the target variable names `ut`.
#' @param loss character, type of loss function. Default is `ml` for maximum likelihood or can be `kl` for KL distance.
#' @param match_moments logical.
#'
#' @export
solarModel.monthly_mixture <- function(data, loss = "ml", match_moments = FALSE){
  i <- 1
  # Normal Mixture Model
  NM_model <- list()
  for(i in 1:12){
    # Monthly data
    eps <- dplyr::filter(data, Month == i)$ut
    # Initial parameters
    p0 <- 0.5
    mu_10 <- mean(eps)
    mu_20 <- -mean(eps)
    sd_10 <- sd_20 <- sd(eps)
    init_params <- c(mu1 = -mu_10, mu2 = mu_20, sd1 = sd_10, sd2 = sd_20, p = p0)
    # Fitted model
    nm <- fit_dnorm_mix(eps, params = init_params, loss = loss)

    # Compute expected value
    e_u <-  nm$par[5]*nm$par[1] + (1 - nm$par[5])*nm$par[2]
    # Compute variance
    v_u <-  nm$par[5]*(nm$par[1]^2 + nm$par[3]^2) + (1 - nm$par[5])*(nm$par[2]^2 + nm$par[4]^2)
    # Compute sample moments
    e_u_hat <- mean(eps, na.rm = TRUE)
    v_u_hat <- var(eps, na.rm = TRUE)
    # Match exactly sample moments
    if (match_moments) {
      nm$par[2] <- (e_u_hat - nm$par[5]*nm$par[1])/(1 - nm$par[5])
      # Update expected value
      e_u <-  nm$par[5]*nm$par[1] + (1 - nm$par[5])*nm$par[2]
      v_2 <- (v_u_hat + e_u^2 - (nm$par[1]^2 + nm$par[3]^2)*nm$par[5] - nm$par[2]^2 + nm$par[5]*nm$par[2]^2)/(1 - nm$par[5])
      nm$par[4] <- sqrt(v_2)
      # Update variance
      v_u <-  nm$par[5]*(nm$par[1]^2 + nm$par[3]^2) + (1 - nm$par[5])*(nm$par[2]^2 + nm$par[4]^2)
    }
    # Update expected value
    e_u <-  nm$par[5]*nm$par[1] + (1 - nm$par[5])*nm$par[2]
    # Update variance
    v_u <-  nm$par[5]*(nm$par[1]^2 + nm$par[3]^2) + (1 - nm$par[5])*(nm$par[2]^2 + nm$par[4]^2)
    # Update log-likelihood
    nm$value <- sum(dnorm_mix(nm$par)(eps, log = TRUE))
    # Fitted parameters
    df_par <- dplyr::bind_cols(dplyr::bind_rows(nm$par), p2 = 1 - nm$par[5])
    colnames(df_par) <- c("mu1", "mu2", "sd1", "sd2", "p1", "p2")
    # Monthly data
    NM_model[[i]] <- dplyr::tibble(Month = i,
                                   df_par,
                                   loss = nm$value,
                                   nobs = length(eps),
                                   mean_loss = loss/nobs,
                                   e_x = e_u,
                                   v_x = v_u,
                                   e_x_hat = e_u_hat,
                                   v_x_hat = v_u_hat)
  }

  NM_model <- dplyr::bind_rows(NM_model) %>%
    dplyr::mutate(
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
      p_dw = 1 -p_up
    ) %>%
    dplyr::select(Month, mu_up:p_dw, loss, nobs, mean_loss,
                  e_x, v_x, e_x_hat, v_x_hat)

  return(NM_model)
}



#' Simulate scenarios of solar model
#'
#' @param from scalar date, starting date for simulations.
#' @param to scalar date, end date for simulations.
#' @param nsim scalar integer, number of simulations.
#' @param lambda scalar numeric, Esscher parameter. When `rf = FALSE`, the input parameter `lambda` will be transformed in negative.
#' @param vol scalar numeric, unconditional mean of GARCH(1,1) standard deviation. If `NA` will be used the estimated one.
#' @param rf logical. When `TRUE` the AR(2) component will be set to zero.
#' @param seed scalar integer, starting random seed.
#' @param quiet logical
#'
#' @rdname simulate.solarModel
#' @name simulate.solarModel
#' @export
simulate.solarModel <- function(object, from = "2010-01-01", to = "2010-12-31", nsim = 1, lambda = 0, vol = NA, rf = FALSE, seed = 1, quiet = FALSE){
  # Initial date
  from <- as.Date(from)
  # End date
  to <- as.Date(to)
  i_start <- object$control$mean.model$arOrder+1
  # Prepare dataset
  max_date_from <- max(object$data$date)
  max_date_to <- max_date_from - i_start
  if (max_date_to >= to) {
    df_emp <- dplyr::filter(object$data, date >= (from - lubridate::days(i_start)) & date <= to)
    df_emp <- dplyr::bind_cols(place = object$place, df_emp)
  } else if (max_date_to >= from & max_date_from >= from) {
    df_emp <- dplyr::filter(object$data, date >= (from - lubridate::days(i_start)))
    df_new_emp <- dplyr::tibble(date = seq.Date(max(df_emp$date) + 1, to, by = "1 day"))
    df_emp <- dplyr::bind_rows(df_emp, df_new_emp)
    df_emp <- dplyr::select(df_emp, -dplyr::any_of(colnames(object$seasonal_data)[-c(1:2)]))
    df_emp <- dplyr::mutate(df_emp,
                            Year = lubridate::year(date),
                            Month = lubridate::month(date),
                            Day = lubridate::day(date))
    df_emp <- dplyr::left_join(dplyr::bind_cols(place = object$place, df_emp),
                               dplyr::select(object$seasonal_data, Month, Day, Ct, GHI_bar, Yt_bar, sigma_bar), by = c("Month", "Day"))
    df_emp$n <- solarr::number_of_day(df_emp$date)
  } else {
    msg <- paste0("The maximum date for starting a simulation is: ", max_date_from)
    if (!quiet) warning(msg)
    return(object)
  }
  # Garch parameters
  omega1 <- object$GARCH$model@fit$coef[2]
  omega2 <- object$GARCH$model@fit$coef[3]
  # "vol" set the level for the unconditional mean
  if (is.na(vol)) {
    vol <- object$GARCH$vol
    omega0 <- (vol^2)*(1 - omega1 - omega2)
  } else {
    omega0 <- object$GARCH$model@fit$coef[1]
  }
  # AR(2) model (GHI)
  AR_model_Yt <- object$AR_model_Yt
  # Initialize the template dataset
  df_sim_init <- dplyr::left_join(df_emp, object$NM_model[,c(1:6)], by = "Month")
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
      df_sim$Xt[i] <- object$Xt(df_sim$Yt[i])
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



#' @rdname simulate.solarModel
#'
#' @export
scenario.solarModel <- function(object, from = "2010-01-01", to = "2010-12-31", by = "1 month", nsim = 1, lambda = 0, vol = NA, rf = FALSE, seed = 1, quiet = FALSE){

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
    sim <- simulate.solarModel(object, from = idx_date[j-1], to = idx_date[j]-1, nsim = nsim,
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



#' Convert a vector of parameter in a structured list
#'
#' @export
from_list_to_parameters.solarModel <- function(params_list){

  params <- update.solarModel(model)
  idx_params <- attributes(params)
  params$Yt <- params_list[idx_params$Yt] # a0, a1, a2
  params$AR <- params_list[idx_params$AR] # alpha1, alpha2
  params$sigma_bar <- params_list[idx_params$sigma_bar] # c0, c1, c2
  params$GARCH <- params_list[idx_params$GARCH] # omega0, omega1, omega2
  params$GARCH[1] <- 1 - params$GARCH[2] - params$GARCH[3]
  params$NM_mu_up <- params_list[idx_params$NM_mu_up] # mu_up(1), ..., mu_up(12)
  params$NM_mu_dw <- params_list[idx_params$NM_mu_dw] # mu_dw(1), ..., mu_dw(12)
  params$NM_sd_up <- params_list[idx_params$NM_sd_up] # sd_up(1), ..., sd_up(12)
  params$NM_sd_dw <- params_list[idx_params$NM_sd_dw] # sd_dw(1), ..., sd_dw(12)
  params$p_up <- params_list[idx_params$NM_p_up]      #  p_up(1), ..., p_up(12)

  return(params)
}


#' Extract and update parameters for Solar Model
#'
#' @export
update.solarModel <- function(object, params = NULL){

  control <- object$control

  update.lm <- function(object, params){
    if (missing(params) | is.null(params)){
      return(object)
    }
    coef_names <- names(object$coefficients)
    object$coefficients <- params
    names(object$coefficients) <- coef_names
    object$fitted.values <- predict.lm(object, new_data = object$model)
    object$residuals <- object$model[,1] -  object$fitted.values
    object
  }

  update.GARCH <- function(object, params){
    if (missing(params) | is.null(params)){
      return(object)
    }
    coef_names <- names(object@fit$coef)
    object@fit$coef <- params
    names(object@fit$coef) <- coef_names
    object
  }

  # Update parameters
  if (!is.null(params)) {
    if (!is.list(params)){
      params <- from_list_to_parameters.solarModel(params)
    }
    object$seasonal_model_Yt <- update.lm(object$seasonal_model_Yt, params = params$Yt)
    object$AR_model_Yt <- update.lm(object$AR_model_Yt, params = params$AR)
    object$seasonal_variance <- update.lm(object$seasonal_variance, params = params$sigma_bar)
    object$GARCH$model <- update.GARCH(object$GARCH$model, params$GARCH)
    object$NM_model$mu_up <- params$NM_mu_up
    object$NM_model$mu_dw <- params$NM_mu_dw
    object$NM_model$mu_up <- params$NM_sd_up
    object$NM_model$sd_dw <- params$NM_sd_dw
    object$NM_model$p_up <- params$NM_p_up
    object$NM_model$p_dw  <- 1 - object$NM_model$p_up
    object$fitted <- FALSE
    return(object)
  }

  # Extract parameters
  # Seasonal Model Yt
  n_seasonal_model_Yt <- control$mean.model$seasonalOrder*2 + 1
  idx_seasonal_model_Yt <- 0:(n_seasonal_model_Yt-1)
  seasonal_model_Yt <- object$seasonal_model_Yt$coefficients
  names(seasonal_model_Yt) <- paste0("a", idx_seasonal_model_Yt)

  # AR model
  n_AR_model_Yt <- control$mean.model$arOrder + ifelse(control$mean.model$include.intercept, 1, 0)
  idx_AR_model_Yt <- ifelse(control$mean.model$include.intercept, list(0:(n_AR_model_Yt-1)), list(1:n_AR_model_Yt))[[1]]
  AR_model_Yt <- object$AR_model_Yt$coefficients
  names(AR_model_Yt) <- paste0("alpha", idx_AR_model_Yt)

  # Unconditional Seasonal variance
  n_seasonal_variance <- control$variance.model$seasonalOrder*2 + 1
  idx_seasonal_variance <- 0:(n_seasonal_variance-1)
  seasonal_variance <- object$seasonal_variance$coefficients
  names(seasonal_variance) <- paste0("c", idx_seasonal_variance)

  # Conditional GARCH variance
  GARCH_model <- object$GARCH$model@fit$coef
  # Fix unconditional mean equal to 1
  GARCH_model[1] <- (1- GARCH_model[2] - GARCH_model[3])
  n_GARCH_model <- length(GARCH_model)
  names(GARCH_model) <- paste0("omega", 0:2)

  # Determine the positions
  idx_seasonal_model_Yt <- 1:n_seasonal_model_Yt
  idx_AR_model_Yt <- (max(idx_seasonal_model_Yt)+1):(max(idx_seasonal_model_Yt)+n_AR_model_Yt)
  idx_seasonal_variance <- (max(idx_AR_model_Yt)+1):(max(idx_AR_model_Yt)+n_seasonal_variance)
  idx_GARCH_model <- (max(idx_seasonal_variance)+1):(max(idx_seasonal_variance)+n_GARCH_model)
  idx_mu_up <- (max(idx_GARCH_model)+1):(max(idx_GARCH_model)+12)
  idx_mu_dw <- (max(idx_mu_up)+1):(max(idx_mu_up)+12)
  idx_sd_up <- (max(idx_mu_dw)+1):(max(idx_mu_dw)+12)
  idx_sd_dw <- (max(idx_sd_up)+1):(max(idx_sd_up)+12)
  idx_p_up <- (max(idx_sd_dw)+1):(max(idx_sd_dw)+12)

  structure(
    list(
      Yt = seasonal_model_Yt,
      AR = AR_model_Yt,
      sigma_bar = seasonal_variance,
      GARCH = GARCH_model,
      NM_mu_up = model$NM_model$mu_up,
      NM_mu_dw = model$NM_model$mu_dw,
      NM_sd_up = model$NM_model$sd_up,
      NM_sd_dw = model$NM_model$sd_dw,
      NM_p_up = model$NM_model$p_up
    ),
    Yt = idx_seasonal_model_Yt,
    AR = idx_AR_model_Yt,
    sigma_bar = idx_seasonal_variance,
    GARCH = idx_GARCH_model,
    NM_mu_up = idx_mu_up,
    NM_mu_dw = idx_mu_dw,
    NM_sd_up = idx_sd_up,
    NM_sd_dw = idx_sd_dw,
    NM_p_up = idx_p_up
  )
}


#' Compute the Log-likelihood for a Solar Model
#'
#' @export
logLik.solarModel <- function(model, params){

  if (missing(params)){
    params <- update.solarModel(model)
  }

  model <- update.solarModel(model, params)
  control <- model$control

  data <- model$data
  params <- unlist(params)
  idx_params <- attributes(update.solarModel(model))
  # Parameters
  seasonal_mean <- params[idx_params$Yt] # a0, a1, a2
  conditional_mean <- params[idx_params$AR] # alpha1, alpha2
  seasonal_variance <- params[idx_params$sigma_bar] # c0, c1, c2
  conditional_variance <- params[idx_params$GARCH] # omega0, omega1, omega2
  conditional_variance[1] <- 1 - conditional_variance[2] - conditional_variance[3]
  mu_up <- params[idx_params$NM_mu_up] # mu_up(1), ..., mu_up(12)
  mu_dw <- params[idx_params$NM_mu_dw] # mu_dw(1), ..., mu_dw(12)
  sd_up <- params[idx_params$NM_sd_up] # sd_up(1), ..., sd_up(12)
  sd_dw <- params[idx_params$NM_sd_dw] # sd_dw(1), ..., sd_dw(12)
  p_up <- params[idx_params$NM_p_up]  #  p_up(1), ..., p_up(12)

  if (any(p_up < 0.01 | p_up > 0.99)){
    print("p_up")
    return(NA)
  }
  if (any(conditional_variance < 0) | any(conditional_variance > 1)){
    print("conditional_variance")
    return(NA)
  }
  if (any(conditional_mean < 0) | sum(conditional_mean) > 1){
    print("conditional_mean")
    return(NA)
  }

  # Update Seasonal mean Yt
  model$data$Yt_bar <- predict.seasonalModel(model$seasonal_model_Yt, n = data$n)
  model$seasonal_data$Yt_bar <- predict.seasonalModel(model$seasonal_model_Yt, n = 1:366)
  # Update Seasonal mean Xt
  data$Xt_bar <- model$Xt(data$Yt_bar)
  model$seasonal_data$Xt_bar <- model$Xt(model$seasonal_data$Yt_bar)
  # Update Seasonal mean GHI
  data$GHI_bar <- model$GHI(data$Ct, data$Xt_bar)
  model$seasonal_data$GHI_bar <- model$GHI(model$seasonal_data$Ct, model$seasonal_data$Xt_bar)
  # Update seasonal variance
  model$seasonal_variance$coefficients <- seasonal_variance
  data$sigma_bar <- sqrt(predict.seasonalModel(model$seasonal_variance, n = data$n))
  model$seasonal_data$sigma_bar <- sqrt(predict.seasonalModel(model$seasonal_variance, n = 1:366))
  # Update AR model
  data$Yt_tilde <- data$Yt - data$Yt_bar
  model$AR_model_Yt$coefficients <- conditional_mean
  data$Yt_tilde_hat <- predict(model$AR_model_Yt, newdata = data)
  data$Yt_tilde_hat[1:(control$mean.model$arOrder)] <- data$Yt_tilde[1:(control$mean.model$arOrder)]
  # Update Fitted values Yt, Xt, GHI
  data$Yt_hat <- data$Yt_bar + data$Yt_tilde_hat
  data$Xt_hat <- model$Xt(data$Yt_hat)
  data$GHI_hat <- model$GHI(data$Ct, data$Xt_hat)
  # Update Residuals AR model
  data$eps <- data$Yt - data$Yt_hat
  # Update Standardized residuals
  data$eps_tilde <- data$eps/data$sigma_bar
  # Initialize log-likelihood
  model$NM_model$loss <- 0
  model$NM_model$type <- "ml"

  i_start <- control$mean.model$arOrder+1
  i <- i_start
  for(i in i_start:nrow(data)){
    # Garch variance
    data$sigma2[i] <- conditional_variance[1] + conditional_variance[2]*data$eps_tilde[i-1]^2 + conditional_variance[3]*data$sigma2[i-1]
    data$sigma[i] <- sqrt(data$sigma2[i])
    # Garch residuals
    data$ut[i] <- data$eps_tilde[i]/data$sigma[i]
    # Normal mixture parameters
    nmonth <- data$Month[i]
    df_nm <- model$NM_model[nmonth,]
    par <- c(mu_up =  df_nm$mu_up, mu_dw = df_nm$mu_dw, sd_up =  df_nm$sd_up, sd_dw = df_nm$sd_dw, p = df_nm$p_up)
    # Log-likelihood for each month
    model$NM_model$loss[nmonth] <- model$NM_model$loss[nmonth] + dnorm_mix(params = par)(data$ut[i], log = TRUE)
    model$NM_model$mean_loss[nmonth] <- model$NM_model$loss[nmonth]/model$NM_model$nobs[nmonth]
  }

  model$data <- data
  model$log_lik <- sum(model$NM_model$loss)
  if(!control$quiet) message(" Log-likelihood: ", model$log_lik, "\r", appendLF = FALSE)
  return(model)
}


#' Log-likelihood optimazation for Solar Model
#'
#' @export
optimize.solarModel <- function(model, maxit = 100, quiet = FALSE){

  init_params <- update.solarModel(model)
  init_params <- unlist(init_params)

  logLik_function <- function(params){
    log_lik <- logLik.solarModel(model, params)$log_lik
    message("Loss: ", log_lik)
    return(log_lik)
  }
  # Optimal parameters
  opt <- optim(par = init_params, fn = logLik_function, control = list(maxit = maxit))
  opt_model <- update.solarModel(model, params = opt$par)
  opt_model
}


