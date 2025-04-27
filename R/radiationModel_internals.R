#' Create a monthly sequence of dates
#'
#' @param t_now Character, today date.
#' @param t_hor Character, horizon date.
#' @param last_day Logical. When `TRUE` the last day will be treated as conditional variance otherwise not.
#' @examples
#' t_now <- "2022-01-01"
#' t_hor <- "2022-03-24"
#' create_monthly_sequence("2022-01-01", "2022-03-24")
#' create_monthly_sequence("2022-01-01", "2022-03-24", TRUE)
#'
#' @keywords internal integral radiationModel
#' @export
create_monthly_sequence <- function(t_now, t_hor, last_day = FALSE){
  # Convert in dates
  t_now <- as.Date(t_now)
  t_hor <- as.Date(t_hor)
  # Compute year and months for reference dates
  # Standard dates
  t_start <- as.Date(paste0(lubridate::year(t_now), "-", lubridate::month(t_now), "-01"))
  t_end <- as.Date(paste0(lubridate::year(t_hor), "-", lubridate::month(t_hor), "-01"))
  # Sequence of dates
  dates_hor <- dates_now <- seq.Date(t_start, t_end, by = "1 month")
  n_months <- length(dates_now)
  dates_now <- c(t_now, dates_now[-1] - lubridate::day(dates_now[-1]))
  dates_hor <- c(dates_hor[-n_months] - lubridate::day(dates_hor[-n_months]) + lubridate::days(lubridate::days_in_month(dates_hor[-n_months])), t_hor)
  n_of_day <- as.numeric(difftime(dates_hor, dates_now, units = "days"))
  df_dates <- data.frame(Year = lubridate::year(dates_now), Month = lubridate::month(dates_hor), n = n_of_day)
  df_dates$N <- df_dates$n + number_of_day(t_now) + c(0, cumsum(lag(df_dates$n, 1)[-1]))
  df_dates$n <- df_dates$N - df_dates$n
  df_dates$tau <- as.numeric(difftime(t_hor, t_now, units = "days")) + df_dates$n[1]
  if (last_day) {
    df_last <- dplyr::mutate(tail(df_dates, 1), n = tau - 1)
    df_dates[nrow(df_dates),]$N <- df_dates[nrow(df_dates),]$N - 1
    df_dates <- dplyr::bind_rows(df_dates, df_last)
  }
  attr(df_dates, "last_day") <- last_day
  return(df_dates)
}
# Martingale estimation for mean-reverson
martingale_method_seasonal <- function(Yt, Yt_bar, e_mu = 0){
  # Dataset
  data <- dplyr::tibble(Yt = Yt, Yt_bar = Yt_bar, e_mu = e_mu)
  data$n <- 1:nrow(data)
  # Quadratic variation
  data$dYt2 <- (lag(data$Yt,1) - lag(data$Yt, 2))^2
  # Seeasonal variance from quadratic variation
  seasonal_variance <- seasonalModel$new()
  seasonal_variance$fit("dYt2 ~ 1", data = data)
  # Estimated seasonal variance
  data$sigma2_bar <- seasonal_variance$predict(1:nrow(data)-1)
  # Martingale estimation
  data$Y_est_L1 <- (lag(data$Yt_bar, 1) - lag(data$Yt, 1)) / data$sigma2_bar
  # Differences from seasonal mean
  data$dYt <-  data$Yt - data$Yt_bar
  data$dYt_L1 <- lag(data$Yt, 1) - lag(data$Yt_bar, 1) - lag(data$e_mu, 1)
  data <- na.omit(data)
  a_n <- sum(data$Y_est_L1 * data$dYt) / sum(data$Y_est_L1 * data$dYt_L1)
  -log(a_n)
}
# Reparametrization for seasonal variance
reparam_seasonal_function <- function(par, theta, omega){
  # Original parameters
  a0 <- par[1]
  a1 <- par[2]
  a2 <- par[3]
  # Correction for long term variance
  c0_long <- a0 * 2 * theta
  c1_long <- a1 * 2 * theta - omega * a2
  c2_long <- a2 * 2 * theta + omega * a1
  # Integral parameters
  gamma0_long <-  c0_long / (2 * theta)
  gamma1_long <- (c1_long * 2 * theta + c2_long * omega) / (4 * theta^2 + omega^2)
  gamma2_long <- (c2_long * 2 * theta - c1_long * omega) / (4 * theta^2 + omega^2)

  # Correction for short term variance
  alpha  <- 1 - exp(-2 * theta) * cos(omega)
  beta  <- exp(-2 * theta) * sin(omega)
  detM  <- alpha^2 + beta^2
  c0 <- (2 * theta * a0) / (1 - exp(-2 * theta))
  c1 <- ((2 * theta * alpha + omega * beta) * a1 + (2 * theta * beta - omega * alpha) * a2) / detM
  c2 <- ((omega * alpha - 2 * theta * beta) * a1 + (omega * beta + 2 * theta * alpha) * a2) / detM
  # Integral parameters
  gamma0 <-  c0 / (2 * theta)
  gamma1 <- (c1 * 2 * theta + c2 * omega) / (4 * theta^2 + omega^2)
  gamma2 <- (c2 * 2 * theta - c1 * omega) / (4 * theta^2 + omega^2)
  structure(
    list(
      alpha = alpha,
      beta = beta,
      detM = detM,
      a_ = c(a0 = a0, a1 = a1, a2 = a2),
      c_ = c(c0 = c0, c1 = c1, c2 = c2),
      c_long = c(c0 = c0_long, c1 = c1_long, c2 = c2_long),
      gamma = c(gamma0 = gamma0, gamma1 = gamma1, gamma2 = gamma2),
      gamma_long = c(gamma0 = gamma0_long, gamma1 = gamma1_long, gamma2 = gamma2_long)
    )
  )
}
# Numeric integral for the square root of the sigma_s * exp(-theta)
integral_sigma_numeric <- function(theta, par, omega = 2*base::pi/365){
  seasonal_function <- function(t) par[1] + par[2] * sin(omega * t) + par[3] * cos(omega * t)
  integrand <- function(tau, T_) sqrt(seasonal_function(tau) * exp(-2*theta*(T_-tau)))
  function(t, s, T_){
    result <- c()
    for(i in 1:length(T_)){
      result[i] <- integrate(integrand, lower = t[i], upper = s[i], T_ = T_[i])$value
    }
    return(result)
  }
}
# Formula integral for sigma_s^2 * exp(-2theta)
integral_sigma2_formula <- function(theta, par, omega = 2*base::pi/365){
  # Functions
  f0 <- function(t, s, T_) exp(-2 * theta * (T_ - s)) - exp(-2 * theta * (T_ - t))
  f1 <- function(t, s, T_) exp(-2 * theta * (T_ - s)) * sin(omega * s) - exp(-2 * theta * (T_ - t)) * sin(omega * t)
  f2 <- function(t, s, T_) exp(-2 * theta * (T_ - s)) * cos(omega * s) - exp(-2 * theta * (T_ - t)) * cos(omega * t)

  function(t, s, T_){
    result <- par[1] * f0(t, s, T_) + par[2] * f1(t, s, T_) + par[3] * f2(t, s, T_)
    return(result)
  }
}
# ===================================================================================================
#' Preprocess the moments for a specific tau to avoid recomputing them each time
#' @export
prepare_dQdP_calibration_tau <- function(tau = 10, model_Rt){
  # Radiation model
  model <- model_Rt$clone(deep = TRUE)
  # Previous filter the dataset
  data_lag <- dplyr::filter(model$model$data, n != 59.5 & weights != 0 & n != 366 & isTrain)
  # Lagged GHI value
  data_lag$L_GHI <- dplyr::lag(data_lag$GHI, tau)
  # Actual date
  data_lag$t_now <- dplyr::lag(data_lag$date, tau)
  # Horizon date
  data_lag$t_hor <- data_lag$date
  # Remove NAs
  data_lag <- na.omit(data_lag)
  data_lag <- data_lag[-nrow(data_lag),]
  # Counter for the number of steps
  nstep <- 1
  # 1) Compute the fixed parameters
  start_time <- Sys.time()
  message("Fitting M_Y (", nstep, "/3) ...", appendLF = FALSE)
  M_Y <- purrr::map(1:nrow(data_lag), ~model$M_Y(data_lag$L_GHI[.x], data_lag$t_now[.x], data_lag$t_hor[.x]))
  # Expectation of Yt under P
  data_lag$M_Y1 <- purrr::map_dbl(M_Y, ~.x[1])
  data_lag$M_Y2 <- purrr::map_dbl(M_Y, ~.x[2])
  end_time <- Sys.time()
  message("Done in ", difftime(end_time, start_time, units = "secs"), " secs \r")
  nstep <- nstep + 1
  # 2) Drifts
  start_time <- Sys.time()
  message("Fitting Mixture drift (", nstep, "/3) ...", appendLF = FALSE)
  e_mix_drift <- purrr::map(1:nrow(data_lag), ~unlist(model$e_mix_drift(data_lag$t_now[.x], data_lag$t_hor[.x])))
  data_lag$M_Y1 <- purrr::map_dbl(e_mix_drift, ~.x[1])
  data_lag$M_Y2 <- purrr::map_dbl(e_mix_drift, ~.x[2])
  data_lag$e_mu1 <- purrr::map_dbl(e_mix_drift, ~.x[3])
  data_lag$e_mu2 <- purrr::map_dbl(e_mix_drift, ~.x[4])
  data_lag$e_sd1 <- purrr::map_dbl(e_mix_drift, ~.x[5])
  data_lag$e_sd2 <- purrr::map_dbl(e_mix_drift, ~.x[6])
  end_time <- Sys.time()
  message("Done in ", difftime(end_time, start_time, units = "secs"), " secs. \r")
  nstep <- nstep + 1
  # 3) Diffusion
  start_time <- Sys.time()
  message("Fitting Mixture diffusions (", nstep, "/3) ...", appendLF = FALSE)
  e_mix_diffusion <- purrr::map(1:nrow(data_lag), ~unlist(model$e_mix_diffusion(data_lag$t_now[.x], data_lag$t_hor[.x])))
  data_lag$S_Y1 <- purrr::map_dbl(e_mix_diffusion, ~.x[1])
  data_lag$S_Y2 <- purrr::map_dbl(e_mix_diffusion, ~.x[2])
  data_lag$v_drift_P <- purrr::map_dbl(e_mix_diffusion, ~.x[3])
  data_lag$v_drift_Q1 <- purrr::map_dbl(e_mix_diffusion, ~.x[4])
  data_lag$v_drift_Q2 <- purrr::map_dbl(e_mix_diffusion, ~.x[5])
  data_lag$v_diffusion <- purrr::map_dbl(e_mix_diffusion, ~.x[6])
  data_lag$last_1 <- purrr::map_dbl(e_mix_diffusion, ~.x[7])
  data_lag$last_2 <- purrr::map_dbl(e_mix_diffusion, ~.x[8])
  end_time <- Sys.time()
  message("Done in ", difftime(end_time, start_time, units = "secs"), " secs. \r")
  nstep <- nstep + 1
  # Add time to maturity
  data_lag$tau <- tau
  # Add bounds parameters
  data_lag$alpha <- model$model$transform$alpha
  data_lag$beta <- model$model$transform$beta

  data_lag <- dplyr::select(data_lag, t_now, t_hor, tau, Rt = "L_GHI", RT = "GHI", M_Y1, M_Y2,
                            e_mu1, e_mu2, e_sd1, e_sd2, v_drift_P, v_drift_Q1, v_drift_Q2, v_diffusion, last_1, last_2,
                            S_Y1, S_Y2, C_T = "Ct", p_T = "p1", alpha, beta)
  return(data_lag)
}

#' Evaluate the loss for a specific time to maturity
#' @export
loss_dQdP_tau <- function(lambda, data, r = 0, nmonths = 1:12, quiet = FALSE){

  data <- dplyr::mutate(data, Month = lubridate::month(t_now))
  data$lambda <- lambda
  # Q-measure
  data$M_Y1_tilde <- data$M_Y1 + data$e_sd1 * data$lambda
  data$M_Y2_tilde <- data$M_Y2 + data$e_sd2 * data$lambda
  common_variance <-  data$v_drift_P +  data$v_drift_Q1 * data$lambda^2 + data$v_drift_Q2 * data$lambda + data$v_diffusion
  data$S_Y1 <- common_variance + data$last_1
  data$S_Y2 <- common_variance + data$last_2
  # Initialize the expected value at maturity
  data$e_RT <- NA
  for(i in 1:nrow(data)){
    pdf_Y = function(x) data$p_T[i] * dnorm(x,  data$M_Y1_tilde[i], sqrt(data$S_Y1[i])) + (1 - data$p_T[i]) * dnorm(x, data$M_Y2_tilde[i], sqrt(data$S_Y2[i]))
    expectation_Xt <- integrate(function(x) exp(-exp(x)) * pdf_Y(x), lower = -Inf, upper = Inf)$value
    data$e_RT[i] <- data$C_T[i] * (1 - data$alpha[i] - data$beta[i] * expectation_Xt)
  }
  # Realized nrownian motion at inception
  data$W_t <- lag(data$z, data$tau[1])
  # Realized brownian motion at maturity
  data$W_T <- data$z
  # Filtered dataset
  data <- dplyr::filter(data, Month %in% nmonths)

  # Clearness index at maturity
  KT <- data$RT / data$C_T
  # Expected value clearness index at maturity under Q
  E_KT_Q <- (data$e_RT / data$C_T) * exp(- r * data$tau)
  # Radon Nickodym derivative
  dQdP <-  exp(-lambda * (data$W_T - data$W_t) - 0.5 * lambda^2 * data$tau)
  # Expected value under Q
  E_KT <- KT * dQdP
  #r <- 0 #mean(E_KT / E_KT_Q, na.rm = TRUE) - 1

  e_loss <- (mean(E_KT / E_KT_Q, na.rm = TRUE) - 1)^2
  if (!quiet) message("Loss: ", e_loss, " Lambda: ", lambda)
  return(e_loss)
}

#' Evaluate the loss for a set of times to maturity
#' @export
loss_dQdP <- function(lambda, data_dQdP, r = 0, nmonths = 1:12, quiet = FALSE){
  loss <- 0
  for(i in 1:length(data_dQdP)){
    loss <- loss + loss_dQdP_tau(lambda, data_dQdP[[i]], r = r, nmonths = nmonths, quiet = TRUE)
  }
  if (!quiet) message("Loss: ", loss, " Lambda: ", lambda)
  return(loss)
}

#' Preprocess the moments for a specific set of time to maturities to avoid recomputing them each time
#' @export
prepare_dQdP_calibration <- function(tau = c(10, 40, 50), model_Rt){
  data_dQdP <- list()
  for(i in 1:length(tau)){
    message("Preparing data for maturity: ", tau[i])
    data_dQdP[[i]] <- prepare_dQdP_calibration_tau(tau = tau[i], model_Rt)
  }
  return(data_dQdP)
}

#' Calibration for the best lambda
#' @export
calibrate_dQdP_solar <- function(data_dQdP, r = 0, nmonths = 1:12){
  optim(par = c(0), loss_dQdP, method = "Brent", lower = -0.5, upper = 0.5, data_dQdP = data_dQdP, r = r, nmonths = nmonths)$par
}
