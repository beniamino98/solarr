#' Radiation model
#'
#' @rdname radiationModel
#' @name radiationModel
#' @note Version 1.0.0
#' @export
radiationModel <- R6::R6Class("radiationModel",
                              # ====================================================================================================== #
                              #                                             Public slots
                              # ====================================================================================================== #
                              public = list(
                                #' @field theta Numeric, mean reversion parameter.
                                theta = NA,
                                #' @field lambda_S Numeric, market risk premium Q-measure.
                                lambda_S = 0,
                                #' @description
                                #' Initialize a `radiationModel` object
                                #' @param model solarModel model fitted
                                initialize = function(model, correction = FALSE){
                                  # Store the discrete model
                                  private$..model <- model$clone(TRUE)
                                  # 1) Estimate mean reversion parameter
                                  df_fit <- filter(private$..model$data, isTrain & weights != 0)
                                  self$theta <-  martingale_method_seasonal(df_fit$Yt, df_fit$Yt_bar)
                                  # Convert AR parameter into mean reversion parameter (alternative)
                                  # self$theta <- -log(self$model$ARMA$coefficients[2])
                                  #
                                  # 2) Reparametrize seasonal function with continuous time parameters
                                  reparam <- reparam_seasonal_function(self$model$seasonal_variance$coefficients, self$theta, omega = 2*base::pi/365)
                                  # Clone seasonal variance to store c parameters
                                  private$seasonal_variance <- self$model$seasonal_variance$clone(TRUE)
                                  private$seasonal_variance$extra_params$reparam <- reparam
                                  names(reparam$c_) <- names(self$model$seasonal_variance$coefficients)
                                  private$seasonal_variance$update(reparam$c_)
                                  private$..integral_variance <- integral_sigma2_formula(self$theta, reparam$gamma, omega = 2*base::pi/365)
                                  private$..integral_expectation <- integral_sigma_numeric(self$theta, reparam$c_, omega = 2*base::pi/365)

                                  # Adjustment for the mean that multiplies I
                                  t <- 1:365
                                  sigma_bar_J <- sqrt(private$..integral_variance(t-1, t, t))
                                  sigma_bar_I <- private$..integral_expectation(t-1, t, t)
                                  private$..k1 <- mean(sigma_bar_I/sigma_bar_J)
                                  # Original variance
                                  v <- self$model$NM_model$moments$variance
                                  # Original parameters
                                  sigma <- self$model$NM_model$sd
                                  probs <- self$model$NM_model$p
                                  # Adjusted means
                                  means <- self$model$NM_model$means * private$..k1
                                  # New expectation
                                  m1 <- means[,1] * probs[,1] + means[,2] * probs[,2]
                                  k2 <- (v + m1 - (means[,1] - means[,2])^2 * probs[,1] * probs[,2]) / (sigma[,1]^2 * probs[,1] + sigma[,2]^2 * probs[,2])
                                  # Adjusted std. dev
                                  private$..k2 <- sqrt(k2)
                                  sigma <- sigma * private$..k2
                                  # Update mixture parameters
                                  if (correction){
                                    private$..model$NM_model$update(means = means)#, sd = sigma)
                                  }
                                  #self$model$filter()
                                },
                                #' @description
                                #' Change the reference probability measure
                                #' @param measure Character, probability measure. Can be `P` or `Q`.
                                change_measure = function(measure){
                                  measure <- match.arg(measure, choices = c("P", "Q"))
                                  private$..measure <- measure
                                  if (measure == "Q"){
                                    private$..lambda <- self$lambda_S
                                  } else {
                                    private$..lambda <- 0
                                  }
                                },
                                #' @description
                                #' Clear sky radiation for a day of the year.
                                #' @param t_now Character, today date.
                                #' @return Clear sky radiation on date t_now.
                                Ct = function(t_now){
                                  self$model$seasonal_model_Ct$predict(number_of_day(t_now))
                                },
                                #' @description
                                #' Seasonal mean for the transformed variable Yt for a given day of the year.
                                #' @param t_now Character, today date.
                                #' @return Seasonal mean for Yt on date t_now.
                                Yt_bar = function(t_now){
                                  self$model$seasonal_model_Yt$predict(number_of_day(t_now))

                                },
                                #' @description
                                #' Compute the seasonal mean for the solar radiation for a given day of the year.
                                #' @param t_now Character, today date.
                                #' @return Seasonal mean for Rt.
                                Rt_bar = function(t_now){
                                  self$model$Y_to_R(self$Yt_bar(t_now), t_now)
                                },
                                #' @description
                                #' Instantaneous seasonal variance for the transformed variable for a given day of the year.
                                #' @param t_now Character, today date.
                                #' @return Seasonal std. deviation for Yt on date t_now.
                                sigma_bar = function(t_now){
                                  sqrt(private$seasonal_variance$predict(number_of_day(t_now)))
                                },
                                #' @description
                                #' Return the mixture drift if B is specified, otherwise it return the average drift.
                                #' @param t_now Character, today date.
                                #' @param B Integer, 1 for the first component, 0 for the second.
                                #' @return Mixture seasonal drift for Yt on date t_now.
                                mu_B = function(t_now, B = 1){
                                  result <- ifelse(B == 1, self$model$NM_model$mu1$predict(t_now), self$model$NM_model$mu2$predict(t_now))
                                  return(result)
                                },
                                #' @description
                                #' Return the mixture diffusion with seasonal jump.
                                #' @param t_now Character, today date.
                                #' @param B Integer, 1 for the first component, 0 for the second.
                                #' @return Mixture seasonal diffusion for Yt.
                                sigma_B = function(t_now, B){
                                  result <- ifelse(B == 1, self$model$NM_model$sd1$predict(t_now), self$model$NM_model$sd2$predict(t_now))
                                  return(result)
                                },
                                #' @description
                                #' Return the drift for the transformed variable Yt.
                                #' @param Yt Numeric, transformed solar radiation.
                                #' @param t_now Character, today date.
                                #' @param B Integer, 1 for the first component, 0 for the second.
                                #' @return Mixture drift for Yt.
                                mu_Y = function(Yt, t_now, B = 1){
                                  # Differential for seasonal function
                                  dYt_bar <- self$model$seasonal_model_Yt$differential(t_now)
                                  # Seasonal variance
                                  sigma_bar <- self$sigma_bar(t_now)
                                  # Drift for Yt
                                  drift_Y <- dYt_bar - self$theta * (Yt - self$Yt_bar(t_now)) + sigma_bar * self$mu_B(t_now, B)
                                  # !!! Change not already implemented for QT measure
                                  drift_Y + sigma_bar * self$sigma_B(t_now, B) * self$lambda * ifelse(self$measure == "Q", 1, 0)
                                },
                                #' @description
                                #' Return the diffusion for solar radiation process
                                #' @param Rt Numeric, solar radiation.
                                #' @param t_now Character, today date.
                                #' @param B Integer, 1 for the first component, 0 for the second.
                                #' @return Diffusion for Rt.
                                sigma_Y = function(t_now, B = 1){
                                  # Seasonal variance
                                  sigma_bar <- self$sigma_bar(t_now)
                                  # Diffusion for Rt process
                                  sigma_bar * self$sigma_B(t_now, B)
                                },
                                #' @description
                                #' Return the drift for solar radiation process
                                #' @param Rt Numeric, solar radiation.
                                #' @param t_now Character, today date.
                                #' @param B Integer, 1 for the first component, 0 for the second.
                                #' @param dt Numeric, time step.
                                #' @return Drift for Rt.
                                mu_R = function(Rt, t_now, B = 1, dt = 1){
                                  # Convert Rt to Yt
                                  Yt <- self$model$R_to_Y(Rt, t_now)
                                  # Approximate the drift for Ct
                                  n <- number_of_day(t_now)
                                  Ct <- self$Ct(n)
                                  dCt <- (self$Ct(n + dt) - Ct)
                                  # Clearness index
                                  Kt <- 1 - self$model$transform$alpha - self$model$transform$beta * exp(-exp(Yt))
                                  Kt * dCt / dt + Ct * self$model$transform$beta * exp(Yt - exp(Yt)) * (self$mu_Y(Yt, t_now, B) + 0.5 * (1 - exp(Yt)) * self$sigma_Y(t_now, B)^2)
                                },
                                #' @description
                                #' Return the diffusion for solar radiation process
                                #' @param Rt Numeric, solar radiation.
                                #' @param t_now Character, today date.
                                #' @param B Integer, 1 for the first component, 0 for the second.
                                #' @return Diffusion for Rt.
                                sigma_R = function(Rt, t_now, B = 1){
                                  # Convert Rt to Yt
                                  Yt <- self$model$R_to_Y(Rt, t_now)
                                  # Diffusion for Rt process
                                  self$Ct(t_now) * self$model$transform$beta * exp(Yt - exp(Yt)) * self$sigma_Y(t_now, B)
                                },
                                #' @description
                                #' Compute the integral for expectation for constant mixture parameters
                                #' @param t_now Character, today date.
                                #' @param t_hor Character, horizon date.
                                #' @param df_date Optional dataframe. See \code{\link{create_monthly_sequence}} for more details.
                                #' @param last_day Logical. When `TRUE` the last day will be treated as conditional variance otherwise not.
                                integral_expectation = function(t_now, t_hor, df_date, last_day = TRUE){
                                  # Create a sequence of dates for constant monthly parameters
                                  if (missing(df_date)) {
                                    df_date <- create_monthly_sequence(t_now, t_hor, last_day = last_day)
                                  }
                                  # Compute the integral
                                  df_date$int_sigma <- private$..integral_expectation(df_date$n, df_date$N, df_date$tau)
                                  return(df_date)
                                },
                                #' @description
                                #' Compute the integral for variance for constant mixture parameters
                                #' @param t_now Character, today date.
                                #' @param t_hor Character, horizon date.
                                #' @param df_date Optional dataframe. See \code{\link{create_monthly_sequence}} for more details.
                                #' @param last_day Logical. When `TRUE` the last day will be treated as conditional variance otherwise not.
                                integral_variance = function(t_now, t_hor, df_date, last_day = TRUE){
                                  # Create a sequence of dates for constant monthly parameters
                                  if (missing(df_date)) {
                                    df_date <- create_monthly_sequence(t_now, t_hor, last_day = last_day)
                                  }
                                  # Compute the integral
                                  df_date$int_sigma2 <- private$..integral_variance(df_date$n, df_date$N, df_date$tau)
                                  return(df_date)
                                },
                                #' @description
                                #' Return the value of the mixture drift of both component of Yt.
                                #' @param t_now Character, today date.
                                #' @param t_hor Character, horizon date.
                                #' @param df_date Optional dataframe. See \code{\link{create_monthly_sequence}} for more details.
                                #' @return Mixture expected value for both component of Yt.
                                e_mix_drift = function(t_now, t_hor, df_date){
                                  if (missing(df_date)){
                                    # Create a sequence of dates
                                    df_date <- create_monthly_sequence(t_now, t_hor, last_day = TRUE)
                                  }
                                  # Adjust parameters for Q-measure
                                  df_params <- self$model$NM_model$coefficients
                                  # Expected drift
                                  df_params$e_mu_B <- df_params$mu1 * df_params$p1 + df_params$mu2 * df_params$p2
                                  # Exprected diffusion drift
                                  df_params$e_sigma_B <- df_params$sd1 * df_params$p1 + df_params$sd2 * df_params$p2
                                  # Combine the datasets
                                  df <- merge(df_date, df_params, by = "Month", all.x = TRUE)
                                  # !! Maybe can be removed !! tocheck
                                  df <- df[order(df$n),]
                                  # Compute the integral
                                  df$int <- private$..integral_expectation(df$n, df$N, df$tau)
                                  # Index for last day
                                  nrows <- nrow(df)
                                  # Total drift from time t up to T-1 for mu
                                  mix_drift_mu <- df[-nrows,]$int * df[-nrows,]$e_mu_B
                                  # Total drift from time T-1 up to T for mu
                                  mix_drift_mu_1 <- sum(mix_drift_mu) + df[nrows,]$int * df[nrows,]$mu1
                                  mix_drift_mu_2 <- sum(mix_drift_mu) + df[nrows,]$int * df[nrows,]$mu2
                                  # Total drift from time t up to T-1 for sd
                                  mix_drift_sd <- df[-nrows,]$int * df[-nrows,]$e_sigma_B
                                  mix_drift_sd_1 <- sum(mix_drift_sd) + df[nrows,]$int * df[nrows,]$sd1
                                  mix_drift_sd_2 <- sum(mix_drift_sd) + df[nrows,]$int * df[nrows,]$sd2
                                  # Realized drift depending on the measure
                                  mix_drift_1 <- mix_drift_mu_1 + mix_drift_sd_1 * self$lambda
                                  mix_drift_2 <- mix_drift_mu_2 + mix_drift_sd_2 * self$lambda
                                  list(e_drift_1 = mix_drift_1, e_drift_2 = mix_drift_2,
                                       e_mu1 = mix_drift_mu_1, e_mu2 = mix_drift_mu_2,
                                       e_sd1 = mix_drift_sd_1, e_sd2 = mix_drift_sd_2)
                                },
                                #' @description
                                #' Return the value of the mixture drift of both component of Yt.
                                #' @param t_now Character, today date.
                                #' @param t_hor Character, horizon date.
                                #' @param df_date Optional dataframe. See \code{\link{create_monthly_sequence}} for more details.
                                #' @return Mixture expected value for both component of Yt.
                                e_mix_diffusion = function(t_now, t_hor, df_date){
                                  if (missing(df_date)){
                                    # Create a sequence of dates
                                    df_date <- create_monthly_sequence(t_now, t_hor, last_day = TRUE)
                                  }
                                  # Adjust parameters for Q-measure
                                  df_params <- self$model$NM_model$coefficients %>%
                                    dplyr::mutate(mu_diff = mu1 - mu2,
                                                  sigma_diff = sd1 - sd2,
                                                  sigma2_1 = sd1^2,
                                                  sigma2_2 = sd2^2)
                                  # Combine the datasets
                                  df <- merge(df_date, df_params, by = "Month", all.x = TRUE)
                                  df <- df[order(df$n),]
                                  # Compute the integral
                                  df$int_sigma2 <- private$..integral_variance(df$n, df$N, df$tau)
                                  # Index for last day
                                  nrows <- nrow(df)
                                  df_t <- df[-nrows,]
                                  df_T <- df[nrows,]
                                  # Drift variance P in t, T-1
                                  variance_drift_P <- sum(df_t$int_sigma2 * df_t$mu_diff^2 * df_t$p1 * df_t$p2)
                                  # Drift variance Q in t, T-1
                                  variance_drift_Q_1 <- sum(df_t$int_sigma2 * df_t$sigma_diff^2 * df_t$p1 * df_t$p2)
                                  variance_drift_Q_2 <- sum(df_t$int_sigma2 * df_t$sigma_diff * df_t$mu_diff * df_t$p1 * df_t$p2)
                                  # Diffusion variance
                                  variance_diffusion <- sum(df_t$int_sigma2 * (df_t$sigma2_1 *df_t$p1 + df_t$sigma2_2 * df_t$p2))
                                  # Total variance
                                  common_variance <- variance_drift_P + variance_drift_Q_1 * self$lambda^2 + variance_drift_Q_2 * self$lambda + variance_diffusion
                                  # Realized variance for each component between T-1 and T
                                  variance_1 <- common_variance + df_T$int_sigma2 * df_T$sigma2_1
                                  variance_2 <- common_variance + df_T$int_sigma2 * df_T$sigma2_2

                                  list(variance_1 = variance_1, variance_2 = variance_2,
                                       variance_drift_P = variance_drift_P, variance_drift_Q_1 = variance_drift_Q_1, variance_drift_Q_2 = variance_drift_Q_2,
                                       common_variance = common_variance + df_T$int_sigma2, variance_diffusion = variance_diffusion,
                                       last_1 = df_T$int_sigma2 * df_T$sigma2_1, last_2 = df_T$int_sigma2 * df_T$sigma2_2)
                                },
                                #' @description
                                #' Return the conditional expectation for Yn for YN.
                                #' @param Rt Numeric, solar radiation.
                                #' @param t_now Character, today date.
                                #' @param t_hor Character, horizon date.
                                #' @param df_date Optional dataframe. See \code{\link{create_monthly_sequence}} for more details.
                                #' @return Conditional mean for Yt
                                M_Y = function(Rt, t_now, t_hor, df_date){
                                  # Create once the sequence of dates
                                  if (missing(df_date)) {
                                    df_date <- create_monthly_sequence(t_now, t_hor, last_day = TRUE)
                                  }

                                  # Convert Rt to Yt
                                  Yt <- self$model$R_to_Y(Rt, t_now)
                                  # Time to maturity in days
                                  tau <- max(df_date$tau) - min(df_date$n)
                                  # Compute the drifts
                                  mix_drift <- self$e_mix_drift(df_date = df_date)
                                  # Forecast for Yt
                                  Y_forecast <- self$Yt_bar(t_hor) + (Yt - self$Yt_bar(t_now)) * exp(-self$theta * tau)
                                  # Expected value for each component
                                  c(up = Y_forecast + mix_drift[[1]], dw = Y_forecast + mix_drift[[2]], e_Yt = Y_forecast)
                                },
                                #' @description
                                #' Return the conditional variance for Yn for YN.
                                #' @param Rt Numeric, solar radiation.
                                #' @param t_now Character, today date.
                                #' @param t_hor Character, horizon date.
                                #' @param df_date Optional dataframe. See \code{\link{create_monthly_sequence}} for more details.
                                #' @return Conditional variance for Yt
                                S_Y = function(t_now, t_hor, df_date){
                                  # Create once the sequence of dates
                                  if (missing(df_date)) {
                                    df_date <- create_monthly_sequence(t_now, t_hor, last_day = TRUE)
                                  }
                                  # Compute the diffusions
                                  mix_diffusion <- self$e_mix_diffusion(df_date = df_date)
                                  # Compute the variances
                                  c(up = mix_diffusion[[1]], dw = mix_diffusion[[2]])
                                },
                                #' @description
                                #' Return the conditional density for Y_N given Yn.
                                #' @param Rt Numeric, solar radiation.
                                #' @param t_now Character, today date.
                                #' @param t_hor Character, horizon date.
                                #' @param B Integer, mixture component, if B is missing will be returned the mixture density, otherwise the component density non weighted.
                                #' @return Conditional density for Y_N
                                pdf_Y = function(Rt, t_now, t_hor, B){
                                  # Create once the sequence of dates
                                  df_date <- create_monthly_sequence(t_now, t_hor, last_day = TRUE)
                                  # Expected values Y at horizon
                                  M_Y <- unlist(self$M_Y(Rt, t_now, t_hor, df_date = df_date))
                                  # Variance of Y at horizon
                                  S_Y <- unlist(self$S_Y(t_now, t_hor, df_date = df_date))
                                  # Mixture probability at horizon time
                                  p <- self$model$NM_model$prob$predict(t_hor)
                                  # Densities
                                  if (missing(B)) {
                                    function(x){p * dnorm(x, M_Y[1], sqrt(S_Y[1])) + (1 - p) * dnorm(x, M_Y[2], sqrt(S_Y[2]))}
                                  } else if (B == 1) {
                                    function(x) {dnorm(x, M_Y[1], sqrt(S_Y[1]))}
                                  } else {
                                    function(x) {dnorm(x, M_Y[2], sqrt(S_Y[2]))}
                                  }
                                },
                                #' @description
                                #' Return the conditional distribution for Y_N given Yn.
                                #' @param Rt Numeric, solar radiation.
                                #' @param t_now Character, today date.
                                #' @param t_hor Character, horizon date.
                                #' @param B Integer, mixture component, if B is missing will be returned the mixture distribution, otherwise the component distribution non weighted.
                                #' @return Conditional distribution for Y_N
                                cdf_Y = function(Rt, t_now, t_hor, B){
                                  # Create once the sequence of dates
                                  df_date <- create_monthly_sequence(t_now, t_hor, last_day = TRUE)
                                  # Expected values Y at horizon
                                  M_Y <- unlist(self$M_Y(Rt, t_now, t_hor, df_date = df_date))
                                  # Variance of Y at horizon
                                  S_Y <- unlist(self$S_Y(t_now, t_hor, df_date = df_date))
                                  # Mixture probability at horizon time
                                  p <- self$model$NM_model$prob$predict(t_hor)
                                  if (missing(B)) {
                                    function(x){p * pnorm(x, M_Y[1], sqrt(S_Y[1])) + (1 - p) * pnorm(x, M_Y[2], sqrt(S_Y[2]))}
                                  } else if (B == 1) {
                                    function(x){pnorm(x, M_Y[1], sqrt(S_Y[1]))}
                                  } else {
                                    function(x){pnorm(x, M_Y[2], sqrt(S_Y[2]))}
                                  }
                                },
                                #' @description
                                #' Return the conditional density for R_N given Rn.
                                #' @param Rt Numeric, solar radiation.
                                #' @param t_now Character, today date.
                                #' @param t_hor Character, horizon date.
                                #' @param B Integer, mixture component, if B is missing will be returned the mixture density, otherwise the component density non weighted.
                                #' @return Conditional density for R_N
                                pdf_R = function(Rt, t_now, t_hor, B){
                                  C_T <- self$Ct(t_hor)
                                  pdf_Y <- self$pdf_Y(Rt, t_now, t_hor, B)
                                  function(x){
                                    z_x <- (1 - x/C_T - self$model$transform$alpha) / self$model$transform$beta
                                    u_x <- suppressWarnings(log(-log(z_x)))
                                    num <- pdf_Y(u_x)
                                    den <- C_T * self$model$transform$beta * log(z_x^z_x)
                                    probs <- -num/den
                                    probs[is.nan(probs)] <- 0
                                    return(probs)
                                  }
                                },
                                #' @description
                                #' Return the conditional distribution for R_N given Rn.
                                #' @param Rt Numeric, solar radiation.
                                #' @param t_now Character, today date.
                                #' @param t_hor Character, horizon date.
                                #' @param B Integer, mixture component, if B is missing will be returned the mixture distribution, otherwise the component distribution non weighted.
                                #' @return Conditional distribution for R_N
                                cdf_R = function(Rt, t_now, t_hor, B){
                                  C_T <- self$Ct(t_hor)
                                  cdf_Y <- self$cdf_Y(Rt, t_now, t_hor, B)
                                  function(x){
                                    z_x <- (1 - x/C_T - self$model$transform$alpha) / self$model$transform$beta
                                    u_x <- suppressWarnings(log(-log(z_x)))
                                    probs <- cdf_Y(u_x)
                                    probs[is.nan(probs)] <- 0
                                    return(probs)
                                  }
                                },
                                #' @description
                                #' Return the conditional expected value for R_N given Rn.
                                #' @param Rt Numeric, solar radiation.
                                #' @param t_now Character, today date.
                                #' @param t_hor Character, horizon date.
                                #' @param B Integer, mixture component, if B is missing will be returned the mixture density, otherwise the component density non weighted.
                                #' @param moment Integer, scalar. Moment order. The default is 1, i.e. the expectation.
                                #' @return Conditional moment for solar radiation
                                e_GHI = function(Rt, t_now, t_hor, B, moment = 1){
                                  pdf_Y <- self$pdf_Y(Rt, t_now, t_hor, B)
                                  C_T <- self$Ct(t_hor)
                                  Rt <- function(y) (C_T * (1 - self$model$transform$alpha - self$model$transform$beta * exp(-exp(y))))^moment * pdf_Y(y)
                                  moment_Rt <- integrate(Rt, lower = -Inf, upper = Inf)$value
                                  return(moment_Rt)
                                },
                                #' @description
                                #' Return the conditional variance value for R_N given Rn.
                                #' @param Rt Numeric, solar radiation.
                                #' @param t_now Character, today date.
                                #' @param t_hor Character, horizon date.
                                #' @param B Integer, mixture component, if B is missing will be returned the mixture density, otherwise the component density non weighted.
                                #' @return Conditional variance for R_N
                                v_GHI = function(Rt, t_now, t_hor, B){
                                  self$e_GHI(Rt, t_now, t_hor, B, 2) - self$e_GHI(Rt, t_now, t_hor, B, 1)^2
                                },
                                #' @description
                                #' Method print
                                print = function(){
                                  self$model$print()
                                  cat("MRP: ", self$lambda, "\n")
                                  cat("Measure: ", self$measure, "\n")
                                  cat("Version: ", private$version, "\n")
                                }
                              ),
                              # ====================================================================================================== #
                              #                                             Private slots
                              # ====================================================================================================== #
                              private = list(
                                version = "1.0.0",
                                ..model = NA,
                                ..measure = "P",
                                ..lambda = 0,
                                ..k1 = 1,
                                ..k2 = 1,
                                ..integral_variance = NA,
                                ..integral_expectation = NA,
                                seasonal_variance = NA
                              ),
                              # ====================================================================================================== #
                              #                                             Active slots
                              # ====================================================================================================== #
                              active = list(
                                #' @field model model
                                model = function(){
                                  private$..model
                                },
                                #' @field measure Character, reference probability measure actually used.
                                measure = function(){
                                  private$..measure
                                },
                                #' @field lambda Numeric, market risk premium actually used.
                                lambda = function(){
                                  private$..lambda
                                }
                              )
)
