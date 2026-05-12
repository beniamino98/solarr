#' Compute the generic conditional moments of a solarModel object
#'
#' @keywords solarMoments
#' @details Version 1.0.1.
#' @examples
#' model <- solarModel$new(spec)
#' model$fit()
#' # Conditional moments
#' mom <- solarMoment$new(model, "2022-01-01", "2025-01-01")
#' mom$filter()
#' mom$data
#' solarOption_model(model, mom$data)
#'
#' @rdname solarMoment
#' @name solarMoment
#' @export
solarMoment <- R6::R6Class("solarMoment",
                           public = list(
                             weights = NA,
                             X0 = NA,
                             H0 = NA,
                             t_now = "",
                             t_hor = "",
                             place = "",
                             step = 1,
                             initialize = function(model, t_now, t_hor){
                               # Reference place
                               place <- model$spec$place
                               # Pricing date
                               t_now <- as.Date(t_now)
                               self$t_now <- t_now
                               # Horizon date
                               t_hor <- as.Date(t_hor)
                               self$t_hor <- t_hor
                               # ARMA model
                               ARMA <- model$spec$mean.model
                               # GARCH model
                               GARCH <- model$spec$variance.model
                               # Maximum order of AR / GARCH
                               lag_max <- max(c(ARMA$order, GARCH$order))
                               # Filter data between (t_now - lag_max + 1) and t_hor
                               data <- dplyr::filter(model$data, date >= (t_now - lag_max + 1) & date <= t_hor)
                               # Check if the data above t_hor are known
                               if (max(data$date) < t_hor) {
                                 # Add unknown variables
                                 new_dates <- tibble(date = seq.Date(max(data$date)+1, t_hor, 1),
                                                     Year = lubridate::year(date),
                                                     Month = lubridate::month(date),
                                                     Day = lubridate::day(date),
                                                     n = number_of_day(date),
                                                     isTrain = FALSE,
                                                     weights = 0,
                                                     Ct = model$spec$seasonal_model_Ct$predict(n),
                                                     Yt_bar = model$spec$seasonal.mean$predict(n),
                                                     sigma_bar = sqrt(model$spec$seasonal.variance$predict(n)),
                                                     GHI_bar = model$spec$transform$iRY(Yt_bar, Ct)) %>%
                                   dplyr::left_join(model$monthly_data, by = "Month")
                                 data <- dplyr::bind_rows(data, new_dates)
                               }
                               # Forecast seasonal variance
                               data$sigma_bar <- data$sigma_bar * data$sigma_uncond
                               # Store bounds parameters
                               Ct <- tail(data$Ct, 1)[[1]]
                               alpha <- model$spec$transform$alpha[[1]]
                               beta <- model$spec$transform$beta[[1]]
                               private$..R_min_max <- list(R_min = Ct * (1 - alpha - beta), R_max = Ct * (1 - alpha))
                               # ****************************************************************
                               # Initialize derived columns
                               weights <- data[-c(1:lag_max), c("date", "Month", "Yt_tilde_hat", "sigma_bar", "Yt_bar", "Yt_tilde_uncond",
                                                                "mu1", "mu2", "sd1", "sd2", "p1", "p2")]
                               # Step ahead
                               weights$step <- 1:nrow(weights)
                               # ARMA intercept
                               weights$ARMA_intercept <- 0
                               # ARMA weights
                               weights$psi_j <- weights$psi2_j <- 0
                               # GARCH moments
                               weights$e_sigma1 <- weights$e_sigma2 <- weights$e_sigma4  <- NA
                               # GARCH variances
                               weights$v_sigma <- weights$v_sigma2 <- 0
                               # GARCH covariances
                               weights$psi_hs <- 0
                               # Mixture moments
                               weights$e_u <- weights$e_u2 <- weights$e_u4 <- NA
                               # Mixture variance
                               weights$v_u <- NA
                               # Cumulated moments of U_{t+h}
                               weights$e_Uh <- 0
                               weights$v_Uh <- 0
                               # Expectation of Yt (unconditional)
                               weights$e_Yt <- 0
                               # Variance of Yt (unconditional)
                               weights$Sigma2 <- 1
                               # Store ARMA weights
                               self$weights <- weights
                               # ****************************************************************
                               # Known data at time t_now (used as vector only in state-space forecast)
                               df_t <- data[c(1:lag_max),]
                               # Must be reversed ordered, i.e. most recent (t, t-1, t-2, ...)
                               df_t <- df_t[order(df_t$date, decreasing = TRUE),]
                               # State vector
                               AR_order <- ARMA$order[1]
                               Y0 <- df_t$Yt_tilde[1:AR_order]
                               names(Y0) <- paste0("Y0_", df_t$date[1:AR_order])
                               eps0 <- c()
                               MA_order <- ARMA$order[2]
                               if (MA_order > 0){
                                 eps0 <- df_t$eps[1:MA_order]
                                 names(eps0) <- paste0("eps0_", df_t$date[1:MA_order])
                               }
                               # Store state vector
                               self$X0 <- c(Y0, eps0)
                               #self$H0 <- c(data$eps_tilde[lag_max]^2, data$sigma[lag_max]^2)
                               #names(self$H0) <- paste0(c("eps2_tilde_0_", "sigma2_0_"), data$date[lag_max])
                               lag.arch  <- max(c(1, GARCH$archOrder))
                               lag.garch <- max(c(1, GARCH$garchOrder))
                               self$H0 <- c(df_t$sigma[1:lag.garch]^2, df_t$eps_tilde[1:lag.arch]^2)
                               lab.arch <- paste0("eps2_tilde_0_", df_t$date[1:lag.arch])
                               lab.garch <- paste0("sigma2_0_", df_t$date[1:lag.garch])
                               names(self$H0) <- c(lab.garch, lab.arch)
                               # ****************************************************************
                               # Initialize moments data
                               df_T <- dplyr::filter(data, date == t_hor)
                               # Numbers of steps ahead
                               h <- nrow(weights)
                               self$step <- h
                               # Remove all the known points except t
                               if (lag_max > 1){
                                 data <- data[-c(1:min(c(1, lag_max-1))),]
                               }
                               # Final object
                               private$..data <- dplyr::tibble(
                                 date = t_hor,
                                 lag_max = lag_max,
                                 h = h,
                                 Year = lubridate::year(t_hor),
                                 Month = lubridate::month(t_hor),
                                 Day = lubridate::day(t_hor),
                                 e_Yt = 0,
                                 sd_Yt = 1,
                                 M_Y1 = 0,
                                 S_Y1 = 1,
                                 M_Y0 = 1,
                                 S_Y0 = 0,
                                 p1 = df_T$p1[1],
                                 GHI_bar = df_T$GHI_bar[1],
                                 Ct = Ct,
                                 alpha = alpha,
                                 beta = beta
                               )
                               # *******************************************************************************
                               #  0) ARMA summations
                               # *******************************************************************************
                               # Store ARMA matrices and vectors
                               private$..ARMA <- list(intercept = ARMA$intercept, A = ARMA$A, b = ARMA$b)
                               # Store GARCH parameters
                               private$..GARCH <- list(A = GARCH$A, b = GARCH$b, d = GARCH$d,
                                                       beta = GARCH$beta, e1 = c(1, GARCH$d[-1]))

                             },
                             filter = function(theta = 0, B = NULL, t_cond = NULL){
                               self$filter_ARMA()
                               self$filter_NM(theta, B, t_cond)
                               self$filter_GARCH()
                               self$filter_weights()
                             },
                             filter_ARMA = function(){
                               # Extract horizon
                               # h <- nrow(self$weights)
                               # ARMA data
                               ARMA <- self$ARMA
                               # ARMA forecast
                               df_tT <- ARMA_forecast(self$step, self$X0, ARMA$A, ARMA$b, intercept = ARMA$intercept)$weights[[1]]
                               # *******************************************************************************
                               # Store ARMA summations
                               self$weights$psi_j <- df_tT$psi_j
                               self$weights$psi2_j <- df_tT$psi2_j
                               self$weights$Yt_tilde_hat <- df_tT$Yt_tilde_hat
                             },
                             filter_NM = function(theta = 0, B = NULL, t_cond = NULL){
                               # Extract weights
                               df_tT <- self$weights
                               # *******************************************************************************
                               #  1) Conditioning variable
                               # *******************************************************************************
                               B1 <- B0 <- rep(1, self$step)
                               if (!is.null(t_cond)) {
                                 if (length(t_cond) == length(B)) {
                                   idx_cond <- df_tT$date %in% t_cond
                                   B1[idx_cond] <- B / df_tT$p1[idx_cond]
                                   B0[idx_cond] <- (1 - B) / df_tT$p2[idx_cond]
                                 }
                               }
                               # *******************************************************************************
                               #  2) Mixture moments
                               # *******************************************************************************
                               if (any(theta != 0)) {
                                 df_tT$mu1 <- df_tT$mu1 + df_tT$sd1 * theta
                                 df_tT$mu2 <- df_tT$mu2 + df_tT$sd2 * theta
                               }
                               # Compute mixture moments
                               df_tT$e_u <- (df_tT$mu1 * df_tT$p1) * B1 + (df_tT$mu2 * df_tT$p2) * B0
                               # Second moment
                               df_tT$e_u2 <- (df_tT$mu1^2 + df_tT$sd1^2) * df_tT$p1 * B1 + (df_tT$mu2^2 + df_tT$sd2^2) * df_tT$p2 * B0
                               # Fourth moment
                               df_tT$e_u4 <- (3 * df_tT$sd1^4 + 6 * df_tT$mu1^2 * df_tT$sd1^2 + df_tT$mu1^4) * df_tT$p1 * B1 + (3 * df_tT$sd2^4 + 6 * df_tT$mu2^2 * df_tT$sd2^2 + df_tT$mu2^4) * df_tT$p2 * B0
                               # Variance
                               df_tT$v_u <- df_tT$e_u2 - df_tT$e_u^2
                               self$weights <- df_tT
                             },
                             filter_GARCH = function(){
                               # Extract horizon
                               h <- self$step # nrow(self$weights)
                               # Extract GARCH parameters
                               GARCH <- self$GARCH
                               # *******************************************************************************
                               #  3) Forecast GARCH moments
                               # *******************************************************************************
                               df_tT <- sGARCH_forecast(h, GARCH$A, GARCH$b, GARCH$d, GARCH$e1, self$H0, self$weights$e_u2, self$weights$e_u4)
                               # Moments GARCH variance (exact)
                               self$weights$e_sigma1 <- df_tT$e_sigma
                               self$weights$e_sigma2 <- df_tT$e_sigma2
                               self$weights$e_sigma4 <- df_tT$e_sigma4
                               self$weights$v_sigma2 <- df_tT$v_sigma2
                               self$weights$v_sigma <- df_tT$v_sigma
                               psi_j <- self$weights$psi_j

                               for(j in 1:h){
                                 self$weights$psi_hs[j] <- sum(df_tT$cv_sigma[[j]] * psi_j[j] * psi_j[-j])
                               }
                             },
                             filter_weights = function(){
                               df_tT <- self$weights
                               h <- self$step # nrow(df_tT)
                               # *******************************************************************************
                               #  3)  Compute the series of psi
                               # *******************************************************************************
                               # Compute weights for expectations
                               df_tT$e_Uh <- df_tT$psi_j * df_tT$sigma_bar * df_tT$e_sigma1 * df_tT$e_u
                               # Variance
                               df_tT$v_Uh <- df_tT$psi2_j * df_tT$sigma_bar^2 * (df_tT$e_sigma2 * df_tT$e_u2 - df_tT$e_sigma1^2 * df_tT$e_u^2)
                               # Forecasted total expectation Yt tilde
                               df_tT$e_Yt_tilde <- df_tT$Yt_tilde_hat + cumsum(df_tT$e_Uh)
                               # Forecasted total expectation Yt
                               df_tT$e_Yt <- df_tT$e_Yt_tilde + df_tT$Yt_tilde_uncond + df_tT$Yt_bar
                               # Forecasted total variance of Yt
                               df_tT$Sigma2 <- cumsum(df_tT$v_Uh)
                               # *******************************************************************************
                               #  4)  Approximate the multinomial mixture with a 2-component GM
                               # *******************************************************************************
                               # Last value
                               df_T <- tail(df_tT, 1)
                               # Next step parameters
                               # Approximated mixture means
                               M_Y <- df_T$Yt_tilde_hat + df_T$Yt_tilde_uncond + df_T$Yt_bar + df_T$sigma_bar * df_T$e_sigma1 * c(df_T$mu1, df_T$mu2)
                               # Approximated mixture variances
                               S2_Y <- c(df_T$sd1^2, df_T$sd2^2) * df_T$sigma_bar^2 * df_T$e_sigma2 + df_T$v_sigma * df_T$sigma_bar^2 * c(df_T$mu1^2, df_T$mu2^2)
                               if (h > 1){
                                 # Add conditional covariances
                                 rho2_U <- sum((df_tT$psi_hs * df_tT$e_u)[1:(h-1)]) + df_T$psi_hs * c(df_T$mu1, df_T$mu2)
                                 # Approximated mixture means
                                 M_Y <- M_Y + sum(df_tT$e_Uh[1:(h-1)])
                                 # Approximated mixture variances
                                 S2_Y <- S2_Y + sum(df_tT$v_Uh[1:(h-1)]) + 2*rho2_U
                                 # Update total variance
                                 df_tT$Sigma2 <- df_tT$Sigma2 + cumsum((df_tT$psi_hs * df_tT$e_u))
                               }
                               S_Y <- sqrt(S2_Y)
                               # *************************************
                               # Update total moments
                               private$..data$e_Yt <- tail(df_tT$e_Yt, 1)
                               private$..data$sd_Yt <- ifelse(h == 1, df_T$sigma_bar * df_T$e_sigma1, sqrt(tail(df_tT$Sigma2, 1)))
                               # Update mixture moment
                               private$..data$M_Y1 <- M_Y[1]
                               private$..data$M_Y0 <- M_Y[2]
                               private$..data$S_Y1 <- S_Y[1]
                               private$..data$S_Y0 <- S_Y[2]
                               # Update psi_j
                               self$weights <- df_tT
                             },
                             print = function(){
                               cat("----------------- Solar Moment  ----------------- \n")
                               cat(paste0("Location: ", self$place, "\n"))
                               cat(paste0("t_now: ", self$t_now, "\n"))
                               cat(paste0("t_hor: ", self$t_hor, "\n"))
                               cat(paste0("Steps ahead: ", self$step, "\n"))
                             }
                           ),
                           private = list(
                             ..data = NA,
                             ..R_min_max = c(),
                             ..ARMA = list(),
                             ..GARCH = list()
                           ),
                           active = list(
                             data = function(){
                               private$..data
                             },
                             ARMA = function(){
                               private$..ARMA
                             },
                             GARCH = function(){
                               private$..GARCH
                             },
                             R_min_max = function(){
                               private$..R_min_max
                             },
                             M_Y = function(){
                               c(private$..data$M_Y1, private$..data$M_Y0)
                             },
                             S_Y = function(){
                               c(private$..data$S_Y1, private$..data$S_Y0)
                             },
                             probs = function(){
                               c(private$..data$p1, private$..data$p1)
                             }
                           ))



#sm <- solarMoment$new(model, "2022-01-01", "2022-03-31")
#sm$filter()
#sm$data

#private <- sm$.__enclos_env__$private
#self <- sm$.__enclos_env__$self

