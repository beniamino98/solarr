# Solar model
control <- control_solarModel(outliers_quantile = 0, mean.model = list(arOrder = 1), garch_variance = TRUE)
spec <- solarModel_spec("Bologna", from="2005-01-01", control_model = control)
Bologna <- solarModel$new(spec)
Bologna$fit()

Bologna$loglik

# Dates for 1 step-ahead
t_now <- as.Date("2020-01-01")
t_hor <- t_now + 1

max_AR_GARCH <- max(c(Bologna$AR_model_Yt$order, Bologna$GARCH$arch_order, Bologna$GARCH$garch_order))
# Filter data for forecast
df_t <- filter(Bologna$data, date %in% c(t_now - max_AR_GARCH + 1, t_now))
# Filter data at horizon
df_T <- filter(Bologna$data, date == t_hor)

# Forecasted seasonal mean
YT_bar <- Bologna$seasonal_model_Yt$predict(t_hor)
# Forecasted seasonal variance
sigmaT_bar <- sqrt(Bologna$seasonal_variance$predict(t_hor))

# Expectation
YT_tilde_hat <- Bologna$AR_model_Yt$next_step(df_t$Yt_tilde[c(2,1)], 1)
YT_bar + YT_tilde_hat
filter(Bologna$moments$conditional, date == t_hor)$e_Yt

# Variance
sigmaT <- Bologna$GARCH$next_step(df_t$eps_tilde[2], df_t$sigma[2], 1)
sigmaT
df_T$sigma

# Mixture parameters
mu1 <- Bologna$NM_model$mu1$predict(t_hor)
mu2 <- Bologna$NM_model$mu2$predict(t_hor)
sd1 <- Bologna$NM_model$sd1$predict(t_hor)
sd2 <- Bologna$NM_model$sd2$predict(t_hor)
p <- Bologna$NM_model$prob$predict(t_hor)

filter(Bologna$moments$conditional, date == t_hor)$e_Yt_up
YT_bar + YT_tilde_hat + sigmaT * sigmaT_bar * mu1

filter(Bologna$moments$conditional, date == t_hor)$e_Yt_dw
YT_bar + YT_tilde_hat + sigmaT * sigmaT_bar * mu2

filter(Bologna$moments$conditional, date == t_hor)$sd_Yt_up
sigmaT * sigmaT_bar * sd1

filter(Bologna$moments$conditional, date == t_hor)$sd_Yt_dw
sigmaT * sigmaT_bar * sd2


pow_matrix <- function(A, pow = 1){
  if (pow == 0){
    A[1,] <- 1
    A[-1,] <- 0
    return(A)
  }
  if(pow == 1) return(A)
  A_pow <- A
  for(i in 2:pow){
    A_pow <- A_pow %*% A
  }
  return(A_pow)
}

A <- Bologna$AR_model_Yt$Phi
d <- Bologna$AR_model_Yt$d

sigmaT * sigmaT_bar * mu2


sigmaT * sigmaT_bar * mu2 +
(pow_matrix(A, 1) %*% c(Bologna$NM_model$moments[1,]$mean, 0))[1] +
(pow_matrix(A, 2) %*% c(Bologna$NM_model$moments[1,]$mean, 0))[1] +
(pow_matrix(A, 3) %*% c(Bologna$NM_model$moments[1,]$mean, 0))[1] +
(pow_matrix(A, 4) %*% c(Bologna$NM_model$moments[1,]$mean, 0))[1] +
(pow_matrix(A, 5) %*% c(Bologna$NM_model$moments[1,]$mean, 0))[1] +
(pow_matrix(A, 6) %*% c(Bologna$NM_model$moments[1,]$mean, 0))[1] +
(pow_matrix(A, 7) %*% c(Bologna$NM_model$moments[1,]$mean, 0))[1] +
(pow_matrix(A, 8) %*% c(Bologna$NM_model$moments[1,]$mean, 0))[1]+
  (pow_matrix(A, 9) %*% c(Bologna$NM_model$moments[1,]$mean, 0))[1]


Bologna$AR_model_Yt$Phi


# Check pdf integrals
df_t <- filter(Bologna$data, date == "2022-07-01")
mom <- filter(Bologna$moments$conditional, date == "2022-07-01")

pdf_Y <- function(x) df_t$p1* dnorm(x, mom$e_Yt_up, mom$sd_Yt_up) + df_t$p2* dnorm(x, mom$e_Yt_dw, mom$sd_Yt_dw)
cdf_Y <- function(x) df_t$p1* pnorm(x, mom$e_Yt_up, mom$sd_Yt_up) + df_t$p2* pnorm(x, mom$e_Yt_dw, mom$sd_Yt_dw)
alpha <- Bologna$transform$alpha
beta <-  Bologna$transform$beta
Ct <- df_t$Ct
Ct * (1- alpha)*cdf_Y(K_Y) - Ct*beta*integrate(function(x) exp(-exp(x))*pdf_Y(x), lower = -Inf, upper = Bologna$R_to_Y(K, df_t$date))$value
K <- df_t$GHI_bar
K_Y <- Bologna$R_to_Y(K, df_t$date)

(K - Ct*(1-alpha))*cdf_Y(K_Y) + Ct* beta*integrate(function(x) exp(-exp(x))*pdf_Y(x), lower = -Inf, upper = K_Y)$value

psolarGHI(K, Ct, alpha, beta, cdf_Y)
cdf_Y(K_Y)

psolarGHI(K, Ct, alpha, beta, cdf_Y) * K -
integrate(function(x) x*dsolarGHI(x, Ct, alpha, beta, pdf_Y), lower = Ct*(1-alpha-beta), upper = K)$value

lm(I(Yt_tilde) ~lag(Yt_tilde, 1) + lag(Yt_tilde, 2)  + lag(Yt_tilde, 3) + lag(Yt_tilde, 4), data = Bologna$data) %>% summary()

mom <- Bologna$moments$conditional
mom$Ct <- Bologna$seasonal_model_Ct$predict(mom$date)
mom$p1 <- Bologna$NM_model$prob$predict(mom$date)
mom$p2 <- 1-mom$p1
mom$GHI <- Bologna$data$GHI
mom$loglik <- NA
for(i in 1:nrow(mom)){
  print(i)
  df_t <- mom[i,]
  pdf_Y <- function(x) df_t$p1* dnorm(x, df_t$e_Yt_up, df_t$sd_Yt_up) + df_t$p2* dnorm(x, df_t$e_Yt_dw, df_t$sd_Yt_dw)
  mom$loglik[i] <- log(dsolarGHI(df_t$GHI, df_t$Ct, alpha, beta, pdf_Y))
}
Bologna$loglik
sum(mom$loglik, na.rm = TRUE)

sum(mom$loglik[!is.infinite(mom$loglik)])

which(is.infinite(mom$loglik))

Bologna$.__enclos_env__$private$outliers$index


model <- Bologna$clone(T)
solarModel_logLik <- function(model){
  # Conditional moments
  mom <- model$moments$conditional
  # Clearksy
  mom$Ct <- model$seasonal_model_Ct$predict(mom$date)
  # Probabilities
  mom$p1 <- model$NM_model$prob$predict(mom$date)
  mom$p2 <- 1 - mom$p1
  # Realized radiation
  mom$Rt <- model$data$GHI
  # Weights
  mom$weights <- model$data$weights
  # Remove non train data
  mom <- dplyr::filter(mom, weights != 0)
  # Compute log likelihood
  mom$loglik <- NA
  for(i in 1:nrow(mom)){
    df_t <- mom[i,]
    pdf_Y <- function(x) df_t$p1* dnorm(x, df_t$e_Yt_up, df_t$sd_Yt_up) + df_t$p2* dnorm(x, df_t$e_Yt_dw, df_t$sd_Yt_dw)
    mom$loglik[i] <- log(dsolarGHI(df_t$Rt, df_t$Ct, alpha, beta, pdf_Y))
  }
  sum(mom$loglik, na.rm = TRUE)
}

solarModel_logLik(model)

