control <- control_solarModel(clearsky = control_seasonalClearsky(order = 1, period = 365, include.intercept = TRUE, delta0 = 1.4,
                                                                  lower = 0, upper = 3, by = 0.001, ntol = 0, quiet = FALSE),
                              stochastic_clearsky = FALSE,
                              seasonal.mean = list(seasonalOrder = 1, include.intercept = TRUE, monthly.mean = TRUE),
                              mean.model = list(arOrder = 2, maOrder = 1, include.intercept = FALSE),
                              seasonal.variance = list(seasonalOrder = 1, correction = TRUE, monthly.mean = TRUE),
                              variance.model = rugarch::ugarchspec(variance.model = list(garchOrder = c(1,1), variance.targeting = 1),
                                                                   mean.model = list(armaOrder = c(0, 0), include.mean = FALSE)),
                              mixture.model = list(abstol = 1e-20, maxit = 100, correction = TRUE), clearsky_threshold = 1.001,
                              threshold = 0.01, outliers_quantile = 0, garch_variance = TRUE, quiet = FALSE)

# Model specification
spec <- solarModel_spec("Parigi", min_date="2005-01-01", from = "2005-01-01", to = "2023-12-31",  control_model = control)
model <- solarModel$new(spec)
# Model fit
model$fit()
# Calibrate the model
model_cal <- solarModel_calibrator(model)

model$loglik
model_cal$loglik

# Check Option prices
control_options <- control_solarOption(nyears = c(2005, 2024))
solarOption_historical(model_cal, control_options = control_options)
# Conditional price
solarOption_model(model_cal, model_cal$moments$conditional, theta = 0, implvol = 1, put = TRUE, control_options = control_options)
solarOption_model(model, model$moments$conditional, theta = 0, implvol = 1, put = TRUE, control_options = control_options)
# Unconditional price
solarOption_model(model_cal, model_cal$moments$unconditional, theta = 0, implvol = 1, put = TRUE, control_options = control_options)
solarOption_model(model, model$moments$unconditional, theta = 0, implvol = 1, put = TRUE, control_options = control_options)
df_m <- filter(model_cal$moments$unconditional, Year == 2006)
solarOption_model(model_cal, df_m, theta = 0.012, implvol = 1, put = TRUE, control_options = control_options)



# Compare parameters
model$ARMA$coefficients
model_cal$ARMA$coefficients
ggplot()+
  geom_line(data = model$seasonal_data, aes(n, Yt_tilde_uncond))+
  geom_line(data = model_cal$seasonal_data, aes(n, Yt_tilde_uncond), color = "blue")
# Seasonal variances
ggplot()+
  geom_line(data = model$seasonal_data, aes(n, sigma_bar))+
  geom_line(data = model_cal$seasonal_data, aes(n, sigma_bar), color = "blue")
ggplot()+
  geom_line(data = model$seasonal_data, aes(n, sigma_uncond))+
  geom_line(data = model_cal$seasonal_data, aes(n, sigma_uncond), color = "blue")
# Mixture parameters
ggplot()+
  geom_line(data = model$seasonal_data, aes(n, mu1))+
  geom_line(data = model_cal$seasonal_data, aes(n, mu1), color = "blue")
ggplot()+
  geom_line(data = model$seasonal_data, aes(n, mu2))+
  geom_line(data = model_cal$seasonal_data, aes(n, mu2), color = "blue")
ggplot()+
  geom_line(data = model$seasonal_data, aes(n, sd1))+
  geom_line(data = model_cal$seasonal_data, aes(n, sd1), color = "blue")
ggplot()+
  geom_line(data = model$seasonal_data, aes(n, sd2))+
  geom_line(data = model_cal$seasonal_data, aes(n, sd2), color = "blue")
ggplot()+
  geom_line(data = model$seasonal_data, aes(n, p1))+
  geom_line(data = model_cal$seasonal_data, aes(n, p1), color = "blue")

# VaR tests
bernoulli_test <- function(x, alpha = 0.05, correct_variance = FALSE){

  # Number of observations
  N <- length(x)
  # Number of violations
  NT <- sum(x)
  # Theoric expectation
  e_NT <- alpha
  # Theoric variance
  v_NT <- alpha * (1 - alpha)
  # Empiric expectation
  e_NT_hat <- NT / N
  # Empiric variance
  v_NT_hat <- e_NT_hat * (1 - e_NT_hat)

  # Statistic test ~ N(0,1)
  stat <- sqrt(N) * (e_NT_hat - e_NT) / sqrt(v_NT)
  # Correction with empiric variance
  if (correct_variance) {
    stat <- sqrt(v_NT) / sqrt(v_NT_hat) * stat
  }
  # Compute pvalue
  pvalue <- 1 - pnorm(abs(stat))

  tibble(
    method = "Bernoulli Test for VaR",
    parameter = alpha,
    stat = stat,
    pvalue = pvalue,
    empiric = e_NT_hat
  )
}

# Confidence intervals
ci <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
VaR_tests <- list()
VaR_tests_cal <- list()
for(i in 1:length(ci)){
  cat(paste0("VaR for ci: ", "\033[1;34m", ci[i], "\033[0m \n"))
  moments <- filter(model$moments$conditional, Year <= 2023)
  # VaR
  VaR <- model$VaR(moments , ci = ci[i])
  # Test
  test <- bernoulli_test(VaR$et, ci[i])
  # Store the results
  VaR_tests[[i]] <- test
  pval <- round(test$pvalue, digits = 4)
  cat(paste0("VaR (model) for ci: ", "\033[1;34m", ci[i], "\033[0m", " pvalue: ", pval, "\n"))
  # VaR
  moments <- filter(model_cal$moments$conditional, Year <= 2023)
  VaR <- model_cal$VaR(moments , ci = ci[i])
  # Test
  test <- bernoulli_test(VaR$et, ci[i])
  # Store the results
  VaR_tests_cal[[i]] <- test
  pval <- round(test$pvalue, digits = 4)
  cat(paste0("VaR (model_cal) for ci: ", "\033[1;34m", ci[i], "\033[0m", " pvalue: ", pval, "\n"))
}



