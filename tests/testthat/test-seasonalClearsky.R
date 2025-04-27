control <- control_solarModel(clearsky = control_seasonalClearsky(order = 1, order_H0 = 4, period = 365, include.intercept = TRUE, include.trend = FALSE,
                                                                  delta0 = 1.4, lower = 0, upper = 3, by = 0.0001, ntol = 0, quiet = FALSE),
                              stochastic_clearsky = FALSE,
                              seasonal.mean = list(seasonalOrder = 1, include.trend = FALSE, include.intercept = TRUE, monthly.mean = TRUE),
                              mean.model = list(arOrder = 2, maOrder = 1, include.intercept = FALSE),
                              seasonal.variance = list(seasonalOrder = 1, correction = FALSE, monthly.mean = TRUE),
                              variance.model = rugarch::ugarchspec(variance.model = list(garchOrder = c(1,1), variance.targeting = 1),
                                                                   mean.model = list(armaOrder = c(0, 0), include.mean = FALSE)),
                              mixture.model = list(abstol = 1e-20, maxit = 200, match.moments = TRUE), clearsky_threshold = 1.00001,
                              threshold = 10e-7, outliers = 0, garch_variance = FALSE, quiet = FALSE)
# Model specification
spec <- solarModel_spec("Parigi", min_date="2005-01-01", from = "2005-01-01", to = "2023-12-31",
                        control_model = control)

# solarModel specification
spec$data
# Extract the required elements
x <- spec$data$GHI
date <- spec$data$date
lat <- spec$coords$lat
clearsky <- spec$data$clearsky

# Initialize the clearsky model
seasonal_model_Ct <- seasonalClearsky$new(control = spec$control$clearsky)
seasonal_model_Ct$control$include.trend <- FALSE
# Fit the model
seasonal_model_Ct$fit(x, date, lat, clearsky)
# Predict the seasonal values
spec$data$t <- max(spec$data$Year)-spec$data$Year
spec$data$Ct <- seasonal_model_Ct$predict(newdata = spec$data)

ggplot(spec$data)+
  geom_line(aes(n, GHI))+
  geom_line(aes(n, Ct), color = "blue")

self <- seasonal_model_Ct$.__enclos_env__$self
private <- seasonal_model_Ct$.__enclos_env__$private
super <- seasonal_model_Ct$.__enclos_env__$super


seasonal_model_Ct <- clearsky_optimizer(seasonal_model_Ct, spec$data)
# Predict the optimized seasonal values
spec$data$Ct <- seasonal_model_Ct$predict(newdata = spec$data)

ggplot(spec$data)+
  geom_line(aes(n, GHI))+
  geom_line(aes(n, Ct), color = "blue")
