library(solarr)
# Control model
control <- control_solarModel(clearsky = control_seasonalClearsky(ntol = 0),
                              stochastic_clearsky = FALSE,
                              seasonal.mean = list(seasonalOrder = 1, include.H0 = FALSE, include.intercept = TRUE, monthly.mean = TRUE),
                              mean.model = list(arOrder = 1, maOrder = 1, include.intercept = FALSE),
                              seasonal.variance = list(seasonalOrder = 1, correction = TRUE, monthly.mean = TRUE),
                              variance.model = rugarch::ugarchspec(variance.model = list(garchOrder = c(1,1), variance.targeting = 1),
                                                                   mean.model = list(armaOrder = c(0, 0), include.mean = FALSE)),
                              mixture.model = list(abstol = 1e-20, maxit = 200), clearsky_threshold = 1.0001,
                              threshold = 0.01, outliers_quantile = 0, garch_variance = TRUE, quiet = FALSE)
# Control Options
control_options <- control_solarOption(nyears = c(2005, 2024))
# Model specification
spec <- solarModel_spec("Oslo", from="2005-01-01", to = "2023-01-01", control_model = control)
model <- solarModel$new(spec)
# Model fit
model$fit()

# Reference dates
t_now <- model$dates$train$to-1
#t_now <- as.Date("2021-12-31")
t_init <- t_now + 1
#t_init <- as.Date("2022-01-01")
t_hor <- as.Date(paste0(lubridate::year(t_init), "-12-31"))
#t_hor <- as.Date("2022-12-31")
cat(paste0("t_now:  ", t_now, "\nt_init: ", t_init, "\nt_hor:  ", t_now))

# SoRad portfolio
portfolio <- SoRadPorfolio(model, t_now, t_init, t_hor)

# TEST 1: conditional prices for one day
object <- portfolio[[1]]
df_n <- dplyr::filter(model$moments$conditional, date == object$t_hor)
solarOption_pricing(object, df_n, theta = , put = TRUE)

solarOption_historical(model, put = TRUE, control = control_solarOption(nyears = c(2005, lubridate::year(t_hor)+1)))

# TEST 2: conditional prices for all the days of the year
moments <- dplyr::filter(model$moments$conditional, date >= tail(portfolio, 1)[[1]]$t_init & date <= tail(portfolio, 1)[[1]]$t_hor)
solarOption_model(model, moments, portfolio, theta = 0, put = TRUE, control = control_solarOption())

# TEST 3: unconditional prices for all the days of the year
moments <- dplyr::filter(model$moments$unconditional, date >= tail(portfolio, 1)[[1]]$t_init & date <= tail(portfolio, 1)[[1]]$t_hor)
solarOption_model(model, moments, portfolio, theta = 0, put = TRUE, control = control_solarOption())

# TEST 4: unconditional prices for all the days of the year
t_seq <- unlist(purrr::map(portfolio, ~as.Date(.x$t_hor)))
moments <- model$Moments(tail(portfolio, 1)[[1]]$t_now-20, as.Date(t_seq))
solarOption_model(model, moments, portfolio, theta = 0, put = TRUE, control = control_solarOption())

moments_2 <- moments
moments <- dplyr::filter(model$moments$unconditional, date >= tail(portfolio, 1)[[1]]$t_init & date <= tail(portfolio, 1)[[1]]$t_hor)

ggplot()+
  geom_line(data = moments_2, aes(date, M_Y1), color = "red")+
  geom_line(data = moments, aes(date, M_Y1))

ggplot()+
  geom_line(data = moments, aes(date, M_Y0))+
  geom_line(data = moments_2, aes(date, M_Y0), color = "red")

ggplot()+
  geom_line(data = moments, aes(date, S_Y1))+
  geom_line(data = moments_2, aes(date, S_Y1), color = "red")

ggplot()+
  geom_line(data = moments, aes(date, S_Y0))+
  geom_line(data = moments_2, aes(date, S_Y0), color = "red")

