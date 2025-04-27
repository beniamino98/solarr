# Test Options prices
test_prices <- function(model){
  # Historical
  df1 <- solarOption_historical(model, put = T, nmonths = 1:12, control_options = control_solarOption(nyears = c(2005, 2024)))
  df1
  # Conditional
  mom <- model$moments$conditional
  df2 <- solarOption_model(model, mom, put = T, control_options = control_solarOption(nyears = c(2005, 2024)))
  df2
  # Unconditional
  mom <- filter(model$moments$unconditional, Year == 2023)
  df3 <- solarOption_model(model, mom, put = T, control_options = control_solarOption(nyears = c(2005, 2024)))
  df3

  t_seq <- seq.Date(as.Date("2018-01-01"), as.Date("2018-12-31"), 1)
  mom <- purrr::map_df(t_seq, ~solarModel_moments(.x-100, .x, model$data, model$ARMA,
                                           model$GARCH, model$NM_model, model$transform))
  #mom <- purrr::map_df(t_seq, ~model$Moments(.x-100, .x))
  mom$GHI_bar <- filter(model$data, Year == 2018)$GHI_bar
  df4 <- solarOption_model(model, mom, put = T, nmonths = 1:12)$payoff
  sum(df4$premium)

  # Cumulated returns
  df_V_hist <- (df1$payoff %>% filter(Year != 2024) %>% group_by(Year) %>% summarise(payoff = sum(payoff)))$payoff
  df_Vt <- (df2$payoff %>% group_by(Year) %>% summarise(payoff = sum(premium)))$payoff
  df_V0 <- rep(df3$payoff_year$premium, length(df_Vt))
  df_V0_m <- rep(sum(df4$premium), length(df_Vt))
  tibble(
    Place = model$place,
    V_hist = df1$payoff_year$premium,
    V0 = df3$payoff_year$premium,
    V0_m = sum(df4$premium),
    Vt = df2$payoff_year$premium,
    V_hist_V0 = 100*(V_hist/V0-1),
    V_hist_Vt = 100*(V_hist/Vt-1),
    Vt_V0 = 100*(Vt/V0-1),
    V0_V0_m = 100*(V0/V0_m-1),
    e_ret_Vt = mean((df_V_hist-df_Vt)/df_Vt)*100,
    e_ret_V0 = mean((df_V_hist-df_V0)/df_V0)*100,
    e_ret_V0_m = mean((df_V_hist-df_V0_m)/df_V0_m)*100
  ) %>%
    mutate(
      V_hist_V0 = paste0(format(V_hist_V0, digits = 3), " %"),
      V_hist_Vt = paste0(format(V_hist_Vt, digits = 3), " %"),
      Vt_V0 = paste0(format(Vt_V0, digits = 3), " %"),
      V0_V0_m = paste0(format(V0_V0_m, digits = 3), " %"),
      e_ret_Vt = paste0(format(e_ret_Vt, digits = 3), " %"),
      e_ret_V0 = paste0(format(e_ret_V0, digits = 3), " %"),
    )
}

# When true the model will be re-fitted otherwise taken from .GlobalEnv
fit_again <- FALSE

#self <- model$.__enclos_env__$self
#private <- model$.__enclos_env__$private

# Model specification
spec <- solarModel_spec$new()
spec$set_params(min_pos = 1, max_pos = 1, transform_delta = 0.05)
spec$set_seasonal.mean(monthly.mean = TRUE)
spec$set_seasonal.variance(monthly.mean = TRUE)
spec$set_mean.model(arOrder = 1, maOrder = 1)
spec$set_mixture.model(match.moments = TRUE)
spec$set_variance.model(garchOrder = 1)
spec$specification("Bologna", min_date="2005-01-01", from = "2005-01-01", to = "2023-12-31", max_date = "2023-12-31",)
# Model fit
if (!exists(spec$place) | fit_again) {
  model <- solarModel$new(spec)
  model$fit()
  assign(model$place, model)
} else {
  model <- get(spec$place)$clone(TRUE)
}
print(test_prices(model))

# Model specification
spec <- solarModel_spec$new()
spec$set_params(min_pos = 1, max_pos = 1, transform_delta = 0.05)
spec$set_mean.model(arOrder = 2, maOrder = 1)
spec$set_mixture.model(match.moments = TRUE)
spec$specification("Milano", min_date="2005-01-01", from = "2005-01-01", to = "2023-12-31", max_date = "2023-12-31",)
# Model fit
if (!exists(spec$place) | fit_again) {
  model <- solarModel$new(spec)
  model$fit()
  assign(model$place, model)
} else {
  model <- get(spec$place)$clone(TRUE)
}
print(test_prices(model))

# Model specification
spec <- solarModel_spec$new()
spec$set_params(min_pos = 1, max_pos = 1, transform_delta = 0.05)
spec$set_mean.model(arOrder = 2, maOrder = 1)
spec$set_variance.model(archOrder = 1, garchOrder = 0)
spec$set_mixture.model(match.moments = TRUE)
spec$specification("Roma", min_date="2005-01-01", from = "2005-01-01", to = "2023-12-31", max_date = "2023-12-31",)
# Model fit
if (!exists(spec$place) | fit_again) {
  model <- solarModel$new(spec)
  model$fit()
  assign(model$place, model)
} else {
  model <- get(spec$place)$clone(TRUE)
}
print(test_prices(model))

# Model specification
spec <- solarModel_spec$new()
spec$set_params(min_pos = 1, max_pos = 1, transform_delta = 0.01)
spec$set_mean.model(arOrder = 2, maOrder = 1)
spec$set_mixture.model(match.moments = TRUE)
spec$specification("Palermo", min_date="2005-01-01", from = "2005-01-01", to = "2023-12-31", max_date = "2023-12-31",)
# Model fit
if (!exists(spec$place) | fit_again) {
  model <- solarModel$new(spec)
  model$fit()
  assign(model$place, model)
} else {
  model <- get(spec$place)$clone(TRUE)
}
print(test_prices(model))

# Model specification
spec <- solarModel_spec$new()
spec$set_params(min_pos = 1, max_pos = 1, transform_delta = 0.01)
spec$set_mean.model(arOrder = 2, maOrder = 1)
spec$set_mixture.model(match.moments = TRUE)
spec$specification("Catania", min_date="2005-01-01", from = "2005-01-01", to = "2023-12-31", max_date = "2023-12-31",)
# Model fit
if (!exists(spec$place) | fit_again) {
  model <- solarModel$new(spec)
  model$fit()
  assign(model$place, model)
} else {
  model <- get(spec$place)$clone(TRUE)
}
print(test_prices(model))

# Model specification
spec <- solarModel_spec$new()
spec$set_params(min_pos = 1, max_pos = 1, transform_delta = 0.01)
spec$set_mean.model(arOrder = 2, maOrder = 1)
spec$set_mixture.model(match.moments = TRUE)
spec$specification("Viterbo", min_date="2005-01-01", from = "2005-01-01", to = "2023-12-31", max_date = "2023-12-31",)
# Model fit
if (!exists(spec$place) | fit_again) {
  model <- solarModel$new(spec)
  model$fit()
  assign(model$place, model)
} else {
  model <- get(spec$place)$clone(TRUE)
}
print(test_prices(model))

# Model specification
spec <- solarModel_spec$new()
spec$set_params(min_pos = 1, max_pos = 1, transform_delta = 0.01)
spec$set_mean.model(arOrder = 2, maOrder = 1)
spec$set_mixture.model(match.moments = TRUE)
spec$specification("Ravenna", min_date="2005-01-01", from = "2005-01-01", to = "2023-12-31", max_date = "2023-12-31",)
# Model fit
if (!exists(spec$place) | fit_again) {
  model <- solarModel$new(spec)
  model$fit()
  assign(model$place, model)
} else {
  model <- get(spec$place)$clone(TRUE)
}
print(test_prices(model))

# Model specification
spec <- solarModel_spec$new()
spec$set_params(min_pos = 1, max_pos = 1, transform_delta = 0.01)
spec$set_mean.model(arOrder = 2, maOrder = 1)
spec$set_mixture.model(match.moments = TRUE)
spec$specification("Rimini", min_date="2005-01-01", from = "2005-01-01", to = "2023-12-31", max_date = "2023-12-31",)
# Model fit
if (!exists(spec$place) | fit_again) {
  model <- solarModel$new(spec)
  model$fit()
  assign(model$place, model)
} else {
  model <- get(spec$place)$clone(TRUE)
}
print(test_prices(model))

# Model specification
spec <- solarModel_spec$new()
spec$set_params(min_pos = 1, max_pos = 1, transform_delta = 0.05)
spec$set_mean.model(arOrder = 1, maOrder = 1)
spec$set_variance.model(archOrder = 1, garchOrder = 1)
spec$set_mixture.model(match.moments = TRUE)
spec$specification("Amsterdam", min_date="2005-01-01", from = "2005-01-01", to = "2023-12-31", max_date = "2023-12-31",)
# Model fit
if (!exists(spec$place) | fit_again) {
  model <- solarModel$new(spec)
  model$fit()
  assign(model$place, model)
} else {
  model <- get(spec$place)$clone(TRUE)
}
print(test_prices(model))

# Model specification
spec <- solarModel_spec$new()
spec$set_params(min_pos = 1, max_pos = 1, transform_delta = 0.05)
spec$set_mean.model(arOrder = 2, maOrder = 1)
spec$set_variance.model(archOrder = 1, garchOrder = 1)
spec$set_mixture.model(match.moments = TRUE)
spec$specification("Berlino", min_date="2005-01-01", from = "2005-01-01", to = "2023-12-31", max_date = "2023-12-31",)
# Model fit
if (!exists(spec$place) | fit_again) {
  model <- solarModel$new(spec)
  model$fit()
  assign(model$place, model)
} else {
  model <- get(spec$place)$clone(TRUE)
}
print(test_prices(model))

# Model specification
spec <- solarModel_spec$new()
spec$set_clearsky(control_seasonalClearsky(order = 1, order_H0 = 1))
spec$set_params(min_pos = 1, max_pos = 1, transform_delta = 0.05)
spec$set_mean.model(arOrder = 1, maOrder = 1)
spec$set_variance.model(archOrder = 1, garchOrder = 1)
spec$set_mixture.model(match.moments = TRUE)
spec$specification("Oslo", min_date="2005-01-01", from = "2005-01-01", to = "2023-12-31", max_date = "2023-12-31",)
# Model fit
if (!exists(spec$place) | fit_again) {
  model <- solarModel$new(spec)
  model$fit()
  assign(model$place, model)
} else {
  model <- get(spec$place)$clone(TRUE)
}
print(test_prices(model))

# Model specification
spec <- solarModel_spec$new()
spec$set_params(min_pos = 1, max_pos = 3, transform_delta = 0.01)
spec$set_mean.model(arOrder = 1, maOrder = 1)
spec$set_mixture.model(match.moments = TRUE)
spec$specification("Parigi", min_date="2005-01-01", from = "2005-01-01", to = "2023-12-31", max_date = "2023-12-31",)
# Model fit
if (!exists(spec$place) | fit_again) {
  model <- solarModel$new(spec)
  model$fit()
  assign(model$place, model)
} else {
  model <- get(spec$place)$clone(TRUE)
}
print(test_prices(model))

