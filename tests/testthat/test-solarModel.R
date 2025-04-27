library(solarr)
library(tidyverse)
# ************************************************
#                    Inputs
# ************************************************
place <- "Bologna"
target <- "GHI"
# ************************************************
# Control list
control <- control_solarModel()
# Model specification
spec <- solarModel_spec(place = place, target = target, min_date = "2005-01-01", to = "2023-01-01", control_model = control)
# Initialize the model
model <- solarModel$new(spec)
# Model fit
model$fit()
# ************************************************
test_that('Test: Check model$target, model$place, model$coords...', {
  # Check response
  expect_true(model$target == target)
  expect_true(model$place == place)
  expect_true(length(model$coords) == 3)
})
# ************************************************
test_that('Test: Check methods $update and $filter...', {
  model_clone <- model$clone(deep = TRUE)
  params <- model_clone$coefficients
  params$NM_mu_up$mu_up_2 <- -0.58

  model_clone$update(params)
  model_clone$filter()

  expect_true(model_clone$coefficients$NM_mu_up$mu_up_2 == params$NM_mu_up$mu_up_2)
  expect_true(model$coefficients$NM_mu_up$mu_up_2 != params$NM_mu_up$mu_up_2)
  expect_true(sum(model_clone$NM_model$fitted$B != model_clone$data$B) == 0)
})
# ************************************************

# Locations info
model$place
model$coords
model$dates
model$target
model$location

# List of parameters
model$coefficients
# Combination for interpolated models
model$combinations
# Moments
model$moments$conditional
model$moments$unconditional

# Models
model$seasonal_model_Ct
model$seasonal_model_Yt
model$AR_model_Yt
model$seasonal_variance
model$GARCH
# Data used for fitting
model$NM_model$data
# Parameters
model$NM_model$coefficients
# Theoric momenta
model$NM_model$moments
# Fitted Bernoulli
model$NM_model$fitted
# Total log-likelihood
model$NM_model$loglik
# List of gaussian mixture models
model$NM_model$model

# Data
model$data
model$seasonal_data
model$monthly_data

# Solar Transform
model$transform
