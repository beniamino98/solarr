library(solarr)
place <- "Bologna"
target <- "GHI"

spec <- solarModel_spec(place = place, target = target)
# ------------  Test for api = "spot" ------------
model <- solarModel$new(spec)
model$fit()

test_that('Test: Check model$target, model$place, model$coords...', {
  # Check response
  expect_true(model$target == target)
  expect_true(model$place == place)
  expect_true(length(model$coords) == 3)
})


# Locations info
model$place
model$coords
model$dates
model$target
model$location

# List of parameters
model$parameters
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
model$NM_model$parameters
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
