library(solarr)
library(dplyr)
library(testthat)
# ***********************************************************************************************
# ***********************************************************************************************
#
#                                 Tests: solarModel class
#
# ***********************************************************************************************
# ***********************************************************************************************
# Model's specification
spec <- solarModel_spec$new()
spec$specification("Bologna")
spec
# ***********************************************************************************************
#                                      solarModel$fit
# ***********************************************************************************************
# Test model's fit
test_that('solarModel$fit: check fit works ...', {
  model <- solarModel$new(spec)
  model$fit()
  # private <- model$.__enclos_env__$private
  # self <- model$.__enclos_env__$self
  expect_true(!is.na(model$loglik))
})
# ***********************************************************************************************
#                                      solarModel$filter
# ***********************************************************************************************
# Test parameters update
test_that('solarModel$update: check params storage ...', {
  spec_train <- model$spec
  model_2 <- solarModel$new(spec_train)
  model_2$filter()
  model_2$update_classification()
  # Update conditional moments
  model_2$update_moments()
  # Update log-likelihood
  model_2$update_logLik()
  # Tests
  expect_true(model_2$loglik == model$loglik)
  expect_true(sum(abs(model$data$GHI-model_2$data$GHI)) == 0)
  expect_true(sum(abs(model$data$Yt-model_2$data$Yt)) == 0)
  expect_true(sum(abs(model$data$u_tilde-model_2$data$u_tilde)) == 0)
})
# ***********************************************************************************************
#                                      solarModel_QMLE
# ***********************************************************************************************
# Test parameters update
test_that('solarModel_QMLE: check params storage ...', {
  model_upd <- solarModel_QMLE(model_2)
  # Tests
  expect_true(model_2$loglik != model_upd$loglik)
  expect_true(sum(abs(model_2$data$GHI-model_upd$data$GHI)) == 0)
  expect_true(sum(abs(model_2$data$Yt-model_upd$data$Yt)) == 0)
  expect_true(sum(abs(model_2$data$u_tilde-model_upd$data$u_tilde)) != 0)
})

# Specification
model$spec

# Data
model$data
model$seasonal_data
model$monthly_data
model$loglik

# List of parameters
model$coefficients
# Moments
model$moments$conditional

# Models
model$spec$seasonal_model_Ct
model$spec$transform
model$spec$seasonal.mean
model$spec$mean.model
model$spec$seasonal.variance
model$spec$variance.model
model$spec$mixture.model
