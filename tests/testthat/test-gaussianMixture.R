library(solarr)
library(tidyverse)
library(testthat)
# ***********************************************************************************************
# ***********************************************************************************************
#
#                                 Tests: ARMA_modelR6 class
#
# ***********************************************************************************************
# ***********************************************************************************************
# True parameters
par_true <- list(
  means = c(0, 2),
  sd = c(0.3, 1),
  p = c(0.3, 0.7)
)
# Simulate a mixture time series
x <- rmixnorm(500, par_true$means, par_true$sd, par_true$p)$X


# Initialize a model
gm <- gaussianMixture$new(components=2, method = "mclust")
gm$fit(x)
gm
# Initialize a model
gm <- gaussianMixture$new(components=2, method = "mixtools")
gm$fit(x)
gm

gm$Hessian()
gm$std.errors

# Fitted moments
gm$moments
gm$fitted
gm$responsabilities
gm$loglik

# Fitted moments
moments_fitted <- gm$moments
fitted <- gm$fitted
loglik_fitted <- gm$loglik
params_fitted <- gm$coefficients

test_that('Test: computation of logLik() function ...', {
  # Check response
  expect_true(gm$logLik(x) == loglik_fitted)
})

# Update the means
new_means <- setNames(c(0, -1, 3), names(gm$means))
gm$update(means = new_means)
test_that('Test: update means parameters ...', {
  # Check response
  expect_true(sum(gm$coefficients$p  == params_fitted$p) == gm$components)
  expect_true(sum(gm$coefficients$sd  == params_fitted$sd) == gm$components)
  expect_true(sum(gm$coefficients$means  == params_fitted$means) == 0)
})

test_that('Test: check that classification do not change after updating the means parameters ...', {
  # Check response
  expect_true(sum(fitted$classification != gm$fitted$classification) == 0)
})

# update the likelihood
test_that('Test: check loglik do not change after updating the means parameters ...', {
  # Check response
  expect_true(loglik_fitted == gm$loglik)
})
gm$update_logLik()
test_that('Test: check loglik change after updating the loglik...', {
  # Check response
  expect_true(loglik_fitted != gm$loglik)
})

gm$filter()
test_that('Test: check classification change after filtering...', {
  # Check response
  expect_true(sum(fitted$classification != gm$fitted$classification) > 0)
})

# Come back to original means
gm$update(means = params_fitted$means)
# Update the likelihood
gm$update_logLik()
test_that('Test: check loglik is equal to initial loglik after coming back to orinal means parameters ...', {
  # Check response
  expect_true(loglik_fitted == gm$loglik)
})

gm$filter()
test_that('Test: check classification is equal to initial classification after coming back to orinal means parameters...', {
  # Check response
  expect_true(sum(fitted$classification != gm$fitted$classification) == 0)
})

gm_new <- gaussianMixture$new()
# Update parameters
gm_new$update(gm$means, gm$sd, gm$p)
gm_new$std.errors
gm_new$Hessian()
gm_new$std.errors

# To update errors a time series is required
gm_new$set_time_series(x)
# Then Hessian
gm_new$Hessian()
gm_new$std.errors
# And filter for fitted data and responsabilities
gm_new$filter()
gm_new$fitted
gm_new$responsabilities


# Checkc
sum(gm$mean == gm_new$mean) == gm$components
sum(gm$sd == gm_new$sd) == gm$components
sum(gm$p == gm_new$p) == gm$components

sum(gm$fitted$uncertanty != gm_new$fitted$uncertanty) == 0
sum(gm$fitted$B1 != gm_new$fitted$B1) == 0

sum(gm$std.errors$means == gm_new$std.errors$means) == gm$components
sum(gm$std.errors$sd == gm_new$std.errors$sd) == gm$components
sum(gm$std.errors$p == gm_new$std.errors$p) == gm$components





