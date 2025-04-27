library(solarr)
library(tidyverse)

means = c(0,0.5,2)
sd = rep(1, 3)
p = c(0.2, 0.3, 0.5)
# Simulated mixture
x <- rmixnorm(500, means, sd, p)

# Define the GM model
gm <- gaussianMixture$new(components=3)
# Fit the model
gm$fit(x$X)
# Fitted moments
moments_fitted <- gm$moments
fitted <- gm$fitted
loglik_fitted <- gm$loglik
params_fitted <- gm$parameters

test_that('Test: computation of logLik() function ...', {
  # Check response
  expect_true(gm$logLik(x$X) == loglik_fitted)
})

# Update the means
new_means <- c(0, -1, 3)
gm$update(means = new_means)
test_that('Test: update means parameters ...', {
  # Check response
  expect_true(sum(gm$parameters$p  == params_fitted$p) == length(p))
  expect_true(sum(gm$parameters$sd  == params_fitted$sd) == length(sd))
  expect_true(sum(gm$parameters$means  == params_fitted$means) == 0)
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

gm$use_empiric_parameters()
empiric_fitted <- gm$parameters
gm$use_empiric_parameters()
test_that('Test: using empirical parameters ...', {
  # Check response
  expect_true(sum(gm$parameters$means != params_fitted$means) == 0)
  expect_true(sum(gm$parameters$sd != params_fitted$sd) == 0)
  expect_true(sum(gm$parameters$p != params_fitted$p) == 0)
})

