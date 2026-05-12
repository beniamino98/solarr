library(solarr)
library(dplyr)
library(testthat)
# ***********************************************************************************************
# ***********************************************************************************************
#
#                                 Tests: ARMA_modelR6 class
#
# ***********************************************************************************************
# ***********************************************************************************************
# Test time series
set.seed(1)
y <- rnorm(100)
# ARMA model
ARMA <- ARMA_modelR6$new(1, 0, include.intercept = TRUE)
#private <- ARMA$.__enclos_env__$private
#self <- ARMA$.__enclos_env__$self

# ***********************************************************************************************
#                                      ARMA$fit
# ***********************************************************************************************
# ARMA model
ARMA <- ARMA_modelR6$new(1, 0, include.intercept = TRUE)
ARMA$fit(y)
test_that('ARMA$fit: check intercept is not zero when include.intercept = TRUE ...', {
  expect_true(sum(ARMA$coefficients != 0) == 2)
})

test_that('ARMA$fit: check consistency with arima residuals when include.intercept = TRUE...', {
  # Fitted residuals from arima
  eps_arima <- ARMA$model$residuals
  # Fitted residuals with ARMA$filter
  eps_hat <- y - ARMA$filter(y)
  expect_true(sum(abs(eps_hat - eps_arima)) < 0.0001)
})

# ***********************************************************************************************
# ARMA model
ARMA <- ARMA_modelR6$new(1, 0, include.intercept = FALSE)
ARMA$fit(y)
test_that('ARMA$fit: check intercept is zero when include.intercept = FALSE ...', {
  expect_true(ARMA$intercept == 0)
})
test_that('ARMA$fit: check consistency with arima residuals when include.intercept = FALSE...', {
  # Fitted residuals from arima
  eps_arima <- ARMA$model$residuals
  # Fitted residuals with ARMA$filter
  eps_hat <- y - ARMA$filter(y)
  expect_true(sum(abs(eps_hat - eps_arima)) < 0.001)
})

# ***********************************************************************************************
test_that('ARMA$fit: check theta is zero when maOrder = 0 ...', {
  expect_true(ARMA$theta == 0)
})
# ***********************************************************************************************
#                                      ARMA$update
# ***********************************************************************************************
test_that('ARMA$update: check theta_1 remains zero after update...', {
  # ARMA$update: check theta_1 remains zero
  x1 <- ARMA$coefficients
  ARMA$update(c(phi_1 = 0.3, theta_1 = 0.1, intercept = 0.4))
  x2 <- ARMA$coefficients
  expect_true(x1["theta_1"] == x2["theta_1"])
})
# ***********************************************************************************************
test_that('ARMA$update: check extra parameters do not update or gives errors...', {
  x1 <- ARMA$coefficients
  ARMA$update(c(phi_2 = 0.3, theta_4 = 0.1, intercepts = 0.4))
  x2 <- ARMA$coefficients
  expect_true(sum(x1 == x2) == 3)
})
# ***********************************************************************************************
test_that('ARMA$update_sigma2: check extra parameters do not update or gives errors...', {
  x1 <- ARMA$sigma2
  ARMA$update_sigma2(0.94)
  x2 <- ARMA$sigma2
  expect_true(x1 != x2)
})
# ***********************************************************************************************
#                                      ARMA$filter
# ***********************************************************************************************
test_that('ARMA$filter: check length filtered series is equal to length of input series', {
  expect_true(length(ARMA$filter(y)) == length(y))
})
# ***********************************************************************************************
#                                      ARMA$next_step
# ***********************************************************************************************
test_that('ARMA$filter: consistency with manual iteration', {
  ## Next step
  x0 <- y[c(10)]
  set.seed(1)
  eps0 <- rnorm(1)
  x_t <- list(`t+0` = c(x0))
  # Simple forecast
  x_t[[2]] <- ARMA$next_step(x_t[[1]], eps = 0)
  # Check
  expect_true(x_t[[2]][1,1] == ARMA$intercept + sum(ARMA$phi * x0) + sum(ARMA$theta * eps0))
  # Two steps ahead
  x_t[[3]] <- ARMA$next_step(x_t[[2]], eps = 0)
  # Check
  expect_true(ARMA$intercept + sum(ARMA$phi * x_t[[2]][1,1]) + sum(ARMA$theta * 0)  == x_t[[3]][1,1])
})
# ***********************************************************************************************
#                                      ARMA$expectation
# ***********************************************************************************************
test_that('ARMA$expectation:', {
  X0 = c(0.2)
  ## Expectation
  expect_true(length(ARMA$expectation(1, X0)) == 1)
  expect_true(length(ARMA$expectation(30, X0)) == 30)
  expect_true(length(ARMA$expectation(100, X0)) == 100)
})
# ***********************************************************************************************
#                                      ARMA$variance
# ***********************************************************************************************
test_that('ARMA$variance:', {
  sigma2 = c(0.9)
  ## Expectation
  expect_true(length(ARMA$variance(1, sigma2)) == 1)
  expect_true(length(ARMA$variance(30, sigma2)) == 30)
  expect_true(length(ARMA$variance(100, sigma2)) == 100)
})
# ***********************************************************************************************
#                                      ARMA Slots
# ***********************************************************************************************
test_that('ARMA$A: check dimensions companion matrix...', {
  expect_true(nrow(ARMA$A) == sum(ARMA$arOrder+ARMA$maOrder))
})
test_that('ARMA$b: check dimensions vector b...', {
  expect_true(length(ARMA$b) == sum(ARMA$arOrder+ARMA$maOrder))
  expect_true(ifelse(ARMA$arOrder != 0 & ARMA$maOrder != 0, ARMA$b[ARMA$arOrder+1] == 1, TRUE))
})

# Slots
## Orders
ARMA$arOrder
ARMA$maOrder
ARMA$order
## Parameters
ARMA$intercept
ARMA$phi
ARMA$theta
ARMA$coefficients
## Parameters std. errors
ARMA$std.errors
# Tidy ouptut
ARMA$tidy
## Matrices
ARMA$A
ARMA$b
## Estimated variance of the residuals
ARMA$sigma2





