library(solarr)
library(dplyr)
library(testthat)
# ***********************************************************************************************
# ***********************************************************************************************
#
#                                 Tests: seasonalModel class
#
# ***********************************************************************************************
# ***********************************************************************************************
# Test time series
set.seed(1)
# True coefs
true_coefs <- c(1, 0.2, 0.4)
# Simulated time series
n <- 1000
eps <- rnorm(n, sd = 0.1)
y <- true_coefs[1] + true_coefs[2] * sin(2*base::pi/365*c(1:n)) + true_coefs[3] * cos(2*base::pi/365*c(1:n)) + eps
plot(y)
# ***********************************************************************************************
#                                         Method: fit
# ***********************************************************************************************
# Fit a clear sky model
sm <- seasonalModel$new("Yt ~ 1")
#self <- sm$.__enclos_env__$self
#private <- sm$.__enclos_env__$private

# Arguments
data <- data.frame(Yt = y, n = 1:n)
sm$fit(data)

ggplot()+
  geom_line(aes(1:n, y))+
  geom_line(aes(1:n, sm$predict(1:n)), color = "red")

sm$coefficients

# ***********************************************************************************************
#                                         Method: differential
# ***********************************************************************************************
n0 <- 5
dt <- 0.0001
test_that('seasonalModel: check consistency between finite difference and differential...', {
  # Differential
  dXt <- sm$differential(n0)
  # Finite difference
  delta_Xt <- (sm$predict(n = n0+dt) - sm$predict(n = n0))/dt
  expect_true(sum(abs(dXt - delta_Xt)) < 0.0001)
})
# ***********************************************************************************************
#                                         Method: update
# ***********************************************************************************************
test_that('seasonalModel: update coefficients...', {
  par_old <- sm$coefficients["intercept"]
  sm$update(c(intercept = -1.3))
  par_new <-  sm$coefficients["intercept"]
  expect_true(par_new == -1.3)
  expect_true(sum(is.na(sm$std.errors)) == 1)
  expect_true(is.na(sm$std.errors["intercept"]))
})
# ***********************************************************************************************
#                                         Method: update_std.errors
# ***********************************************************************************************
test_that('seasonalModel: update std.errors', {
  par_old <- sm$std.errors["intercept"]
  sm$update_std.errors(c(intercept = 0.3))
  par_new <-  sm$std.errors["intercept"]
  expect_true(par_new == 0.3)
  expect_true(sum(is.na(sm$std.errors)) == 0)
})

# ***********************************************************************************************
#                                         Method: update_coefs_names
# ***********************************************************************************************
test_that('seasonalModel: update coefficients names...', {
  par_old <- sm$coefs_names
  sm$update_coefs_names(c(intercept = "delta_intercept"))
  par_new <-  names(sm$coefficients)[1]
  expect_true(par_new == "delta_intercept")
  sm$update_coefs_names(c(delta_intercept = "intercept"))
})























# ***********************************************************************************************
# New seasonal model
sm <- seasonalModel$new("Yt ~ 1")
self <- sm$.__enclos_env__$self
private <- sm$.__enclos_env__$private

# Fit the model
data <- data.frame(Yt = y, n = 1:n)

sm$fit(data)
sm

sm$coefficients
sm$dcoefficients
# ***********************************************************************************************
# Update coefficients names
coefs_names <- c(intercept = "a_0")
sm$update_coefs_names(coefs_names)
sm$coefficients
sm$std.errors
# ***********************************************************************************************
# Update coefficients values
coefficients <- c(a_0 = 0.3)
sm$update(coefficients)
sm$coefficients
sm$dcoefficients

# No update if names do not matches
coefficients <- c(a_45 = 0.3)
sm$update(coefficients)
sm$coefficients

# ***********************************************************************************************
# Check forecasts
sm$predict(3)

# Differential
dt <- 0.00001
sm$differential(3)
(sm$predict(3+dt)-sm$predict(3))/dt


