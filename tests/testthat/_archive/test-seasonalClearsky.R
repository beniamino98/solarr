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
# Reference latitude
lat <- spec$coords$lat
# Test time series
set.seed(1)
# True coefs
true_coefs <- c(1, 1.3, 0.2, 0.4)
# Simulated time series
n <- 3000
ssf <- seasonalSolarFunctions$new()
H0 <- ssf$Hon(1:n, lat)
eps <- rnorm(n, sd = 1.2)
y <- true_coefs[1] + H0 * true_coefs[2] + true_coefs[3] * sin(2*base::pi/365*c(1:n)) + true_coefs[4] * cos(2*base::pi/365*c(1:n)) + eps
plot(y)
# ***********************************************************************************************
#                                         Method: fit
# ***********************************************************************************************
# Fit a clear sky model
sc <- seasonalClearsky$new()
self <- sc$.__enclos_env__$self
super <- sc$.__enclos_env__$super
private <- sc$.__enclos_env__$private

# Arguments
x <- y
dates <- seq.Date(from =  as.Date("2010-01-01"), length.out = n, by = 1)
lat <- 44.4949
clearsky <- y*1.3
sc$fit(x, dates, lat, clearsky)

ggplot()+
  geom_line(aes(1:n, y))+
  geom_line(aes(1:n, sc$predict(1:n)), color = "blue")

sc$coefficients
sc$control

# ***********************************************************************************************
#                                         Method: differential
# ***********************************************************************************************
n0 <- 5
dt <- 0.0001
test_that('seasonalClearsky: check consistency between finite difference and differential...', {
  # Differential
  dCt <- sc$differential(n0)
  # Finite difference
  delta_Ct <- (sc$predict(n = n0+dt)-sc$predict(n = n0))/dt
  expect_true(sum(abs(dCt - delta_Ct)) < 0.0001)
})
# ***********************************************************************************************
#                                         Method: update
# ***********************************************************************************************
test_that('seasonalClearsky: update coefficients...', {
  par_old <- sc$coefficients["delta_0"]
  sc$update(c(delta_0 = -1.3))
  par_new <-  sc$coefficients["delta_0"]
  expect_true(par_new == -1.3)
  expect_true(sum(is.na(sc$std.errors)) == 1)
  expect_true(is.na(sc$std.errors["delta_0"]))
})
# ***********************************************************************************************
#                                         Method: update_std.errors
# ***********************************************************************************************
test_that('seasonalClearsky: update std.errors', {
  par_old <- sc$std.errors["delta_0"]
  sc$update_std.errors(c(delta_0 = 0.3))
  par_new <-  sc$std.errors["delta_0"]
  expect_true(par_new == 0.3)
  expect_true(sum(is.na(sc$std.errors)) == 0)
})

# ***********************************************************************************************
#                                         Method: update_coefs_names
# ***********************************************************************************************
test_that('seasonalClearsky: update coefficients names...', {
  par_old <- sc$coefs_names
  sc$update_coefs_names(c(delta_0 = "delta_intercept"))
  par_new <-  names(sc$coefficients)[1]
  expect_true(par_new == "delta_intercept")
  sc$update_coefs_names(c(delta_intercept = "delta_0"))
})
