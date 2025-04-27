# ***********************************************************************************************
#                                  Test for class "ARMA_modelR6"
# ***********************************************************************************************
# Time series to use for fitting
x <- Bologna$data$Yt_tilde

# Define the model
arma <- ARMA_modelR6$new(maOrder = 2, arOrder = 1)
# Fit the model
arma$fit(x)

#private <- arma$.__enclos_env__$private
#self <- arma$.__enclos_env__$self


# Check residuals computation
# Fitted residuals (true)
eps_true <- arma$model$residuals

# Method `filter`
# Initial residuals
eps0 <- eps_true[1:arma$order[2]]
# Fitted residuals
eps_hat <- arma$filter(x, eps0)$residuals
# Average Percentage difference
mean(abs((eps_hat - eps_true) / eps_true))

eps_true[1:5]
eps_hat[1:5]

# Check autocorrelation
solarrfigs::fig_acf(eps_hat)
solarrfigs::fig_acf(eps_true)

ggplot()+
  geom_line(aes(1:length(eps_hat), eps_hat-eps_true))

Box.test(eps_hat)
Box.test(eps_true)


# Check next step function
x0 <- x[c(10)]
eps0 <- eps_hat[c(10, 9)]
x_t <- c(x0, eps0)
# Simple forecast
arma$next_step(x_t, eps = 0)[1,1]
sum(arma$phi * x0) + sum(arma$theta * eps0)
# Next step
arma$next_step(x_t, eps = eps_hat[11])[1,1]
sum(arma$phi * x0) + sum(arma$theta * eps0) + eps_hat[11]

# Simulations
nsim <- 50
eps_sim <- rnorm(nsim)
x_sim <- purrr::map_dbl(1:nsim, ~arma$next_step(x_t, n.ahead = .x, eps = eps_sim[1:.x])[1,1])

ggplot()+
  geom_line(aes(1:nsim, x_sim))

sum(arma$phi * x[c(10,9)] + arma$theta * e[c(10,9)])

x[c(11)] - e[11]

# Update the coefficients
new_coef <- arma$coefficients * 1.01
arma$update(new_coef)
arma$coefficients == new_coef



