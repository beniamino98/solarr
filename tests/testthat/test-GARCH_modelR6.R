library(solarr)
library(rugarch)
# Reference solar model for time series
model <- solarModel$new(spec)
model$fit()
# Reference time series to fit
x <- model$data$eps_tilde

# Specification
GARCH_spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1), variance.targeting = 1),
                         mean.model = list(armaOrder = c(0,0), include.mean = FALSE))

# Initialize the model
GARCH <- GARCH_modelR6$new(GARCH_spec, x, sigma20 = 1)
# Fit the model
GARCH$fit()
# self <- GARCH$.__enclos_env__$self
# private <- GARCH$.__enclos_env__$private
GARCH

# Check computation of the hessian are corrects
prev_H <- GARCH$.__enclos_env__$private$..hessian
GARCH$update_hessian()
next_H <- GARCH$.__enclos_env__$private$..hessian
prev_H - next_H

# Check computation of the std.errors are corrects
prev_std.error <- GARCH$tidy$std.error
GARCH$update_std.errors()
next_std.error <- GARCH$tidy$std.error
next_std.error - prev_std.error

# Check outputs parameters
GARCH$omega
GARCH$alpha
GARCH$beta
GARCH$coefficients
# Check unconditional vol
GARCH$vol
# Check order
GARCH$order
# Tidy output
GARCH$tidy

# Check next step prediction
eps0 <- c(0.3)
sigma0 <- c(0.9)
sigmat <- sigma0
# One step ahead
sigmat <- c(sqrt(GARCH$omega + GARCH$alpha * eps0^2 + + GARCH$beta * sigma0^2), sigma0)
sigmat[1]
GARCH$next_step(eps0, sigma0)
# Two steps ahead
sigmat <- c(sqrt(GARCH$omega + GARCH$alpha * sigmat[1]^2 + + GARCH$beta * sigmat[1]^2), sigmat)
sigmat[1]
GARCH$next_step(eps0, sigma0, n.ahead = 2)
# Three steps ahead
sigmat <- c(sqrt(GARCH$omega + GARCH$alpha * sigmat[1]^2 + + GARCH$beta * sigmat[1]^2), sigmat)
sigmat[1]
GARCH$next_step(eps0, sigma0, n.ahead = 3)

# All sample prediction
GARCH$filter(x = x[1:2] )$sigma


