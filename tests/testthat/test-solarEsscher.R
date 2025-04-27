library(solarr)
# ********************* Settings *********************
place <- "Berlino"
n_key_points = 3
control_model = control_solarModel(clearsky = control_seasonalClearsky(),
                                   seasonal.mean = list(seasonalOrder = 1, include.H0 = FALSE, include.intercept = TRUE, monthly.mean = FALSE),
                                   seasonal.variance = list(seasonalOrder = 1, correction = FALSE,  monthly.mean = FALSE),
                                   mean.model = list(arOrder = 2, correction = TRUE, include.intercept = FALSE),
                                   mixture.model = list(abstol = 1e-20, maxit = 100),
                                   threshold = 0.01, outliers_quantile = 0.0001, quiet = FALSE)
control_options = control_solarOption()
# *****************************************************

# 1) Model specification and fit
spec <- solarModel_spec(place, control_model = control_model)
model <- solarModel$new(spec)
model$fit()

# Calibrate GM parameters around the seasonal mean
#model <- solarOption_calibrator(model)
# 2) Initialize a payoff object
payoffs <- solarOptionPayoffs(model, control_options = control_options)
# Add model price (P-measure)
payoffs$call$model$P <- solarOption_model(model, theta = 0, put = FALSE, control_options = control_options)
payoffs$put$model$P <- solarOption_model(model, theta = 0, put = TRUE, control_options = control_options)
# Add bootstrapped price (P-measure)
payoffs$call$model$boot <- solarOption_historical_bootstrap(model, put = FALSE,  control_options = control_options)
payoffs$put$model$boot <- solarOption_historical_bootstrap(model, put = TRUE, control_options = control_options)

# 3) Calibrate Esscher bounds and optimal parameters
esscher <- solarEsscher$new(n_key_points = n_key_points, control_options = control_options)
esscher$calibrate_bounds(model, payoffs)
# Add model price (Q-measure)
payoffs$call$model$Q <- solarOption_model(model, theta = esscher$bounds$Q, put = FALSE, control_options = control_options)
payoffs$put$model$Q  <- solarOption_model(model, theta = esscher$bounds$Q, put = TRUE, control_options = control_options)
# Add model price (Q-measure)
payoffs$call$model$Qdw <- solarOption_model(model, theta = esscher$bounds$Qdw, put = FALSE, control_options = control_options)
payoffs$put$model$Qdw  <- solarOption_model(model, theta = esscher$bounds$Qdw,  put = TRUE, control_options = control_options)
# Add model price (Q-measure)
payoffs$call$model$Qup <- solarOption_model(model, theta = esscher$bounds$Qup, put = FALSE, control_options = control_options)
payoffs$put$model$Qup  <- solarOption_model(model, theta = esscher$bounds$Qup,  put = TRUE, control_options = control_options)

# 4) Parametrize Esscher Theta
pay <- payoffs$put$model
# Benchmark price with r = 0
benchmark_price <- pay$P$payoff_year$premium
# Benchmark for best case prices with r > 0
lower <- pay$Qdw$payoff_year$premium
# Benchmark for worste case prices with r < 0
upper <- pay$Qup$payoff_year$premium
# Create the grid of esscher parameters for Yt
esscher$create_grid(model, benchmark_price, lower_price = lower, upper_price = upper)
