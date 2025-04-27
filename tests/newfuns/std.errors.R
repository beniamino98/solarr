
model_cal <- model$clone(TRUE)

params <- model_cal$GARCH$coefficients[c(2)]


loss <- function(params){
  model_cal$update(params)
  model_cal$filter()
  model_cal$fit_NM_model()
  model_cal$update_moments()
  model_cal$update_logLik()
  print(model_cal$loglik)
  -model_cal$loglik
}

opt <- optim(params, loss)
model_cal$GARCH
model_cal$loglik
model_cal$fit_NM_model()
model_cal$update_moments()
model_cal$update_logLik()

print(test_prices(model_cal))




solarModel_std.errors <- function(params, model) {
  # Log-likelihood Yt_tilde given the mixture parameters
  logLik <- function(params, model){
    # Clone the model
    model_cal <- model$clone(TRUE)
    # Update the parameters
    model_cal$update(params)
    # Filter the time series
    model_cal$filter()
    # Update the moments
    model_cal$update_moments()
    # Update log-likelihood
    model_cal$update_logLik()
    # Negative Log-likelihood
    loglik <- model_cal$loglik
    return(loglik)
  }
  # Numerically compute the std. errors
  H <- numDeriv::hessian(logLik, x = params, model = model)
  std.errors <- sqrt(diag(solve(-H)))
  names(std.errors) <- names(params)
  return(std.errors)
}
model_cal$GARCH

model_cal$GARCH$update_hessian(params, logLik, model = model_cal)
model_cal$GARCH$update_std.errors()
params <- model_cal$GARCH$coefficients


solarModel_std.errors(params, model)
