#' Seasonal model for solar radiation radiation
#'
#' Fit a seasonal model for solar radiation
#'
#' @param spec an object with class `solarModelSpec`. See the function \code{\link{solarModel_spec}} for details.
#' @examples
#' library(ggplot2)
#' # Seasonal model for GHI
#' spec <- solarModel_spec("Oslo", target = "GHI")
#' model <- seasonalRadiation(spec)
#' spec$data$GHI_bar <- model$predict(spec$data$n)
#' ggplot(spec$data)+
#'  geom_line(aes(n, GHI))+
#'  geom_line(aes(n, GHI_bar), color = "blue")
#'
#' # Seasonal model for clear sky
#' spec <- solarModel_spec("Oslo", target = "clearsky")
#' model <- seasonalRadiation(spec)
#' spec$data$Ct_bar <- model$predict(spec$data$n)
#' ggplot(spec$data)+
#'  geom_line(aes(n, clearsky))+
#'  geom_line(aes(n, Ct_bar), color = "blue")
#' @rdname seasonalRadiation
#' @name seasonalRadiation
#' @export
seasonalRadiation <- function(spec){
  # Control parameters
  control <- spec$control$seasonal.mean
  target <- spec$target
  # Complete data
  data <- spec$data[spec$data$isTrain & spec$data$weights != 0,]
  # Custom formula
  formula_ <- ifelse(control$include.H0, paste0(target, " ~ 1 + H0"), paste0(target, " ~ 1"))
  # Seasonal model
  seasonal_model <- seasonalModel$new(order = control$seasonalOrder)
  # Seasonal model
  seasonal_model$fit(formula_, data = data)
  # Check that the predictions are positive for all n's
  negative_condition <- any(seasonal_model$predict(1:366) < 0)
  if (negative_condition) {
    if (!spec$control$quiet) message("Seasonal radiation is negarive! Optimizing the parameters...")
    # Initial parameters
    params <- seasonal_model$coefficients
    if (control$include.intercept) {
      params[1] <- params[1]*1.5
    }
    # Loss function to ensure positivity
    loss_function <- function(params, model, data, target){
      model$update(params)
      pred <- model$predict(data$n)
      if (any(pred < 0)) {
        return(NA)
      }
      sum((pred - data[[target]])^2)
    }
    # Optimal parameters
    opt <- optim(params, loss_function, model = seasonal_model, data = data, target = target)
    # Update parameters
    seasonal_model$update(opt$par)
  }
  # Update class
  class(seasonal_model) <- c("seasonalRadiation", class(seasonal_model))
  return(seasonal_model)
}

