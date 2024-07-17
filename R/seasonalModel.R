#' Fit a seasonal model
#'
#' Fit a seasonal model as a linear combination of sine and cosine functions.
#'
#' @param formula formula, an object of class `formula` (or one that can be coerced to that class).
#' It is a symbolic description of the model to be fitted and can be used to include or exclude the intercept or external regressors in `data`.
#' @param order numeric, of sine and cosine elements.
#' @param period numeric, periodicity. The default is `2*base::pi/365`.
#' @param data 	an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
#' If not found in data, the variables are taken from environment(formula), typically the environment from which `lm` is called.
#'
#' @rdname seasonalModel_fit
#' @name seasonalModel_fit
#' @export
seasonalModel_fit <- function(formula = "Yt ~ 1", order = 1, period = 365, data){

  df_fit <- data
  df_fit$n <- 1:nrow(df_fit)
  if (order > 0) {
    for (i in 1:order){
      formula <- paste0(formula, " + ", "cos((2*base::pi)/", eval(period), "*n*", i, ")", " + ",
                        "sin((2*base::pi)/", eval(period), "*n*", i, ")")
    }
  }
  # Fit seasonal model
  seasonal_model <- lm(as.formula(formula), data = df_fit)
  seasonal_model$period <- period
  seasonal_model$order <- order
  class(seasonal_model) <- c("seasonalModel", class(seasonal_model))
  return(seasonal_model)
}

#' Predict method for the class `seasonalModel`.
#'
#' @param object object of the class `seasonalModel`.
#' @param n integer, number of day of the year.
#'
#' @rdname seasonalModel_predict
#' @name seasonalModel_predict
#' @export
seasonalModel_predict <- function(object, n = 1){
  predict.lm(object, newdata = data.frame(n = n %% object$period))
}
