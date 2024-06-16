#' Seasonal Model
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param order number of sine/cosine expansions.
#' @param period periodicity period 2pi/365.
#' @param data 	an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
#' If not found in data, the variables are taken from environment(formula), typically the environment from which `lm` is called.
#'
#' @rdname seasonalModel
#' @name seasonalModel
#' @export
seasonalModel <- function(formula = "Yt ~ 1", order = 1, period = 365, data){

  #' @examples
  #' formula = "GHI ~ 1"
  #' order = 1
  #' period = 365
  #' data = model$data
  df_fit <- data
  df_fit$n <- 1:nrow(df_fit)
  if (order > 0) {
    for (i in 1:order){
      formula <- paste0(formula, " + ", "cos((2*base::pi)/", eval(period), "*n*", i, ")", " + ",
                        "sin((2*base::pi)/", eval(period), "*n*", i, ")")
    }
  }
  # Fit clearsky max model
  seasonal_model <- lm(as.formula(formula), data = df_fit)
  seasonal_model$period <- period
  seasonal_model$order <- order
  seasonal_model$n <- 1:period
  class(seasonal_model) <- c("seasonalModel", class(seasonal_model))
  return(seasonal_model)
}

#' Optimizer for Solar Clear sky
#'
#' Find the best parameter delta for fitting clear sky radiation.
#'
#' @param GHI vector of solar radiation
#' @param G0 vector of extraterrestrial radiation
#'
#' @return a numeric vector containing the fitted clear sky radiation.
#'
#' @name solar_clearsky_optimizer
#' @rdname solar_clearsky_optimizer
#' @export
solar_clearsky_optimizer <- function(GHI, G0, control = control.seasonalModel()){
  delta_grid <- seq(control$delta0[1], control$delta0[2], by = control$delta0[3])
  opt <- data.frame(delta = delta_grid, loss = purrr::map_dbl(delta_grid, ~sum(.x*G0 - GHI < 0)))
  delta_opt <- opt[which(opt$loss <= control$tol)[1],]$delta
  Ct <- delta_opt*G0
  attr(Ct, "delta") <- delta_opt
  return(Ct)
}


#' Seasonal Model Clear sky Radiation
#'
#' @rdname clearsky.seasonalModel
#' @name clearsky.seasonalModel
#' @export
clearsky.seasonalModel <- function(object, control = control.seasonalModel()){

  #' @examples
  #' object <- CAMS("Amsterdam")
  #' control = control.seasonalModel()

  # Control parameters
  delta_init = control$delta_init
  quiet = control$quiet
  seed = control$seed
  method = control$method
  order = control$order
  include.intercept <- control$include.intercept
  # Dataset
  data <- object$data
  # Method I: Estimate only with Extraterrestrial radiation
  if (method == "I") {
    df_fit <- data
    formula <- paste0("GHI ~ H0", ifelse(include.intercept, " + 1", " - 1"))
    # Fit clearsky max model
    seasonal_model_Ct <- seasonalModel(formula = formula, order = order, period = 365, data = df_fit)
    # Fitted seasonal mean
    data$Ct <- predict.lm(seasonal_model_Ct, newdata = df_fit)
    data$Ct <- solar_clearsky_optimizer(data$GHI, data$Ct*delta_init, control = control)
    # Update parameter
    seasonal_model_Ct$coefficients_orig <- seasonal_model_Ct$coefficients
    seasonal_model_Ct$delta <- attr(data$Ct, "delta")
    seasonal_model_Ct$coefficients <- seasonal_model_Ct$coefficients*seasonal_model_Ct$delta*delta_init
    # Fitted clearsky
    data$Ct <- predict.lm(seasonal_model_Ct, newdata = df_fit)
    # Detect violations
    outliers <- data$GHI > data$Ct
    # Compute the Risk Driver
    data$Xt <- 1 - data$GHI/data$Ct
    # Method II: Estimate with Extraterrestrial and clearsky radiation
  } else if (method == "II") {
    # Compute maximum clearsky for each day
    df_fit <- dplyr::group_by(data, Month, Day)
    df_fit <- dplyr::reframe(df_fit, Ct_max = max(clearsky), H0 = mean(H0))
    df_fit$n <- 1:nrow(df_fit)
    formula <- paste0("Ct_max ~ H0", ifelse(include.intercept, " + 1", " - 1"))
    # Fit clearsky max model
    seasonal_model_Ct <- seasonalModel(formula = formula, order = order, period = 365, data = df_fit)
    # Predict clearsky max
    df_fit$Ct_max_fit <- predict.lm(seasonal_model_Ct, newdata = df_fit)
    df_fit <- dplyr::select(df_fit, Month, Day, Ct_max, Ct_max_fit)
    df_fit <- dplyr::left_join(data, df_fit, by = c("Month", "Day"))
    data$Ct <- solar_clearsky_optimizer(df_fit$GHI, df_fit$Ct_max_fit*delta_init, control = control)
    # Update parameters
    seasonal_model_Ct$coefficients_orig <- seasonal_model_Ct$coefficients
    seasonal_model_Ct$delta <- attr(data$Ct, "delta")
    seasonal_model_Ct$coefficients <- seasonal_model_Ct$coefficients*seasonal_model_Ct$delta*delta_init
    # Detect violations
    outliers <- data$GHI > data$Ct
    # Compute the Risk Driver
    data$Xt <- 1 - data$GHI/data$Ct
  }
  if (!quiet) {
    message("Outliers: ", sum(outliers), " (",
            format(sum(outliers)/nrow(data)*100, digits = 2), "%)",
            " - sum(GHI - Ct): ", sum((data$GHI - data$Ct)[outliers]))
  }
  # Seed for imputed violations
  set.seed(seed)
  syntetic_Xt <- quantile(data$Xt[!outliers], probs = runif(length(outliers), min = 0.001, max = sum(outliers)/nrow(data)))
  # Add extraterrestrial radiation for forecasting
  seasonal_model_Ct$H0 <- dplyr::select(object$seasonal_data, n, H0)
  # Add control
  seasonal_model_Ct$control <- append(control, c(order = order))
  # Add outliers informations
  seasonal_model_Ct$outliers <- list(index = which(outliers),
                                     n = sum(outliers),
                                     orig_Xt = data$Xt[outliers],
                                     syntetic_Xt = syntetic_Xt[outliers])
  # Add model to object
  object$seasonal_model_Ct <- seasonal_model_Ct
  # Impute outliers
  data$Xt <- ifelse(outliers, syntetic_Xt, data$Xt)
  data$GHI <- ifelse(outliers, data$Ct*(1 - data$Xt), data$GHI)
  object$data <- data
  class(seasonal_model_Ct) <- c("seasonalModel", class(seasonal_model_Ct)[-1])
  return(object)
}

#' predict method for seasonalModel
#'
#' @rdname predict.seasonalModel
#' @name predict.seasonalModel
#' @export
predict.seasonalModel <- function(object, n = 1){
  if (is.null(object$H0)){
    predict.lm(object, newdata = data.frame(n = n))
  } else {
    df_n <- dplyr::left_join(data.frame(n = n), object$H0, by = "n")
    predict.lm(object, newdata = df_n)
  }
}




