#' Control for `clearskyModel`
#'
#' Control parameters for a `clearskyModel` object.
#'
#' @param method character, method for clearsky estimate, can be `I` or `II`.
#' @param include.intercept logical. When `TRUE`, the default, the intercept will be included in the model.
#' @param order numeric, of sine and cosine elements.
#' @param period numeric, periodicity. The default is `2*base::pi/365`.
#' @param seed numeric, random seed for reproducibility. It is used to random impute the violations.
#' @param delta_init Value for delta init in the clearsky model.
#' @param tol integer, numeric tolerance for clearsky > GHI condition. Maximum number of violations admitted. Used in `clearskyModel_optimize()`.
#' @param lower numeric, lower bound for delta grid, see `clearskyModel_optimize()`.
#' @param upper numeric, upper bound for delta grid, see `clearskyModel_optimize()`.
#' @param by numeric, step for delta grid, see `clearskyModel_optimize()`.
#' @param quiet logical. When `FALSE`, the default, the functions displays warning or messages.
#'
#' @rdname clearskyModel_control
#' @export
clearskyModel_control <- function(method = "II", include.intercept = TRUE, order = 1, period = 365, seed = 1,
                                  delta_init = 1.1, tol = 30, lower = 0, upper = 1, by = 0.001, quiet = FALSE){
  structure(
    list(
      method = method,
      include.intercept = include.intercept,
      order = order,
      period = period,
      seed = seed,
      delta_init = delta_init,
      tol = tol,
      delta0 = c(lower = lower, upper = upper, by = by),
      quiet = quiet
    ),
    class = c("control", "list")
  )
}

#' Optimizer for Solar Clear sky
#'
#' Find the best parameter delta for fitting clear sky radiation.
#'
#' @param GHI vector of solar radiation
#' @param G0 vector of extraterrestrial radiation
#' @param control list of control parameters. See `clearskyModel_control()` for details.
#'
#' @return a numeric vector containing the fitted clear sky radiation.
#'
#' @name clearskyModel_optimize
#' @rdname clearskyModel_optimize
#' @export
clearskyModel_optimize <- function(GHI, G0, control = clearskyModel_control()){

  # Control parameters used
  tol = control$tol
  lower = control$delta0[1]
  upper = control$delta0[2]
  by = control$delta0[3]

  # Grid
  delta_grid <- seq(lower, upper, by = by)
  # Optimization
  opt <- data.frame(delta = delta_grid, loss = purrr::map_dbl(delta_grid, ~sum(.x*G0 - GHI < 0)))
  delta_opt <- opt[which(opt$loss <= tol)[1],]$delta
  # Fitted clear sky
  Ct <- delta_opt*G0
  # Output
  attr(Ct, "delta") <- delta_opt
  return(Ct)
}

#' Seasonal model for clear sky radiation
#'
#' Fit a seasonal model for clear sky radiation in a location.
#'
#' @param data dataset
#' @param seasonal_data dataset with two columns: `n` with the number of the day in 1:365
#' and `H0` with the extraterrestrial radiation.
#' @param control list of control parameters. See `clearskyModel_control()` for details.
#'
#' @rdname clearskyModel_fit
#' @name clearskyModel_fit
#' @export
clearskyModel_fit <- function(data, seasonal_data, control = clearskyModel_control()){

  # Control parameters
  method = control$method
  include.intercept = control$include.intercept
  order = control$order
  period = control$period
  seed = control$seed
  delta_init = control$delta_init
  quiet = control$quiet

  # Method I: Estimate only with Extraterrestrial radiation
  if (method == "I") {
    df_fit <- data
    formula <- paste0("GHI ~ H0", ifelse(include.intercept, " + 1", " - 1"))
    # Fit the parameters of the clear sky max model
    seasonal_model_Ct <- seasonalModel_fit(formula = formula, order = order, period = period, data = df_fit)
    # Fitted clear sky max
    data$Ct <- predict.lm(seasonal_model_Ct, newdata = df_fit)
    # Optimize the fit
    data$Ct <- clearskyModel_optimize(data$GHI, data$Ct*delta_init, control = control)
    # Update the parameters inside the `lm` object
    seasonal_model_Ct$coefficients_orig <- seasonal_model_Ct$coefficients
    seasonal_model_Ct$delta <- attr(data$Ct, "delta")
    seasonal_model_Ct$coefficients <- seasonal_model_Ct$coefficients*seasonal_model_Ct$delta*delta_init
    # Update the fitted clearsky
    data$Ct <- predict.lm(seasonal_model_Ct, newdata = df_fit)
    # Detect violations
    outliers <- data$GHI > data$Ct
    # Method II: Estimate with Extraterrestrial and clear sky radiation
  } else if (method == "II") {
    # Compute maximum clear sky for each day
    df_fit <- dplyr::group_by(data, Month, Day)
    df_fit <- dplyr::reframe(df_fit, Ct_max = max(clearsky), H0 = mean(H0))
    df_fit$n <- 1:nrow(df_fit)
    formula <- paste0("Ct_max ~ H0", ifelse(include.intercept, " + 1", " - 1"))
    # Fit the parameters of the clear sky max model
    seasonal_model_Ct <- seasonalModel_fit(formula = formula, order = order, period = period, data = df_fit)
    # Fitted clear sky max
    df_fit$Ct_max_fit <- predict.lm(seasonal_model_Ct, newdata = df_fit)
    # Optimize the fit
    df_fit <- dplyr::select(df_fit, Month, Day, Ct_max, Ct_max_fit)
    df_fit <- dplyr::left_join(data, df_fit, by = c("Month", "Day"))
    data$Ct <- clearskyModel_optimize(df_fit$GHI, df_fit$Ct_max_fit*delta_init, control = control)
    # Update the parameters inside the `lm` object
    seasonal_model_Ct$coefficients_orig <- seasonal_model_Ct$coefficients
    seasonal_model_Ct$delta <- attr(data$Ct, "delta")
    seasonal_model_Ct$coefficients <- seasonal_model_Ct$coefficients*seasonal_model_Ct$delta*delta_init
    # Detect violations
    outliers <- data$GHI > data$Ct
  }
  # Add period for predictions
  seasonal_model_Ct$period <- control$period
  # Add extraterrestrial radiation for predictions
  seasonal_model_Ct$H0 <- seasonal_data

  # Verbose message
  if (!quiet) {
    message("Outliers: ", sum(outliers), " (",
            format(sum(outliers)/nrow(data)*100, digits = 2), "%)",
            " - sum(GHI - Ct): ", sum((data$GHI - data$Ct)[outliers]))
  }

  # Compute the Risk Driver
  data$Xt <- 1 - data$GHI/data$Ct
  # Seed for reproducibility
  set.seed(seed)
  # Randomly impute the violations
  rand_probs <- runif(length(outliers), min = 0.001, max = sum(outliers)/nrow(data))
  syntetic_Xt <- quantile(data$Xt[!outliers], probs = rand_probs)
  # Add control object
  seasonal_model_Ct$control <- append(control, c(order = order))
  # Add outliers information
  seasonal_model_Ct$outliers <- list(index = which(outliers),
                                     n = sum(outliers),
                                     orig_Xt = data$Xt[outliers],
                                     syntetic_Xt = syntetic_Xt[outliers])
  class(seasonal_model_Ct) <- c("clearskyModel", class(seasonal_model_Ct)[-1])

  # Impute outliers
  data$Xt <- ifelse(outliers, syntetic_Xt, data$Xt)
  data$GHI <- ifelse(outliers, data$Ct*(1 - data$Xt), data$GHI)

  structure(
    list(
      data = data,
      seasonal_model_Ct = seasonal_model_Ct,
      control = control
    )
  )
}

#' Predict method
#'
#' Predict method for the class `cleaskyModel`.
#'
#' @param object an object of the class `cleaskyModel`
#' @param n number of day of the year.
#'
#' @rdname clearskyModel_predict
#' @name clearskyModel_predict
#' @export
clearskyModel_predict <- function(object, n = 1){
  df_n <- data.frame(n = n %% object$period)
  df_n <- dplyr::left_join(df_n, object$H0, by = "n")
  predict.lm(object, newdata = df_n)
}
