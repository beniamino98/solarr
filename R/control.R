#' Control function
#'
#' @export
control <- function(object){
  UseMethod("control", object)
}

#' seasonalModel control parameters
#'
#' @export
control.seasonalModel <- function(object, method = "II", include.intercept = TRUE, order = 1, seed = 1,
                                  delta_init = 1.1, tol = 30, lower = 0, upper = 1, by = 0.001, quiet = FALSE){
  if (missing(object)) {
    structure(
      list(
        method = method,
        include.intercept = include.intercept,
        order = order,
        seed = seed,
        delta_init = delta_init,
        tol = tol,
        delta0 = c(lower = lower, upper = upper, by = by),
        quiet = quiet
      ),
      class = c("control", "list")
    )
  } else {
    object$control
  }
}

#' solarModel control parameters
#'
#' @param loss type of loss function for mixture model, `ml` stands for maximum-likelihood, while
#' `kl` for KL-distance.
#' @param mean.model a list of parameters
#' @param variance.model a list of parameters
#' @param threshold Threshold for the estimation of alpha and beta
#' @param quiet logical, when `TRUE` the function will not display any message.
#'
#' @export
control.solarModel <- function(object, loss = "ml",
                               clearsky.model = control.seasonalModel(),
                               mean.model = list(seasonalOrder = 1, arOrder = 2, include.intercept = FALSE),
                               variance.model = list(seasonalOrder = 1, match_moments = FALSE),
                               threshold = 0.001, quiet = FALSE){
  if (missing(object)) {
    structure(
      list(
        loss = loss,
        threshold = threshold,
        clearsky.model = clearsky.model,
        mean.model = mean.model,
        variance.model = variance.model,
        quiet = quiet
      ),
      class = c("control", "list")
    )
  } else {
    object$control
  }
}

#' solarOption control parameters
#'
#' @export
control.solarOption <- function(nyears = c(2010, 2022), K = 0, put = TRUE, B = discount_factor()){
  structure(
    list(
      nyears = nyears,
      from = as.Date(paste0(nyears[1], "-01-01")),
      to = as.Date(paste0(nyears[2], "-01-01")),
      K = K,
      put = put,
      B = B
    )
  )
}

#' solarEsscher control parameters
#'
#' @export
control.solarEsscher <- function(nsim = 200, ci = 0.05, seed = 1, quiet = FALSE, n_key_points = 15,
                                 init_lambda = 0, lower_lambda = -1, upper_lambda = 1){
  list(
    nsim = nsim,
    ci = ci,
    seed = seed,
    quiet = quiet,
    n_key_points = n_key_points,
    init_lambda = 0,
    lower_lambda = -1,
    upper_lambda = 1
  )
}


