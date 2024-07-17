#' solarOption control parameters
#' @rdname control_solarOption
#' @export
control_solarOption <- function(nyears = c(2010, 2022), K = 0, put = TRUE, B = discount_factor()){
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
#' @rdname control_solarEsscher
#' @export
control_solarEsscher <- function(nsim = 200, ci = 0.05, seed = 1, quiet = FALSE, n_key_points = 15,
                                 init_lambda = 0, lower_lambda = -1, upper_lambda = 1){
  list(
    nsim = nsim,
    ci = ci,
    seed = seed,
    quiet = quiet,
    n_key_points = n_key_points,
    init_lambda = init_lambda,
    lower_lambda = lower_lambda,
    upper_lambda = upper_lambda
  )
}
