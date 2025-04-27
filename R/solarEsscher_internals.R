#' Function to establish up and down parameters of the Esscher transform
#' A positive theta identify a `down` parameter, a negative theta identify an `up` parameter.
#'
#' @param theta Esscher parameter.
#' @return A `list` with first element named `up` with the positive parameter and second element named `dw` with the negative one.
#'
#' @keywords internals
#' @noRd
#' @export
solarEsscher_bounds <- function(theta){
  # Identify the up and down parameters of Esscher transform
  # dw parameter (positive): theta > 0
  # up parameter (negative): theta < 0
  par <- list(up = 0, bar = 0, dw = 0)
  if (theta >= 0) {
    par[["dw"]] <- theta
    par[["up"]] <- -theta/(1+theta)
  } else {
    par[["dw"]] <- -theta/(1-theta)
    par[["up"]] <- theta
  }
  par[["bar"]] <- 0.5*(par[["dw"]] + par[["up"]])
  return(par)
}

#' Change probability according to Esscher parameters
#'
#' @rdname solarEsscher_probability
#' @name solarEsscher_probability
#' @export
solarEsscher_probability <- function(params = c(0,0,1,1,0.5), df_n, theta = 0){
  params <- list(
    mu_up = df_n$Yt_bar + df_n$Yt_tilde_uncond + df_n$Yt_tilde_hat + df_n$sigma*df_n$sigma_bar*params[1],
    mu_dw = df_n$Yt_bar + df_n$Yt_tilde_uncond + df_n$Yt_tilde_hat + df_n$sigma*df_n$sigma_bar*params[2],
    sd_up = params[3]*df_n$sigma_bar*df_n$sigma,
    sd_dw = params[4]*df_n$sigma_bar*df_n$sigma,
    p_up = params[5]
  )
  params <- unlist(params)
  num <- params[5]*exp(theta*params[1] + 0.5*(theta^2*params[3])^2)
  den <- (1-params[5])*exp(theta*params[2] + 0.5*(theta^2*params[4])^2)
  num/(num + den)
}




