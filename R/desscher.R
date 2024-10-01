#' Esscher transform of a density
#'
#' Given a function of `x`, i.e. \eqn{f_X(x)}, compute its Esscher transform and return again a function of `x`.
#'
#' @param pdf density function.
#' @param theta Esscher parameter.
#' @param lower numeric, lower bound for integration, i.e. the lower bound for the pdf.
#' @param upper numeric, lower bound for integration, i.e. the upper bound for the pdf.
#'
#' @examples
#' # Grid of points
#' grid <- seq(-3, 3, 0.1)
#' # Density function of x
#' pdf <- function(x) dnorm(x, mean = 0)
#' # Esscher density (no transform)
#' esscher_pdf <- desscher(pdf, theta = 0)
#' pdf(grid) - esscher_pdf(grid)
#' # Esscher density (transform)
#' esscher_pdf_1 <- function(x) dnorm(x, mean = -0.1)
#' esscher_pdf_2 <- desscher(pdf, theta = -0.1)
#' esscher_pdf_1(grid) - esscher_pdf_2(grid)
#' # Log-probabilities
#' esscher_pdf(grid, log = TRUE)
#' esscher_pdf_2(grid, log = TRUE)
#'
#' @details Given a pdf \eqn{f_X(x)} the function computes its Esscher transform, i.e.
#'
#' \deqn{\mathcal{E}_{\theta}\{f_X(x)\} = \frac{e^{\theta x} f_X(x)}{\int_{-\infty}^{\infty} e^{\theta x} f_X(x) dx}}
#'
#' @rdname desscher
#' @name desscher
#' @export
desscher <- function(pdf, theta = 0, lower = -Inf, upper = Inf){
  # Esscher Numerator
  esscher_num <- function(x, theta = 0) ifelse(is.infinite(exp(theta*x)), 0, exp(theta*x))
  # Esscher Denominator
  esscher_den <- function(x, theta = 0) esscher_num(x, theta)*pdf(x)
  # Normalization factor
  den <- integrate(esscher_den, lower = lower, upper = upper, theta = theta)$value
  # Esscher density
  esscher_pdf <- function(x) (esscher_num(x, theta)*pdf(x))/den
  # Esscher pdf depending on `x`
  function(x, log = FALSE){
    # Probabilities
    probs <- purrr::map_dbl(x, esscher_pdf)
    # Log-probabilities
    if (log) {
      probs <- base::log(probs)
    }
    return(probs)
  }
}

#' Change probability according to esscher parameters
#' @rdname esscher_probability
#' @name esscher_probability
#' @export
esscher_probability <- function(params = c(0,0,1,1,0.5), df_n, theta = 0){
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


