#' Control parameters for a solar option
#'
#' @param nyears numeric vector. Interval of years considered. The first element will be the minimum and the second the maximum years used in
#' the computation of the fair payoff.
#' @param K numeric, level for the strike with respect to the seasonal mean. The seasonal mean is multiplied by `exp(K)`.
#' @param leap_year logical, when `FALSE`, the default, the year will be considered of 365 days, otherwise 366.
#' @param nsim integer, number of simulations used to bootstrap the premium's bounds. See \code{\link{solarOption_historical_bootstrap}}.
#' @param ci numeric, confidence interval for bootstrapping. See \code{\link{solarOption_historical_bootstrap}}.
#' @param seed integer, random seed for reproducibility. See \code{\link{solarOption_historical_bootstrap}}.
#' @examples
#' control_options <- control_solarOption()
#'
#' @rdname control_solarOption
#' @name control_solarOption
#' @keywords control
#' @note Version 1.0.0.
#' @export
control_solarOption <- function(nyears = c(2005, 2025), K = 0, leap_year = FALSE, nsim = 200, ci = 0.05, seed = 1){
  structure(
    list(
      nyears = nyears,
      from = as.Date(paste0(nyears[1], "-01-01")),
      to = as.Date(paste0(nyears[2], "-01-01")),
      K = K,
      leap_year = leap_year,
      nsim = nsim,
      ci = ci,
      seed = seed
    ),
    class = c("controlSolarOption", "list")
  )
}

#' Payoff function of a Solar Option
#'
#' @param R Numeric, vector of values of solar radiation at maturity.
#' @param K Numeric, scalar or vector of strikes.
#' @param put Logical, when `TRUE`, the default, the function will return the output of a put payoff
#' otherwise a call payoff. See the details.
#' @details When `put = TRUE`, the function evaluate a put payoff:
#' \deqn{(K - R)^{+} = (K - R) 1_{K > R}}
#' When `put = TRUE`, the function evaluate a call payoff:
#' \deqn{(R-K)^{+} = (R - K) 1_{R \ge K}}
#' @examples
#' # Put
#' solarOption_payoff(10, 9, put = TRUE)
#' solarOption_payoff(8, 9, put = TRUE)
#' # Call
#' solarOption_payoff(10, 9, put = FALSE)
#' solarOption_payoff(8, 9, put = FALSE)
#' @rdname solarOption_payoff
#' @name solarOption_payoff
#' @keywords solarOption
#' @note Version 1.0.0.
#' @export
solarOption_payoff <- function(R, K = 0, put = TRUE){
  put <- ifelse(put, -1, 1)
  payoff <- (R - K) * put
  payoff[payoff < 0] <- 0
  return(payoff)
}

#' Discount factor of a Solar Option
#'
#' @param tau Numeric, time to maturity in days.
#' @param P Numeric, price of the contract.
#' @param Gamma_h Numeric, hedged payoff.
#' @param r Numeric, daily risk-free rate.
#' @details The discount factor reads:
#' \deqn{B(\tau, P, \Gamma^h, r) = e^{-r \tau} + \frac{\Gamma^h}{P} (1 - e^{-r \tau})}
#' @examples
#' solarOption_discount(365. 0.6, 2, 0.000008)
#' solarOption_discount(365, 0.3, 2, 0.00008)
#'
#' @rdname solarOption_discount
#' @name solarOption_discount
#' @keywords solarOption
#' @note Version 1.0.0.
#' @export
solarOption_discount <- function(tau, P = 1, Gamma_h = 0, r = 0.03/365){
  exp(-r * tau) + (Gamma_h / P) * (1 - exp(-r * tau))
}

#' Compute the price of a `solarOption`
#'
#' @param K Numeric, strike price.
#' @param pdf_R Function, density function of solar radiation at maturity.
#' @param R_min Numeric, minimum value of solar radiation at maturity.
#' @param R_max Numeric, maximum value of solar radiation at maturity.
#' @details When `put = TRUE`, the function evaluate a put payoff:
#' \deqn{\mathbb{E}\{(K - R)^{+}\} = K \mathbb{P}\{R_T < K\} - \mathbb{E}\{R 1_{K > R_T}\}}
#' When `put = TRUE`, the function evaluate a call payoff:
#' \deqn{\mathbb{E}\{(K - R)^{+}\} = \mathbb{E}\{R 1_{K < R_T}\} - K \mathbb{P}\{R_T > K\}}
#' @examples
#' # Parameters for the pdf of Rt
#' Ct <- 7
#' alpha <- 0.01
#' beta <- 0.92
#' Rt_min_max <- Ct * c(1-alpha-beta, 1-alpha)
#' pdf_R <- function(x) dsolarGHI(x, Ct, alpha, beta, function(x) dnorm(x))
#' cdf_R <- function(x) psolarGHI(x, Ct, alpha, beta, function(x) pnorm(x))
#' # Put price
#' solarOption_price(5, R_min = Rt_min_max[1], R_max = Rt_min_max[2], TRUE)
#' # Call price
#' solarOption_price(5, pdf_R, Rt_min_max[1], Rt_min_max[2], FALSE)
#' @rdname solarOption_price
#' @name solarOption_price
#' @keywords solarOption
#' @note Version 1.0.0.
#' @export
solarOption_price <- function(K, pdf_R, R_min, R_max, put = TRUE){
  if(!is.function(pdf_R)) stop("pdf_R is not a function!")
  # Bounds
  lower <- ifelse(put, R_min, K)
  upper <- ifelse(put, K, R_max)
  # Payoff function
  Gamma <- function(x) solarOption_payoff(x, K, put = put) * pdf_R(x)
  # Premium
  price <- integrate(function(x) solarOption_payoff(x, K, put = put) * pdf_R(x),
                     stop.on.error = FALSE, lower = lower, upper = upper)$value
  # Return the option price
  return(price)
}


#' Structure the outputs of solarOption functions
#'
#' @param data df_payoff
#' @param leap_year control params
#'
#' @return An object of the class `solarOptionPayoff`.
#'
#' @rdname solarOptionPayoff
#' @name solarOptionPayoff
#' @keywords solarOption
#' @note Version 1.0.0.
#' @export
solarOptionPayoff <- function(data, leap_year = FALSE){
  # Reorder variables
  df_payoff <- dplyr::select(data, side, date, Year, Month, Day, Rt, strike, premium, payoff, exercise, max_payoff)
  # Include or not 29-th of February from computation
  if (!leap_year) {
    df_payoff <- dplyr::filter(df_payoff, !(Month == 2 & Day == 29))
  }
  # Aggregation by Year and Month
  df_year_month <- df_payoff %>%
    dplyr::group_by(Year, Month, side) %>%
    dplyr::reframe(ndays = dplyr::n(),
                   payoff = sum(payoff),
                   premium = sum(premium),
                   exercise = mean(exercise),
                   Rt = sum(Rt),
                   strike = sum(strike),
                   max_payoff = sum(max_payoff))

  # Aggregation by Month and Day
  df_month_day <- df_payoff %>%
    dplyr::group_by(Month, Day, side) %>%
    dplyr::reframe(premium = mean(premium),
                   payoff = mean(payoff),
                   exercise = mean(exercise),
                   Rt = mean(Rt),
                   strike = mean(strike),
                   max_payoff = mean(max_payoff))

  # Aggregation by Month
  df_month <- df_year_month %>%
    dplyr::group_by(Month, side)  %>%
    dplyr::reframe(ndays = dplyr::n(),
                   payoff = mean(payoff),
                   premium = mean(premium),
                   exercise = mean(exercise),
                   Rt = mean(Rt),
                   strike = mean(strike),
                   max_payoff = mean(max_payoff))

  # Aggregation by Year
  df_year <- df_month %>%
    dplyr::group_by(side)  %>%
    dplyr::reframe(ndays = sum(ndays),
                   payoff = sum(payoff),
                   premium = sum(premium),
                   exercise = mean(exercise),
                   Rt = sum(Rt),
                   strike = sum(strike),
                   max_payoff = sum(max_payoff))
  # Structured output
  structure(
    list(
      payoff = df_payoff,
      payoff_year_month = df_year_month,
      payoff_month_day = df_month_day,
      payoff_month = df_month,
      payoff_year = df_year
    ),
    class = c("solarOptionPayoff", "list")
  )
}


#' Calibrate the implied Choquet parameter
#'
#' @param P_target Numeric, target price to match.
#' @param r_imp Numeric, implied return from measure under P.
#' @inheritParams solarOption_model
#'
#' @details The optimization function will find the best parameter \eqn{\lambda} such that
#' \deqn{\underset{-0.5 < \lambda < 1}{\text{argmin}} \left(P^{\mathbb{P}(\lambda)} (1 + r^{\text{imp}}) - P_{\text{target}}\right)^2}
#'
#' @rdname solarOption_lambda
#' @name solarOption_lambda
#' @keywords solarOption
#' @note Version 1.0.0.
#' @export
solarOption_lambda <- function(model, moments, P_target, r_imp = 0, put = TRUE, quiet = FALSE, control_options = control_solarOption()){
  # Yearly premium to establish the calibration
  # price_P <- solarOption_model(model, moments, lambda = 0, put = put, control_options = control_options)$payoff_year$premium
  # price_P <- price_P * (1 + r_imp)
  # Loss function
  loss_function <- function(par, r_imp, P_target){
    # Price under P(lambda)
    price_P_lambda <- solarOption_model(model, moments, lambda = par, put = put, control_options = control_options)$payoff_year$premium
    # Fitted target price
    P_target_fit <- price_P_lambda * (1 + r_imp)
    # Square loss function
    loss <- (P_target_fit - P_target)^2
    if(!quiet) cli::cli_alert_info(paste0("Implied  parameter: ", round(par, 4), " loss: ", round(loss, 7)))
    return(loss)
  }
  # Optimization
  lambda <- optim(0, loss_function, method = "Brent", lower = -0.5, upper = 1, r_imp = r_imp, P_target = P_target)$par
  return(lambda)
}
