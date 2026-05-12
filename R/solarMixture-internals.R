#' Correct the moments to ensure moments matching
#'
#' @keywords solarMixture internal
#' @export
#' @noRd
solarMixture_moments_match <- function(coefficients, e_target, v_target, sk_target, kt_target, x){
  # Loss function
  loss_function <- function(params, prob, e_target, v_target, sk_target, kt_target, x){
    # Extract the parameters
    means <- params[stringr::str_detect(names(params), "mu")]
    sigma <- params[stringr::str_detect(names(params), "sd")]
    probs <- c(p1 = prob[[1]], p2 = 1 - prob[[1]])
    # Compute the moments
    mom <- GM_moments(means, sigma, probs)
    loss <- 0
    # Compute the loss from target moments
    if (!is.na(e_target)){
      loss <- loss + abs(mom$mean - e_target)
    }
    if (!is.na(v_target)){
      loss <- loss + abs(mom$variance - v_target)
    }
    if (!is.na(sk_target)){
      loss <- loss + abs(mom$skewness - sk_target)
    }
    if (!is.na(kt_target)){
      loss <- loss + abs(mom$kurtosis - kt_target)
    }
    # Compute log-likelihood
    loss <- 10000*loss^2 #- GM_loglik(means, sigma, probs, x)
    return(loss)
  }
  # Monthly optimization
  opt_coefficients <- coefficients
  for(nmonth in 1:nrow(coefficients)) {
    # Initial parameters
    params <- unlist(purrr::flatten(coefficients[nmonth,]))
    # Optimal parameters
    opt <- optim(par = params[-5], fn = loss_function, x = x[[nmonth]], prob = params[5],
                 e_target = e_target[nmonth], v_target = v_target[nmonth], sk_target = sk_target[nmonth],
                 kt_target = kt_target[nmonth])
    opt_coefficients[nmonth, 1:4] <- dplyr::bind_rows(opt$par)
  }
  opt_coefficients$p2 <- 1 - opt_coefficients$p1
  return(opt_coefficients)
}

#' Compute the Value at Risk and Expected Shortfall of a SolarMixture
#'
#' @param model solarMixture
#' @param alpha Numeric vector of confidence levels. Allows for more than one alpha.
#' @param ci Numeric scalar, confidence levels used to state if the Null is rejected or not on VaR tests.
#' @param ES Logical, when `TRUE` the expected shortfall will be also computed for each alpha.
#' @param type Numeric, type of evaluation, `full` on the complete data, `train` on the train data, `test` on the test data.
#'
#' @rdname solarMixture_VaR
#' @name solarMixture_VaR
#' @export
solarMixture_VaR <- function(solarMix, date, x, alpha = 0.05, ci = 0.05, ES = FALSE){
  # Create a base dataset
  data <- dplyr::tibble(date = date, x = x)
  # Number of observations
  N <- length(date)
  # Number of VaRs
  k <- length(alpha)
  # 1) Compute VaR, Violations and ES
  VaR_alpha <- solarMix$VaR(date, alpha)
  # Violations of the VaR
  Viol_VaR_alpha <- VaR_viol(x, VaR_alpha)
  colnames(Viol_VaR_alpha) <- paste0("e_", alpha)
  # 2) Perform the tests on the violations
  VaR_tests <- VaR_test(Viol_VaR_alpha, alpha, ci = ci)
  # Expected shortfall
  if (ES) {
    ES_alpha <- solarMix$ES(date, alpha)
  }
  # 3) Summarise the results
  # Initialize a matrix to store the VaR
  VaR_alpha_emp <- dplyr::bind_rows(VaR_alpha[1,])
  VaR_alpha_emp[1,] <- NA
  if (ES) {
    # Initialize a matrix to store the model's Expected Shortfall
    ES_alpha_mod <- dplyr::bind_rows(ES_alpha[1,])
    ES_alpha_mod[1,] <- NA
    # Initialize a matrix to store the empiric Expected Shortfall
    ES_alpha_emp <- dplyr::bind_rows(ES_alpha[1,])
    ES_alpha_emp[1,] <- NA
  }
  # Evaluate Empiric VaR and Expected Shortfalls
  for(i in 1:k){
    # Index of the violations
    idx_violations <- which(Viol_VaR_alpha[,i] == 1)
    if (!purrr::is_empty(idx_violations)){
      VaR_alpha_emp[1,i] <- sum(Viol_VaR_alpha[,i] == 1)/nrow(Viol_VaR_alpha)
      if (ES) {
        ES_alpha_mod[1,i] <- mean(ES_alpha[idx_violations,i])
        ES_alpha_emp[1,i] <- mean(x[idx_violations])
      }
    }
  }
  # Standard output
  output <- structure(
    list(
      VaR = dplyr::bind_cols(data, VaR_alpha),
      Viol = dplyr::bind_cols(data, Viol_VaR_alpha),
      VaR_emp = VaR_alpha_emp,
      VaR_test = VaR_tests
    )
  )
  if (ES) {
    output$ES <- dplyr::bind_cols(data, ES_alpha)
    output$ES_tot <- dplyr::bind_cols(Type = c("Model", "Empiric"),
                                      dplyr::bind_rows(ES_alpha_mod, ES_alpha_emp))
  }
  return(output)
}

#' Create a function of time for monthly parameters
#'
#' @param params Vector of length 12 with the monthly parameters.
#'
#' @examples
#' set.seed(1)
#' params <- runif(12)
#' mp <- monthlyParams$new(params)
#' t_now <- as.Date("2022-01-01")
#' t_hor <- as.Date("2024-12-31")
#' dates <- seq.Date(t_now, t_hor, by = "1 day")
#' plot(mp$predict(dates), type = "l")
#' @note Version 1.0.0
#' @export
monthlyParams <- R6::R6Class("monthlyParams",
                             public = list(
                               #' @description
                               #' Initialize a `monthlyParams` object
                               #' @param params numeric vector of parameters with length 12.
                               initialize = function(params){
                                 if (length(params) != 12) {
                                   stop("The length of the vector of parameters must be 12!")
                                 }
                                 private$..parameters <- params
                               },
                               #' @description
                               #' Predict the monthly paramete
                               #' @param x date as character or month as numeric.
                               predict = function(x){
                                 nmonths <- lubridate::month(x)
                                 par <- c()
                                 for(i in 1:length(nmonths)){
                                   par[i] <- self$parameters[nmonths[i]]
                                 }
                                 return(par)
                               },
                               #' @description
                               #' Update the monthly parameters
                               #' @param params numeric vector of parameters with length 12.
                               update = function(params){
                                 if (length(params) != 12) {
                                   stop("The length of the vector of parameters must be 12!")
                                 }
                                 private$..parameters <- params
                               }
                             ),
                             private = list(
                               version = "1.0.0",
                               ..parameters = NA
                             ),
                             active = list(
                               #' @field parameters vector of parameters with length 12.
                               parameters = function(){
                                 private$..parameters
                               }
                             ))
