#' Distribution test
#'
#' Evaluate a Kolmogorov-Smirnov test on the residuals of a `solarModel` model
#' object against the estimated Gaussian mixture distribution.
#'
#' @param model An object of the class `solarModel`
#' @param H0 Character, null hypothesis for the residuals distribution. Can be:
#' \describe{
#'  \item{`"gm"`}{for a test against a gaussian mixture distribution;}
#'  \item{`"norm"`}{for a test against a gaussian distribution.}
#'}
#' @param type character(1), type of test. Can be:
#' \describe{
#'  \item{`"train"`}{the test is performed only on train data;}
#'  \item{`"test"`}{the test is performed only on test data;}
#'  \item{`"full"`}{the test is performed only on the full data.}
#'}
#' @inheritParams ks_test
#'
#' @examples
#' model = solarModel$new(spec)
#' model$fit()
#' solarModel_test_distribution(model)
#'
#' @rdname solarModel_test_distribution
#' @name solarModel_test_distribution
#' @keywords solarModel
#' @note Version 1.0.0.
#' @export
solarModel_test_distribution <- function(model, H0 = "gm", type = "full", ci = 0.05, min_quantile = 0.025, max_quantile = 0.985){

  # Match H0 for the target distribution
  H0 <- match.arg(H0, choices = c("gm", "norm"))
  # Type of data used
  type = match.arg(type, choices = c("train", "test", "full"))

  # Extract data
  data <- model$data
  # Filter to exclude or not test data
  if (type == "train") {
    data <- data[data$isTrain & data$weights != 0,]
  } else if (type == "test") {
    data <- data[!data$isTrain,]
  }
  # Monthly tests
  tests <- list()
  for(nmonth in 1:12){
    # Residuals
    x <- dplyr::filter(data, Month == nmonth)$u_tilde
    if (H0 == "norm"){
      mu_x <- mean(x)
      sd_x <- sd(x)
      # Normal CDF
      cdf_Yt <- function(x) pnorm(x, mu_x, sd_x)
    } else {
      # Gaussian mixture parameters
      gm <- model$spec$mixture.model$model[[nmonth]]
      mu_x <- gm$means
      sd_x <- gm$sd
      prob <- gm$p
      # Mixture CDF
      cdf_Yt <- function(x) pmixnorm(x, mu_x, sd_x, prob)
    }
    # Distribution test
    tests[[nmonth]] <- ks_test(x, cdf_Yt, ci = ci, min_quantile = min_quantile, max_quantile = max_quantile)
  }
  # Unique dataset
  tests <- dplyr::bind_rows(tests)
  # Output
  tests <- dplyr::bind_cols(Month = 1:12,
                            H0 = H0,
                            type = type,
                            tests[,-6],
                            result = ifelse(tests$H0 == "Non-Rejected", "Passed", "Not-passed"))
  return(tests)
}

#' Autocorrelation test
#'
#' Evaluate the autocorrelation in the components of a `solarModel` object.
#'
#' @param model An object of the class `solarModel`
#' @param lag.max integer(1), maximum lag to consider for the test.
#' @inheritParams solarModel_test_distribution
#' @param method character(1), method of test performed. Can be:
#' \describe{
#'  \item{`"bg"`}{for Breush-Godfrey test;}
#'  \item{`"bp"`}{for Box-pierce test;}
#'  \item{`"lb"`}{for Ljung-Box test.}
#'}
#'
#' @examples
#' model = solarModel$new(spec)
#' model$fit()
#' solarModel_test_autocorr(model, method = "lb")
#'
#' @rdname solarModel_test_autocorr
#' @name solarModel_test_autocorr
#' @keywords solarModel
#' @note Version 1.0.0.
#' @export
solarModel_test_autocorr <- function(model, nyears, method = "bg", type = "full", lag.max = 3, ci = 0.05){
  # Type of test
  type <- match.arg(type, choices = c("train", "test", "full"))
  # Method
  method <- match.arg(method, choices = c("bg", "bp", "lb"))

  # 1) Extract the data
  data <- model$data
  if (type == "train") {
    data <- data[data$isTrain & data$weights != 0,]
  } else if (type == "test") {
    data <- data[!data$isTrain,]
  }

  if (!missing(nyears)){
    data <- dplyr::filter(data, Year %in% nyears)
  }

  # Mixture monthly moments
  mix.mom <- dplyr::select(model$spec$mixture.model$moments, Month, mean, variance)
  data <- dplyr::left_join(data, mix.mom, by = c("Month"))
  # Mixture standardization
  data$z_tilde <- (data$u_tilde - data$mean) / sqrt(data$variance)
  # Squared random variables
  data$Yt_tilde2 <- data$Yt_tilde^2
  data$eps2 <- data$eps^2
  data$eps_tilde2 <- data$eps_tilde^2
  data$u_tilde2 <- data$u_tilde^2
  data$z_tilde2 <- data$z_tilde^2

  # Autocorrelation test Breusch-Godfrey
  bg_test <- function(target, data, lag.max, ci, expected = "Not-rejected"){
    formula <- as.formula(paste0(target, "~ 1"))
    test <- lmtest::bgtest(formula, order = lag.max, data = data)
    test <- broom::tidy(test)
    test$target <- target
    test <- dplyr::select(test, target, statistic, p.value, H0 = "method", lags = "parameter")
    test$H0 <- ifelse(test$p.value > ci, "Not-rejected", "Rejected")
    test$p.value <- round(test$p.value, digits = 5)
    test$result <- ifelse(test$H0 == expected, "passed", "Not-passed")
    return(test)
  }
  # Autocorrelation test Ljung-Box
  lb_test <- function(target, data, lag.max, ci, expected = "Not-rejected"){
    test <- Box.test(data[[target]], lag = lag.max, type = "Ljung-Box")
    test <- broom::tidy(test)
    test$target <- target
    test <- dplyr::select(test, target, statistic, p.value, H0 = "method", lags = "parameter")
    test$H0 <- ifelse(test$p.value > ci, "Not-rejected", "Rejected")
    test$p.value <- round(test$p.value, digits = 5)
    test$result <- ifelse(test$H0 == expected, "passed", "Not-passed")
    return(test)
  }
  # Autocorrelation test Box-pierce
  bp_test <- function(target, data, lag.max, ci, expected = "Not-rejected"){
    test <- Box.test(data[[target]], lag = lag.max, type = "Box-Pierce")
    test <- broom::tidy(test)
    test$target <- target
    test <- dplyr::select(test, target, statistic, p.value, H0 = "method", lags = "parameter")
    test$H0 <- ifelse(test$p.value > ci, "Not-rejected", "Rejected")
    test$p.value <- round(test$p.value, digits = 5)
    test$result <- ifelse(test$H0 == expected, "passed", "Not-passed")
    return(test)
  }

  # Variables to test
  ## Residuals
  targets <- c("Yt_tilde", "eps", "eps_tilde", "u_tilde", "z_tilde")
  ## Squared residuals
  targets <- c(targets, "Yt_tilde2", "eps2", "eps_tilde2", "u_tilde2", "z_tilde2")

  # Expected results
  expected <- c("Rejected", "Not-rejected", "Not-rejected", "Not-rejected", "Not-rejected")
  expected <- c(expected, "Rejected", "Rejected", "Rejected", "Not-rejected", "Not-rejected")
  names(expected) <- targets

  # Iterate the test for all the variables
  tests <- list()
  for(target in targets){
    if (method == "bg"){
      tests[[target]] <- bg_test(target, data, lag.max, ci, expected = expected[target])
    } else if (method == "bp") {
      tests[[target]] <- bp_test(target, data, lag.max, ci, expected = expected[target])
    } else if (method == "lb") {
      tests[[target]] <- lb_test(target, data, lag.max, ci, expected = expected[target])
    }
  }
  tests <- dplyr::bind_rows(tests)
  return(tests)
}


#' Autocorrelation and distribution tests
#'
#' Evaluate a Kolmogorov-Smirnov test on the residuals of a `solarModel` model
#' object against the estimated Gaussian mixture distribution and a Breush-pagan or Box-pierce
#' test on the residuals.
#' @param lags Numeric vector. Lags on which perform the autocorrelation tests.
#' @inheritParams solarModel_test_distribution
#' @inheritParams solarModel_test_autocorr
#'
#' @examples
#' model = solarModel$new(spec)
#' model$fit()
#' solarModel_tests(model, train_data = TRUE)
#'
#' @rdname solarModel_tests
#' @name solarModel_tests
#' @keywords solarModel
#' @note Version 1.0.0.
#' @export
solarModel_tests <- function(model, lags = c(7),  method = "bg", type = "full",
                             ci = 0.05, min_quantile = 0.025, max_quantile = 0.985){

  # Number of lags to test
  n.lags <- length(lags)


  # Tests for absence of autocorrelation
  autocorr_tests <- list()
  for(lag in lags){
    test <- solarModel_test_autocorr(model, method = method, type = type, lag.max = lag, ci = ci)
    autocorr_tests <- append(autocorr_tests, setNames(list(test), paste0("lag_", lag)))
  }

  # Tests for normal distribution
  normality_test <- solarModel_test_distribution(model, H0 = "norm", type = type,
                                                 ci = ci, min_quantile = min_quantile, max_quantile = max_quantile)
  # Tests for gaussian mixture distribution
  mixture_test <- solarModel_test_distribution(model, H0 = "gm", type = type,
                                               ci = ci, min_quantile = min_quantile, max_quantile = max_quantile)
  structure(
    list(
      autocorr = autocorr_tests,
      normality = normality_test,
      mixture = mixture_test
    )
  )
}

#' Evaluate a KS test on the PIT.
#'
#' @inheritParams solarModel_test_distribution
#' @rdname solarModel_PIT_test
#' @name solarModel_PIT_test
#' @details See the functions \code{\link{solarModel_PIT}} and \code{\link{ks_test}} for more details.
#' @keywords solarModel
#' @note Version 1.0.0.
#' @export
solarModel_PIT_test <- function(model, type = "full", nyears, ci = 0.05, min_quantile = 0.015, max_quantile = 0.985){
  # Type of computation
  type = match.arg(type, choices = c("train", "test", "full"))
  # Extract conditional moments
  moments <- model$moments$conditional
  if (type == "train") {
    moments <- moments[model$data$isTrain & model$data$weights != 0,]
  } else if (type == "test") {
    moments <- moments[!model$data$isTrain,]
  }
  if (!missing(nyears)){
    moments <- dplyr::filter(moments, Year %in% nyears)
  }

  # Compute the grades
  moments <- solarModel_PIT(model, moments)
  # Perform a KS test on PIT
  dplyr::bind_cols(link = model$spec$transform$link,
                   ks_test(moments$grade, function(x) punif(x), ci = ci, min_quantile = min_quantile, max_quantile = max_quantile))
}


#' Compute the Log-predictive density of a solarModel
#'
#' @rdname solarModel_test_LPD
#' @name solarModel_test_LPD
#' @keywords solarModel
#' @note Version 1.0.0.
#' @export
solarModel_test_LPD <- function(model, type = "full"){
  # Match argument
  type = match.arg(type, choices = c("train", "test", "full"))
  # Extract conditional moments
  moments <- model$moments$conditional
  if (type == "train") {
    moments <- moments[model$data$isTrain & model$data$weights != 0,]
  } else if (type == "test") {
    moments <- moments[!model$data$isTrain,]
  }
  LPD <- model$logLik(moments, target = "GHI")
  LPD <- LPD[!is.infinite(LPD)]

  dplyr::tibble(
    type = type,
    link = model$spec$transform$link,
    LPD =  mean(LPD, na.rm = TRUE)
  )
}

#' Compute metrics to test forecasts
#'
#' @rdname solarModel_test_forecast
#' @name solarModel_test_forecast
#' @keywords solarModel
#' @note Version 1.0.0.
#' @export
solarModel_test_forecast <- function(model, ci = 0.1, type = c("train", "test", "full") ){
  # Type of dataset
  type = match.arg(type, choices = c("train", "test", "full"))
  # Extract conditional moments
  moments <- model$moments$conditional
  if (type == "train") {
    moments <- moments[model$data$isTrain & model$data$weights != 0,]
  } else if (type == "test") {
    moments <- moments[!model$data$isTrain,]
  }
  # Compute forecasts
  forecast <- solarModel_forecast(model, moments, ci = ci)

  # Evaluate the goodness of fit with different metrics
  errors <- forecast$Rt - forecast$e_Rt
  # SSE: sum of the squared errors
  SSE <- sum(errors^2)
  # RMSE: root mean squared error
  RMSE <- sqrt(mean(errors^2))
  # MAE: mean absolute error
  MAE <- mean(abs(errors))
  # MAPE: mean absolute percentage error
  MAPE <- mean(abs(errors) / forecast$Rt * 100)
  # Violations of Upper VaR
  viol_ci_hi <- mean(forecast$Rt > forecast$ci_mix_hi)
  # Violations of Lower VaR
  viol_ci_lo <- mean(forecast$Rt < forecast$ci_mix_lo)

  dplyr::tibble(type = type,
                SSE = SSE,
                RMSE = RMSE,
                MAE = MAE,
                MAPE = MAPE,
                ci = ci*2,
                viol_ci_hi = viol_ci_hi,
                viol_ci_lo = viol_ci_lo)
}


#' Compute metrics to test option pricing
#'
#' @rdname solarModel_test_pricing
#' @name solarModel_test_pricing
#' @keywords solarModel_test
#' @note Version 1.0.0.
#' @export
solarModel_test_pricing <- function(model, type = c("train", "test", "full"), control = control_solarOption()){
  # Type of dataset
  type = match.arg(type, choices = c("train", "test", "full"))
  # Extract conditional moments
  moments <- model$moments$conditional
  if (type == "train") {
    moments <- moments[model$data$isTrain & model$data$weights != 0,]
  } else if (type == "test") {
    moments <- moments[!model$data$isTrain,]
  }
  # Option prices
  put_prices <- solarOption_model(model, moments, put = TRUE, control_options = control)
  call_prices <- solarOption_model(model, moments, put = FALSE, control_options = control)
  # PUT / Call Stats
  put_call_stats <- function(price, exercise, Gamma){
    errors <- price - Gamma
    dplyr::tibble(
      SSE = sum(errors^2),
      premium = sum(price),
      payoff = sum(Gamma),
      diff = premium - payoff,
      cPt = mean(price[exercise != 0]),
      cGamma = mean(Gamma[exercise != 0]),
      cdiff = cGamma - cPt
    )
  }
  # PUT Stats
  put_stats <- put_call_stats(put_prices$payoff$premium, put_prices$payoff$exercise, put_prices$payoff$payoff)
  # Call Stats
  call_stats <- put_call_stats(call_prices$payoff$premium, call_prices$payoff$exercise, call_prices$payoff$payoff)

  # SoRadIDX Put stats
  # Average payoff
  Pt_hist = solarOption_historical(model, put = TRUE, control_options = control)$payoff_year$premium
  # Average premium
  Pt <- put_prices$payoff_year$premium
  # Realized payoff
  Gamma_p <- put_prices$payoff_year$payoff
  # Stats
  put_idx_stats <- bind_cols(Pt_hist = Pt_hist, Pt = Pt, Gamma = Gamma_p) %>%
    mutate(diff_Gamma_Pt = Gamma - Pt, diff_Gamma_Pt_hist = Gamma - Pt_hist, diff_Pt_hist_Pt = Pt_hist-Pt)

  # SoRadIDX Call stats
  Ct_hist = solarOption_historical(model, put = FALSE, control_options = control)$payoff_year$premium
  # Average premium
  Ct <- call_prices$payoff_year$premium
  # Realized payoff
  Gamma_c <- call_prices$payoff_year$payoff
  # Stats
  call_idx_stats <- bind_cols(Ct_hist = Ct_hist, Ct = Ct, Gamma = Gamma_c)%>%
    mutate(diff_Gamma_Ct = Gamma - Ct, diff_Gamma_Ct_hist = Gamma - Ct_hist, diff_Ct_hist_Ct = Ct_hist-Ct)

  structure(
    list(
      put = bind_cols(type = type, K = control$K, n = nrow(moments), put_stats),
      call = bind_cols(type = type, K = control$K, n = nrow(moments), call_stats),
      put_idx = bind_cols(type = type, K = control$K, n = nrow(moments), put_idx_stats),
      call_idx = bind_cols(type = type, K = control$K, n = nrow(moments), call_idx_stats)
    )
  )
}

