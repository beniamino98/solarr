#' PIT computation
#'
#' @param Rt vector of solar radiation.
#' @param Ct vector of clear-sky radiation.
#' @keywords diagnostic
#' @rdname PIT
#' @name PIT
#' @export
PIT_psolarGHI <- function(Rt, Ct, alpha, beta, M_Y, S_Y, probs, link = "invgumbel"){
  # Number of observations
  N <- length(Rt)
  # PIT grades
  grade <- vector("numeric", length = N)
  for(n in 1:N){
    # Distribution of Yt
    cdf_Yt <- function(x) pmixnorm(x, mean = M_Y[n,], sd = S_Y[n,], alpha = probs[n,])
    # Grades on Rt
    grade[n] <- psolarGHI(Rt[n], Ct[n], alpha, beta, cdf_Yt, link = link)
  }
  grade
}
#' Perform KS and LB test on PIT grades
#'
#' @param u vector, PIT.
#' @param type type of autocorrelation test ("bp" or "lb")
#' @param lag vector of lags at which perform autocorrelation test.
#' @keywords diagnostic
#' @rdname PIT_test
#' @name PIT_test
#' @export
PIT_test <- function(u, type = "bp", lag = 5){
  # 1) KS-tests vs uniform
  KS.test <- ks_test(u, punif, ci = 0.01, min_quantile = 0.015, max_quantile = 0.985, k = 1000, plot = FALSE)
  KS.test$Test <- "KS"
  KS.test$parameter <- NA
  KS.test$statistic <- KS.test$KS
  # 2) Autocorrelations tests
  type <- names(match.arg(type, choices = c(`Box-Pierce` = "bp", `Ljung-Box` = "lb")))
  LB.test <- purrr::map_df(lag, ~broom::tidy(Box.test(u, lag = .x, type = type)))
  LB.test$Test <- toupper(type)
  LB.test$H0 <- ifelse(LB.test$p.value <= 0.01, "Rejected", "Non-Rejected")
  dplyr::bind_rows(
    KS.test[,c("Test", "parameter", "statistic", "p.value", "H0")],
    LB.test[,c("Test", "parameter", "statistic", "p.value", "H0")]
  )
}
