#' Perform normality tests
#'
#' @param x numeric, a vector of observation.
#' @param pvalue numeric, the desiderd level of `p.value` at which the null hypothesis will be rejected.
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(1000, 0, 1) + rchisq(1000, 1)
#' test_normality(x)
#' x <- rnorm(1000, 0, 1)
#' test_normality(x)
#'
#' @return a tibble with the results of the normality tests.
#'
#' @rdname test_normality
#' @name test_normality
#'
#' @export
test_normality <- function(x = NULL, pvalue = 0.05){

  safe_tests <- list(
    ad.test = purrr::safely(nortest::ad.test),
    cvm.test = purrr::safely(nortest::cvm.test),
    lillie.test = purrr::safely(nortest::lillie.test),
    pearson.test = purrr::safely(nortest::pearson.test),
    sf.test = purrr::safely(nortest::sf.test),
    shapiro.test = purrr::safely(stats::shapiro.test)
  )
  tests <- list()
  for(i in 1:length(safe_tests)){
    test <- suppressWarnings(safe_tests[[i]](x)$result)
    if (!is.null(test)) {
      tests[[i]] <- broom::tidy(test)
    }
  }
  tests <- dplyr::bind_rows(tests)
  tests <- dplyr::mutate(tests, H0 = ifelse(p.value <= pvalue, "Rejected", "Non Rejected"))
  tests <- dplyr::select(tests, method, statistic, p.value, H0)

  return(tests)
}
