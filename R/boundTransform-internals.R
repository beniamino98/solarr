#' List all the available link functions
#'
#' @examples
#' boundTransform_av_links()
#'
#' @keywords utils
#' @note Version 1.0.0
#'
#' @export
boundTransform_av_links <- function(){
  av_links <- ls(pattern = "boundTransform_link\\.", envir = .GlobalEnv)
  av_links <- c(av_links, ls(getNamespace("solarr"), pattern = "boundTransform_link"))
  av_links <- unique(av_links)
  stringr::str_remove_all(av_links, "boundTransform_link.")
}

#' Link function (invGumbel)
#'
#' @keywords utils
#' @note Version 1.0.0
boundTransform_link.invgumbel <- function(){
  structure(
    list(
      g = function(x) log(-log(x)),
      ig = function(x) exp(-exp(x)),
      # eqivalent to. 1/(x * log(x))
      g_prime = function(x) 1/(log(x^x)),
      monotonicity = "decreasing",
      link = "invgumbel"
    ),
    class = c("boundTransform_link", "list")
  )
}

#' Link function (Gumbel)
#'
#' @keywords utils
#' @note Version 1.0.0
boundTransform_link.gumbel <- function(){
  structure(
    list(
      g = function(x) -log(-log(x)),
      ig = function(x) exp(-exp(-x)),
      # eqivalent to. -1/(x * log(x))
      g_prime = function(x) -1/(log(x^x)),
      monotonicity = "increasing",
      link = "gumbel"
    ),
    class = c("boundTransform_link", "list")
  )
}

#' Link function (logit)
#'
#' @keywords utils
#' @note Version 1.0.0
boundTransform_link.logis <- function(){
  structure(
    list(
      g = function(x) log(x/(1-x)),
      ig = function(x) 1/(1+exp(-x)),
      g_prime = function(x) 1 / (x * (1 - x)),
      monotonicity = "increasing",
      link = "logis"
    ),
    class = c("boundTransform_link", "list")
  )
}
#' Link function (probit)
#'
#' @keywords utils
#' @note Version 1.0.0
boundTransform_link.norm <- function(){
  structure(
    list(
      g = function(x) qnorm(x),
      ig = function(x) pnorm(x),
      g_prime = function(x) 1 / dnorm(qnorm(x)),
      monotonicity = "increasing",
      link = "norm"
    ),
    class = c("boundTransform_link", "list")
  )
}

#' Link function (identity)
#'
#' @keywords utils
#' @note Version 1.0.0
boundTransform_link.identity <- function(){
  structure(
    list(
      g = function(x) x,
      ig = function(x) x,
      g_prime = function(x) 1,
      monotonicity = "increasing",
      link = "identity"
    ),
    class = c("boundTransform_link", "list")
  )
}



