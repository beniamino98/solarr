#' Kullback–Leibler divergence
#'
#' Compute the Kullback–Leibler distance between two probability measure.
#'
#' @param p Numeric, probability vector. Usually, the empiric probabilities.
#' @param q Numeric, probability vector. Usually, the model probabilities.
#' @param quiet Boolean, default is `TRUE`.
#' When set to `FALSE` the function will not display warnings.
#'
#' @details The function implements:
#' \deqn{\sum_{i} p_i \log(\frac{p_i}{q_i}) \quad p_i, q_i > 0 \; \forall i}
#'
#'
#' @references https://en.wikipedia.org/wiki/Kullback–Leibler_divergence
#'
#' @examples
#'
#' p <- dnorm(rnorm(100))
#' q <- dnorm(rnorm(100))
#' kl_dist(p, q)
#'
#' pdf_1 <- function(x) dnorm(x, mean = 2, sd = 1)
#' pdf_2 <- function(x) dnorm(x, mean = -2, sd = 3)
#' kl_dist_cont(pdf_1, pdf_2, lower = -Inf, upper = Inf)
#'
#' @export
#'
#' @rdname kl_dist
#' @name kl_dist
kl_dist <- function(p, q, quiet = FALSE){

  # Check equal length
  len_p <- length(p); len_q <- length(q)
  if (len_p != len_q) {
    msg <- paste0("The length of `p`",
                  " (", len_p, ") ",
                  "differ from the lenght of `q`",
                  " (", len_q, ")")
    if (!quiet) warning(msg)
    return(NA)
  }

  # Remove NA-values
  if (any(is.na(p)) | any(is.na(q))) {
    idx_non_na <- which(!is.na(q) & !is.na(p))
    q <- q[idx_non_na]
    p <- p[idx_non_na]
    msg <- paste0("Some elements in `q` or `p` are NA. Removed ", len_p - length(idx_non_na), " elements.")
    if (!quiet) warning(msg)
  }

  # Remove non-zero probabilities
  len_p <- length(p)
  if (any(p == 0) | any(q == 0)) {
    idx_nonzero <- which(q != 0 & p != 0)
    q <- q[idx_nonzero]
    p <- p[idx_nonzero]
    msg <- paste0("Some elements in `p` or `q` are zero. Removed ", len_p - length(idx_nonzero), " elements.")
    if (!quiet) warning(msg)
  }
  p <- p/sum(p)
  q <- q/sum(q)
  sum(p*log(p/q))
}

#' @export
#' @rdname kl_dist
kl_dist_cont <- function(pdf_1, pdf_2, lower = -Inf, upper = Inf, quiet = FALSE){
  dist <- function(x){
    ifelse(pdf_1(x) == 0 | pdf_2(x) == 0, 0, pdf_1(x)*log(pdf_1(x)/pdf_2(x)))
  }
  integrate(dist, lower = lower, upper = upper)$value
}

