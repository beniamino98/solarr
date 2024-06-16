#' Riccati Square Root
#'
#' Square root of a symmetric matrix.
#'
#' @param x symmetric matrix.
#'
#' @examples
#' x <- matrix(c(1, 0.3, 0.3, 1), nrow = 2, byrow = TRUE)
#' riccati_root(x)

#' @rdname riccati_root
#' @name riccati_root
#' @export

riccati_root <- function(x){
  dec <- eigen(x)
  e <- dec$vectors
  lam <- dec$values
  e %*% diag(sqrt(lam)) %*% t(e)
}
