#' Make a matrix positive semidefined
#'
#' The matrix is decomposed using a spectral decomposition.
#' Then the negative eigenvalues are imputed with `neg_values` and the original matrix is constructed again.
#'
#' @param x symmetric matric
#' @param neg_values numeric
#' @keywords spatialCorrelations internals
#' @noRd
#' @export
makeSemiPositive <- function(x, neg_values = 1e-5){

  mat <- x
  dec <- eigen(x)
  e <- dec$vectors
  lam <- dec$values
  if (any(lam < 0)){
    lam[lam < 0] <- neg_values
    mat <- e %*% diag(lam) %*% t(e)
  }
  attr(mat, "index_neg_values") <- which(dec$values < 0)
  attr(mat, "original_values") <- dec$values[dec$values < 0]
  return(mat)
}

#' Check admissibility of common probabilities matrix
#'
#' @param commonprob symmetric matrix of joint probabilities.
#' @param check.commonprob logical
#' @param simulvals `bindata::SimulVals`
#' @keywords spatialCorrelations internals
#' @noRd
#' @export
interpolate_commonprob <- function(commonprob, check.commonprob = TRUE, simulvals = NULL){

  if (is.null(simulvals)) {
    simulvals <- bindata::SimulVals
  }

  # Check marginal probabilities
  if (check.commonprob) {
    check <- bindata::check.commonprob(commonprob)
    if (!check) {
      cat(attr(check, "message"), sep = "\n")
      message("Matrix commonprob not admissible. Run check.commonprob(commonprob) for more details.")
    }
  }
  margprob <- diag(commonprob)
  for (m in 1:(ncol(commonprob) - 1)) {
    for (n in (m + 1):nrow(commonprob)) {
      x <- cbind(margprob[m], margprob[n], as.numeric(dimnames(simulvals)[[3]]))
      y <- e1071::interpolate(x, simulvals)
      # Check admissibility and impute errors
      if (commonprob[m, n] > max(y)) {
        print(paste0("Error: (", n, " ", m, ")", " substitued from ", commonprob[m, n], " to ", max(y)))
        commonprob[m, n] <- max(y)
      } else if (commonprob[m, n] < min(y)){
        print(paste0("Error: (", n, " ", m, ")", " substitued from ", commonprob[m, n], " to ", min(y)))
        commonprob[m, n] <- min(y)
      }
    }
    message(m, "/", (ncol(commonprob) - 1), "\r", appendLF = FALSE)
  }
  # Check if any elements is NAs
  if (any(is.na(commonprob))) {
    stop("Some elements in commonprob are NAs... margprob and commonprob are not compatible?")
  }
  return(commonprob)
}
