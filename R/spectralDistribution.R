#' Compute the spectral distribution for a black body
#'
#' @param x Numeric vector. Wavelength values. Values are interpreted in
#'   nanometers when `measure = "nanometer"` and in micrometers when
#'   `measure = "micrometer"`.
#' @param measure Character. Unit used for the wavelength scale and returned
#'   spectral density. One of `"nanometer"` or `"micrometer"`.
#'
#' @return A numeric vector of spectral density values.
#'
#' @examples
#' spectralDistribution(c(500, 600))
#' spectralDistribution(c(0.5, 0.6), measure = "micrometer")
#'
#' @rdname spectralDistribution
#' @name spectralDistribution
#'
#' @export
spectralDistribution <- function(x, measure = "nanometer"){

  # choose the measure: nanometer or mircrometer
  measure <- match.arg(measure, choices = c("nanometer", "micrometer"))

  Ts <- 5777          # constant, black-body temperature in (kelvin)
  Rs <- 6.955*10^(8)  # constant, radius of the sun in (meter)
  R  <- 1.495*10^(11) # constant, sun-earth distance in (meter)
  C1 <- 3.742*10^(8)  # constant, in (W micro-meter^4)/meter^2
  C2 <- 1.439*10^(4)  # constant, in (micro-meter*Kelvin)

  # if "lambda" is "nanometer", the final energy will be in (W/m2 x nanometer)
  if (measure == "nanometer") {
    # from micrometer to nanometer
    lambda <- x/1000
    spectra <- (Rs/R)^(2)*(C1/(lambda^(5)*(exp(C2/lambda*Ts) - 1))) # Plank's Law
    spectra <- spectra/1000
  } else if (measure == "micrometer") {
    # if "lambda" is "micrometer", the final energy will be in (W/m2 x micrometer)
    spectra <- (Rs/R)^(2)*(C1/(x^(5)*(exp(C2/x*Ts) - 1))) # Plank's Law
  }
  return(spectra)
}
