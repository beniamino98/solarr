# Compute the spectral distribution for a blackbody (sun)
spectralDistribution <- function(lambda = NULL, measure = c("nanometer", "micrometer")){

  # choose the measure: nanometer or mircrometer
  measure <- match.arg(measure, choices = c("nanometer", "micrometer"))

  Ts <- 5777          # constant, black-body temperature in (kelvin)
  Rs <- 6.955*10^(8)  # constant, radius of the sun in (meter)
  R  <- 1.495*10^(11) # constant, sun-earth distance in (meter)
  C1 <- 3.742*10^(8)  # constant, in (W micro-meter^4)/meter^2
  C2 <- 1.439*10^(4)  # constant, in (micro-meter*Kelvin)

  # if "lambda" is "nanometer", the final energy will be in (W/m2 x nano-meter)
  if (measure == "nanometer") {
    lambda <- lambda/1000  # from micrometer to nanometer
    spectra <- (Rs/R)^(2)*(C1/(lambda^(5)*(exp(C2/lambda*Ts) - 1))) # Plank's Law
    return(spectra/1000)
    # if "lambda" is "micrometer", the final energy will be in (W/m2 x micro-meter)
  } else if(measure == "micrometer") {
    spectra <- (Rs/R)^(2)*(C1/(lambda^(5)*(exp(C2/lambda*Ts) - 1))) # Plank's Law
    return(spectra)
  }
}
