% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solarOption.R
\name{solarOption_pricing}
\alias{solarOption_pricing}
\title{Compute the price of a \code{solarOption}}
\usage{
solarOption_pricing(
  moments,
  sorad,
  theta = 0,
  put = TRUE,
  control_options = control_solarOption()
)
}
\arguments{
\item{moments}{description}

\item{sorad}{An object of the class \code{solarOption}.}

\item{put}{Logical. When \code{TRUE}, the default, will be computed the price for a \code{put} contract, otherwise for a \code{call} contract.}

\item{control_options}{Named list. Control parameters, see \code{\link{control_solarOption}} for more details.}
}
\description{
Compute the price of a \code{solarOption}
}
\examples{
model <- Bologna$clone(TRUE)
moments <- filter(model$moments$conditional, Year == 2022)
# Pricing without contracts
solarOption_pricing(moments[1,])
# Pricing with contracts specification
sorad <- solarOption$new()
sorad$set_contract("2021-12-31", "2022-01-01", "2022-04-20", moments$GHI_bar[1])
solarOption_pricing(moments[1,], sorad)
solarOption_pricing(moments[1,], sorad, theta = 0.02)
solarOption_pricing(moments[1,], sorad, theta = -0.02)

}
