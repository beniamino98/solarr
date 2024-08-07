% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solarFunctions.R
\name{from_radiant_to_degree}
\alias{from_radiant_to_degree}
\alias{from_degree_to_radiant}
\title{Conversion in Radiant or Degrees}
\usage{
from_radiant_to_degree(x)

from_degree_to_radiant(x)
}
\arguments{
\item{x}{numeric vector, angles in radiant or degrees.}
}
\value{
numeric vector.
}
\description{
Convert an angle in radiant into an angle in degrees.
}
\examples{
# convert 0.34 radiant in degrees
from_radiant_to_degree(0.34)
# convert 19.48 degree in radiant
from_degree_to_radiant(19.48)

}
