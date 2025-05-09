% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatialGrid.R
\name{havDistance}
\alias{havDistance}
\title{Haversine distance}
\usage{
havDistance(lat_1, lon_1, lat_2, lon_2)
}
\arguments{
\item{lat_1}{numeric, latitude of first point.}

\item{lon_1}{numeric, longitude of first point.}

\item{lat_2}{numeric, latitude of second point.}

\item{lon_2}{numeric, longitude of second point.}
}
\value{
Numeric vector the distance in kilometers.
}
\description{
Compute the Haversine distance between two points.
}
\examples{
havDistance(43.3, 12.1, 43.4, 12.2)
havDistance(43.35, 12.15, 43.4, 12.2)
}
