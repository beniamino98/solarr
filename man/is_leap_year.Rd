% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\name{is_leap_year}
\alias{is_leap_year}
\title{Is leap year?}
\usage{
is_leap_year(x)
}
\arguments{
\item{x}{numeric value or dates vector in the format \code{YYYY-MM-DD}.}
}
\value{
Boolean. \code{TRUE} if it is a leap year, \code{FALSE} otherwise.
}
\description{
Check if a given year is leap (366 days) or not (365 days).
}
\examples{
is_leap_year("2024-02-01")
is_leap_year(c(2023:2030))
is_leap_year(c("2024-10-01", "2025-10-01"))
is_leap_year("2029-02-01")
}
