% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solarFunctions.R
\name{number_of_day}
\alias{number_of_day}
\title{Number of Day}
\usage{
number_of_day(day_date = NULL)
}
\arguments{
\item{day_date}{dates vector in the format `%YYYY-%MM-%DD`.}
}
\value{
Numeric vector with the number of the day during the year. Can vary from `1` up to `365` or `366`.
}
\description{
Compute the number of day of the year given a vector of dates.
}
\examples{
# detect the number of the day in 2040-01-31
number_of_day("2040-01-31")

# detect the number of the day for a vector of dates
number_of_day(c("2040-01-31", "2023-04-01", "2015-09-02"))

}
