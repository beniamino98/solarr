% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solarModel_internals.R
\name{solarModel_forecaster_plot}
\alias{solarModel_forecaster_plot}
\title{Plot a forecast from a solarModel object}
\usage{
solarModel_forecaster_plot(
  model,
  date = "2021-05-29",
  ci = 0.1,
  type = "mix",
  unconditional = FALSE
)
}
\description{
Plot a forecast from a solarModel object
}
\examples{
model <- Bologna
day_date <- "2013-01-13"
solarModel_forecaster_plot(model, date = day_date)
solarModel_forecaster_plot(model, date = day_date, unconditional = TRUE)
solarModel_forecaster_plot(model, date = day_date, type = "dw")
solarModel_forecaster_plot(model, date = day_date, type = "dw", unconditional = TRUE)
solarModel_forecaster_plot(model, date = day_date, type = "up")
solarModel_forecaster_plot(model, date = day_date, type = "up", unconditional = TRUE)
}
