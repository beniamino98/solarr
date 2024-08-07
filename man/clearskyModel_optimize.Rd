% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clearskyModel.R
\name{clearskyModel_optimize}
\alias{clearskyModel_optimize}
\title{Optimizer for Solar Clear sky}
\usage{
clearskyModel_optimize(GHI, G0, control = clearskyModel_control())
}
\arguments{
\item{GHI}{vector of solar radiation}

\item{G0}{vector of extraterrestrial radiation}

\item{control}{list of control parameters. See `clearskyModel_control()` for details.}
}
\value{
a numeric vector containing the fitted clear sky radiation.
}
\description{
Find the best parameter delta for fitting clear sky radiation.
}
