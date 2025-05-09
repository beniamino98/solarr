% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PDF.R
\name{PDF}
\alias{PDF}
\alias{CDF}
\alias{Quantile}
\title{Density, distribution and quantile function}
\usage{
PDF(.f, ...)

CDF(.f, lower = -Inf, ...)

Quantile(cdf, interval = c(-100, 100))
}
\arguments{
\item{.f}{density function}

\item{...}{other parameters to be passed to \code{.f}.}

\item{lower}{lower bound for integration (CDF).}

\item{cdf}{cumulative distribution function.}

\item{interval}{lower and upper bounds for unit root (Quantile).}
}
\description{
Return a function of \code{x} given the specification of a function of \code{x}.
}
\examples{
# Density
pdf <- PDF(dnorm, mean = 0.3, sd = 1.3)
pdf(3)
dnorm(3, mean = 0.3, sd = 1.3)
# Distribution
cdf <- CDF(dnorm, mean = 0.3, sd = 1.3)
cdf(3)
pnorm(3, mean = 0.3, sd = 1.3)
# Numeric quantile function
pnorm(Quantile(pnorm)(0.9))
}
