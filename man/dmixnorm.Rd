% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dmixnorm.R
\name{dmixnorm}
\alias{dmixnorm}
\alias{pmixnorm}
\alias{qmixnorm}
\alias{rmixnorm}
\title{Gaussian mixture random variable}
\usage{
dmixnorm(x, mean = rep(0, 2), sd = rep(1, 2), alpha = rep(1/2, 2), log = FALSE)

pmixnorm(
  q,
  mean = rep(0, 2),
  sd = rep(1, 2),
  alpha = rep(1/2, 2),
  lower.tail = TRUE,
  log.p = FALSE
)

qmixnorm(
  p,
  mean = rep(0, 2),
  sd = rep(1, 2),
  alpha = rep(1/2, 2),
  lower.tail = TRUE,
  log.p = FALSE
)

rmixnorm(n, mean = rep(0, 3), sd = rep(1, 3), alpha = rep(1/3, 3))
}
\arguments{
\item{x, q}{vector of quantiles.}

\item{mean}{vector of means parameters.}

\item{sd}{vector of std. deviation parameters.}

\item{alpha}{vector of probability parameters for each component.}

\item{log}{logical; if \code{TRUE}, probabilities are returned as \code{log(p)}.}

\item{lower.tail}{logical; if TRUE (default), probabilities are \code{P[X < x]} otherwise, \code{P[X > x]}.}

\item{log.p}{logical; if \code{TRUE}, probabilities p are given as \code{log(p)}.}

\item{p}{vector of probabilities.}

\item{n}{number of observations. If \code{length(n) > 1}, the length is taken to be the number required.}
}
\description{
Gaussian mixture density, distribution, quantile and random generator.
}
\examples{
# Parameters
mean = c(-3,0,3)
sd = rep(1, 3)
alpha = c(0.2, 0.3, 0.5)
# Density function
dmixnorm(3, mean, sd, alpha)
# Distribution function
dmixnorm(c(1.2, -3), mean, sd, alpha)
# Quantile function
qmixnorm(0.2, mean, sd, alpha)
# Random generator
rmixnorm(1000, mean, sd, alpha)

}
