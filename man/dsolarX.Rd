% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dsolarX.R
\name{dsolarX}
\alias{dsolarX}
\alias{psolarX}
\alias{qsolarX}
\alias{rsolarX}
\title{Solar risk driver random variable}
\usage{
dsolarX(x, alpha, beta, pdf_Y, log = FALSE)

psolarX(x, alpha, beta, cdf_Y, log.p = FALSE, lower.tail = TRUE)

qsolarX(p, alpha, beta, cdf_Y, log.p = FALSE, lower.tail = TRUE)

rsolarX(n, alpha, beta, cdf_Y)
}
\arguments{
\item{x}{vector of quantiles.}

\item{alpha}{parameter \code{alpha > 0}.}

\item{beta}{parameter \code{beta > 0} and \code{alpha + beta < 1}.}

\item{pdf_Y}{density of Y.}

\item{log}{logical; if \code{TRUE}, probabilities are returned as \code{log(p)}.}

\item{cdf_Y}{distribution function of Y.}

\item{log.p}{logical; if \code{TRUE}, probabilities p are given as \code{log(p)}.}

\item{lower.tail}{logical; if \code{TRUE}, the default, the computed probabilities are \code{P[X < x]}. Otherwise, \code{P[X > x]}.}

\item{p}{vector of probabilities.}
}
\description{
Solar risk driver density, distribution, quantile and random generator.
}
\details{
Consider a random variable \eqn{Y \in [-\infty, \infty]} with a known density function \code{pdf_Y}. Then
the funtion \code{dsolarX} compute the density function of the following transformed random variable, i.e.
\deqn{X(Y) = \alpha+\beta \exp(-\exp(Y))}
where \eqn{X(Y) \in [\alpha, \alpha+\beta]}.
}
\examples{
# Parameters
alpha = 0.001
beta = 0.9
# Grid of points
grid <- seq(alpha, alpha+beta, length.out = 50)[-50]

# Density
dsolarX(0.4, alpha, beta, function(x) dnorm(x))
dsolarX(0.4, alpha, beta, function(x) dnorm(x, sd = 2))
plot(grid, dsolarX(grid, alpha, beta, function(x) dnorm(x, sd = 0.2)), type="l")

# Distribution
psolarX(0.493, alpha, beta, function(x) pnorm(x))
dsolarX(0.493, alpha, beta, function(x) pnorm(x, sd = 2))
plot(grid, psolarX(grid, alpha, beta, function(x) pnorm(x, sd = 0.2)), type="l")

# Quantile
qsolarX(c(0.05, 0.95), alpha, beta, function(x) pnorm(x))
qsolarX(c(0.05, 0.95), alpha, beta, function(x) pnorm(x, sd = 1.3))

# Random generator (I)
set.seed(1)
Kt <- rsolarX(366, alpha, beta, function(x) pnorm(x, sd = 0.8))
plot(1:366, Kt, type="l")

# Random generator (II)
cdf <- function(x) pmixnorm(x, c(-1.8, 0.9), c(0.5, 0.7), c(0.6, 0.4))
Kt <- rsolarX(366, alpha, beta, cdf)
plot(1:366, Kt, type="l")
}
