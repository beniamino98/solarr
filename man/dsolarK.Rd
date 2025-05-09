% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dsolarK.R
\name{dsolarK}
\alias{dsolarK}
\alias{psolarK}
\alias{qsolarK}
\alias{rsolarK}
\title{Clearness index random variable}
\usage{
dsolarK(x, alpha, beta, pdf_Y, log = FALSE)

psolarK(x, alpha, beta, cdf_Y, log.p = FALSE, lower.tail = TRUE)

qsolarK(p, alpha, beta, cdf_Y, log.p = FALSE, lower.tail = TRUE)

rsolarK(n, alpha, beta, cdf_Y)
}
\arguments{
\item{x}{vector of quantiles.}

\item{alpha}{parameter \code{alpha > 0}.}

\item{beta}{parameter \code{beta > 0} and \code{alpha + beta < 1}.}

\item{pdf_Y}{density function of Y.}

\item{log}{logical; if \code{TRUE}, probabilities are returned as \code{log(p)}.}

\item{cdf_Y}{distribution function of Y.}

\item{log.p}{logical; if \code{TRUE}, probabilities p are given as \code{log(p)}.}

\item{lower.tail}{logical; if \code{TRUE}, the default, the computed probabilities are \code{P[X < x]}. Otherwise, \code{P[X > x]}.}

\item{p}{vector of probabilities.}
}
\description{
Clearness index density, distribution, quantile and random generator.
}
\details{
Consider a random variable \eqn{Y \in [-\infty, \infty]} with a known density function \code{pdf_Y}. Then
the funtion \code{dsolarK} compute the density function of the following transformed random variable, i.e.
\deqn{K(Y) = 1-\alpha-\beta \exp(-\exp(Y))}
where \eqn{K(Y) \in [1-\alpha-\beta, 1-\alpha]}.
}
\examples{
# Parameters
alpha = 0.001
beta = 0.9
# Grid of points
grid <- seq(1-alpha-beta, 1-alpha, length.out = 50)[-50]

# Density
dsolarK(0.4, alpha, beta, function(x) dnorm(x))
dsolarK(0.4, alpha, beta, function(x) dnorm(x, sd = 2))
plot(grid, dsolarK(grid, alpha, beta, function(x) dnorm(x, sd = 0.2)), type="l")

# Distribution
psolarK(0.493, alpha, beta, function(x) pnorm(x))
psolarK(0.493, alpha, beta, function(x) pnorm(x, sd = 2))
plot(grid, psolarK(grid, alpha, beta, function(x) pt(0.2*x, 3)), type="l")
plot(grid, psolarK(grid, alpha, beta, function(x) pnorm(x, sd = 0.2)), type="l")

# Quantile
qsolarK(c(0.05, 0.95), alpha, beta, function(x) pnorm(x))
qsolarK(c(0.05, 0.95), alpha, beta, function(x) pnorm(x, sd = 2))

# Random generator (I)
Kt <- rsolarK(366, alpha, beta, function(x) pnorm(x, sd = 1.3))
plot(1:366, Kt, type="l")

# Random generator (II)
pdf <- function(x) pmixnorm(x, c(-1.8, 0.8), c(0.5, 0.7), c(0.6, 0.4))
Kt <- rsolarK(36, alpha, beta, pdf)
plot(1:36, Kt, type="l")
}
