% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dsn.R
\name{dsnorm}
\alias{dsnorm}
\alias{psnorm}
\alias{qsnorm}
\alias{rsnorm}
\title{Skewed Normal random variable}
\usage{
dsnorm(x, location = 0, scale = 1, shape = 0, log = FALSE)

psnorm(q, location = 0, scale = 1, shape = 0, log.p = FALSE, lower.tail = TRUE)

qsnorm(p, location = 0, scale = 1, shape = 0, log.p = FALSE, lower.tail = TRUE)

rsnorm(n, location = 0, scale = 1, shape = 0)
}
\arguments{
\item{x, q}{vector of quantiles.}

\item{location}{location parameter.}

\item{scale}{scale parameter.}

\item{shape}{skewness parameter.}

\item{log}{logical; if \code{TRUE}, probabilities are returned as \code{log(p)}.}

\item{log.p}{logical; if \code{TRUE}, probabilities p are given as \code{log(p)}.}

\item{lower.tail}{logical; if TRUE (default), probabilities are \code{P[X < x]} otherwise, \code{P[X > x]}.}

\item{p}{vector of probabilities.}

\item{n}{number of observations. If \code{length(n) > 1}, the length is taken to be the number required.}
}
\description{
Skewed Normal density, distribution, quantile and random generator.
}
\examples{
# Grid of points
x <- seq(-5, 5, 0.01)

# Density function
# right tailed
plot(x, dsnorm(x, shape = 1.9), type = "l")
# left tailed
plot(x, dsnorm(x, shape = -1.9), type = "l")

# Distribution function
plot(x, psnorm(x, shape = 4.9), type = "l")
plot(x, psnorm(x, shape = -4.9), type = "l")

# Quantile function
dsnorm(0.1, shape = 4.9)
dsnorm(0.1, shape = -4.9)
psnorm(qsnorm(0.9, shape = 3), shape = 3)

# Random generator
set.seed(1)
plot(rsnorm(100, shape = 5), type = "l")

}
