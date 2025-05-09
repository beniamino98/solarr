% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dgumbel.R
\name{dgumbel}
\alias{dgumbel}
\alias{pgumbel}
\alias{qgumbel}
\alias{rgumbel}
\title{Gumbel random variable}
\usage{
dgumbel(x, location = 0, scale = 1, log = FALSE)

pgumbel(q, location = 0, scale = 1, log.p = FALSE, lower.tail = TRUE)

qgumbel(p, location = 0, scale = 1, log.p = FALSE, lower.tail = TRUE)

rgumbel(n, location = 0, scale = 1)
}
\arguments{
\item{x, q}{vector of quantiles.}

\item{location}{location parameter.}

\item{scale}{scale parameter.}

\item{log}{logical; if \code{TRUE}, probabilities are returned as \code{log(p)}.}

\item{log.p}{logical; if \code{TRUE}, probabilities p are given as \code{log(p)}.}

\item{lower.tail}{logical; if TRUE (default), probabilities are \code{P[X < x]} otherwise, \code{P[X > x]}.}

\item{p}{vector of probabilities.}

\item{n}{number of observations. If \code{length(n) > 1}, the length is taken to be the number required.}
}
\description{
Gumbel density, distribution, quantile and random generator.
}
\examples{
# Grid
x <- seq(-5, 5, 0.01)

# Density function
p <- dgumbel(x, location = 0, scale = 1)
plot(x, p, type = "l")

# Distribution function
p <- pgumbel(x, location = 0, scale = 1)
plot(x, p, type = "l")

# Quantile function
qgumbel(0.1)
pgumbel(qgumbel(0.1))

# Random Numbers
rgumbel(1000)
plot(rgumbel(1000), type = "l")

}
