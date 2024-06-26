% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dgumbel.R
\name{gumbel}
\alias{gumbel}
\alias{dgumbel}
\alias{pgumbel}
\alias{qgumbel}
\alias{rgumbel}
\title{Gumbel Random Variable}
\usage{
dgumbel(x, mean = 0, scale = 1, log.p = FALSE, invert = FALSE)

pgumbel(
  x,
  mean = 0,
  scale = 1,
  log.p = FALSE,
  lower.tail = TRUE,
  invert = FALSE
)

qgumbel(
  p,
  mean = 0,
  scale = 1,
  log.p = FALSE,
  lower.tail = TRUE,
  invert = FALSE
)

rgumbel(n, mean = 0, scale = 1, invert = FALSE)
}
\arguments{
\item{x}{vector of quantiles.}

\item{mean}{vector of means.}

\item{scale}{vector of scale parameter.}

\item{log.p}{logical; if `TRUE`, probabilities p are given as `log(p)`.}

\item{invert}{logical, use the inverted Gumbel distribution}

\item{lower.tail}{logical; if TRUE (default), probabilities are `P[X < x]` otherwise, `P[X > x]`.}

\item{p}{vector of probabilities.}

\item{n}{number of observations. If `length(n) > 1`, the length is taken to be the number required.}
}
\description{
Probability density function for a gumbel random variable
}
\examples{
x <- seq(-5, 5, 0.01)

# Density function
p <- dgumbel(x, mean = 0, scale = 1)
plot(x, p, type = "l")

# Distribution function
p <- pgumbel(x, mean = 0, scale = 1)
plot(x, p, type = "l")

# Quantile function
qgumbel(0.1)
pgumbel(qgumbel(0.1))

# Random Numbers
rgumbel(1000)
plot(rgumbel(1000), type = "l")

}
