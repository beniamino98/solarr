% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kl_dist.R
\name{kl_dist}
\alias{kl_dist}
\alias{kl_dist_cont}
\title{Kullback–Leibler divergence}
\usage{
kl_dist(p, q, quiet = FALSE)

kl_dist_cont(pdf_1, pdf_2, lower = -Inf, upper = Inf, quiet = FALSE)
}
\arguments{
\item{p}{Numeric, probability vector. Usually, the empiric probabilities.}

\item{q}{Numeric, probability vector. Usually, the model probabilities.}

\item{quiet}{Boolean, default is `TRUE`.
When set to `FALSE` the function will not display warnings.}
}
\description{
Compute the Kullback–Leibler distance between two probability measure.
}
\details{
The function implements:
\deqn{\sum_{i} p_i \log(\frac{p_i}{q_i}) \quad p_i, q_i > 0 \; \forall i}
}
\examples{

p <- dnorm(rnorm(100))
q <- dnorm(rnorm(100))
kl_dist(p, q)

pdf_1 <- function(x) dnorm(x, mean = 2, sd = 1)
pdf_2 <- function(x) dnorm(x, mean = -2, sd = 3)
kl_dist_cont(pdf_1, pdf_2, lower = -Inf, upper = Inf)

}
\references{
https://en.wikipedia.org/wiki/Kullback–Leibler_divergence
}
