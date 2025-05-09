% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ks_test.R
\name{ks_test}
\alias{ks_test}
\alias{ks_ts_test}
\title{Kolmogorov Smirnov test for a distribution}
\usage{
ks_test(
  x,
  cdf,
  ci = 0.05,
  min_quantile = 0.015,
  max_quantile = 0.985,
  k = 1000,
  plot = FALSE
)
}
\arguments{
\item{x}{a vector.}

\item{ci}{p.value for rejection.}

\item{min_quantile}{minimum quantile for the grid of values.}

\item{max_quantile}{maximum quantile for the grid of values.}

\item{k}{finite value for approximation of infinite sum.}

\item{plot}{when \code{TRUE} a plot is returned, otherwise a \code{tibble}.}

\item{pdf}{a function. The theoric density to use for comparison.}

\item{seed}{random seed for two sample test.}
}
\value{
when \code{plot = TRUE} a plot is returned, otherwise a \code{tibble}.
}
\description{
Test against a specific distribution
}
