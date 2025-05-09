% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solarScenario.R
\name{solarScenario_filter}
\alias{solarScenario_filter}
\title{Simulate trajectories from a a \code{solarScenario_spec}}
\usage{
solarScenario_filter(simSpec)
}
\arguments{
\item{simSpec}{object with the class \code{solarScenario_spec}. See the function \code{\link{solarScenario_spec}} for details.}
}
\description{
Simulate trajectories from a a \code{solarScenario_spec}
}
\examples{
model <- Bologna
simSpec <- solarScenario_spec(model, from = "2023-01-01", to = "2023-12-31")
simSpec <- solarScenario_residuals(simSpec, nsim = 1, seed = 1)
simSpec <- solarScenario_filter(simSpec)
# Empiric data
df_emp <- simSpec$emp
# First simulation
df_sim <- simSpec$simulations[[1]]
ggplot()+
geom_line(data = df_emp, aes(date, GHI))+
geom_line(data = df_sim, aes(date, GHI), color = "red")

}
