% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatialGrid.R
\name{spatialGrid}
\alias{spatialGrid}
\title{Spatial Grid}
\value{
a tibble with two columns \code{lat} and \code{lon}.
}
\description{
Spatial Grid

Spatial Grid
}
\details{
Create a grid from a range of latitudes and longitudes.
}
\examples{
grid <- spatialGrid$new()
grid$make_grid()
grid$grid
grid$weights <- IDW(beta = 2)
grid$is_known_point(49.95, 12.15)
grid$known_point(c(44.85, 44.9), c(12.15, 11.2))
grid$is_inside_bounds(44.8, 10.9)
grid$neighborhoods(44.9, 12.1)
grid$neighborhoods(44.95, 12.15)
filter(grid$grid, lat == 44.95 & lon == 12.15)

}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-spatialGrid-new}{\code{spatialGrid$new()}}
\item \href{#method-spatialGrid-make_grid}{\code{spatialGrid$make_grid()}}
\item \href{#method-spatialGrid-is_inside_bounds}{\code{spatialGrid$is_inside_bounds()}}
\item \href{#method-spatialGrid-is_known_point}{\code{spatialGrid$is_known_point()}}
\item \href{#method-spatialGrid-known_point}{\code{spatialGrid$known_point()}}
\item \href{#method-spatialGrid-neighborhoods}{\code{spatialGrid$neighborhoods()}}
\item \href{#method-spatialGrid-clone}{\code{spatialGrid$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-spatialGrid-new"></a>}}
\if{latex}{\out{\hypertarget{method-spatialGrid-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{spatialGrid$new(
  lat_min = 43.7,
  lat_max = 45.1,
  lon_min = 9.2,
  lon_max = 12.7,
  lat_by = 0.1,
  lon_by = 0.1,
  weights = IDW(2)
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-spatialGrid-make_grid"></a>}}
\if{latex}{\out{\hypertarget{method-spatialGrid-make_grid}{}}}
\subsection{Method \code{make_grid()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{spatialGrid$make_grid(labels)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-spatialGrid-is_inside_bounds"></a>}}
\if{latex}{\out{\hypertarget{method-spatialGrid-is_inside_bounds}{}}}
\subsection{Method \code{is_inside_bounds()}}{
Check if a point is inside the bounds of the spatial grid.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{spatialGrid$is_inside_bounds(lat, lon)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{lat}}{numeric, latitude of a location.}

\item{\code{lon}}{numeric, longitude of a location.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
\code{TRUE} when the point is inside the limits and \code{FALSE} otherwise.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-spatialGrid-is_known_point"></a>}}
\if{latex}{\out{\hypertarget{method-spatialGrid-is_known_point}{}}}
\subsection{Method \code{is_known_point()}}{
Check if a point is already in the spatial grid
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{spatialGrid$is_known_point(lat, lon)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{lat}}{numeric, latitude of a location.}

\item{\code{lon}}{numeric, longitude of a location.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
\code{TRUE} when the point is a known point and \code{FALSE} otherwise.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-spatialGrid-known_point"></a>}}
\if{latex}{\out{\hypertarget{method-spatialGrid-known_point}{}}}
\subsection{Method \code{known_point()}}{
Return the ID and coordinates of a point that is already in the spatial grid
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{spatialGrid$known_point(lat, lon)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{lat}}{numeric, latitude of a location.}

\item{\code{lon}}{numeric, longitude of a location.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-spatialGrid-neighborhoods"></a>}}
\if{latex}{\out{\hypertarget{method-spatialGrid-neighborhoods}{}}}
\subsection{Method \code{neighborhoods()}}{
Find the n-closest neighborhoods of a point
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{spatialGrid$neighborhoods(lat, lon, n = 4)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{lat}}{numeric, latitude of a point in the grid.}

\item{\code{lon}}{numeric, longitude of a point in the grid.}

\item{\code{n}}{number of neighborhoods}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-spatialGrid-clone"></a>}}
\if{latex}{\out{\hypertarget{method-spatialGrid-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{spatialGrid$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
