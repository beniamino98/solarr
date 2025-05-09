% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seasonalClearsky.R
\name{seasonalClearsky}
\alias{seasonalClearsky}
\title{R6 implementation for a clear sky seasonal model}
\description{
R6 implementation for a clear sky seasonal model

R6 implementation for a clear sky seasonal model
}
\section{Super class}{
\code{\link[solarr:seasonalModel]{solarr::seasonalModel}} -> \code{seasonalClearsky}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{control}}{Named list. Control parameters. See the function \code{\link{control_seasonalClearsky}} for details.}

\item{\code{lat}}{Numeric, scalar. Latitude of the location considered.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-seasonalClearsky-new}{\code{seasonalClearsky$new()}}
\item \href{#method-seasonalClearsky-fit}{\code{seasonalClearsky$fit()}}
\item \href{#method-seasonalClearsky-H0}{\code{seasonalClearsky$H0()}}
\item \href{#method-seasonalClearsky-predict}{\code{seasonalClearsky$predict()}}
\item \href{#method-seasonalClearsky-print}{\code{seasonalClearsky$print()}}
\item \href{#method-seasonalClearsky-clone}{\code{seasonalClearsky$clone()}}
}
}
\if{html}{\out{
<details open><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="solarr" data-topic="seasonalModel" data-id="differential"><a href='../../solarr/html/seasonalModel.html#method-seasonalModel-differential'><code>solarr::seasonalModel$differential()</code></a></span></li>
<li><span class="pkg-link" data-pkg="solarr" data-topic="seasonalModel" data-id="fit_differential"><a href='../../solarr/html/seasonalModel.html#method-seasonalModel-fit_differential'><code>solarr::seasonalModel$fit_differential()</code></a></span></li>
<li><span class="pkg-link" data-pkg="solarr" data-topic="seasonalModel" data-id="update"><a href='../../solarr/html/seasonalModel.html#method-seasonalModel-update'><code>solarr::seasonalModel$update()</code></a></span></li>
<li><span class="pkg-link" data-pkg="solarr" data-topic="seasonalModel" data-id="update_std.errors"><a href='../../solarr/html/seasonalModel.html#method-seasonalModel-update_std.errors'><code>solarr::seasonalModel$update_std.errors()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-seasonalClearsky-new"></a>}}
\if{latex}{\out{\hypertarget{method-seasonalClearsky-new}{}}}
\subsection{Method \code{new()}}{
Initialize a \code{seasonalClearsky} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{seasonalClearsky$new(control = control_seasonalClearsky())}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{control}}{Named list. Control parameters. See the function \code{\link{control_seasonalClearsky}} for details.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-seasonalClearsky-fit"></a>}}
\if{latex}{\out{\hypertarget{method-seasonalClearsky-fit}{}}}
\subsection{Method \code{fit()}}{
Fit the seasonal model for clear sky radiation.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{seasonalClearsky$fit(x, date, lat, clearsky)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{Numeric vector. Time series of solar radiation.}

\item{\code{date}}{Character or Date vector. Time series of dates.}

\item{\code{lat}}{Numeric scalar. Reference latitude.}

\item{\code{clearsky}}{Numeric vector. Time series of CAMS clear sky radiation.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-seasonalClearsky-H0"></a>}}
\if{latex}{\out{\hypertarget{method-seasonalClearsky-H0}{}}}
\subsection{Method \code{H0()}}{
Compute the extraterrestrial radiation at a given location.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{seasonalClearsky$H0(n)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{Integer, scalar or vector. Number of day of the year.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-seasonalClearsky-predict"></a>}}
\if{latex}{\out{\hypertarget{method-seasonalClearsky-predict}{}}}
\subsection{Method \code{predict()}}{
Predict method for \code{seasonalClearsky} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{seasonalClearsky$predict(n, newdata)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{Integer, scalar or vector. number of day of the year.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-seasonalClearsky-print"></a>}}
\if{latex}{\out{\hypertarget{method-seasonalClearsky-print}{}}}
\subsection{Method \code{print()}}{
Print method for \code{seasonalClearsky} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{seasonalClearsky$print()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-seasonalClearsky-clone"></a>}}
\if{latex}{\out{\hypertarget{method-seasonalClearsky-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{seasonalClearsky$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
