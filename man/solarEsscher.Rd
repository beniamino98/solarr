% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solarEsscher.R
\name{solarEsscher}
\alias{solarEsscher}
\title{Calibrate Esscher Bounds and parameters}
\description{
Calibrate Esscher Bounds and parameters

Calibrate Esscher Bounds and parameters
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{control}}{list containing the control parameters}

\item{\code{grid}}{list containing the grids}
}
\if{html}{\out{</div>}}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{bounds}}{calibrated bounds with respect to bootstrapped payoff.}

\item{\code{models}}{models to predict the optimal theta given the expected return.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-solarEsscher-new}{\code{solarEsscher$new()}}
\item \href{#method-solarEsscher-calibrator}{\code{solarEsscher$calibrator()}}
\item \href{#method-solarEsscher-calibrate_bounds}{\code{solarEsscher$calibrate_bounds()}}
\item \href{#method-solarEsscher-create_grid}{\code{solarEsscher$create_grid()}}
\item \href{#method-solarEsscher-fit_theta}{\code{solarEsscher$fit_theta()}}
\item \href{#method-solarEsscher-predict}{\code{solarEsscher$predict()}}
\item \href{#method-solarEsscher-clone}{\code{solarEsscher$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-solarEsscher-new"></a>}}
\if{latex}{\out{\hypertarget{method-solarEsscher-new}{}}}
\subsection{Method \code{new()}}{
Initialize the settings for calibration of Esscher parameter.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{solarEsscher$new(
  n_key_points = 15,
  init_lambda = 0,
  lower_lambda = -1,
  upper_lambda = 1,
  put = TRUE,
  target.Yt = TRUE,
  quiet = FALSE,
  control_options = control_solarOption()
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n_key_points}}{integer, number of key points for interpolation.}

\item{\code{init_lambda}}{numeric, initial value for the Esscher parameter.}

\item{\code{lower_lambda}}{numeric, lower value for the Esscher parameter.}

\item{\code{upper_lambda}}{numeric, upper value for the Esscher parameter.}

\item{\code{put}}{logical, when \code{TRUE} will be considered a put contract otherwise a call contract.}

\item{\code{target.Yt}}{logical, when \code{TRUE} will be distorted with esscher parameter the pdf of Yt otherwise the pdf of the GHI.}

\item{\code{quiet}}{logical}

\item{\code{control_options}}{control function. See \code{\link{control_solarOption}} for details.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-solarEsscher-calibrator"></a>}}
\if{latex}{\out{\hypertarget{method-solarEsscher-calibrator}{}}}
\subsection{Method \code{calibrator()}}{
Calibrate the optimal Esscher parameter given a target price
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{solarEsscher$calibrator(model, target_price, nmonths = 1:12, target.Yt)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{model}}{solar model}

\item{\code{target_price}}{the \code{target_price} represent the model price under the target Q-measure.}

\item{\code{nmonths}}{month or months}

\item{\code{target.Yt}}{logical, when \code{TRUE} will be distorted with esscher parameter the pdf of Yt otherwise the pdf of the GHI.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-solarEsscher-calibrate_bounds"></a>}}
\if{latex}{\out{\hypertarget{method-solarEsscher-calibrate_bounds}{}}}
\subsection{Method \code{calibrate_bounds()}}{
Calibrate Esscher upper and lower bounds
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{solarEsscher$calibrate_bounds(model, payoffs, target.Yt)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{model}}{object with the class \code{solarModel}. See the function \code{\link{solarModel}} for details.}

\item{\code{payoffs}}{object with the class \code{solarOptionPayoffs}. See the function \code{\link{solarOptionPayoffs}} for details.}

\item{\code{target.Yt}}{logical, when \code{TRUE} will be distorted with esscher parameter the pdf of Yt otherwise the pdf of the GHI.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-solarEsscher-create_grid"></a>}}
\if{latex}{\out{\hypertarget{method-solarEsscher-create_grid}{}}}
\subsection{Method \code{create_grid()}}{
Create a grid of optimal theta and expected returns with respect of the benchmark price.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{solarEsscher$create_grid(
  model,
  benchmark_price,
  lower_price,
  upper_price,
  target.Yt
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{model}}{object with the class \code{solarModel}. See the function \code{\link{solarModel}} for details.}

\item{\code{benchmark_price}}{benchmark price for an expected return equal to zero.}

\item{\code{lower_price}}{lower price in the grid.}

\item{\code{upper_price}}{upper price in the grid.}

\item{\code{target.Yt}}{logical, when \code{TRUE} will be distorted with esscher parameter the pdf of Yt otherwise the pdf of the GHI.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-solarEsscher-fit_theta"></a>}}
\if{latex}{\out{\hypertarget{method-solarEsscher-fit_theta}{}}}
\subsection{Method \code{fit_theta()}}{
Fit the models to predict the optimal Esscher parameters given the grid.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{solarEsscher$fit_theta()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-solarEsscher-predict"></a>}}
\if{latex}{\out{\hypertarget{method-solarEsscher-predict}{}}}
\subsection{Method \code{predict()}}{
Predict the optimal Esscher parameters given a certain level of expected return.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{solarEsscher$predict(r, target.Yt = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{r}}{expected return}

\item{\code{target.Yt}}{logical, when \code{TRUE} will be distorted with esscher parameter the pdf of Yt otherwise the pdf of the GHI.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-solarEsscher-clone"></a>}}
\if{latex}{\out{\hypertarget{method-solarEsscher-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{solarEsscher$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
