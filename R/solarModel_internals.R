#' Control parameters for a `solarModel` object
#'
#' Control function for a solarModel
#'
#' @param clearsky list with control parameters for the clear sky seasonal model. See \code{\link{control_seasonalClearsky}} for details.
#' @param seasonal.mean a list of parameters. Available choices are:
#'\describe{
#'  \item{`seasonalOrder`}{An integer specifying the order of the seasonal component. The default is `1`.}
#'  \item{`include.intercept`}{Logical, when `TRUE` the intercept will be included in the seasonal model, otherwise will be omitted. The dafault if `TRUE`.}
#'  \item{`include.H0`}{Logical, when `TRUE` the extraterrestrial radiation will be included in the seasonal model, otherwise will be omitted. The dafault if `FALSE`.}
#'  \item{`monthly.mean`}{Logical, when `TRUE` a vector of 12 monthly means is computed on the deseasonalized series and it is subtracted to ensure that it is perfectly centered around zero.}
#'}
#' @param mean.model a list of parameters for the AR model.
#' \describe{
#'  \item{`arOrder`}{An integer specifying the order of the AR component. The default is `2`.}
#'  \item{`include.intercept`}{When `TRUE` an intercept will be included in the AR equation, otherwise is omitted. The dafault if `FALSE`.}
#'}
#' @param seasonal.variance a list of parameters for the seasonal variance. Available choices are:
#'\describe{
#'  \item{`seasonalOrder`}{Integer, it specify the order of the seasonal component in the model. The default is `1`.}
#'  \item{`correction`}{Logical, when `TRUE` the parameters of seasonal variance are corrected to ensure that the standardize the residuals have exactly a unitary variance.}
#'  \item{`monthly.mean`}{Logical, when `TRUE` a vector of 12 monthly std. deviations is computed on the GARCH residuals. Then, they are divided by this quantity to ensure that the monthly variance is one.}
#'}
#' @param variance.model an `ugarchspec` object for GARCH variance. Default is `GARCH(1,1)` specification.
#' @param mixture.model a list of parameters for the monthly Gaussian mixture model. Available choices are:
#'\describe{
#'  \item{`abstol`}{Numeric, absolute level for convergence. The default is `1e-3`.}
#'  \item{`maxit`}{Integer, maximum number of iterations. The default is `200`.}
#'}
#' @param threshold numeric, threshold used to estimate the transformation parameters alpha and beta. See \code{\link{solarTransform}} for details.
#' @param outliers_quantile quantile for outliers detection. If different from 0, the observations that are below the quantile at confidence levels `outliers_quantile` and
#' the observation above the quantile at confidence level 1-`outliers_quantile` will have a weight equal to zero and will be excluded from estimations.
#' @param quiet logical, when `TRUE` the function will not display any message.
#'
#' @examples
#' control <- control_solarModel()
#'
#' @rdname control_solarModel
#' @name control_solarModel
#'
#' @export
control_solarModel <- function(clearsky = control_seasonalClearsky(),
                               stochastic_clearsky = FALSE,
                               seasonal.mean = list(seasonalOrder = 1, include.H0 = FALSE, include.intercept = TRUE, monthly.mean = TRUE),
                               mean.model = list(arOrder = 2, include.intercept = FALSE),
                               seasonal.variance = list(seasonalOrder = 1, correction = TRUE, monthly.mean = TRUE),
                               variance.model = rugarch::ugarchspec(variance.model = list(garchOrder = c(1,1)),
                                                                    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE)),
                               mixture.model = list(abstol = 1e-3, maxit = 150),
                               threshold = 0.01, outliers_quantile = 0, garch_variance = TRUE, quiet = FALSE){

  # Seasonal mean model default parameters
  seasonal_mean = list(seasonalOrder = 1, include.H0 = FALSE, include.intercept = TRUE, monthly.mean = TRUE)
  names_seasonal_mean <- names(seasonal_mean)
  for(name in names_seasonal_mean){
    arg <- seasonal.mean[[name]]
    if (!is.null(arg)) {
      seasonal_mean[[name]] <- seasonal.mean[[name]]
    }
  }

  # Seasonal variance model default parameters
  seasonal_variance = list(seasonalOrder = 1, correction = TRUE,  monthly.mean = TRUE)
  names_seasonal_variance <- names(seasonal_variance)
  for(name in names_seasonal_variance){
    arg <- seasonal.variance[[name]]
    if (!is.null(arg)) {
      seasonal_variance[[name]] <- seasonal.variance[[name]]
    }
  }

  # Mean model default parameters
  mean_model = list(arOrder = 2, include.intercept = FALSE)
  names_mean_model <- names(mean_model)
  for(name in names_mean_model){
    arg <- mean.model[[name]]
    if (!is.null(arg)) {
      mean_model[[name]] <- mean.model[[name]]
    }
  }

  # Variance model default parameters
  mixture_model = list(abstol = 1e-20, maxit = 100, EM = TRUE)
  names_mixture_model <- names(mixture_model)
  for(name in names_mixture_model){
    arg <- mixture.model[[name]]
    if (!is.null(arg)) {
      mixture_model[[name]] <- mixture.model[[name]]
    }
  }

  structure(
    list(
      clearsky = clearsky,
      stochastic_clearsky = stochastic_clearsky,
      seasonal.mean = seasonal_mean,
      seasonal.variance = seasonal_variance,
      mean.model = mean_model,
      variance.model = variance.model,
      mixture.model = mixture_model,
      threshold = threshold,
      outliers_quantile = outliers_quantile,
      garch_variance = garch_variance,
      quiet = quiet
    ),
    class = c("control", "list")
  )
}


#' Specification function for a `solarModel`
#'
#' @param place Character, name of an element in the `CAMS_data` list.
#' @param target Character, target variable to model. Can be `GHI` or `clearsky`.
#' @param min_date Character. Date in the format `YYYY-MM-DD`. Minimum date for the complete data. If `missing` will be used the minimum data available.
#' @param max_date Character. Date in the format `YYYY-MM-DD`. Maximum date for the complete data. If `missing` will be used the maximum data available.
#' @param from Character. Date in the format `YYYY-MM-DD`. Starting date to use for training data.
#' If `missing` will be used the minimum data available after filtering for `min_date`.
#' @param to character. Date in the format `YYYY-MM-DD`. Ending date to use for training data.
#' If `missing` will be used the maximum data available after filtering for `max_date`.
#' @param CAMS_data named list with radiation data for different locations.
#' @param control_model list with control parameters, see \code{\link{control_solarModel}} for details.
#'
#' @examples
#' control <- control_solarModel(outliers_quantile = 0)
#' spec <- solarModel_spec("Bologna", from="2005-01-01", to="2022-01-01", control_model = control)
#'
#' @rdname solarModel_spec
#' @name solarModel_spec
#'
#' @export
solarModel_spec <- function(place, target = "GHI", min_date, max_date, from, to, CAMS_data = solarr::CAMS_data, control_model = control_solarModel()){

  # Match the target variable to model
  target <- match.arg(target, choices = c("GHI", "clearsky"))

  # Match a location in the dataset
  place <- match.arg(place, choices = names(CAMS_data), several.ok = FALSE)
  # Extract CAMS data for the selected location
  data <- CAMS_data[[place]]

  # Minimum date for the complete data
  if (missing(min_date) || is.null(min_date) || is.na(min_date)) {
    min_date <- min(data$date, na.rm = TRUE)
  } else {
    min_date <- as.Date(min_date)
  }
  # Maximum date for the complete data
  if (missing(max_date) || is.null(max_date) || is.na(max_date)) {
    max_date <- max(data$date, na.rm = TRUE)
  } else {
    max_date <- as.Date(max_date)
  }
  # Minimum date for train data
  if (missing(from) || is.null(from) || is.na(from)) {
    from <- min(data$date, na.rm = TRUE)
  } else {
    from <- as.Date(from)
  }
  # Maximum date for train data
  if (missing(to) || is.null(to) || is.na(to)) {
    to <- max(data$date, na.rm = TRUE)
  } else {
    to <- as.Date(to)
  }
  # Filter for min and max dates the complete dataset
  data <- dplyr::filter(data, date >= min_date & date <= max_date)
  # Label for data used for estimation
  data <- dplyr::mutate(data,
                        isTrain = ifelse(date >= from & date <= to, TRUE, FALSE),
                        weights = ifelse(isTrain, 1, 0))
  # Normalize the weights
  data$weights <-  data$weights/sum(data$weights)
  # Move clearsky max radiation
  data$clearsky <- data$clearsky*1.005

  # Train observations and percentage
  nobs_train <- length(data$isTrain[data$isTrain])
  perc_train <- nobs_train/nrow(data)
  # Train observations and percentage
  nobs_test <- length(data$isTrain[!data$isTrain])
  perc_test <- nobs_test/nrow(data)
  # Model dates
  model_dates = list(data = list(from = min_date, to = max_date, nobs = nrow(data), perc = 1),
                     train = list(from = from, to = to, nobs = nobs_train, perc = perc_train),
                     test = list(from = to, to = max_date, nobs = nobs_test, perc = perc_test))

  structure(
    list(
      place = attr(data, "place"),
      coords = attr(data, "coords"),
      data = data,
      dates = model_dates,
      target = target,
      control = control_model
    ),
    class = c("solarModelSpec", "list")
  )
}

#' Compute conditional moments from a `solarModel` object
#'
#' @examples
#' model <- Bologna
#' solarModel_conditional_moments(model)
#' solarModel_conditional_moments(model, date = "2022-01-01")
#'
#' @rdname solarModel_conditional_moments
#' @name solarModel_conditional_moments
#'
#' @export
solarModel_conditional_moments <- function(model, date){

  # Extract complete data
  data <- model$data
  # Filter for a set of dates
  if (!missing(date)) {
    data <- data[data$date %in% as.Date(date),]
  }
  # Compute conditional moments
  data <- dplyr::mutate(data,
                        # Conditional expectation Yt
                        e_Yt = Yt_bar + Yt_tilde_hat + Yt_tilde_uncond,
                        # Conditional std. deviation Yt
                        sd_Yt = sigma*sigma_bar*sigma_m,
                        # Conditional moments Yt (state up)
                        e_Yt_up = e_Yt + sd_Yt*mu1,
                        sd_Yt_up = sd_Yt*sd1,
                        # Conditional moments Yt (state dw)
                        e_Yt_dw = e_Yt + sd_Yt*mu2,
                        sd_Yt_dw = sd_Yt*sd2,
                        # Fitted value depending on B
                        e_Yt_mix = e_Yt +  e_Yt_up*B + e_Yt_dw*(1-B),
                        # Values for target variable
                        e_Rt = model$transform$GHI_y(e_Yt, Ct),
                        e_Rt_mix = model$transform$GHI_y(e_Yt_mix, Ct),
                        e_Rt_up = model$transform$GHI_y(e_Yt_up, Ct),
                        e_Rt_dw = model$transform$GHI_y(e_Yt_dw, Ct),
                        # Compute bounds for Risk driver
                        lower_Xt = model$transform$alpha,
                        upper_Xt = model$transform$alpha + model$transform$beta,
                        # Compute bounds for Clearness index
                        lower_Kt = 1-upper_Xt,
                        upper_Kt = 1-lower_Xt,
                        # Compute bounds for GHI
                        lower_Rt = Ct*lower_Kt,
                        upper_Rt = Ct*upper_Kt)
  # Recoded variables
  new_vars_names <- paste0(paste0("e_", model$target), c("", "_mix", "_up", "_dw"))
  new_vars_names <- c(new_vars_names, paste0("lower_", model$target), paste0("upper_", model$target))
  names(new_vars_names) <- c("e_Rt", "e_Rt_mix", "e_Rt_up", "e_Rt_dw", "lower_Rt", "upper_Rt")
  for(i in 1:length(new_vars_names)){
    data[[new_vars_names[i]]] <- data[[names(new_vars_names[i])]]
  }
  # Add target and seasonal mean
  vars_names <- c(model$target, paste0(model$target, "_bar"), new_vars_names)
  names(vars_names) <- NULL
  # Extract only relevant variables
  data <- dplyr::select(data, date, Month, Day, Ct, Yt, e_Yt, sd_Yt,
                        e_Yt_up, sd_Yt_up, e_Yt_dw, sd_Yt_dw, p1, B,
                        tidyr::any_of(vars_names),
                        lower_Kt, upper_Kt, lower_Xt, upper_Xt)
  return(data)
}


#' Compute conditional moments from a `solarModel` object
#'
#' @examples
#' model <- Bologna
#' solarModel_unconditional_moments(model)
#' solarModel_unconditional_moments(model, nmonths = 1)
#' solarModel_unconditional_moments(model, nmonths = 1, ndays = 1)
#' solarModel_unconditional_moments(model, date = "2022-01-01")
#'
#' @rdname solarModel_unconditional_moments
#' @name solarModel_unconditional_moments
#'
#' @export
solarModel_unconditional_moments <- function(model, nmonths, ndays, date){

  # Seasonal data
  data <- model$seasonal_data
  # Filter for a set of dates
  if (missing(nmonths)){
    nmonths <- 1:12
  }
  # Filter for a set of dates
  if (missing(ndays)){
    ndays <- 1:31
  }
  # Filter for a set of dates
  if (!missing(date)){
    nmonths <- lubridate::month(date)
    ndays <- lubridate::day(date)
  }
  # Filter for a range of dates
  data <- dplyr::filter(data, Month %in% nmonths & Day %in% ndays)
  # Compute conditional moments
  data <- dplyr::mutate(data,
                        # Unconditional expectation Yt
                        e_Yt = Yt_bar + Yt_tilde_uncond,
                        # Unconditional std. deviation Yt
                        sd_Yt = sigma_bar*sigma_m,
                        # Unconditional moments Yt (state up)
                        e_Yt_up = e_Yt + sd_Yt*mu1,
                        sd_Yt_up = sd_Yt*sd1,
                        # Unconditional moments Yt (state dw)
                        e_Yt_dw = e_Yt + sd_Yt*mu2,
                        sd_Yt_dw = sd_Yt*sd2,
                        # Fitted value mixture
                        e_Yt_mix = e_Yt +  e_Yt_up*p1 + e_Yt_dw*(1-p1),
                        # Values for target variable
                        e_Rt = model$transform$GHI_y(e_Yt, Ct),
                        e_Rt_mix = model$transform$GHI_y(e_Yt_mix, Ct),
                        e_Rt_up = model$transform$GHI_y(e_Yt_up, Ct),
                        e_Rt_dw = model$transform$GHI_y(e_Yt_dw, Ct),
                        # Compute bounds for Risk driver
                        lower_Xt = model$transform$alpha,
                        upper_Xt = model$transform$alpha + model$transform$beta,
                        # Compute bounds for Clearness index
                        lower_Kt = 1-upper_Xt,
                        upper_Kt = 1-lower_Xt,
                        # Compute bounds for GHI
                        lower_Rt = Ct*lower_Kt,
                        upper_Rt = Ct*upper_Kt)
  # Recoded variables
  new_vars_names <- paste0(paste0("e_", model$target), c("", "_mix", "_up", "_dw"))
  new_vars_names <- c(new_vars_names, paste0("lower_", model$target), paste0("upper_", model$target))
  names(new_vars_names) <- c("e_Rt", "e_Rt_mix", "e_Rt_up", "e_Rt_dw", "lower_Rt", "upper_Rt")
  for(i in 1:length(new_vars_names)){
    data[[new_vars_names[i]]] <- data[[names(new_vars_names[i])]]
  }
  # Add seasonal mean
  vars_names <- c(paste0(model$target, "_bar"), new_vars_names)
  names(vars_names) <- NULL
  # Extract only relevant variables
  data <- dplyr::select(data, Month, Day, Ct, e_Yt, sd_Yt,
                        e_Yt_up, sd_Yt_up, e_Yt_dw, sd_Yt_dw, p1,
                        tidyr::any_of(vars_names),
                        lower_Kt, upper_Kt, lower_Xt, upper_Xt)
  return(data)
}


#' Produce a forecast from a solarModel object
#'
#' @examples
#' model <- Bologna
#' solarModel_forecaster(model, date = "2010-04-01")
#' object <- solarModel_forecaster(model, date = "2020-04-01", unconditional = TRUE)
#' object
#' @rdname solarModel_forecaster
#' @name solarModel_forecaster
#' @export
solarModel_forecaster <- function(model, date = "2020-01-01", ci = 0.1, unconditional = FALSE){

  message("Forecasting solar radiation for ", date, "\r", appendLF = FALSE)

  if (unconditional) {
    df_n <- solarModel_unconditional_moments(model, date = date)
  } else {
    df_n <- solarModel_conditional_moments(model, date = date)
  }
  # Number of points for the grid
  n_points <- 100
  # Upper and lower bounds
  lower_bound <- df_n[[paste0("lower_", model$target)]]
  upper_bound <- df_n[[paste0("upper_", model$target)]]
  # Grid for PDF plot
  grid_x <- seq(lower_bound, upper_bound, length.out = n_points+2)[-c(1,n_points+2)]
  grid <- dplyr::tibble(x = grid_x, label = as.character(date))
  # Density of Yt
  pdf_Yt <- function(x) dmixnorm(x, mean = c(df_n$e_Yt_up, df_n$e_Yt_dw), sd = c(df_n$sd_Yt_up, df_n$sd_Yt_dw), alpha = c(df_n$p1, 1-df_n$p1))
  # Distribution of Yt
  cdf_Yt <- function(x) pmixnorm(x, mean = c(df_n$e_Yt_up, df_n$e_Yt_dw), sd = c(df_n$sd_Yt_up, df_n$sd_Yt_dw), alpha = c(df_n$p1, 1-df_n$p1))
  # Density of GHI
  pdf_Rt <- function(x, pdf_Yt) dsolarGHI(x, df_n$Ct, model$transform$alpha, model$transform$beta, pdf_Yt)
  # Add grid of points
  grid[[paste0("pdf_", model$target)]] <- pdf_Rt(grid_x, pdf_Yt)
  grid[[paste0("pdf_", model$target, "_up")]] <- pdf_Rt(grid_x, function(x) df_n$p1*dnorm(x, df_n$e_Yt_up, df_n$sd_Yt_up))
  grid[[paste0("pdf_", model$target, "_dw")]] <- pdf_Rt(grid_x, function(x) (1-df_n$p1)*dnorm(x, df_n$e_Yt_dw, df_n$sd_Yt_dw))

  # Empiric variables for the day to be extracted
  emp <- model$data
  emp <- emp[emp$date == date,]
  emp <- dplyr::select(emp, date, Month, Day, tidyr::any_of(model$target), B)
  # Add fitted expected values
  vars_names <- paste0(paste0("e_", model$target), c("", "_mix", "_up", "_dw"))
  emp <- dplyr::bind_cols(emp,  dplyr::select(df_n, any_of(vars_names)))
  # ==========================  Probabilistic Forecasts ==============================
  # Up and down density non weighted
  pdf_Yt_up <- function(x) dnorm(x, df_n$e_Yt_up, df_n$sd_Yt_up)
  pdf_Yt_dw <- function(x) dnorm(x, df_n$e_Yt_dw, df_n$sd_Yt_dw)
  # Up and down distributions non weighted
  cdf_Yt_up <- function(x) pnorm(x, df_n$e_Yt_up, df_n$sd_Yt_up)
  cdf_Yt_dw <- function(x) pnorm(x, df_n$e_Yt_dw, df_n$sd_Yt_dw)
  # Expected value function
  e_Rt <- function(pdf_Yt, .f = function(x) x) integrate(function(x) .f(x)*pdf_Rt(x, pdf_Yt), lower = lower_bound, upper = upper_bound)$value
  # Confidence intervals function
  conf_Rt <- function(p, cdf_Yt) qsolarGHI(p, df_n$Ct, model$transform$alpha, model$transform$beta, cdf_Yt)

  # Expected value and variance
  # Mixture
  emp$e_Rt <- e_Rt(pdf_Yt)
  emp$v_Rt <- e_Rt(pdf_Yt, .f = function(x) x^2) - emp$e_Rt^2
  emp$pdf_e_Rt <- pdf_Rt(emp$e_Rt, pdf_Yt)
  # Sunny state
  emp$e_Rt_up <- e_Rt(pdf_Yt_up)
  emp$v_Rt_up <- e_Rt(pdf_Yt_up, .f = function(x) x^2) - emp$e_Rt_up^2
  # Cloudy state
  emp$e_Rt_dw <- e_Rt(pdf_Yt_dw)
  emp$v_Rt_dw <- e_Rt(pdf_Yt_dw, .f = function(x) x^2) - emp$e_Rt_dw^2

  # Confidence intervals
  # Mixture
  emp$ci_mix_lo <- conf_Rt(ci, cdf_Yt)
  emp$ci_mix_hi <- conf_Rt(1-ci, cdf_Yt)
  emp$pdf_ci_mix_lo <- pdf_Rt(emp$ci_mix_lo, pdf_Yt)
  emp$pdf_ci_mix_hi <- pdf_Rt(emp$ci_mix_hi, pdf_Yt)
  # Sunny state
  emp$ci_up_lo <- conf_Rt(ci, cdf_Yt_up)
  emp$ci_up_hi <- conf_Rt(1-ci, cdf_Yt_up)
  emp$pdf_ci_up_lo <- pdf_Rt(emp$ci_up_lo, pdf_Yt)
  emp$pdf_ci_up_hi <- pdf_Rt(emp$ci_up_hi, pdf_Yt)
  # Cloudy state
  emp$ci_dw_lo <- conf_Rt(ci, cdf_Yt_dw)
  emp$ci_dw_hi <- conf_Rt(1-ci, cdf_Yt_dw)
  emp$pdf_ci_dw_lo <- pdf_Rt(emp$ci_dw_lo, pdf_Yt)
  emp$pdf_ci_dw_hi <- pdf_Rt(emp$ci_dw_hi, pdf_Yt)
  emp$ci <- ci


  structure(
    list(
      grid = grid,
      emp = emp
    ),
    class = c("solarModelForecast", "list")
  )
}


#' Iterate the forecast on multiple dates
#' @examples
#' model <- Bologna
#' dates <- seq.Date(as.Date("2020-01-01"), as.Date("2020-01-31"), 1)
#' solarModel_forecast(model, date = dates)
#' @rdname solarModel_forecast
#' @name solarModel_forecast
#' @export
solarModel_forecast <- function(model, date, ci = 0.1, unconditional = FALSE){
  fun <- function(date){
    safe_forecaster <- purrr::safely(solarModel_forecaster)
    smf <- safe_forecaster(model, date = date, ci = ci, unconditional = unconditional)$result
    if (is.null(smf)){
      return(NULL)
    }
    return(smf[-1][[1]])
  }
  out <- purrr::map_df(date, fun)
  return(out)
}

#' Plot a forecast from a solarModel object
#'
#' @examples
#' model <- Bologna
#' day_date <- "2013-01-13"
#' solarModel_forecaster_plot(model, date = day_date)
#' solarModel_forecaster_plot(model, date = day_date, unconditional = TRUE)
#' solarModel_forecaster_plot(model, date = day_date, type = "dw")
#' solarModel_forecaster_plot(model, date = day_date, type = "dw", unconditional = TRUE)
#' solarModel_forecaster_plot(model, date = day_date, type = "up")
#' solarModel_forecaster_plot(model, date = day_date, type = "up", unconditional = TRUE)
#' @export
solarModel_forecaster_plot <- function(model, date = "2021-05-29", ci = 0.1, type = "mix", unconditional = FALSE){

  object <- solarModel_forecaster(model, date = date, ci = ci, unconditional = unconditional)

  grid <- object$grid
  emp <- object$emp
  type <- match.arg(type, choices = c("mix", "up", "dw"))
  vals_colors <- c(realized = "black", expected_value = "orange",seasonal = "red", bounds = "purple", fitted = "magenta")
  vals_labels <- c(realized = "Realized", expected_value = "Expectation", seasonal = "Seasonal", bounds = "Conf int", fitted = "Fitted")

  pdf_plot <-  ggplot()+
    # Mixture density
    geom_line(data = grid, aes(x, pdf_GHI))+
    # Components densities (weighted by prior probs)
    geom_line(data = grid, aes(x, pdf_GHI_dw), color = "green")+
    geom_line(data = grid, aes(x, pdf_GHI_up), color = "red")

  if (type == "mix") {
    pdf_plot <- pdf_plot +
      geom_segment(data = emp, aes(x = ci_mix_lo, xend = ci_mix_lo, y = 0, yend = pdf_ci_mix_lo, color = "bounds"), linewidth = 0.5)+
      geom_segment(data = emp, aes(x = ci_mix_hi, xend = ci_mix_hi, y = 0, yend = pdf_ci_mix_hi, color = "bounds"), linewidth = 0.5)+
      geom_segment(data = emp, aes(x = ci_mix_lo, xend = ci_mix_hi, y = 0, yend = 0, color = "bounds"))+
      geom_point(data = emp, aes(ci_mix_lo, 0, color = "bounds"), size = 3)+
      geom_point(data = emp, aes(ci_mix_hi, 0, color = "bounds"), size = 3)+
      geom_point(data = emp, aes(e_GHI, 0, color = "fitted"), shape = 17, size = 2)+
      geom_point(data = emp, aes(e_Rt, 0, color = "expected_value"), size = 3)


  } else if (type == "up") {
    pdf_plot <- pdf_plot +
      geom_segment(data = emp, aes(x = ci_up_lo, xend = ci_up_lo, y = 0, yend = pdf_ci_up_lo, color = "bounds"), linewidth = 0.5)+
      geom_segment(data = emp, aes(x = ci_up_hi, xend = ci_up_hi, y = 0, yend = pdf_ci_up_hi, color = "bounds"), linewidth = 0.5)+
      geom_segment(data = emp, aes(x = ci_up_lo, xend = ci_up_hi, y = 0, yend = 0, color = "bounds"))+
      geom_point(data = emp, aes(ci_up_lo, 0, color = "bounds"), size = 3)+
      geom_point(data = emp, aes(ci_up_hi, 0, color = "bounds"), size = 3)+
      geom_point(data = emp, aes(e_GHI_up, 0, color = "fitted"), shape = 17, size = 2)+
      geom_point(data = emp, aes(e_Rt_up, 0, color = "expected_value"), size = 3)

  } else if (type == "dw") {
    pdf_plot <- pdf_plot +
      geom_segment(data = emp, aes(x = ci_dw_lo, xend = ci_dw_lo, y = 0, yend = pdf_ci_dw_lo, color = "bounds"), linewidth = 0.5)+
      geom_segment(data = emp, aes(x = ci_dw_hi, xend = ci_dw_hi, y = 0, yend = pdf_ci_dw_hi, color = "bounds"), linewidth = 0.5)+
      geom_segment(data = emp, aes(x = ci_dw_lo, xend = ci_dw_hi, y = 0, yend = 0, color = "bounds"))+
      geom_point(data = emp, aes(ci_dw_lo, 0, color = "bounds"), size = 3)+
      geom_point(data = emp, aes(ci_dw_hi, 0, color = "bounds"), size = 3)+
      geom_point(data = emp, aes(e_GHI_dw, 0, color = "fitted"), shape = 17, size = 2)+
      geom_point(data = emp, aes(e_Rt_dw, 0, color = "expected_value"), size = 3)
  }
  type <- ifelse(type == "mix", "Mixture", ifelse(type == "up", "Sunny", "Cloudy"))
  plot_title <- ifelse(unconditional, "Unconditional", "Conditional")
  plot_title <- paste0(plot_title, " forecast: ", date, " (", type, ")")
  pdf_plot+
    geom_point(data = emp, aes(GHI, 0, color = "realized"), size = 3)+
    theme_bw()+
    theme(legend.position = "top")+
    labs(color = NULL, y = NULL, x = NULL, title = plot_title)+
    scale_color_manual(values = vals_colors, labels = vals_labels)
}


#' Stationarity and distribution test (Gaussian mixture) for a `solarModel`
#'
#' @examples
#' model <- Bologna
#' solarModel_test_residuals(model)
#' @rdname solarModel_test_residuals
#' @name solarModel_test_residuals
#' @export
solarModel_test_residuals <- function(model, nrep = 50, ci = 0.05, min_quantile = 0.015, max_quantile = 0.985, seed = 1){

  distribution_test <- list()
  stationary_test <- list()
  for(nmonth in 1:12){
    # Residuals
    x <- dplyr::filter(model$data, Month == nmonth)$u_tilde
    # Gaussian mixture parameters
    means = model$NM_model$means[nmonth,]
    sd = model$NM_model$sd[nmonth,]
    p = model$NM_model$p[nmonth,]

    # Mixture CDF
    cdf_GM <- function(x) pmixnorm(x, means = means, sd = sd, p = p)
    # Distribution test
    distribution_test[[nmonth]] <- dplyr::bind_cols(Month = nmonth,
                                                    ks_test(x, cdf_GM, ci = ci, min_quantile = min_quantile, max_quantile = max_quantile))
    # Stationary test
    iid_test <- function(seed) ks_ts_test(x, ci = ci, min_quantile = min_quantile, max_quantile = max_quantile, seed = seed)
    iid_tests <- purrr::map_df(1:nrep, ~iid_test(.x))
    Rejected <- mean(iid_tests$H0 == "Non-Rejected")
    stationary_test[[nmonth]] <- dplyr::bind_cols(Month = nmonth, Rejected = 1-Rejected, Non_Rej = Rejected)
    stationary_test[[nmonth]] <- dplyr::bind_cols(stationary_test[[nmonth]], iid_test(seed))
  }

  structure(
    list(
      distribution = dplyr::bind_rows(distribution_test),
      stationary = dplyr::bind_rows(stationary_test)
    )
  )
}


#' Long term variance AR process
#' @param phi AR parameters
#' @param sigma2 variance
#' @examples
#'AR_variance(0.4)
#'AR_variance(c(0.4, 0.2))
#'AR_variance(c(0.4, 0.2, 0.1))
#'AR_variance(c(0.3, 0.2, 0.1, 0.05))
#' @export
AR_variance <- function(phi, sigma2 = 1){
  var <- NA
  if (length(phi) == 1) {
    var <- sigma2/(1-phi[1]^2)
  } else if (length(phi) == 2) {
    var <- sigma2*(1-phi[2])/((1 - phi[2])*(1 - phi[1]^2 - phi[2]^2) - 2*phi[1]^2*phi[2])
  } else if (length(phi) == 3) {
    phi_tilde_1 <- (phi[1] + phi[2]*phi[3])/(1 - phi[2] - phi[3]^2 - phi[1]*phi[3])
    phi_tilde_2 <- (phi[1] + phi[3])*phi_tilde_1 + phi[2]
    phi_tilde_3 <- (phi[1]*phi_tilde_2 + phi[2]*phi_tilde_1 + phi[3])
    phi_tilde_0 <- 1/(1 - phi[1]*phi_tilde_1 - phi[2]*phi_tilde_2 - phi[3]*phi_tilde_3)
    var <- phi_tilde_0*sigma2
  } else if (length(phi) == 4) {
    phi_1 <- (phi[1] + phi[3])/(1 - phi[4])
    phi_0 <- phi[2]/(1 - phi[4])
    psi_1 <- (phi_1*phi[3] + phi[1]*phi_0*phi[4]+ phi[1] + phi[3]*phi[4])
    psi_1 <- psi_1/(1-phi_1*(phi[3] + phi[1]*phi[4]) - phi[2]*(1 + phi[4]) - phi[4]^2)
    psi_2 <- phi_1*psi_1 + phi_0
    psi_3 <- phi[1]*psi_2 + phi[2]*psi_1 + phi[4]*psi_1 + phi[3]
    psi_4 <- phi[1]*psi_3 + phi[2]*psi_2 + phi[3]*psi_1 + phi[4]
    var <- sigma2/(1 - phi[1]*psi_1 - phi[2]*psi_2 - phi[3]*psi_3 - phi[4]*psi_4)
  } else {
    warning("Variance for AR(", length(phi), ") not implemented!")
  }
  return(var)
}


