#' Control parameters for a `solarModel` object
#'
#' Control function for a solarModel
#'
#' @param clearsky list with control parameters, see \code{\link{control_seasonalClearsky}} for details.
#' @param seasonal.mean a list of parameters. Available choices are:
#'\describe{
#'  \item{`seasonalOrder`}{An integer specifying the order of the seasonal component in the model. The default is `1`.}
#'  \item{`include.intercept`}{When `TRUE` the intercept will be included in the seasonal model. The dafault if `TRUE`.}
#'  \item{`monthly.mean`}{When `TRUE` a set of 12 monthly means parameters will be computed from the deseasonalized time series to center it perfectly around zero.}
#'}
#' @param seasonal.variance a list of parameters. Available choices are:
#'\describe{
#'  \item{`seasonalOrder`}{An integer specifying the order of the seasonal component in the model. The default is `1`.}
#'  \item{`correction`}{When true the seasonal variance is corrected to ensure that the standardize the residuals with a unitary variance.  }
#'  \item{`monthly.mean`}{When `TRUE` a set of 12 monthly variances parameters will be computed from the deseasonalized time series to center it perfectly around zero.}
#'}
#' @param mean.model a list of parameters.
#' \describe{
#'  \item{`arOrder`}{An integer specifying the order of the AR component in the model. The default is `2`.}
#'  \item{`include.intercept`}{When `TRUE` the intercept will be included in the AR model. The dafault if `FALSE`.}
#'}
#' @param variance.model an `ugarchspec` object for GARCH variance. Default is `GARCH(1,1)`.
#' @param mixture.model a list of parameters.
#' @param threshold numeric, threshold for the estimation of alpha and beta.
#' @param outliers_quantile quantile for outliers detection. If different from 0, the observations that are below the quantile at confidence levels `outliers_quantile` and
#' the observation above the quantile at confidence level 1-`outliers_quantile` will have a weight equal to zero and will be excluded from estimations.
#' @param quiet logical, when `TRUE` the function will not display any message.
#' @examples
#' control <- control_solarModel()
#' @rdname control_solarModel
#' @export
control_solarModel <- function(clearsky = control_seasonalClearsky(),
                               stochastic_clearsky = FALSE,
                               seasonal.mean = list(seasonalOrder = 1, include.H0 = FALSE, include.intercept = TRUE, monthly.mean = TRUE),
                               mean.model = list(arOrder = 2, include.intercept = FALSE),
                               seasonal.variance = list(seasonalOrder = 1, correction = TRUE, monthly.mean = TRUE),
                               variance.model = rugarch::ugarchspec(variance.model = list(garchOrder = c(1,1)),
                                                                    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE)),
                               mixture.model = list(match_moments = FALSE, abstol = 1e-3, maxit = 150),
                               threshold = 0.01, outliers_quantile = 0, quiet = FALSE){

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
  mixture_model = list(match_moments = FALSE, abstol = 1e-20, maxit = 100, prior_p = NA)
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
      quiet = quiet
    ),
    class = c("control", "list")
  )
}


#' Specification function for a `solarModel`
#'
#' @param place character, name of an element in the `CAMS_data` list.
#' @param target character, target variable to model. Can be `GHI` or `clearsky`.
#' @param min_date character. Date in the format `YYYY-MM-DD`. Minimum date for the complete data. If `missing` will be used the minimum data available.
#' @param max_date character. Date in the format `YYYY-MM-DD`. Maximum date for the complete data. If `missing` will be used the maximum data available.
#' @param from character. Date in the format `YYYY-MM-DD`. Starting date to use for training data.
#' If `missing` will be used the minimum data available after filtering for `min_date`.
#' @param to character. Date in the format `YYYY-MM-DD`. Ending date to use for training data.
#' If `missing` will be used the maximum data available after filtering for `max_date`.
#' @param CAMS_data named list with radiation data for different locations.
#' @param control_model list with control parameters, see \code{\link{control_solarModel}} for details.
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

  # Initialize seasonal dataset
  # Seasonal data by month and day for an year with 366 days
  seasonal_data <- dplyr::tibble(date = seq.Date(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "1 day"))
  seasonal_data <- dplyr::mutate(seasonal_data,
                                 Month = lubridate::month(date),
                                 Day = lubridate::day(date),
                                 n = number_of_day(date))
  seasonal_data <- dplyr::select(seasonal_data, -date)
  seasonal_data <- dplyr::arrange(seasonal_data, Month, Day)

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

#' Monthly Gaussian mixture with two components
#'
#' @param x arg
#' @param date arg
#' @param match_moments arg
#' @inheritParams gaussianMixture
#'
#' @rdname solarModel_mixture
#' @name solarModel_mixture
#' @export
solarModel_mixture <- function(x, date, weights, match_moments = FALSE, maxit = 100, abstol = 10e-15){

  data <- dplyr::tibble(date = date, Month = lubridate::month(date), eps = x, w = weights)
  # Gaussian Mixture parameters
  GM_model <- list()
  # Monthly data
  data_months <- list()
  for(m in 1:12){
    data_months[[m]] <- dplyr::filter(data, Month == m)
    w <- data_months[[m]]$w
    # Monthly data
    eps <- data_months[[m]]$eps
    # Initial parameters
    mu_0 <- quantile(eps, c(0.2, 0.8), na.rm = TRUE)
    sd_0 <- sd(eps, na.rm = TRUE)
    init_means <- c(mu1 = mu_0[1], mu2 = mu_0[2])
    init_sd <- c(sd1 = sd_0, sd2 = sd_0)
    init_p <- c(0.5, 0.5)
    # Fitted model
    gm <- gaussianMixture(eps, means = init_means, components = 2, sd = init_sd, p = init_p,
                          weights = w, maxit = maxit, abstol = abstol, na.rm = TRUE)
    # Compute the sample mean
    e_u_hat <- mean(eps[w!=0], na.rm = TRUE)
    # Compute the sample variance
    v_u_hat <- var(eps[w!=0], na.rm = TRUE)
    # Match moments approach
    if (match_moments) {
      if (gm$par$p[1] > gm$par$p[2]) {
        gm$par$mean[2] <- (e_u_hat - gm$par$mean[1]*gm$par$p[1])/gm$par$p[2]
        gm$par$sd[2] <- sqrt((v_u_hat - gm$par$sd[1]^2)/gm$par$p[2] - gm$par$p[1]*(gm$par$mean[1]-gm$par$mean[2])^2)
      } else {
        gm$par$mean[1] <- (e_u_hat-gm$par$mean[2]*gm$par$p[2])/gm$par$p[1]
        gm$par$sd[1] <- sqrt((v_u_hat - gm$par$sd[2]^2)/gm$par$p[1] - gm$par$p[2]*(gm$par$mean[1]-gm$par$mean[2])^2)
      }
    }
    # Compute the theoretical expected value
    e_u <- sum(gm$par$mean*gm$par$p)
    # Compute the theoric variance
    v_u <- sum((gm$par$mean^2 + gm$par$sd^2)*gm$par$p) - e_u^2
    # Fitted parameters
    df_par <- dplyr::bind_cols(purrr::map(gm$par, ~dplyr::bind_rows(.x)))
    # Fitted expected value
    data_months[[m]]$e_x <- gm$fitted$B1*gm$par$mean[1] + gm$fitted$B2*gm$par$mean[2]
    # Add fitted Bernoulli series
    data_months[[m]]$B <- gm$fitted$B1
    # Reorder parameters
    if (df_par$mu1 < df_par$mu2) {
      data_months[[m]]$B <- gm$fitted$B2
    }
    # Monthly data
    GM_model[[m]] <- dplyr::tibble(Month = m,
                                   df_par,
                                   loss = gm$log_lik,
                                   nobs = length(eps[w!=0]),
                                   e_x = e_u,
                                   v_x = v_u,
                                   e_x_hat = e_u_hat,
                                   v_x_hat = v_u_hat)
  }

  # Reorder the variables
  GM_model <- dplyr::bind_rows(GM_model)
  GM_model <- dplyr::mutate(GM_model,
                            mu_up = dplyr::case_when(
                              mu1 > mu2 ~ mu1,
                              TRUE ~ mu2),
                            mu_dw = dplyr::case_when(
                              mu1 > mu2 ~ mu2,
                              TRUE ~ mu1),
                            sd_up = dplyr::case_when(
                              mu1 > mu2 ~ sd1,
                              TRUE ~ sd2),
                            sd_dw = dplyr::case_when(
                              mu1 > mu2 ~ sd2,
                              TRUE ~ sd1),
                            p_up = dplyr::case_when(
                              mu1 > mu2 ~ p1,
                              TRUE ~ p2),
                            p_dw = 1 - p_up)
  # Reorder variables
  GM_model <- dplyr::select(GM_model, Month, mu_up, mu_dw, sd_up, sd_dw, p_up, p_dw,
                            loss, nobs, e_x, v_x, e_x_hat, v_x_hat)
  # Fitted series
  data_months <- dplyr::bind_rows(data_months)

  structure(
    list(
      model = GM_model,
      data = data_months
    ),
    class = c("solarModelMixture", "list")
  )
}


#' Monthly multivariate Gaussian mixture with two components
#'
#' @param model_Ct arg
#' @param model_GHI arg
#'
#' @rdname solarModel_mvmixture
#' @name solarModel_mvmixture
#' @export
solarModel_mvmixture <- function(model_Ct, model_GHI){

  data_GHI <- dplyr::select(model_GHI$data, date, Month, GHI = "u_tilde")
  data_Ct <- dplyr::select(model_Ct$data, date, Month, Ct = "u_tilde")
  data <- dplyr::inner_join(data_GHI, data_Ct, by = c("date", "Month"))
  data <- dplyr::filter(data, !(date %in% model_Ct$outliers$date) & !(date %in% model_GHI$outliers$date))

  # Gaussian Mixture parameters
  m <- 1
  # Monthly data
  data_months <- list()
  for(m in 1:12){
    data_months[[m]] <- dplyr::filter(data, Month == m)
    # Monthly data
    eps <- data_months[[m]][,c(3,4)]
    # Fitted model
    gm <- mvgaussianMixture(eps, components = 2, na.rm = TRUE)

    model_GHI$NM_model$mu_up[m] <- gm$params$means[1,2]
    model_GHI$NM_model$mu_dw[m] <- gm$params$means[2,2]
    model_GHI$NM_model$sd_up[m] <- sqrt(gm$params$sigma2[1,2])
    model_GHI$NM_model$sd_dw[m] <- sqrt(gm$params$sigma2[2,2])

    model_Ct$NM_model$mu_up[m] <- gm$params$means[1,1]
    model_Ct$NM_model$mu_dw[m] <- gm$params$means[2,1]
    model_Ct$NM_model$sd_up[m] <- sqrt(gm$params$sigma2[1,1])
    model_Ct$NM_model$sd_dw[m] <- sqrt(gm$params$sigma2[2,1])
    model_Ct$NM_model$p_up[m] <- model_GHI$NM_model$p_up[m] <- gm$params$p[1]
    model_Ct$NM_model$p_dw[m] <- model_GHI$NM_model$p_dw[m] <- gm$params$p[2]
    model_Ct$NM_model$rho_up[m] <- model_GHI$NM_model$rho_up[m] <- gm$params$rho[1]
    model_Ct$NM_model$rho_dw[m] <- model_GHI$NM_model$rho_dw[m] <- gm$params$rho[2]

    # Add fitted Bernoulli series
    data_months[[m]]$B <- gm$B_hat$B1
  }

  # Fitted series
  data_months <- dplyr::select(dplyr::bind_rows(data_months), date, B)

  model_Ct$data <- dplyr::left_join(dplyr::select(model_Ct$data, -B), data_months, by = "date")
  model_GHI$data <- dplyr::left_join(dplyr::select(model_GHI$data, -B), data_months, by = "date")
  model_GHI$data$B[is.na(model_GHI$data$B)] <- 0
  model_Ct$data$B[is.na(model_Ct$data$B)] <- 0

  structure(
    list(
      model_Ct = model_Ct,
      model_GHI = model_GHI
    ),
    class = c("solarModelMixture", "list")
  )
}


#' Compute conditional moments from a `solarModel` object
#'
#' @examples
#' model <- Bologna
#' solarModel_conditional_moments(model)
#' solarModel_conditional_moments(model, date = "2022-01-01")
#' @export
solarModel_conditional_moments <- function(model, date){

  # Extract complete data
  data <- model$data
  # Filter for a set of dates
  if (!missing(date)){
    data <- data[data$date %in% as.Date(date),]
  }
  # Compute conditional moments
  data <- dplyr::mutate(data,
                        # Conditional expectation Yt
                        e_Yt = Yt_bar + Yt_tilde_hat + Yt_tilde_uncond,
                        # Conditional std. deviation Yt
                        sd_Yt = sigma*sigma_bar*sigma_m,
                        # Conditional moments Yt (state up)
                        e_Yt_up = e_Yt + sd_Yt*mu_up,
                        sd_Yt_up = sd_Yt*sd_up,
                        # Conditional moments Yt (state dw)
                        e_Yt_dw = e_Yt + sd_Yt*mu_dw,
                        sd_Yt_dw = sd_Yt*sd_dw,
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
                        e_Yt_up, sd_Yt_up, e_Yt_dw, sd_Yt_dw, p_up, B,
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
                        e_Yt_up = e_Yt + sd_Yt*mu_up,
                        sd_Yt_up = sd_Yt*sd_up,
                        # Unconditional moments Yt (state dw)
                        e_Yt_dw = e_Yt + sd_Yt*mu_dw,
                        sd_Yt_dw = sd_Yt*sd_dw,
                        # Fitted value mixture
                        e_Yt_mix = e_Yt +  e_Yt_up*p_up + e_Yt_dw*(1-p_up),
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
                        e_Yt_up, sd_Yt_up, e_Yt_dw, sd_Yt_dw, p_up,
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
  pdf_Yt <- function(x) dmixnorm(x, means = c(df_n$e_Yt_up, df_n$e_Yt_dw), sd = c(df_n$sd_Yt_up, df_n$sd_Yt_dw), p = c(df_n$p_up, 1-df_n$p_up))
  # Density of GHI
  pdf_Rt <- function(x, pdf_Yt) dsolarGHI(x, df_n$Ct, model$transform$alpha, model$transform$beta, pdf_Yt)
  # Add grid of points
  grid[[paste0("pdf_", model$target)]] <- pdf_Rt(grid_x, pdf_Yt)
  grid[[paste0("pdf_", model$target, "_up")]] <- pdf_Rt(grid_x, function(x) df_n$p_up*dnorm(x, df_n$e_Yt_up, df_n$sd_Yt_up))
  grid[[paste0("pdf_", model$target, "_dw")]] <- pdf_Rt(grid_x, function(x) (1-df_n$p_up)*dnorm(x, df_n$e_Yt_dw, df_n$sd_Yt_dw))

  # Empiric variables for the day to be extracted
  emp <- solarModel_data(model, FALSE, FALSE)
  emp <- emp[emp$date == date,]
  emp <- dplyr::select(emp, date, Month, Day, tidyr::any_of(model$target), B)
  # Add fitted expected values
  vars_names <- paste0(paste0("e_", model$target), c("", "_mix", "_up", "_dw"))
  emp <- dplyr::bind_cols(emp,  dplyr::select(df_n, any_of(vars_names)))
  # ==========================  Probabilistic Forecasts ==============================
  # Up and down density non weighted
  pdf_Yt_up <- function(x) dnorm(x, df_n$e_Yt_up, df_n$sd_Yt_up)
  pdf_Yt_dw <- function(x) dnorm(x, df_n$e_Yt_dw, df_n$sd_Yt_dw)
  # Expected value function
  e_Rt <- function(pdf_Yt, .f = function(x) x) integrate(function(x) .f(x)*pdf_Rt(x, pdf_Yt), lower = lower_bound, upper = upper_bound)$value
  # Confidence intervals function
  conf_Rt <- function(p, pdf_Yt) qsolarGHI(p, df_n$Ct, model$transform$alpha, model$transform$beta, pdf_Yt)

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
  emp$ci_mix_lo <- conf_Rt(ci, pdf_Yt)
  emp$ci_mix_hi <- conf_Rt(1-ci, pdf_Yt)
  emp$pdf_ci_mix_lo <- pdf_Rt(emp$ci_mix_lo, pdf_Yt)
  emp$pdf_ci_mix_hi <- pdf_Rt(emp$ci_mix_hi, pdf_Yt)
  # Sunny state
  emp$ci_up_lo <- conf_Rt(ci, pdf_Yt_up)
  emp$ci_up_hi <- conf_Rt(1-ci, pdf_Yt_up)
  emp$pdf_ci_up_lo <- pdf_Rt(emp$ci_up_lo, pdf_Yt)
  emp$pdf_ci_up_hi <- pdf_Rt(emp$ci_up_hi, pdf_Yt)
  # Cloudy state
  emp$ci_dw_lo <- conf_Rt(ci, pdf_Yt_dw)
  emp$ci_dw_hi <- conf_Rt(1-ci, pdf_Yt_dw)
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
#' @export
solarModel_forecast <- function(model, date, ci = 0.1, unconditional = FALSE){
  fun <- function(date){
    safe_forecaster <- purrr::safely(solarModel_forecaster)
    smf <- safe_forecaster(model, date = date, ci = ci, unconditional = unconditional)$result
    if (is.null(smf)){
      return(NULL)
    }
    smf[-1][[1]]
  }
  out <- purrr::map_df(date, fun)
  return(out)
}

#' Plot a forecast from a solarModel object
#'
#' @examples
#' model <- Bologna
#' day_date <- "2024-01-01"
#' solarModel_forecaster_plot(model, date = day_date)
#' solarModel_forecaster_plot(model, date = day_date, unconditional = TRUE)
#' solarModel_forecaster_plot(model, date = day_date, type = "dw")
#' solarModel_forecaster_plot(model, date = day_date, type = "dw", unconditional = TRUE)
#' solarModel_forecaster_plot(model, date = day_date, type = "up")
#' solarModel_forecaster_plot(model, date = day_date, type = "up", unconditional = TRUE)
#'
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
    geom_line(data = grid, aes(x, pdf_GHI_up), color = "green")+
    geom_line(data = grid, aes(x, pdf_GHI_dw), color = "red")

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
#' @export
solarModel_test_residuals <- function(model, nrep = 50, ci = 0.05, min_quantile = 0.015, max_quantile = 0.985, seed = 1){

  distribution_test <- list()
  stationary_test <- list()
  for(nmonth in 1:12){
    # Residuals
    x <- dplyr::filter(model$data, Month == nmonth)$u_tilde
    # Gaussian mixture parameters
    nm <- model$NM_model[nmonth,]
    means = c(nm$mu_up,  nm$mu_dw)
    sd = c(nm$sd_up,  nm$sd_dw)
    p = c(nm$p_up,  1-nm$p_up)

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


#' Empiric Gaussian Mixture parameters
#'
#' @keywords OLD
solarModel_empiric_GM <- function(model, match_moments = FALSE){

  # Select only train data
  data <- model$data
  data$weights <- model$outliers$weights
  data <- dplyr::filter(data, isTrain & weights != 0)
  data <- dplyr::select(data, Month, u_tilde, B)
  for(m in 1:12){
    # Monthly data
    df_m <- dplyr::filter(data, Month == m)
    u_up <- dplyr::filter(df_m, B == 1)$u_tilde
    u_dw <- dplyr::filter(df_m, B == 0)$u_tilde
    nm <- model$NM_model[m,]
    # Up moments
    n1 <- length(u_up)
    nm$mu_up <- mean(u_up, na.rm = TRUE)
    nm$sd_up <- sd(u_up, na.rm = TRUE)*sqrt(n1/(n1-1))
    nm$p_up <- mean(df_m$B, na.rm = TRUE)
    # Down moments
    n2 <- length(u_dw)
    nm$mu_dw <- mean(u_dw, na.rm = TRUE)
    nm$sd_dw <- sd(u_dw, na.rm = TRUE)*sqrt(n2/(n2-1))
    nm$p_dw <- 1 - nm$p_up
    df_m
    if (match_moments) {
      if (n1 >= n2) {
        nm$mu_dw <- (nm$e_x_hat - nm$mu_up*nm$p_up)/(1-nm$p_up)
        nm$sd_dw <- sqrt((nm$v_x_hat - nm$sd_up^2*nm$p_up)/(1-nm$p_up) - nm$p_up*(nm$mu_up - nm$mu_dw)^2)
      } else {
        nm$mu_up <- (nm$e_x_hat - nm$mu_dw*nm$p_dw)/(1-nm$p_dw)
        nm$sd_up <- sqrt((nm$v_x_hat - nm$sd_dw^2*nm$p_dw)/(1-nm$p_dw) - nm$p_dw*(nm$mu_up - nm$mu_dw)^2)
      }
    }
    # Extract parameters
    means = c(nm$mu_up, nm$mu_dw)
    sd = c(nm$sd_up, nm$sd_dw)
    p = c(nm$p_up, nm$p_dw)
    # Update log-likelihood
    pdf <- function(x) dmixnorm(x, means = means, sd = sd, p = p)
    nm$loss <- sum(log(pdf(df_m$u_tilde)), na.rm = TRUE)
    nm$nobs <- nrow(df_m)
    # Update moments
    nm$e_x <- sum(means*p)
    nm$v_x <- sum((means^2 + sd^2)*p) - nm$e_x^2
    # Update data
    params <- unlist(nm[1,c(2:7)])
    model <- solarModel_update_GM(model, params, m)
  }
  model$log_lik <- solarModel_loglik(model)
  model
}

