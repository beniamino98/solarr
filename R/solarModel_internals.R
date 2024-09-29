#' Extract a dataset from a solarModel object
#'
#' @examples
#' model <- Bologna
#' solarModel_data(model)
#' # Do not add monthly data
#' solarModel_data(model, monthly = FALSE)
#' # Do not add seasonal data
#' solarModel_data(model, monthly = FALSE, GM = FALSE)
#' @export
solarModel_data <- function(model, monthly = TRUE, GM = TRUE){
  # Seasonal data
  seasonal_data <- solarModel_seasonal_data(model, monthly = monthly, GM = GM)
  # Dataset
  data <- dplyr::left_join(model$data, seasonal_data, by = c("Month", "Day"))
  return(data)
}

#' Extract a dataset with seasonal data from `solarModel`
#'
#' @examples
#' model <- Bologna
#' # All seasonal and monthly data
#' solarModel_seasonal_data(model)
#' # Do not add monthly data
#' solarModel_seasonal_data(model, monthly = FALSE)
#' # Do not add Gaussian mixture data
#' solarModel_seasonal_data(model, monthly = FALSE, GM = FALSE)
#' # Extract onlu a particular day
#' solarModel_seasonal_data(model, date="2022-01-01")
#' @export
solarModel_seasonal_data <- function(model, monthly = TRUE, GM = TRUE, nmonths, ndays, date){
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
  # Add monthly data
  if (monthly) {
    data <- dplyr::left_join(data, model$monthly_data, by = "Month")
  }
  # Add Gaussian mixture parameters
  if (GM) {
    data <- dplyr::left_join(data, model$NM_model[,c(1:6)], by = "Month")
  }
  return(data)
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
  data <- solarModel_data(model)
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

  # Extract complete data
  data <- solarModel_seasonal_data(model, nmonths = nmonths, ndays = ndays, date = date)
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


#' Compute the log-likelihood of a `solarModel` object
#'
#' @param model `solarModel` object
#' @param nmonths months to consider
#' @examples
#' model <- Bologna
#' solarModel_loglik(model)
#' @export
solarModel_loglik <- function(model, nmonths = 1:12){
  # Conditional moments
  moments <- solarModel_conditional_moments(model)
  # Add weights
  moments$weights <- model$outliers$weights
  # Filter for a set of months
  moments <- dplyr::filter(moments, Month %in% nmonths)
  # Compute log likelihood
  moments$log_lik <- 0
  for(i in 1:nrow(moments)){
    df_n <- moments[i,]
    if (df_n$weights == 0){
      moments$log_lik[i] <- 0
      next
    }
    # Conditional mixture Pdf
    pdf_Yt <- function(x) dmixnorm(x, means = c(df_n$e_Yt_up, df_n$e_Yt_dw), sd = c(df_n$sd_Yt_up, df_n$sd_Yt_dw), p = c(df_n$p_up, 1-df_n$p_up))
    moments$log_lik[i] <- log(pdf_Yt(df_n$Yt))
  }
  # model$log_lik <- sum(moments$log_lik, na.rm = TRUE)
  return(sum(moments$log_lik, na.rm = TRUE))
}


#' Produce a forecast from a solarModel object
#'
#' @examples
#' model <- Bologna
#' solarModel_forecaster(model, date = "2020-04-01")
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
#' dates <- seq.Date(as.Date("2020-01-01"), as.Date("2020-12-31"), 1)
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
      geom_segment(data = emp, aes(x = ci_mix_lo, xend = ci_mix_lo, y = 0, yend = pdf_ci_mix_lo, color = "bounds"), size = 0.5)+
      geom_segment(data = emp, aes(x = ci_mix_hi, xend = ci_mix_hi, y = 0, yend = pdf_ci_mix_hi, color = "bounds"), size = 0.5)+
      geom_segment(data = emp, aes(x = ci_mix_lo, xend = ci_mix_hi, y = 0, yend = 0, color = "bounds"))+
      geom_point(data = emp, aes(ci_mix_lo, 0, color = "bounds"), size = 3)+
      geom_point(data = emp, aes(ci_mix_hi, 0, color = "bounds"), size = 3)+
      geom_point(data = emp, aes(e_GHI, 0, color = "fitted"), shape = 17, size = 2)+
      geom_point(data = emp, aes(e_Rt, 0, color = "expected_value"), size = 3)


  } else if (type == "up") {
    pdf_plot <- pdf_plot +
      geom_segment(data = emp, aes(x = ci_up_lo, xend = ci_up_lo, y = 0, yend = pdf_ci_up_lo, color = "bounds"), size = 0.5)+
      geom_segment(data = emp, aes(x = ci_up_hi, xend = ci_up_hi, y = 0, yend = pdf_ci_up_hi, color = "bounds"), size = 0.5)+
      geom_segment(data = emp, aes(x = ci_up_lo, xend = ci_up_hi, y = 0, yend = 0, color = "bounds"))+
      geom_point(data = emp, aes(ci_up_lo, 0, color = "bounds"), size = 3)+
      geom_point(data = emp, aes(ci_up_hi, 0, color = "bounds"), size = 3)+
      geom_point(data = emp, aes(e_GHI_up, 0, color = "fitted"), shape = 17, size = 2)+
      geom_point(data = emp, aes(e_Rt_up, 0, color = "expected_value"), size = 3)

  } else if (type == "dw") {
    pdf_plot <- pdf_plot +
      geom_segment(data = emp, aes(x = ci_dw_lo, xend = ci_dw_lo, y = 0, yend = pdf_ci_dw_lo, color = "bounds"), size = 0.5)+
      geom_segment(data = emp, aes(x = ci_dw_hi, xend = ci_dw_hi, y = 0, yend = pdf_ci_dw_hi, color = "bounds"), size = 0.5)+
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





