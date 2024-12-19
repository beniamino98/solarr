#' solarModelCont
#' @export
solarModelCont <- R6::R6Class("solarModelCont",
                              inherit = solarModel,
                              public = list(
                                theta = NA,
                                initialize = function(place, ...){
                                  control <- control_solarModel(clearsky = control_seasonalClearsky(),
                                                                stochastic_clearsky = FALSE,
                                                                seasonal.mean = list(seasonalOrder = 1, include.H0 = FALSE, include.intercept = TRUE, monthly.mean = FALSE),
                                                                mean.model = list(arOrder = 1, include.intercept = FALSE),
                                                                seasonal.variance = list(seasonalOrder = 1, correction = TRUE, monthly.mean = TRUE),
                                                                variance.model = rugarch::ugarchspec(variance.model = list(garchOrder = c(1,1)),
                                                                                                     mean.model = list(armaOrder = c(0, 0), include.mean = FALSE)),
                                                                mixture.model = list(abstol = 1e-3, maxit = 150),
                                                                threshold = 0.01, outliers_quantile = 0, garch_variance = FALSE, quiet = FALSE)
                                  # Model specification
                                  spec <- solarModel_spec(place, ..., control_model = control)
                                  super$initialize(spec)
                                  super$fit()
                                  self$theta <- -log(1-self$AR_model_Yt$coefficients[1])
                                },
                                # Clear sky radiation
                                Ct = function(n){
                                  self$seasonal_model_Ct$predict(n)
                                },
                                # Seasonal radiation
                                Yt_bar = function(n){
                                  self$seasonal_model_Yt$predict(n)
                                },
                                # Seasonal radiation
                                Rt_bar = function(n){
                                  self$Y_to_R(self$Yt_bar(n))
                                },
                                R_to_Y = function(x, n){
                                  Xt <- self$transform$iGHI(x, self$Ct(n))
                                  self$transform$Y(Xt)
                                },
                                Y_to_R = function(x, n){
                                  self$transform$GHI_y(x, self$Ct(n))
                                },
                                # Conditional expectation of Y
                                M_Y = function(Rt, n){
                                  Yt <- self$R_to_Y(Rt, n)
                                  NM <- self$seasonal_data
                                  Y_forecast <- self$Yt_bar(n+1) + Yt + (self$Yt_bar(n) - Yt)*exp(-self$theta)
                                  M_Y_dw <- Y_forecast + sqrt(self$seasonal_variance$predict(n+1))*NM$mu1[n+1]
                                  M_Y_up <- Y_forecast + sqrt(self$seasonal_variance$predict(n+1))*NM$mu2[n+1]
                                  dplyr::tibble(dw = M_Y_dw, up = M_Y_up, p1 = NM$p1[n+1], p2 = NM$p2[n+1])
                                },
                                S_Y = function(n){
                                  NM <- self$seasonal_data
                                  S_Y_dw <- sqrt(self$seasonal_variance$predict(n+1)*NM$sd1[n+1]*(1-exp(-2*self$theta))/(2*self$theta))
                                  S_Y_up <- sqrt(self$seasonal_variance$predict(n+1)*NM$sd2[n+1]*(1-exp(-2*self$theta))/(2*self$theta))
                                  dplyr::tibble(dw = S_Y_dw, up = S_Y_up, p1 = NM$p1[n+1], p2 = NM$p2[n+1])
                                },
                                pdf_Y = function(x, Rt, n){
                                  mean <- unlist(self$M_Y(Rt, n)[,2:1])
                                  sd <- unlist(self$S_Y(n)[,2:1])
                                  alpha <- unlist(self$M_Y(Rt, n)[,3:4])
                                  extraDistr::dmixnorm(x, mean, sd, alpha)
                                }
                                )
                              )

#' solarScenarioCont_filter
#' @export
solarScenarioCont_filter <- function(simSpec, lam = 1){

  if (purrr::is_empty(simSpec$residuals)) {
    stop("The slot `simSpec$residuals` is empty! Consider running `simSpec <- solarScenario_residuals(simSpec)` before!")
  }
  # Number of lags to consider
  i_start <- simSpec$i_start
  # Number of simulations
  nsim <- nrow(simSpec$residuals)
  simulations <- list()
  theta <- 1-simSpec$AR_model_Yt$coefficients
  j <- 1
  for(j in 1:nsim){
    # Initialize dataset for storing the simulation
    df_sim <- simSpec$sim
    # Add simulated residuals
    df_sim$z <- 0
    df_sim$z[(i_start):nrow(df_sim)] <- simSpec$residuals[j,]$X[[1]][, simSpec$place][[1]]
    df_sim$B[(i_start):nrow(df_sim)] <- simSpec$residuals[j,]$B[[1]][, simSpec$place][[1]]
    # Verbose message
    if (!simSpec$quiet) message("Simulation: ", j, "/", nsim, " (", round(j/nsim*100, 4), " %) \r", appendLF = FALSE)
    # Routine
    i <- i_start
    df_sim$lam <- 0
    df_sim$dQdP <- 1
    for(i in i_start:nrow(df_sim)){
      # Simulated GARCH standard deviation
      df_sim$sigma[i] <- 1 # simSpec$GARCH_next_step(df_sim$eps_tilde[(i-simSpec$archOrder):(i-1)], df_sim$sigma[(i-simSpec$garchOrder):(i-1)])
      # Simulated standardized monthly normal mixture
      df_sim$u_tilde[i] <- (df_sim$mu1[i] + df_sim$sd1[i]*df_sim$z[i])*df_sim$B[i] + (df_sim$mu2[i] + df_sim$sd2[i]*df_sim$z[i])*(1-df_sim$B[i])
      # Simulated GARCH residuals
      df_sim$u[i] <- df_sim$u_tilde[i] #df_sim$sigma_m[i]
      # Simulated deseasonalized residuals
      df_sim$eps_tilde[i] <- df_sim$sigma[i]*df_sim$u[i]
      # Simulated AR residuals
      df_sim$eps[i] <- df_sim$eps_tilde[i]*df_sim$sigma_bar[i]
      # Simulated Yt_tilde
      df_sim$Yt_tilde[i] <- df_sim$Yt_tilde[i-1] -theta * (df_sim$Yt[i-1] - df_sim$Yt_bar[i-1]) + df_sim$eps[i]
      # Simulated Yt
      df_sim$Yt[i] <- df_sim$Yt_bar[i] + df_sim$Yt_tilde[i]
      # Simulated Xt
      df_sim$Xt[i] <- simSpec$transform$iY(df_sim$Yt[i])
      # Simulated GHI
      df_sim[[simSpec$target]][i] <- simSpec$transform$GHI(df_sim$Xt[i], df_sim$Ct[i])
      # Drift for the change of measure
      df_sim$lam[i] <- (df_sim$mu1[i]*df_sim$B[i] + df_sim$mu2[i]*(1-df_sim$B[i]))
      df_sim$lam[i] <- df_sim$lam[i] + 0.5*df_sim$sigma_bar[i]*(df_sim$B[i]*df_sim$sd1[i]^2 + (1 - df_sim$B[i])*df_sim$sd2[i]^2)*(1 - exp(df_sim$Yt[i]))
      df_sim$lam[i] <- lam*df_sim$lam[i]/(df_sim$B[i]*df_sim$sd1[i] + (1 - df_sim$B[i])*df_sim$sd2[i])
      # Doléans-Dade exponential
      df_sim$dQdP[i] <- exp(-df_sim$lam[i]*df_sim$z[i] - 0.5*df_sim$lam[i]^2)
      df_sim$Yt[i] <- df_sim$Yt[i] - df_sim$sigma_bar[i]*df_sim$lam[i]*(df_sim$sd1[i]*df_sim$B[i] +  df_sim$sd2[i])*(1-df_sim$B[i])
    }
    #df_sim$dQdP <- df_sim$dQdP/mean(df_sim$dQdP)
    # Remove redundant variables
    df_sim <- dplyr::select(df_sim, -mu1, -mu2, -sd1, -sd2, -p1, -Ct, -Yt_bar, -sigma_bar, -sigma_m, -Yt_tilde_hat, -Yt_tilde_uncond)
    # Remove initial values
    if (simSpec$exclude_known) {
      df_sim <- dplyr::filter(df_sim, date >= simSpec$from & date <= simSpec$to)
    }
    # Store simulations
    simulations[[j]] <- dplyr::bind_cols(nsim = j, df_sim)
  }
  # Add simulations
  simSpec$simulations <- simulations
  return(simSpec)
}

#' solarScenarioCont
#' @export
solarScenarioCont <- function(model, from = "2010-01-01", to = "2011-01-01", by = "1 year", theta = 0, nsim = 1, seed = 1, quiet = FALSE, ...){

  idx_date <- seq.Date(as.Date(from), as.Date(to), by = by)
  scenarios <- list()
  df_emp <- dplyr::tibble()
  n_scenario <- length(idx_date)
  j <- 2
  for(j in 2:n_scenario){
    if (!quiet) {
      # To report progress
      pb <- txtProgressBar(min = 1,            # Minimum value of the progress bar
                           max = n_scenario,   # Maximum value of the progress bar
                           style = 3,          # Progress bar style (also available style = 1 and style = 2)
                           width = 50,         # Progress bar width. Defaults to getOption("width")
                           char = "#")
      setTxtProgressBar(pb, j)
    }
    simSpec <- solarScenario_spec(model, from = idx_date[j-1], to = idx_date[j]-1, theta = theta, exclude_known = TRUE, quiet = TRUE)
    simSpec <- solarScenario_residuals(simSpec, nsim = nsim, seed = seed)
    simSpec <- solarScenarioCont_filter(simSpec, ...)
    df_emp <- dplyr::bind_rows(df_emp, simSpec$emp)
    scenarios[[j]] <- dplyr::bind_cols(seed = seed, dplyr::bind_rows(simSpec$simulations))
    seed <- seed + j - 1
  }
  if (!quiet) close(pb)

  simSpec$emp <- df_emp[!duplicated(df_emp),]
  simSpec$simulations <- scenarios
  return(as_solarScenario(simSpec))
}

#' solarOption_modelCont
#' @export
solarOption_modelCont <- function(model, nmonths = 1:12, theta = 0, combinations = NA, implvol = 1, put = TRUE, target.Yt = TRUE, control_options = control_solarOption()){

  # Options control
  K <- control_options$K
  leap_year = control_options$leap_year
  # Target and seasonal mean
  target <- model$target
  target_bar <- paste0(model$target, "_bar")
  target_plus <- paste0(model$target, "_plus")
  # AR(2) stationary variance
  ar_variance <- model$AR_model_Yt$variance
  # Extract seasonal data
  seasonal_data <- model$seasonal_data
  # Filter for the selected months
  seasonal_data <- dplyr::filter(seasonal_data, Month %in% nmonths)
  # Remove 29 of February from computations
  if (!leap_year) {
    seasonal_data <- dplyr::filter(seasonal_data, !(Month == 2 & Day == 29))
  }
  # Option Strike
  seasonal_data$strike <- seasonal_data[[target_bar]]*exp(K)
  # Add upper and lower bounds and unconditional moments
  seasonal_data <- dplyr::mutate(seasonal_data,
                                 e_Yt = Yt_bar + Yt_tilde_uncond,
                                 sd_Yt = implvol*sigma_bar*sigma_m*sqrt(ar_variance))
  # Add transform parameters
  seasonal_data$alpha <- model$transform$alpha
  seasonal_data$beta <- model$transform$beta
  if (target.Yt) {
    # Bounds for integration
    seasonal_data <- dplyr::mutate(seasonal_data,
                                   lower = -Inf,
                                   upper = Inf,
                                   z_x = model$transform$Y(1-strike/Ct),
                                   f_x = list(function(x, Ct) model$transform$GHI_y(x, Ct)))
  } else {
    # Bounds for integration
    seasonal_data <- dplyr::mutate(seasonal_data,
                                   lower = Ct*(1-alpha-beta),
                                   upper = Ct*(1-alpha),
                                   z_x = strike,
                                   f_x = list(function(x, Ct) x))
  }

  # Add Esscher parameter
  if (is.function(theta)) {
    seasonal_data$theta <- purrr::map_dbl(seasonal_data$Month, ~theta(.x))
  } else {
    seasonal_data$theta <- theta
  }


  # Select only relevant variables
  seasonal_data <- dplyr::select(seasonal_data, n, Month, Day, strike, Ct, theta, Yt_bar, z_x, e_Yt, sd_Yt, alpha, beta,
                                 mu1, mu2, sd1, sd2, p1, sigma_bar, lower, upper, f_x)

  # Pricing function
  pricing_month_day <- function(seasonal_data, combinations = NA, nmonth = 1, nday = 1, put = TRUE){

    # Seasonal Data
    df_n <- dplyr::filter(seasonal_data, Month == nmonth & Day == nday)
    if (df_n$n == 1){
      df_L1 <- dplyr::filter(seasonal_data, n == 365)
    } else {
      df_L1 <- dplyr::filter(seasonal_data, n == df_n$n - 1)
    }
    df_L1 <- dplyr::mutate(df_L1,
                           e_Yt_up = e_Yt + sd_Yt*mu1,
                           e_Yt_dw = e_Yt + sd_Yt*mu2,
                           sd_Yt_up = sd_Yt*sd1,
                           sd_Yt_dw = sd_Yt*sd2)

    E_expYt <- 0.5*df_L1$p1*(1 - exp(df_L1$e_Yt_up + 0.5*df_L1$sd_Yt_up^2)) + 0.5*(1-df_L1$p1)*(1 - exp(df_L1$e_Yt_dw + 0.5*df_L1$sd_Yt_dw^2))
    df_n$lam <- 0#df_L1$sigma_bar*(df_L1$mu1*df_L1$p1 + df_L1$mu2*(1-df_L1$p1))
    #df_n$lam <- df_n$lam + E_expYt*df_n$sigma_bar^2*(df_L1$sd1^2*df_L1$p1 + df_L1$sd2^2*(1-df_L1$p1))
    df_n$lam1 <- 0.5*(1 - exp(df_L1$e_Yt_up + 0.5*df_L1$sd_Yt_up^2))*df_L1$sigma_bar^2*df_L1$sd1^2 + df_L1$mu1*df_L1$sigma_bar #+ model$seasonal_model_Yt$differential(df_L1$n)/df_L1$sigma_bar
    df_n$lam2 <- 0.5*(1 - exp(df_L1$e_Yt_dw + 0.5*df_L1$sd_Yt_dw^2))*df_L1$sigma_bar^2*df_L1$sd2^2 + df_L1$mu2*df_L1$sigma_bar #+ model$seasonal_model_Yt$differential(df_L1$n)/df_L1$sigma_bar


    # Mixture combinations
    if (is.na(combinations) && length(combinations) == 1) {
      df_n <- dplyr::mutate(df_n,
                            e_Yt_up = e_Yt + sd_Yt*mu1 - lam1,
                            e_Yt_dw = e_Yt + sd_Yt*mu2 - lam2,
                            sd_Yt_up = sd_Yt*sd1,
                            sd_Yt_dw = sd_Yt*sd2)
      # Create combinations table
      comb <- dplyr::tibble(mean = c(df_n$e_Yt_up, df_n$e_Yt_dw), sd = c(df_n$sd_Yt_up, df_n$sd_Yt_dw), probs = c(df_n$p1, 1-df_n$p1))
    } else {
      comb <- dplyr::filter(combinations, Month == nmonth)
      # Update mixture parameters
      comb$mean <- df_n$e_Yt + df_n$sd_Yt*comb$mean
      comb$sd <- df_n$sd_Yt*comb$sd
    }

    if (target.Yt) {
      if (df_n$theta == 0) {
        # Mixture Pdf
        pdf_Yt <- function(x) dmixnorm(x, comb$mean, comb$sd, comb$probs)
        # Distribution Yt
        cdf_Rt <- function(x) pmixnorm(x, comb$mean, comb$sd, comb$probs)
      } else {
        # Esscher Mixture Pdf
        pdf_Yt <- desscherMixture(comb$mean, comb$sd, comb$probs, df_n$theta)
        # Esscher Distribution Yt
        cdf_Rt <- function(x) integrate(pdf_Yt, lower = df_n$lower, upper = x)$value
      }
      # First moment for solar radiation
      e_Rt <- function(lower, upper) integrate(function(x) df_n$f_x[[1]](x, df_n$Ct)*pdf_Yt(x), lower = lower, upper = upper)$value
    } else {
      # Mixture Pdf
      pdf_Yt <- function(x) dmixnorm(x, comb$mean, comb$sd, comb$probs)
      # Density for solar radiation
      pdf_Rt <- function(x) dsolarGHI(x, df_n$Ct, df_n$alpha, df_n$beta, pdf_Yt)

      # Esscher transform
      if (df_n$theta != 0) {
        pdf_Rt <- desscher(pdf_Rt, theta = df_n$theta, lower = df_n$lower, upper = df_n$upper)
        # Distribution for solar radiation
        cdf_Rt <- function(x) integrate(pdf_Rt, lower = df_n$lower, upper = x)$value
      } else {
        # Distribution for solar radiation
        cdf_Rt <- function(x) psolarGHI(x, df_n$Ct, model$transform$alpha, model$transform$beta, pdf_Yt)
      }
      # First moment for solar radiation
      e_Rt <- function(lower, upper) integrate(function(x) x*pdf_Rt(x), lower = lower, upper = upper)$value
    }

    # Option pricing
    if (put) {
      # Expected value (Put)
      df_n$Rt_plus <- e_Rt(df_n$lower, df_n$z_x)
      # Probability of exercise
      df_n$exercise <- cdf_Rt(df_n$z_x)
      # Option expected value
      df_n$premium <- df_n$strike*df_n$exercise - df_n$Rt_plus
      # Option type
      df_n$side <- "put"
    } else {
      # Expected value (Call)
      df_n$Rt_plus <-  e_Rt(df_n$z_x, df_n$upper)
      # Probability of exercise
      df_n$exercise <- 1 - cdf_Rt(df_n$z_x)
      # Option expected value
      df_n$premium <- df_n$Rt_plus - df_n$strike*df_n$exercise
      # Option type
      df_n$side <- "call"
    }
    # Expected value (GHI)
    df_n$Rt <- e_Rt(df_n$lower, df_n$upper)
    # Select only relevant variables
    df_n <- dplyr::select(df_n, Month, Day, side, premium, exercise, Rt_plus, Rt, strike)
    return(df_n)
  }

  # Model premium for each day and month
  df_month_day <- purrr::map2_df(seasonal_data$Month, seasonal_data$Day,
                                 ~pricing_month_day(seasonal_data, combinations, nmonth = .x, nday = .y, put = put))

  day_date <- paste0(ifelse(leap_year, "2020-", "2019-"), df_month_day$Month, "-", df_month_day$Day)
  df_month_day$n <- number_of_day(day_date)
  # Model premium aggregated for each Month
  df_month <- df_month_day %>%
    dplyr::group_by(Month, side) %>%
    dplyr::reframe(ndays = dplyr::n(),
                   premium = sum(premium),
                   daily_premium = premium/ndays,
                   exercise = mean(exercise),
                   Rt = sum(Rt),
                   Rt_plus = sum(Rt_plus),
                   strike = sum(strike))

  # Model premium aggregated by Year
  df_year <- df_month %>%
    dplyr::group_by(side) %>%
    dplyr::reframe(ndays = sum(ndays),
                   premium = sum(premium),
                   exercise = mean(exercise),
                   Rt = sum(Rt),
                   Rt_plus = sum(Rt_plus),
                   strike = sum(strike))

  # Rename columns to match target name
  rename_columns <- function(data, old_names, new_names){
    col_names <- colnames(data)
    for(i in 1:length(new_names)) {
      idx_old_names <- which(col_names == old_names[i])
      col_names[idx_old_names] <- new_names[i]
    }
    colnames(data) <- col_names
    return(data)
  }
  df_month_day <- rename_columns(df_month_day, c("Rt", "Rt_plus"), c(target, target_plus))
  df_month <- rename_columns(df_month, c("Rt", "Rt_plus"), c(target, target_plus))
  df_year <- rename_columns(df_year, c("Rt", "Rt_plus"), c(target, target_plus))

  structure(
    list(
      payoff_month_day = df_month_day,
      payoff_month = df_month,
      payoff_year = df_year
    ),
    class = c("solarOptionPayoff", "list")
  )
}
