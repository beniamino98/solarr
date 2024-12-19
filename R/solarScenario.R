#' Simulate multiple scenarios
#'
#' Simulate multiple scenarios of solar radiation with a `solarModel` object.
#'
#' @param model object with the class `solarModel`. See the function \code{\link{solarModel}} for details.
#' @param from character, start Date for simulations in the format `YYYY-MM-DD`.
#' @param to character, end Date for simulations in the format `YYYY-MM-DD`.
#' @param by character, steps for multiple scenarios, e.g. `1 day` (day-ahead simulations), `15 days`, `1 month`, `3 months`, ecc.
#' For each step are simulated `nsim` scenarios.
#' @param nsim integer, number of simulations.
#' @param theta numeric, Esscher parameter.
#' @param seed scalar integer, starting random seed.
#' @param quiet logical
#'
#' @examples
#' model <- Bologna
#' scen <- solarScenario(model, "2024-01-01", "2025-12-31")
#'
#' @rdname solarScenario
#' @name solarScenario
#' @export
solarScenario <- function(model, from = "2010-01-01", to = "2011-01-01", by = "1 year", theta = 0, nsim = 1, seed = 1, quiet = FALSE){

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
    simSpec <- solarScenario_filter(simSpec)
    df_emp <- dplyr::bind_rows(df_emp, simSpec$emp)
    scenarios[[j]] <- dplyr::bind_cols(seed = seed, dplyr::bind_rows(simSpec$simulations))
    seed <- seed + j - 1
  }
  if (!quiet) close(pb)

  simSpec$emp <- df_emp[!duplicated(df_emp),]
  simSpec$simulations <- scenarios
  return(as_solarScenario(simSpec))
}

#' Specification of a solar scenario
#'
#' @inheritParams solarScenario
#' @param exclude_known when true the two starting points (equals for all the simulations) will be excluded from the output.
#'
#' @examples
#' model <- Bologna
#' simSpec <- solarScenario_spec(model)
#'
#' @rdname solarScenario_spec
#' @name solarScenario_spec
#' @export
solarScenario_spec <- function(model, from = "2010-01-01", to = "2010-12-31", theta = 0, exclude_known = FALSE, quiet = FALSE){
  # Extract informations
  data <- model$data
  place <- model$place
  # Number of lags to consider
  i_start <- model$control$mean.model$arOrder+1
  i_start <- max(c(model$control$variance.model@model$maxOrder, i_start))
  # Initial date
  from <- as.Date(from)
  # End date
  to <- as.Date(to)

  # Selected columns
  cols_emp <- c("date", "n", "Year", "Month", "Day", "GHI", "clearsky",
                 "Xt", "Yt", "Yt_tilde", "Yt_tilde_hat", "eps", "eps_tilde", "sigma", "u", "u_tilde", "B", "z1", "z2")
  cols_sim <- c(cols_emp, "Ct", "Yt_bar", "GHI_bar", "sigma_bar", "Yt_tilde_uncond", "sigma_m", "mu1", "mu2", "sd1", "sd2", "p1")

  # Initialize a dataset
  max_date_from <- max(data$date)
  max_date_to <- max_date_from - i_start
  if (max_date_to >= to) {
    df_emp <- dplyr::filter(data, date >= (from - lubridate::days(i_start-1)) & date <= to)
    df_emp <- dplyr::bind_cols(place = place, df_emp)
  } else if (max_date_to >= from & max_date_from >= from) {
    df_emp <- dplyr::filter(data, date >= (from - lubridate::days(i_start-1)))
    df_new_emp <- dplyr::tibble(date = seq.Date(max(df_emp$date) + 1, to, by = "1 day"))
    df_emp <- dplyr::bind_rows(df_emp, df_new_emp)
    df_emp <- dplyr::mutate(df_emp,
                            Year = lubridate::year(date),
                            Month = lubridate::month(date),
                            Day = lubridate::day(date))
    df_emp$n <- solarr::number_of_day(df_emp$date)
    # Add seasonal variables
    df_emp <- df_emp[, cols_emp]
    df_emp <- dplyr::left_join(df_emp, model$seasonal_data, by = c("Month", "Day", "n"))
  } else {
    msg <- paste0("The maximum date for starting a simulation is: ", max_date_from)
    if (!quiet) warning(msg)
    return(NULL)
  }

  # Initialize simulation dataset
  df_sim <- df_emp[, cols_sim]
  # Initialize lambda
  df_sim$theta <- theta
  # Arch parameters
  archOrder <- model$GARCH$p
  # Garch parameters
  garchOrder <- model$GARCH$q
  # Garch next step function
  GARCH_next_step <- GARCH_pq_next_step(model$GARCH$omega, model$GARCH$alpha, model$GARCH$beta)
  # Filter df_emp to be in [from - to] dates
  df_emp <- df_emp[, cols_emp]
  if (exclude_known) {
    df_emp <- dplyr::filter(df_emp, date >= from & date <= to)
  }

  # Output structure
  structure(
    list(
      # Base dataset for simulations
      sim = df_sim,
      # Dataset with empirical data
      emp = df_emp,
      # Model place
      place = model$place,
      # Coordinates
      coords = model$coords,
      # Model target
      target = model$target,
      # Model transform
      transform = model$transform,
      # Monthly probability
      p_up = function(nmonth){model$NM_model$parameters$p1[nmonth]},
      # AR model
      AR_model_Yt = model$AR_model_Yt,
      # GARCH model
      arOrder = model$control$mean.model$arOrder,
      archOrder = archOrder,
      garchOrder = garchOrder,
      GARCH_next_step = GARCH_next_step,
      # Number of lags for simulations
      i_start = i_start,
      # Other info
      exclude_known = exclude_known,
      from = from,
      to = to,
      seed = NA,
      nsim = NA,
      quiet = quiet,
      # List to store residuals and simulations
      residuals = list(),
      simulations = list()
    ),
    class = c("solarScenarioSpec", "list")
  )
}


#' Simulate residuals for a `solarScenario_spec`
#'
#' @inheritParams solarScenario
#' @param simSpec object with the class `solarScenario_spec`. See the function \code{\link{solarScenario_spec}} for details.
#'
#' @examples
#' model <- Bologna
#' simSpec <- solarScenario_spec(model)
#' simSpec <- solarScenario_residuals(simSpec, nsim = 10)
#' @rdname solarScenario_residuals
#' @name solarScenario_residuals
#' @export
solarScenario_residuals <- function(simSpec, nsim = 1, seed = 1){
  # Random seed
  set.seed(seed)
  # Initialize the sequence of dates
  # Initialize a dataset for storing simulations
  df_sim <- dplyr::tibble(date = simSpec$sim$date[-c(1:(simSpec$i_start-1))]) %>%
    dplyr::mutate(Year = lubridate::year(date), Month = lubridate::month(date)) %>%
    dplyr::group_by(Year, Month) %>%
    dplyr::mutate(n = dplyr::n()) %>%
    tidyr::nest(data_B = date, data_X = date)

  m <- 1
  # Simulating scenarios
  for(m in 1:nrow(df_sim)) {
    # Number of days of the month
    n <- df_sim$n[m]
    # Complete Normal simulation
    sim_x12 <- mvtnorm::rmvnorm(n*nsim, mean = rep(0, 2), sigma = diag(c(1,1)))
    colnames(sim_x12) <- c("x1_1", "x2_1")
    sim_x12 <- dplyr::as_tibble(sim_x12)
    # 2) Bernoulli simulation
    # sim_B <- rbinom(n*nsim, 1, prob = simSpec$p_up(m))
    sim_B <- bindata::rmvbin(n*nsim, margprob = simSpec$p_up(df_sim$Month[m]))
    colnames(sim_B) <- simSpec$place
    sim_B <- dplyr::as_tibble(sim_B)
    # Slit the mixture simulations
    sim_x1 <- dplyr::select(sim_x12, dplyr::contains("x1_"))
    sim_x2 <- dplyr::select(sim_x12, dplyr::contains("x2_"))
    sim_X <- sim_x1
    for(i in 1:nrow(sim_B)){
      for(j in 1:ncol(sim_B)){
        sim_X[i,j] <- ifelse(sim_B[i,j] == 1, sim_x1[i,j], sim_x2[i,j])
      }
    }
    colnames(sim_B) <- simSpec$place
    colnames(sim_X) <- simSpec$place
    # Structure binomial simulations
    sim_B <- dplyr::bind_cols(scenario = rep(1:n, nsim), j = rep(1:nsim, n),  sim_B)
    sim_B <- tidyr::nest(dplyr::group_by(sim_B, scenario))
    sim_B <- dplyr::mutate(sim_B, data = purrr::map(data, ~dplyr::mutate(.x, j = 1:nsim)))
    # Structure normal simulations
    sim_X <- dplyr::bind_cols(scenario = rep(1:n, nsim), j = rep(1:nsim, n),  sim_X)
    sim_X <- dplyr::group_by(sim_X, scenario) %>% tidyr::nest()
    sim_X <- dplyr::mutate(sim_X, data = purrr::map(data, ~dplyr::mutate(.x, j = 1:nsim)))
    # Store data
    df_sim$data_B[[m]] <- dplyr::bind_cols(df_sim$data_B[[m]], sim_B)
    df_sim$data_X[[m]] <- dplyr::bind_cols(df_sim$data_X[[m]], sim_X)
  }
  # Structure binomial simulations
  df_sim_X <- dplyr::bind_rows(df_sim$data_X) %>%
    tidyr::unnest(cols = "data") %>%
    dplyr::select(-scenario) %>%
    dplyr::group_by(j) %>%
    dplyr::rename(nsim = "j") %>%
    tidyr::nest()%>%
    dplyr::rename(X = "data")
  # Structure normal simulations
  df_sim_B <- dplyr::bind_rows(df_sim$data_B) %>%
    tidyr::unnest(cols = "data") %>%
    dplyr::select(-scenario) %>%
    dplyr::group_by(j) %>%
    dplyr::rename(nsim = "j") %>%
    tidyr::nest() %>%
    dplyr::rename(B = "data")
  # Create a unique dataset
  df_sim <- dplyr::left_join( df_sim_X, df_sim_B, by = "nsim")
  df_sim <- dplyr::bind_cols(seed = seed, df_sim)
  simSpec$residuals <- df_sim
  simSpec$nsim <- nsim
  simSpec$seed <- seed
  return(simSpec)
}

#' Simulate trajectories from a a `solarScenario_spec`
#'
#' @inheritParams solarScenario
#' @inheritParams solarScenario_residuals
#'
#' @examples
#' model <- Bologna
#' simSpec <- solarScenario_spec(model, from = "2023-01-01", to = "2023-12-31")
#' simSpec <- solarScenario_residuals(simSpec, nsim = 1, seed = 3)
#' simSpec <- solarScenario_filter(simSpec)
#' # Empiric data
#' df_emp <- simSpec$emp
#' # First simulation
#' df_sim <- simSpec$simulations[[1]]
#' ggplot()+
#' geom_line(data = df_emp, aes(date, GHI))+
#' geom_line(data = df_sim, aes(date, GHI), color = "red")
#'
#' @rdname solarScenario_filter
#' @name solarScenario_filter
#' @export
solarScenario_filter <- function(simSpec){

  if (purrr::is_empty(simSpec$residuals)) {
    stop("The slot `simSpec$residuals` is empty! Consider running `simSpec <- solarScenario_residuals(simSpec)` before!")
  }
  # Number of lags to consider
  i_start <- simSpec$i_start
  # Number of simulations
  nsim <- nrow(simSpec$residuals)
  simulations <- list()
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
    for(i in i_start:nrow(df_sim)){
      # Simulated GARCH standard deviation
      df_sim$sigma[i] <- simSpec$GARCH_next_step(df_sim$eps_tilde[(i-simSpec$archOrder):(i-1)], df_sim$sigma[(i-simSpec$garchOrder):(i-1)])
      # Simulated Yt_tilde
      df_sim$Yt_tilde_hat[i] <- predict(simSpec$AR_model_Yt, newdata = df_sim[(i-(simSpec$arOrder+1)):i,])[i_start]
      # Simulated standardized monthly normal mixture
      df_sim$u_tilde[i] <- (df_sim$mu1[i] + df_sim$sd1[i]*df_sim$z[i])*df_sim$B[i] + (df_sim$mu2[i] + df_sim$sd2[i]*df_sim$z[i])*(1-df_sim$B[i])
      # IF theta != 0 the probabilities has to be distorted and B and u_tilde re-simulated
      if (df_sim$theta[1] != 0) {
        # Distort probability according to Esscher parameter
        params <- c(df_sim$mu_up[i], df_sim$mu_dw[i], df_sim$sd_up[i], df_sim$sd_dw[i], df_sim$p_up[i])
        df_sim$p_up[i] <- solarEsscher_probability(params, df_n = df_sim[i,], df_sim$theta[i])
        # Simulated bernoulli jump
        df_sim$B[i] <- rbinom(1, 1, df_sim$p_up[i])
        # Simulated standardized monthly normal mixture
        df_sim$u_tilde[i] <- (df_sim$mu1[i] + df_sim$sd1[i]*df_sim$z[i])*df_sim$B[i] + (df_sim$mu2[i] + df_sim$sd2[i]*df_sim$z[i])*(1-df_sim$B[i])
        # Simulated Esscher parameter
        df_sim$theta[i] <- df_sim$theta[i]*(df_sim$sigma[i]*df_sim$sigma_bar[i])^2*(df_sim$sd1[i]^2*df_sim$B[i] + df_sim$sd2[i]^2*(1-df_sim$B[i]))
      }
      # Simulated GARCH residuals
      df_sim$u[i] <- df_sim$sigma_m[i]*df_sim$u_tilde[i]
      # Simulated deseasonalized residuals
      df_sim$eps_tilde[i] <- df_sim$sigma[i]*df_sim$u[i]
      # Simulated AR residuals
      df_sim$eps[i] <- df_sim$eps_tilde[i]*df_sim$sigma_bar[i]
      # Simulated Yt_tilde
      df_sim$Yt_tilde[i] <- df_sim$Yt_tilde_hat[i] + df_sim$eps[i]
      # Simulated Yt
      df_sim$Yt[i] <- df_sim$Yt_bar[i] + df_sim$Yt_tilde[i] + df_sim$Yt_tilde_uncond[i] + df_sim$theta[i]
      # Simulated Xt
      df_sim$Xt[i] <- simSpec$transform$iY(df_sim$Yt[i])
      # Simulated GHI
      df_sim[[simSpec$target]][i] <- simSpec$transform$GHI(df_sim$Xt[i], df_sim$Ct[i])
    }
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


#' Extract and structure simulations from a `solarScenarioSpec`
#'
#' @inheritParams solarScenario_residuals
#' @rdname as_solarScenario
#' @name as_solarScenario
#' @keywords solarScenario internal
#' @export
as_solarScenario <- function(simSpec) {

  if (purrr::is_empty(simSpec$simulations)) {
    stop("The slot `simSpec$simulations` is empty! Consider running `simSpec <- solarScenario_filter(simSpec)` before!")
  }

  df_sim <- dplyr::bind_rows(simSpec$simulations) %>%
    dplyr::group_by(date, Year, Month, Day) %>%
    tidyr::nest() %>%
    dplyr::ungroup()

  structure(
    list(
      emp = simSpec$emp,
      sim = df_sim,
      target = simSpec$target
    ),
    class = c("solarScenario", "list")
  )
}



