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
#' @param theta numeric, shift parameter for the mixture.
#' @param seed scalar integer, starting random seed.
#' @param quiet logical
#'
#' @examples
#' # Solar Model
#' model <- solarModel$new(spec)
#' model$fit()
#' scenario <- solarScenario_by(model, "2005-01-10", "2020-01-01", theta = 0, nsim = 4, by = "1 year")
#' # Plot
#' solarScenario_plot(scenario, nsim = 2)
#' # Solar Option
#' solarOption_scenario(model, scenario)
#' solarOption_historical(model)
#'
#' @rdname solarScenario_by
#' @name solarScenario_by
#' @keywords solarScenario
#' @note Version 1.0.1
#' @export
solarScenario_by <- function(model, from = "2010-01-01", to = "2011-01-01", by = "1 year", theta = 0, nsim = 1, seed = 1, quiet = FALSE){
  idx_date <- seq.Date(as.Date(from), as.Date(to), by = by)
  # Number of scenarios
  n.scenarios <- length(idx_date) - 1
  # Initialization
  df_sim <- list()
  df_res <- list()
  if (n.scenarios == 1) {
    # Specification
    simSpec <- solarScenario_spec$new(model, from = idx_date[1], to = idx_date[2], theta = theta, exclude_known = TRUE, quiet = TRUE)
    scenario_j <- solarScenario$new(simSpec, seed = seed)
    scenario_j$simulate_residuals(nsim = nsim)
    scenario_j$filter()
    df_sim[[1]] <- scenario_j$spec$simulations
    df_res[[1]] <- scenario_j$residuals
  } else {
    for(j in 1:n.scenarios){
      if (!quiet) {
        # To report progress
        pb <- txtProgressBar(min = 0,            # Minimum value of the progress bar
                             max = n.scenarios,  # Maximum value of the progress bar
                             style = 3,          # Progress bar style (also available style = 1 and style = 2)
                             width = 50,         # Progress bar width. Defaults to getOption("width")
                             char = "#")
        setTxtProgressBar(pb, j)
      }
      # Specification
      simSpec <- solarScenario_spec$new(model, from = idx_date[1], to = idx_date[2], theta = theta, exclude_known = TRUE, quiet = TRUE)
      scenario_j <- solarScenario$new(simSpec, seed = seed)
      scenario_j$simulate_residuals(nsim = nsim)
      scenario_j$filter()
      df_sim[[j]] <- scenario_j$spec$simulations
      df_res[[j]] <- scenario_j$residuals
      # Update seed
      seed <- seed + j - 1
    }
    if (!quiet) close(pb)
  }
  # Create an object to store all the simulations
  simSpec <- solarScenario_spec$new(model, from = from, to = to, theta = theta, exclude_known = TRUE, quiet = TRUE)
  simSpec$simulations <- purrr::flatten(df_sim)
  # Final object
  scenarios <- solarScenario$new(simSpec, seed = seed)
  scenarios$residuals <- dplyr::mutate(dplyr::bind_rows(df_res), filter = FALSE)
  return(scenarios)
}

#' Simulate
#'
#' Specification of scenarios of a `solarModel`
#' @examples
#' model <- solarModel$new(spec)
#' model$fit()
#' # Specification
#' spec <- solarScenario_spec$new(model)
#' solarScenario_residuals(spec, 4, 1)
#'
#' @rdname solarScenario_residuals
#' @name solarScenario_residuals
#' @keywords solarScenario
#' @note Version 1.0.2
#' @export
solarScenario_residuals <- function(spec, nsim = 1, seed = 1){
  # Random seed
  set.seed(seed)
  # Simulatioon dates
  dates <- spec$data_sim$date
  i_start <- spec$control$i_start

  # Initialize a dataset for storing simulations
  df_sim <- dplyr::tibble(date = dates[-c(1:(i_start-1))]) %>%
    dplyr::mutate(Year = lubridate::year(date), Month = lubridate::month(date)) %>%
    dplyr::group_by(Year, Month) %>%
    dplyr::mutate(n = dplyr::n()) %>%
    tidyr::nest(data_B = date, data_X = date)
  # Compute marginal probabilities once
  margprob <- spec$prob(df_sim$Month)
  # Reference location
  place <- spec$place

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
    sim_B <- bindata::rmvbin(n*nsim, margprob = margprob[m])
    colnames(sim_B) <- place
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
    colnames(sim_B) <- place
    colnames(sim_X) <- place
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
  df_sim <- dplyr::left_join(df_sim_X, df_sim_B, by = "nsim")
  df_sim <- dplyr::bind_cols(seed = seed, df_sim)

  return(df_sim)
}


#' Filter
#'
#' Specification of scenarios of a `solarModel`
#'
#' @examples
#' model <- solarModel$new(spec)
#' model$fit()
#' # Specification
#' spec <- solarScenario_spec$new(model)
#' residuals <- solarScenario_residuals(spec, 4, 1)
#' solarScenario_filter(spec, residuals)
#'
#' @rdname solarScenario_filter
#' @name solarScenario_filter
#' @keywords solarScenario
#' @note Version 1.0.2
#' @export
solarScenario_filter <- function(spec, residuals) {

  if (missing(residuals) || purrr::is_empty(residuals)) {
    stop("The slot `scenario$residuals` is empty! Run `simulate()` first.")
  }

  # Reference location
  place <- spec$place
  # Reference target variable
  target <- spec$target
  # Control
  i_start <- spec$control$i_start
  # ARMA model
  ARMA <- spec$ARMA
  intercept <- ARMA$intercept
  # AR order
  arOrder <- ARMA$arOrder
  phi <- ARMA$phi
  # ARMA order
  maOrder <- ARMA$maOrder
  theta <- ARMA$theta
  # GARCH model
  GARCH <- spec$GARCH
  omega <- GARCH$omega
  # ARCH order
  archOrder <- GARCH$archOrder
  alpha <- GARCH$alpha
  # GARCH order
  garchOrder <- GARCH$garchOrder
  beta <- GARCH$beta
  # Number of simulations
  nsim <- nrow(residuals)
  # Transform function
  transform <- spec$transform
  # Remove initial values
  dates <- spec$data_sim$date
  if (spec$control$exclude_known) {
    idx_exclude_known <- dates >= spec$control$from & dates <= spec$control$to
  } else {
    idx_exclude_known <- rep(TRUE, length(dates))
  }

  simulations <- vector("list", nsim)
  for (j in seq_len(nsim)) {
    df_sim <- spec$data_sim
    # Inject residuals for this scenario
    df_sim$z <- 0
    df_sim$z[i_start:nrow(df_sim)] <- residuals[j, ]$X[[1]][, place][[1]]
    df_sim$B[i_start:nrow(df_sim)] <- residuals[j, ]$B[[1]][, place][[1]]
    if (!spec$quiet) {
      message("Simulation: ", j, "/", nsim, " (",
              round(j / nsim * 100, 4), " %) \r", appendLF = FALSE)
    }

    # ----------------------------------------------------------------
    # Call C core for the ARMA-GARCH + mixture recursion
    # ----------------------------------------------------------------
    res <- .Call("solarScenario_filter_c",
                 as.numeric(df_sim$z),
                 as.numeric(df_sim$B),
                 as.numeric(df_sim$theta),          # mixture shift
                 as.numeric(df_sim$mu1),
                 as.numeric(df_sim$mu2),
                 as.numeric(df_sim$sd1),
                 as.numeric(df_sim$sd2),
                 as.numeric(df_sim$sigma_bar),
                 as.numeric(df_sim$sigma_uncond),
                 as.numeric(df_sim$Yt_tilde),
                 as.numeric(df_sim$eps),
                 as.numeric(df_sim$eps_tilde),
                 as.numeric(df_sim$sigma),
                 as.numeric(df_sim$Yt_tilde_hat),
                 as.numeric(omega),
                 as.numeric(alpha),
                 as.numeric(beta),
                 as.numeric(intercept),
                 as.numeric(phi),
                 as.numeric(theta),
                 as.integer(i_start))

    # res is a list: z, sigma, eps_tilde, eps, Yt_tilde, Yt_tilde_hat, B
    df_sim$z            <- res$z
    df_sim$sigma        <- res$sigma
    df_sim$eps_tilde    <- res$eps_tilde
    df_sim$eps          <- res$eps
    df_sim$Yt_tilde     <- res$Yt_tilde
    df_sim$Yt_tilde_hat <- res$Yt_tilde_hat
    df_sim$B            <- res$B

    # ----------------------------------------------------------------
    # Back-transform (still in R – cheap compared to the recursion)
    # ----------------------------------------------------------------
    df_sim$Yt <- df_sim$Yt_bar + df_sim$Yt_tilde + df_sim$Yt_tilde_uncond
    df_sim$Xt <- transform$iX_prime(transform$iY(df_sim$Yt))
    df_sim[[target]] <- transform$iX(df_sim$Xt, df_sim$Ct)

    df_sim <- dplyr::select(
      df_sim,
      -mu1, -mu2, -sd1, -sd2, -p1, -Ct,
      -Yt_bar, -sigma_bar, -Yt_tilde_hat,
      -Yt_tilde_uncond, -sigma_uncond
    )
    simulations[[j]] <- dplyr::bind_cols(nsim = j, df_sim[idx_exclude_known, ])
  }

  simulations
}


#' Compute the Value at Risk from simulated values
#'
#' @param scenarios An object of the class `solarScenario`
#' @param alpha Confidence level for the VaR
#'
#' @examples
#' # solarModel
#' model <- solarModel$new(spec)
#' model$fit()
#' # Scenarios
#' scenario <- solarScenario(model, "2016-01-01", "2017-01-01", nsim = 10, by = "1 month")
#' # VaR
#' solarScenario_VaR(scenario, 0.05)
#' @rdname solarScenario_VaR
#' @name solarScenario_VaR
#' @keywords solarScenario
#' @note Version 1.0.0.
#' @export
solarScenario_VaR <- function(scenario, alpha = 0.05){
  # Extract simulated data
  simulated <- scenario$spec$scenarios
  # Compute the VaR
  simulated <- simulated %>%
    tidyr::unnest(cols = c("data")) %>%
    dplyr::group_by(date) %>%
    dplyr::mutate(VaR_alpha = quantile(GHI, probs = alpha)) %>%
    dplyr::select(date, VaR_alpha) %>%
    dplyr::ungroup() %>%
    unique()
  # Empiric values
  data <- dplyr::left_join(simulated, scenario$spec$emp, by = "date")
  # Violations of VaR
  data$et <- ifelse(data$GHI <= data$VaR_alpha, 1, 0)
  # Select only relevant variables
  data <- dplyr::select(data, date, Year, Month, Day, VaR_alpha, et)
  # Output structure
  structure(
    list(
      data = data,
      alpha = alpha
    )
  )
}

#' Plot scenarios from a solarScenario object
#'
#' @examples
#' # solarModel
#' model <- solarModel$new(spec)
#' model$fit()
#' # Scenarios
#' scenario <- solarScenario(model, "2016-01-01", "2017-01-01", nsim = 10, by = "1 month")
#' # Plot
#' solarScenario_plot(scenario)
#' @export
solarScenario_plot <- function(scenario, target = "GHI", nsim = 10, empiric = TRUE, ci = 0.05){
  # Empirical data
  df_emp <- scenario$spec$emp
  # Simulated data
  sim <- scenario$spec$scenarios
  df_sim <- sim
  df_sim <- dplyr::mutate(df_sim, data = purrr::map(data, ~.x[1:nsim,]))
  df_sim <- tidyr::unnest(df_sim, cols = c("data"))
  # Rename target columm
  df_sim$target <- df_sim[[target]]
  df_emp$target <- df_emp[[target]]

  # Plot scenarios
  plot_scenario <- ggplot()+
    geom_line(data = df_sim, aes(date, target, group = nsim, color = "simulated"), alpha = 0.1)

  # Add simulated value at risk
  if (!is.na(ci) & ci > 0) {
    VaR <- tidyr::unnest(sim, cols = c("data"))
    # Rename target columm
    VaR$target <- VaR[[target]]
    VaR <- VaR %>%
      group_by(date) %>%
      mutate(VaR_up = quantile(target, probs = 1-ci),
             VaR_dw = quantile(target, probs = ci)) %>%
      select(date, VaR_up, VaR_dw)

    plot_scenario <- plot_scenario +
      geom_line(data = VaR, aes(date, VaR_up, color = "VaR_up"), linewidth = 0.5)+
      geom_line(data = VaR, aes(date, VaR_dw, color = "VaR_dw"), linewidth = 0.5)
  }

  # Add realized GHI
  if (empiric) {
    plot_scenario <- plot_scenario +
      geom_point(data = df_emp, aes(date, target), color = "black", size = 1.5)+
      geom_point(data = df_emp, aes(date, target, color = "empiric"), size = 0.8)
  }
  vals_colors <- c(empiric = "orange", VaR_up = "red",VaR_dw = "red", simulated = "black")
  vals_labels <- c(empiric = "Realized", VaR_up = paste0("VaR (", (ci)*100, " %)"),
                   VaR_dw = paste0("VaR (", (1-ci)*100, " %)"),  simulated = paste0(nsim, " simulations"))
  plot_scenario+
    theme_bw()+
    theme(legend.position = "top")+
    labs(color = NULL, y = NULL, x = NULL)+
    scale_color_manual(values = vals_colors, labels = vals_labels)

}
