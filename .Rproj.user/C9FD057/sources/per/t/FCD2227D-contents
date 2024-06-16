#' Normal Mixture Pdf
#'
#' Probability density function for a normal mixture with two components
#'
#' @param params parameters of the two components, (mu1, mu2, sd1, sd2, p)
#'
#' @return a function of x.
#'
#' @examples
#' params <- c(mu1 = -2, mu2 = 2, sd1 = 3, sd2 = 1, p = 0.5)
#' # Density function
#' pdf <- dnorm_mix(params)
#' # Distribution function
#' cdf <- dnorm_mix(params)
#' # Quantile function
#' q <- qnorm_mix(params)
#' # Random numbers
#' rnorm_mix(100, params)
#'
#' @rdname norm_mix
#' @name norm_mix
#' @aliases dnorm_mix
#' @aliases pnorm_mix
#' @aliases qnorm_mix
#' @aliases rnorm_mix
#' @export
dnorm_mix <- function(params) {
  # parameters
  mu1 = params[1]
  mu2 = params[2]
  sd1 = params[3]
  sd2 = params[4]
  p = params[5]
  # Density function
  function(x, log = FALSE){
    probs <- p*stats::dnorm(x, mean = mu1, sd = sd1) + (1-p)*stats::dnorm(x, mean = mu2, sd = sd2)
    if (log) {
      probs <- base::log(probs)
    }
    return(probs)
  }
}


#' @export
#' @rdname norm_mix
pnorm_mix <- function(params) {
  # parameters
  mu1 = params[1]
  mu2 = params[2]
  sd1 = params[3]
  sd2 = params[4]
  p = params[5]
  # Distribution function
  function(x, lower.tail = TRUE, log.p = FALSE){
    probs <- p*stats::pnorm(x, mean = mu1, sd = sd1, lower.tail = lower.tail)
    probs <- probs + (1-p)*stats::pnorm(x, mean = mu2, sd = sd2, lower.tail = lower.tail)
    if (log.p) {
      probs <- base::log(probs)
    }
    return(probs)
  }
}


#' @export
#' @rdname norm_mix
qnorm_mix <- function(params) {
  cdf <- pnorm_mix(params)
  loss_function <- function(p, lower.tail = lower.tail, log.p = log.p){
    function(x){
      (cdf(x, lower.tail = lower.tail, log.p = log.p) - p)^2
    }
  }
  function(p, lower.tail = TRUE, log.p = FALSE, interval = c(-100, 100)){
    x <- c()
    for(i in 1:length(p)){
      loss <- loss_function(p[i], lower.tail = lower.tail, log.p = log.p)
      x[i] <- suppressWarnings(optim(par = 0, loss)$par)
    }
    return(x)
  }
}


#' @export
#' @rdname norm_mix
rnorm_mix <- function(n, params){
  # parameters
  mu1 = params[1]
  mu2 = params[2]
  sd1 = params[3]
  sd2 = params[4]
  p = params[5]

  B <- rbinom(n, 1, prob = p)
  Z1 <- rnorm(n, mean = mu1, sd = sd1)
  Z2 <- rnorm(n, mean = mu2, sd = sd2)
  Z <- B*Z1 + (1-B)*Z2
  return(Z)
}


#' Fit Gaussian Mixture Pdf
#'
#' Fit the parameters for the density function of a Gaussian mixture with two components.
#'
#' @param x vector
#' @param params initial parameters
#' @param loss loss type. Can be `ml` for maximum likelihood or `kl` for kl_dist.
#'
#' @examples
#' params <- c(mu1 = -2, mu2 = 2, sd1 = 3, sd2 = 1, p = 0.5)
#' n1 <- rnorm(t_bar, mean = params[1], sd = params[3])
#' n2 <- rnorm(t_bar, mean = params[2], sd = params[4])
#' Z <- rbinom(t_bar, 1, params[5])
#' x <- Z*n1 + (1-Z)*n2
#' fit_dnorm_mix(x, params = init_params, loss = "ml")$par
#' fit_dnorm_mix(x, params = init_params, loss = "kl")$par
#'
#' @rdname fit_dnorm_mix
#' @name fit_dnorm_mix
#' @export
fit_dnorm_mix <- function(x, params, loss = "ml", na.rm = TRUE){
  # Loss type
  loss <- match.arg(loss, choices = c("ml", "kl"))
  # Loss function
  log_lik <- function(params){
    # Parameters
    mu1 = params[1]
    mu2 = params[2]
    sd1 = params[3]
    sd2 = params[4]
    p = params[5]
    # Ensure that probability is in (0,1)
    if(p > 0.99 | p < 0.01 | sd1 < 0 | sd2 < 0){ s
      return(NA_integer_)
    }
    # Mixture density
    pdf_mix <- dnorm_mix(params)
    if (loss == "kl") {
      # Empirical min and max
      min_x <- min(x, na.rm = na.rm)
      max_x <- max(x, na.rm = na.rm)
      # Empirical density
      ker <- density(x, from = min_x, to = max_x)
      # kl-distance
      loss <- kl_dist(ker$y, pdf_mix(ker$x))
    } else if (loss == "ml") {
      # Log-likelihood
      loss <- sum(pdf_mix(x, log = TRUE), na.rm = na.rm)
    }
    return(loss)
  }
  # Optimal parameters
  # fnscale = -1 to maximize (or use negative likelihood)
  opt <- optim(par = params, log_lik, control = list(maxit = 500000, fnscale = -1))
  return(opt)
}

#' Fit a monthly Gaussian Mixture Pdf
#'
#' Fit the monthly parameters for the density function of a Gaussian mixture with two components.
#'
#' @param x vector
#' @param date vector of dates
#' @param params initial parameters
#' @param loss loss type. Can be `ml` for maximum likelihood or `kl` for kl_dist.
#' @param match_moments when true the theoric moments will match the empirical ones
#'
#' @rdname fit_dnorm_mix.monthly
#' @name fit_dnorm_mix.monthly
#' @export

fit_dnorm_mix.monthly <- function(x, date, loss = "ml", match_moments = FALSE){
  i <- 1
  data <- dplyr::tibble(Month = lubridate::month(date), ut = x)
  # Normal Mixture Model
  NM_model <- list()
  for(i in 1:12){
    # Monthly data
    eps <- dplyr::filter(data, Month == i)$ut
    # Initial parameters
    p0 <- 0.5
    mu_10 <- mean(eps)
    mu_20 <- -mean(eps)
    sd_10 <- sd_20 <- sd(eps)
    init_params <- c(mu1 = -mu_10, mu2 = mu_20, sd1 = sd_10, sd2 = sd_20, p = p0)
    # Fitted model
    nm <- fit_dnorm_mix(eps, params = init_params, loss = loss)
    # Compute expected value
    e_u <-  nm$par[5]*nm$par[1] + (1 - nm$par[5])*nm$par[2]
    # Compute variance
    v_u <-  nm$par[5]*(nm$par[1]^2 + nm$par[3]^2) + (1 - nm$par[5])*(nm$par[2]^2 + nm$par[4]^2)
    # Compute sample moments
    e_u_hat <- mean(eps, na.rm = TRUE)
    v_u_hat <- var(eps, na.rm = TRUE)
    # Match exactly sample moments
    if (match_moments) {
      nm$par[2] <- (e_u_hat - nm$par[5]*nm$par[1])/(1 - nm$par[5])
      # Update expected value
      e_u <-  nm$par[5]*nm$par[1] + (1 - nm$par[5])*nm$par[2]
      v_2 <- (v_u_hat + e_u^2 - (nm$par[1]^2 + nm$par[3]^2)*nm$par[5] - nm$par[2]^2 + nm$par[5]*nm$par[2]^2)/(1 - nm$par[5])
      nm$par[4] <- sqrt(v_2)
      # Update variance
      v_u <-  nm$par[5]*(nm$par[1]^2 + nm$par[3]^2) + (1 - nm$par[5])*(nm$par[2]^2 + nm$par[4]^2)
    }
    # Update expected value
    e_u <-  nm$par[5]*nm$par[1] + (1 - nm$par[5])*nm$par[2]
    # Update variance
    v_u <-  nm$par[5]*(nm$par[1]^2 + nm$par[3]^2) + (1 - nm$par[5])*(nm$par[2]^2 + nm$par[4]^2)
    # Update log-likelihood
    nm$value <- sum(dnorm_mix(nm$par)(eps, log = TRUE))
    # Fitted parameters
    df_par <- dplyr::bind_cols(dplyr::bind_rows(nm$par), p2 = 1 - nm$par[5])
    colnames(df_par) <- c("mu1", "mu2", "sd1", "sd2", "p1", "p2")
    # Monthly data
    NM_model[[i]] <- dplyr::tibble(Month = i,
                                   df_par,
                                   loss = nm$value,
                                   nobs = length(eps),
                                   mean_loss = loss/nobs,
                                   e_x = e_u,
                                   v_x = v_u,
                                   e_x_hat = e_u_hat,
                                   v_x_hat = v_u_hat)
  }

  NM_model <- dplyr::bind_rows(NM_model) %>%
    dplyr::mutate(
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
      p_dw = 1 -p_up
    ) %>%
    dplyr::select(Month, mu_up:p_dw, loss, nobs, mean_loss,
                  e_x, v_x, e_x_hat, v_x_hat)

  return(NM_model)
}


#' Bootstrap Parameters of Gaussian Mixture
#'
#' (experimental) Bootstrap the parameters for the density function for a Gaussian mixture with two components.
#'
#' @param x vector
#' @param params initial parameters
#' @param B number of bootstraps
#' @param ci confidence interval for empirical quantiles
#' @param seed random seed
#' @param loss loss type. Can be `ml` for maximum likelihood or `kl` for kl_dist.
#' @param na.rm logical.
#' @examples
#' params <- c(mu1 = -2, mu2 = 2, sd1 = 3, sd2 = 1, p = 0.5)
#' n1 <- rnorm(t_bar, mean = params[1], sd = params[3])
#' n2 <- rnorm(t_bar, mean = params[2], sd = params[4])
#' Z <- rbinom(t_bar, 1, params[5])
#' x <- Z * n1 + (1-Z)*n2
#' boot_dnorm_mix(x, params = init_params,  B = 50, ci = 0.95, sample_perc = 0.8, loss = "ml")$par
#' boot_dnorm_mix(x, params = init_params,  B = 50, ci = 0.95, sample_perc = 0.8, loss = "kl")$par
#'
#' @rdname boot_dnorm_mix
#' @name boot_dnorm_mix
#' @export
boot_dnorm_mix <- function(x, params, B = 50, ci = 0.95, sample_perc = 0.8, loss = "kl", seed = 1, na.rm = TRUE){

  t_bar <- length(x)
  nobs_sample <- round(t_bar*sample_perc)
  par <- matrix(0, nrow = B, ncol = 5)
  opt_par <- params
  set.seed(seed)
  for(i in 1:B){
    x_sample <- x[sample(t_bar, nobs_sample)]
    opt_par <- fit_dnorm_mix(x_sample, params = params, loss = loss, na.rm = na.rm)$par
    par[i,] <- opt_par
    #opt_par <- opt_par*runif(5, 0.5, 1.5)
  }
  tibble(
    theta = c("mu1", "mu2", "sd1", "sd2", "p"),
    B = B,
    nobs = nobs_sample,
    loss = loss,
    e_theta = colMeans(par),
    sd_theta = apply(par, 2, sd, na.rm = na.rm),
    e_theta_lo = apply(par, 2, quantile, probs = 1 - ci),
    e_theta_up = apply(par, 2, quantile, probs = ci)
  )
}
