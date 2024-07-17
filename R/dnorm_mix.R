#' Gaussian mixture density
#'
#' Probability density function for a Gaussian mixture with two components.
#'
#' @param params parameters of the two components, (mu1,mu2,sd1,sd2,p)
#'
#' @return a density function.
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

#' Fit the Gaussian mixture parameters with Maximum likelihood.
#'
#' Fit the parameters for the density function of a Gaussian mixture with two components.
#'
#' @param x vector
#' @param params initial parameters.
#' @param match_moments logical. When `TRUE`, the parameters of the second distribution will be estimated such that
#' the empirical first two moments of `x` matches the theoretical Gaussian mixture moments.
#' @param na.rm logical. When `TRUE`, the default, `NA` values will be excluded from the computations.
#' @examples
#' t_bar <- 1000
#' params <- c(mu1 = -2, mu2 = 2, sd1 = 3, sd2 = 1, p = 0.5)
#' n1 <- rnorm(t_bar, mean = params[1], sd = params[3])
#' n2 <- rnorm(t_bar, mean = params[2], sd = params[4])
#' Z <- rbinom(t_bar, 1, params[5])
#' x <- Z*n1 + (1-Z)*n2
#' fit_dnorm_mix_ml(x, params = params)$par
#' fit_dnorm_mix_ml(x, params = params)$par
#'
#' @rdname fit_dnorm_mix_ML
#' @name fit_dnorm_mix_ML
#' @export
fit_dnorm_mix_ml <- function(x, params, match_moments = FALSE, na.rm = TRUE){
  # Loss function
  log_lik <- function(params){
    # Parameters
    mu1 = params[1]
    mu2 = params[2]
    sd1 = params[3]
    sd2 = params[4]
    p = params[5]
    # Ensure that probability is in (0,1)
    if(p > 0.99 | p < 0.01 | sd1 < 0 | sd2 < 0){
      return(NA_integer_)
    }
    # Mixture density
    pdf_mix <- dnorm_mix(params)
    # Log-likelihood
    loss <- sum(pdf_mix(x, log = TRUE), na.rm = na.rm)
    return(loss)
  }
  # Optimal parameters
  # fnscale = -1 to maximize (or use negative likelihood)
  opt <- optim(par = params, log_lik, control = list(maxit = 500000, fnscale = -1))

  # Match sample moments
  if (match_moments) {
    # Compute sample moments
    e_x_hat <- mean(x, na.rm = na.rm)
    v_x_hat <- var(x, na.rm = na.rm)
    # Match exactly sample moments
    if (match_moments) {
      opt$par[2] <- (e_x_hat - opt$par[5]*opt$par[1])/(1 - opt$par[5])
      # Compute theoric expected value
      e_x <-  opt$par[5]*opt$par[1] + (1 - opt$par[5])*opt$par[2]
      opt$par[4] <- sqrt((v_x_hat + e_x^2 - (opt$par[1]^2 + opt$par[3]^2)*opt$par[5] - opt$par[2]^2 + opt$par[5]*opt$par[2]^2)/(1 - opt$par[5]))
      # Update log-likelihood
      opt$value <- log_lik(opt$par)
    }
  }

  return(opt)
}

#' Fit the Gaussian mixture parameters with EM algorithm.
#'
#' Fit the parameters for the density function of a Gaussian mixture with two components.
#'
#' @param x vector
#' @param params initial parameters
#' @param absotol absolute level for convergence.
#' @param maxit maximum number of iterations.
#' @param match_moments logical. When `TRUE`, the parameters of the second distribution will be estimated such that
#' the empirical first two moments of `x` matches the theoretical Gaussian mixture moments.
#' @param na.rm logical. When `TRUE`, the default, `NA` values will be excluded from the computations.
#' @examples
#' t_bar <- 1000
#' params <- c(mu1 = -2, mu2 = 2, sd1 = 3, sd2 = 1, p = 0.5)
#' n1 <- rnorm(t_bar, mean = params[1], sd = params[3])
#' n2 <- rnorm(t_bar, mean = params[2], sd = params[4])
#' Z <- rbinom(t_bar, 1, params[5])
#' x <- Z*n1 + (1-Z)*n2
#' fit_dnorm_mix_em(x, params = params)$par
#' fit_dnorm_mix_em(x, params = params)$par
#'
#' @rdname fit_dnorm_mix_em
#' @name fit_dnorm_mix_em
#' @export
fit_dnorm_mix_em <- function(x, params = NULL, abstol = 1e-30, maxit = 50000, match_moments = FALSE, na.rm = FALSE){

  n <- length(x)
  # Empirical moments
  e_x_emp <- mean(x, na.rm = na.rm)
  v_x_emp <- var(x, na.rm = na.rm)
  sd_x_emp <- sqrt(v_x_emp)

  # Initialization
  log_likelihood <- 0
  previous_log_likelihood <- -Inf
  responsibilities <- matrix(0, nrow = n, ncol = 2)
  previous_params <- params
  # EM Algorithm
  for (iteration in 1:maxit) {
    # 1. E-step: Calculate the responsibilities
    for (i in 1:n) {
      responsibilities[i, 1] <- previous_params[5] * dnorm(x[i], previous_params[1], previous_params[3])
      responsibilities[i, 2] <- (1 - previous_params[5]) * dnorm(x[i], previous_params[2], previous_params[4])
    }

    responsibilities <- responsibilities/rowSums(responsibilities)

    # 2. M-step: Update the parameters
    n1 <- sum(responsibilities[, 1], na.rm = na.rm)
    n2 <- sum(responsibilities[, 2], na.rm = na.rm)
    ## Moments second distribution
    mu1 <- sum(responsibilities[, 1] * x)/n1 # means
    ## Std. deviations
    sigma1 <- sqrt(sum(responsibilities[, 1] * (x - mu1)^2, na.rm = na.rm)/(n1-1)) # std. deviation
    if (match_moments) {
      ## Match moments approach for (mu2, sigma2)
      mu2 <- (e_x_emp - p*mu1)/(1-p)
      sigma2 <- sqrt((v_x_emp - p*sigma1^2)/(1-p) - p*(mu1 - mu2)^2)
    } else {
      ## Moments second distribution
      mu2 <- sum(responsibilities[, 2] * x, na.rm = na.rm)/n2 # means
      sigma2 <- sqrt(sum(responsibilities[, 2] * (x - mu2)^2, na.rm = na.rm)/(n2-1)) # std. deviation
    }
    ## Bernoulli probability
    p <- n1/n

    # 3. Calculate the log-likelihood
    log_likelihood <- sum(log(rowSums(responsibilities*cbind(p * dnorm(x, mu1, sigma1), (1 - p) * dnorm(x, mu2, sigma2)))), na.rm = na.rm)

    # 4. Check for convergence
    check_convergence <- abs(log_likelihood - previous_log_likelihood) < abstol
    if (check_convergence) {
      break
    }
    # Update log-likelihood
    previous_log_likelihood <- log_likelihood
    # Update parameters
    previous_params <- c(mu1, mu2, sigma1, sigma2, p)
  }

  # Final classification
  B_hat <- ifelse(responsibilities[, 1] > responsibilities[, 2], 1, 0)
  x1_hat <- x[B_hat == 1]
  x2_hat <- x[B_hat == 0]
  # Fitted sequence
  df_fitted <- dplyr::tibble(B = B_hat, x = x, x1 = x*B, x2 = x*(1-B))

  # parameters with EM algorithm
  par <- previous_params
  names(par) <- c("mu1","mu2","sd1","sd2","p")

  # Parameters on estimated sub-samples
  par_hat <- c(mu1 = mean(x1_hat, na.rm = na.rm),
               mu2 = mean(x2_hat, na.rm = na.rm),
               sd1 = sd(x1_hat, na.rm = na.rm),
               sd2 = sd(x2_hat, na.rm = na.rm),
               p = mean(B_hat, na.rm = na.rm))

  # Theoric moments
  e_x_th = par[1]*par[5] + par[2]*(1-par[5])
  sd_x_th <- sqrt(par[5]*(1-par[5])*(par[1] - par[2])^2 + (par[3]^2)*par[5] + (par[4]^2)*(1-par[5]))
  # Moments on the classified series
  e_x_hat = par_hat[1]*par_hat[5] + par_hat[2]*(1-par_hat[5])
  sd_x_hat <- sqrt(par_hat[5]*(1-par_hat[5])*(par_hat[1] - par_hat[2])^2 + (par_hat[3]^2)*par_hat[5] + (par_hat[4]^2)*(1-par_hat[5]))
  # Moments
  df_moments <- dplyr::tibble(e_x = e_x_emp, sd_x = sd_x_emp,
                              e_x_th = e_x_th, sd_x_th = sd_x_th,
                              e_x_hat = e_x_hat, sd_x_hat = sd_x_hat)
  structure(
    list(
      par = par,
      par_hat = par_hat,
      value = previous_log_likelihood,
      mom =  df_moments,
      fitted = df_fitted
    )
  )
}


#' Fit a monthly Gaussian Mixture Pdf
#'
#' Fit the monthly parameters for the density function of a Gaussian mixture with two components.
#'
#' @param x vector
#' @param date vector of dates
#' @param loss loss type. Can be `ml` for maximum likelihood or `kl` for kl_dist.
#' @param algo character, type of algorithm. Can be `em` for Expectation-maximization or `ml` for maximum likelihood.
#' @param match_moments when true the theoric moments will match the empirical ones.
#'
#' @rdname fit_dnorm_mix_monthly
#' @name fit_dnorm_mix_monthly
#' @export
fit_dnorm_mix_monthly <- function(x, date, algo = c("em", "ml"), ...){

  algo <- match.arg(algo, choices = c("em", "ml"))
  data <- dplyr::tibble(Month = lubridate::month(date), eps = x)
  # Normal Mixture Model
  NM_model <- list()
  for(m in 1:12){
    # Monthly data
    eps <- dplyr::filter(data, Month == m)$eps
    # Initial parameters
    mu_0 <- mean(eps)
    sd_0 <- sd(eps)
    init_params <- c(mu1 = -mu_0, mu2 = mu_0, sd1 = sd_0, sd2 = sd_0, p = 0.5)
    if (algo == "em") {
      # Fitted model
      nm <- fit_dnorm_mix_em(eps, params = init_params, ...)
    } else {
      nm <- fit_dnorm_mix_ml(eps, params = init_params, ...)
    }

    # Compute the theoretical expected value
    e_u <-  nm$par[5]*nm$par[1] + (1 - nm$par[5])*nm$par[2]
    # Compute the theoric variance
    v_u <-  nm$par[5]*(nm$par[1]^2 + nm$par[3]^2) + (1 - nm$par[5])*(nm$par[2]^2 + nm$par[4]^2)
    # Compute the sample mean
    e_u_hat <- mean(eps, na.rm = TRUE)
    # Compute the sample variance
    v_u_hat <- var(eps, na.rm = TRUE)
    # Fitted parameters
    df_par <- dplyr::bind_cols(dplyr::bind_rows(nm$par), p2 = 1 - nm$par[5])
    colnames(df_par) <- c("mu1", "mu2", "sd1", "sd2", "p1", "p2")
    # Monthly data
    NM_model[[m]] <- dplyr::tibble(Month = m,
                                   df_par,
                                   loss = nm$value,
                                   nobs = length(eps),
                                   e_x = e_u,
                                   v_x = v_u,
                                   e_x_hat = e_u_hat,
                                   v_x_hat = v_u_hat)
  }

  # Reorder the variables
  NM_model <- dplyr::bind_rows(NM_model)
  NM_model <- dplyr::mutate(NM_model,
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
  NM_model <- dplyr::select(NM_model, Month, mu_up:p_dw, loss, nobs, e_x, v_x, e_x_hat, v_x_hat)
  return(NM_model)
}
