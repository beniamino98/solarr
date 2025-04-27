#' Gaussian mixture
#'
#' Fit the parameters of a gaussian mixture with k-components.
#'
#' @examples
#' means = c(0,0.5,2)
#' sd = rep(1, 3)
#' p = c(0.2, 0.3, 0.5)
#' # Grid
#' grid <- seq(-4, 4, 0.01)
#' plot(dmixnorm(grid, means, sd, p))
#' # Simulated sample
#' x <- rmixnorm(5000, means, sd, p)
#' # Gaussian mixture model
#' gm <- gaussianMixture$new(components=3)
#' # Fit the model
#' gm$fit(x$X)
#' # EM-algo
#' gm$EM(x$X)
#' # Model parameters
#' gm$coefficients
#' # Fitted series
#' gm$fitted
#' # Theoric moments
#' gm$moments
#' gm$update(means = c(-2, 0, 2))
#' @rdname gaussianMixture
#' @name gaussianMixture
#' @note Version 1.0.0
#' @export
#self <- gm$.__enclos_env__$self
#private <- gm$.__enclos_env__$private
gaussianMixture <- R6::R6Class("gaussianMixture",
                               public = list(
                                 #' @field maxit Integer, maximum number of iterations.
                                 maxit = 5000,
                                 #' @field abstol Numeric, absolute level for convergence.
                                 abstol = 1e-08,
                                 #' @field components Integer, number of mixture components.
                                 components = 2,
                                 #' @description
                                 #' Initialize a gaussianMixture object
                                 #' @param components Integer, number of components.
                                 #' @param maxit (`integer(1)`)\cr
                                 #'   Numeric, maximum number of iterations.
                                 #' @param abstol (`numeric(1)`) Numeric, absolute level for convergence.
                                 initialize = function(components = 2, maxit = 5000, abstol = 1e-08){
                                   # Control parameters
                                   self$maxit <- maxit
                                   self$abstol <- abstol
                                   # Mixture components
                                   self$components <- components

                                   # Initialize the means parameters
                                   init_means <- seq(-3, 3, length.out = components)
                                   names(init_means) <- paste0("mu_", 1:components)
                                   # Initialize the std. deviations parameters
                                   init_sd <- rep(1, components)
                                   names(init_sd) <- paste0("sd_", 1:components)
                                   # Initialize the slots with the parameters
                                   init_p <- rep(1/components, components)
                                   names(init_p) <- paste0("p_", 1:components)
                                   # Update private slots
                                   private$..means = init_means
                                   private$..sd = init_sd
                                   private$..p = init_p/sum(init_p)
                                 },
                                 #' @description
                                 #' Compute the log-likelihood
                                 #' @param x vector
                                 #' @param params Optional. Named list with mixture parameters.
                                 logLik = function(x, params){
                                   if (missing(params)) {
                                     params <- self$coefficients
                                   }
                                   # Calculate the likelihoods
                                   likelihoods <- matrix(NA, nrow = length(x), ncol = self$components)
                                   for(k in 1:self$components){
                                     likelihoods[,k] <- params$p[k] * dnorm(x, params$means[k], params$sd[k])
                                   }
                                   # Calculate the total log-likelihood
                                   log_likelihood <- sum(log(rowSums(likelihoods)), na.rm = TRUE)
                                   return(log_likelihood)
                                 },
                                 #' @description
                                 #' Compute the posterior probabilities (E-step)
                                 #' @param x vector
                                 #' @param params a list of mixture parameters
                                 E_step = function(x, params){
                                   if (missing(x)) {
                                     x <- private$x
                                   }
                                   if (missing(params)) {
                                     params <- self$coefficients
                                   }
                                   responsabilities <- matrix(0, nrow = length(x), ncol = self$components)
                                   for(k in 1:self$components){
                                     responsabilities[, k] <- params$p[k] * dnorm(x, params$means[k], params$sd[k])
                                   }
                                   # Normalize the posterior probabilities
                                   responsabilities <- apply(responsabilities, 2, function(x) ifelse(is.na(x)|is.nan(x), 0, x))
                                   responsabilities <- responsabilities / rowSums(responsabilities)
                                   colnames(responsabilities) <- paste0("B", 1:self$components)
                                   return(dplyr::as_tibble(responsabilities))
                                 },
                                 #' @description
                                 #' Classify the time series in its components
                                 #' @param x vector
                                 classify = function(x){
                                   if (missing(x)) {
                                     x <- private$x
                                   }
                                   # Optimal parameters
                                   params <- self$coefficients
                                   # Number of observations
                                   n <- length(x)
                                   # E-step: posterior probabilities
                                   responsabilities <- self$E_step(x)
                                   # Classification of each component
                                   x_hat <- matrix(0, nrow = n, ncol = self$components)
                                   B_hat <- matrix(0, nrow = n, ncol = self$components)
                                   classification <- c()
                                   uncertanty <- c()
                                   for(i in 1:n) {
                                     classification[i] <- which.max(responsabilities[i,])
                                     uncertanty[i] <- 1-max(responsabilities[i,])
                                     B_hat[i, classification[i]] <- 1
                                     x_hat[i,] <- B_hat[i,]*x[i]
                                   }
                                   colnames(B_hat) <- paste0("B", 1:self$components)
                                   colnames(x_hat) <- paste0("x", 1:self$components)
                                   B_hat <- dplyr::as_tibble(B_hat)
                                   x_hat <- dplyr::as_tibble(x_hat)
                                   z_hat <- dplyr::as_tibble((x_hat-params$means)/params$sd*B_hat)
                                   colnames(z_hat) <- paste0("z", 1:self$components)
                                   # Output
                                   dplyr::bind_cols(classification = classification, B_hat, x_hat, z_hat,
                                                    uncertanty = uncertanty)
                                 },
                                 #' @description
                                 #' Fit the parameters with mclust package
                                 #' @param x vector
                                 #' @param weights observations weights, if a weight is equal to zero the observation is excluded, otherwise is included with unitary weight.
                                 #' When `missing` all the available observations will be used.
                                 fit = function(x, weights, B = 50, method = "mixtools"){
                                   # Number of observations
                                   n <- length(x)
                                   # Weights
                                   if (missing(weights)){
                                     w <- rep(1, n)
                                   } else {
                                     w <- ifelse(weights == 0, 0, 1)
                                   }
                                   if (!purrr::is_empty(w[is.na(x)])){
                                     w[is.na(x)] <- 0
                                   }
                                   # Add time series
                                   private$x <- x
                                   private$w <- w
                                   if (B > 0){
                                     clust <- GM_fit_rob(x[w!=0], B = B, method = method, components = self$components, maxit = self$maxit)
                                   } else {
                                     clust <- GM_fit(x[w!=0], method = method, components = self$components, maxit = self$maxit)
                                   }
                                   # Update the parameters
                                   private$..means <- clust$means
                                   private$..sd <- clust$sd
                                   private$..p <- clust$p
                                   # Assign a name to the parameters
                                   names(private$..means) <- paste0("mu", 1:self$components)
                                   names(private$..sd) <- paste0("sd", 1:self$components)
                                   names(private$..p) <- paste0("p", 1:self$components)

                                   # E-step: posterior probabilities
                                   private$responsabilities <- self$E_step(x)
                                   # Calculate the log-likelihood
                                   self$update_logLik()
                                   # Final classification of each component
                                   private$..fitted <- self$classify(x)
                                   # Compute empiric parameters
                                   self$update_empiric_parameters()
                                   # Compute hessian matrix
                                   self$Hessian()
                                 },
                                 #' @description
                                 #' Fit the parameters with EM-algorithm
                                 #' @param x vector
                                 #' @param weights observations weights, if a weight is equal to zero the observation is excluded, otherwise is included with unitary weight.
                                 #' When `missing` all the available observations will be used.
                                 EM = function(x, weights){
                                   if (missing(x)){
                                     x <- private$x
                                   }
                                   if (missing(weights)){
                                     w <- private$w
                                   } else {
                                     w <- ifelse(weights == 0, 0, 1)
                                   }
                                   w[is.na(x)] <- 0
                                   # Initialization
                                   n <- length(x)
                                   # Update n parameter
                                   n_w <- length(x[w != 0])
                                   # Initialization
                                   log_likelihood <- 0
                                   previous_log_likelihood <- -Inf
                                   prev_responsabilities <- matrix(0, nrow = n, ncol = self$components)
                                   previous_params <- list(means = self$means, sd = self$sd, p = self$p)
                                   # Initialize parameters
                                   if (any(c(purrr::is_empty(previous_params$means), is.na(previous_params$means)))) {
                                     previous_params <- list()
                                     previous_params$p <- rep(1/self$components, self$components)
                                     previous_params$means <- quantile(x[w != 0], seq(0.2, 0.8, length.out = self$components))
                                     previous_params$sd <- rep(sd(x[w != 0]), self$components)
                                   }
                                   previous_params
                                   # EM Algorithm
                                   for (iteration in 1:self$maxit){
                                     # E-step: posterior probabilities
                                     responsabilities <- self$E_step(x, previous_params)
                                     # Optimal parameters
                                     params <- previous_params
                                     # M-step: Update the parameters
                                     for(k in 1:self$components){
                                       # Normalizing factor for each group
                                       n_k <- sum(w*responsabilities[, k], na.rm = TRUE)
                                       # Mean parameter k-component
                                       params$means[k] <- sum(responsabilities[, k]*w*x, na.rm = TRUE)/n_k
                                       # Std. deviation k-component
                                       params$sd[k] <- sqrt(sum(responsabilities[, k]*w*(x - params$means[k])^2, na.rm = TRUE)/n_k)
                                       # Probability k-component
                                       params$p[k] <- n_k/n_w
                                     }
                                     # Check divergence
                                     if (any(params$p > 0.98)) {
                                       warning("Probs > 0.98 break")
                                       params <- previous_params
                                       break
                                     }
                                     # Calculate the log-likelihood
                                     log_likelihood <- self$logLik(x[w!=0], params)
                                     # Check for convergence
                                     stop_condition <- abs(log_likelihood - previous_log_likelihood) < self$abstol
                                     if (stop_condition) {
                                       break
                                     } else {
                                       # Update log-likelihood
                                       previous_log_likelihood <- log_likelihood
                                       # Update parameters
                                       previous_params <- params
                                     }
                                   }
                                   # Update the parameters
                                   private$..means <- params$means
                                   private$..sd <- params$sd
                                   private$..p <- params$p
                                   # Assign a name to the parameters
                                   names(private$..means) <- paste0("mu", 1:self$components)
                                   names(private$..sd) <- paste0("sd", 1:self$components)
                                   names(private$..p) <- paste0("p", 1:self$components)
                                   # Add time series
                                   private$x <- x
                                   private$w <- w
                                   # E-step: posterior probabilities
                                   private$responsabilities <- self$E_step(x)
                                   # Calculate the log-likelihood
                                   self$update_logLik()
                                   # Final classification of each component
                                   private$..fitted <- self$classify(x)
                                   # Compute hessian matrix
                                   self$Hessian()
                                 },
                                 #' @description
                                 #' Update only the parameters (means, sd and p) inside the object.
                                 #' @param means Numeric vector of means parameters.
                                 #' @param sd Numeric vector of std. deviation parameters.
                                 #' @param p Numeric vector of probability parameters.
                                 update = function(means, sd, p){
                                   # Update mean parameters
                                   if (!missing(means)) {
                                     private$..means <- unlist(means)
                                   }
                                   # Update Std. deviations parameters
                                   if (!missing(sd)) {
                                     private$..sd <- unlist(sd)
                                   }
                                   # Update probability parameters
                                   if (!missing(p)) {
                                     private$..p <- unlist(p)
                                   }
                                   # Reorder the parameters
                                   private$..means <- private$..means
                                   private$..sd <- private$..sd
                                   private$..p <- private$..p
                                   # Assign a unique name to the parameters
                                   names(private$..means) <- paste0("mu", 1:self$components)
                                   names(private$..sd) <- paste0("sd", 1:self$components)
                                   names(private$..p) <- paste0("p", 1:self$components)
                                 },
                                 #' @description
                                 #' Update the log-likelihood with the current parameters
                                 update_logLik = function(){
                                   # Default x and weights
                                   x <- private$x
                                   w <- private$w
                                   # Set NA weights equal to zero
                                   w[is.na(x)] <- 0
                                   # Update log-likelihood
                                   private$..loglik <- self$logLik(x[w!=0])
                                 },
                                 #' @description
                                 #' Compute the parameters on the classified time series.
                                 #' @details Applied after updating the parameters
                                 update_empiric_parameters = function(){
                                   # Compute empiric parameters
                                   df_emp <- private$..fitted %>%
                                     dplyr::mutate(x = private$x) %>%
                                     dplyr::group_by(classification) %>%
                                     dplyr::arrange(classification) %>%
                                     dplyr::summarise(mu = mean(x, na.rm = TRUE),
                                                      sd = sd(x, na.rm = TRUE),
                                                      p = dplyr::n()/nrow(private$..fitted))
                                   if (length(df_emp$mu) < self$components){
                                     cli::cli_alert_warning("Classification with ML gives only one class!")
                                     # Store the parameters
                                     private$empiric_params$means <- private$..means
                                     private$empiric_params$sd <- private$..sd
                                     private$empiric_params$p <- private$..p
                                   } else {
                                     # Store the parameters
                                     private$empiric_params$means <- df_emp$mu
                                     private$empiric_params$sd <- df_emp$sd
                                     private$empiric_params$p <- df_emp$p
                                   }
                                   # Standard names
                                   names(private$empiric_params$means) <- names(private$..means)
                                   names(private$empiric_params$sd) <- names(private$..sd)
                                   names(private$empiric_params$p) <- names(private$..p)
                                 },
                                 #' @description
                                 #' Update the responsibilities, the log-likelihood, classify again the points and recompute empiric parameters.
                                 #' @details Applied after updating the parameters
                                 filter = function(){
                                   # Update posterior probabilities
                                   private$responsabilities <- self$E_step()
                                   # Update log-likelihood
                                   self$update_logLik()
                                   # Final classification of each component
                                   private$..fitted <- self$classify()
                                   # Update empiric parameters
                                   self$update_empiric_parameters()
                                   # Compute hessian matrix
                                   self$Hessian()
                                 },
                                 #' @description
                                 #' Hessian matrix `gaussianMixture` class.
                                 Hessian = function(){
                                   # Estimated parameters
                                   params <- self$coefficients
                                   # Remove last probability
                                   params$p <- params$p[-c(2)]
                                   params <- unlist(c(params$means, params$sd, params$p))
                                   # Log likelihood
                                   logLik <- function(params){
                                     par <- list()
                                     par$means <- params[stringr::str_detect(names(params), "mu")]
                                     par$sd <- params[stringr::str_detect(names(params), "sd")]
                                     par$p <- params[stringr::str_detect(names(params), "p")]
                                     par$p <- c(par$p, 1-sum(par$p))
                                     -self$logLik(x = private$x[private$w != 0], params = par)
                                   }
                                   # Numeric computation of the hessian matrix
                                   private[["..hessian"]] <- numDeriv::hessian(logLik, x = params)
                                   # Std. errors
                                   std.errors <- sqrt(diag(solve(private[["..hessian"]])))
                                   names(std.errors) <- names(params)
                                   # Store the std. errors
                                   private$..std.means <- std.errors[stringr::str_detect(names(std.errors), "mu")]
                                   private$..std.sd <- std.errors[stringr::str_detect(names(std.errors), "sd")]
                                   private$..std.p <- c(NA_integer_, std.errors[stringr::str_detect(names(std.errors), "p")])
                                   names(private$..std.p) <- names(self$p)
                                 },
                                 #' @description
                                 #' Substitute the empiric parameters with EM parameters. If evaluated again
                                 #' the EM parameters will be substituted back.
                                 use_empiric_parameters = function(){
                                   private$..use_empiric <- !private$..use_empiric
                                 },
                                 #' @description
                                 #' Print method for `gaussianMixture` class.
                                 #' @param label Character, optional label.
                                 print = function(label){
                                   # Format parameters
                                   means <- purrr::map2_chr(self$means, self$std.errors$means, ~paste0(format(.x, digits = 3), " (\033[1;31m", format(.y, digits = 3), "\033[0m)"))
                                   sd <- purrr::map2_chr(self$sd, self$std.errors$sd, ~paste0(format(.x, digits = 3), " (\033[1;31m", format(.y, digits = 3), "\033[0m)"))
                                   p <- purrr::map2_chr(self$p, self$std.errors$p, ~paste0(format(.x, digits = 3), " (\033[1;31m", format(.y, digits = 3), "\033[0m)"))
                                   lbl <- ifelse(missing(label), paste0(self$components, " components"), paste0("\033[1;35m", label, "\033[0m"))
                                   msg_0 <- paste0("#################### ", "gaussianMixture", " (", lbl, ") ", "####################\n")
                                   msg_1 <- paste0("Means: ", paste0(means, collapse = " "), " \n")
                                   msg_2 <- paste0("Stddev: ", paste0(sd, collapse = " "), " \n")
                                   msg_3 <- paste0("Probs: ", paste0(p, collapse = " "), " \n")
                                   msg_4 <- paste0("Nobs: ", length(private$x), "\n")
                                   msg_5 <- paste0("Use empirical moments: ", self$use_empiric, "\n")
                                   msg_6 <- paste0("Log-Likelihood: ", format(self$loglik, digits = 3), "\n")
                                   msg_7 <- paste0("Version: ", private$version, "\n")
                                   line <- paste0(c(rep("-", ifelse(missing(label), 0, -11) + length(strsplit(msg_0, "")[[1]])), "\n"), collapse = "")
                                   cat(paste0(line, msg_0, line, msg_1, msg_2, msg_3, line,
                                              msg_4, msg_5, line, msg_6, msg_7))
                                 }
                               ),
                               private  = list(
                                 version = "1.0.0",
                                 x = NA,
                                 w = NA,
                                 ..loglik = NA,
                                 ..fitted = NA,
                                 ..means = NA,
                                 ..sd = NA,
                                 ..p = NA,
                                 ..hessian = NA,
                                 ..std.means = NA,
                                 ..std.sd = NA,
                                 ..std.p = NA,
                                 ..use_empiric = FALSE,
                                 empiric_params = list(),
                                 responsabilities = NA,
                                 ..moments = list(m1 = NA, m2 = NA, m3 = NA, m4 = NA,
                                                  mean = NA, variance = NA, skewness = NA, kurtosis = NA, nobs = NA)
                               ),
                               active = list(
                                 #' @field means Numeric vector containing the location parameter for each component.
                                 means = function(){
                                   if (self$use_empiric) {
                                     private$empiric_params$means
                                   } else {
                                     private$..means
                                   }
                                 },
                                 #' @field sd Numeric vector containing the scale parameter for each component.
                                 sd = function(){
                                   if (self$use_empiric) {
                                     private$empiric_params$sd
                                   } else {
                                     private$..sd
                                   }
                                 },
                                 #' @field p Numeric vector containing the probability for each component.
                                 p = function(){
                                   if (self$use_empiric) {
                                     private$empiric_params$p
                                   } else {
                                     private$..p
                                   }
                                 },
                                 #' @field coefficients named list with mixture coefficients.
                                 coefficients = function(){
                                   # Output data
                                   list(means = self$means, sd = self$sd, p = self$p)
                                 },
                                 #' @field use_empiric logical to denote if empiric parameters are currently used
                                 use_empiric = function(){
                                   private$..use_empiric
                                 },
                                 #' @field std.errors named list with mixture parameters.
                                 std.errors = function(){
                                   # Output data
                                   list(means = private$..std.means, sd = private$..std.sd, p = private$..std.p)
                                 },
                                 #' @field model Tibble with mixture parameters, in order means, sd, p.
                                 model = function(){
                                   dplyr::bind_cols(as.list(self$means), as.list(self$sd), as.list(self$p))
                                 },
                                 #' @field loglik log-likelihood of the fitted series.
                                 loglik = function(){
                                   private$..loglik
                                 },
                                 #' @field fitted fitted series
                                 fitted = function(){
                                   df <- private$..fitted
                                   limit_uncertanty <- median(df$uncertanty)
                                   dplyr::mutate(df, uncertanty = ifelse(uncertanty < limit_uncertanty, uncertanty, 0.5))
                                 },
                                 #' @field moments Tibble with the theoric moments and the number of observations used for fit.
                                 moments = function(){
                                   # Compute the mixture moments
                                   private$..moments <- GM_moments(self$means, self$sd, self$p)
                                   # Add number of observations
                                   private$..moments$nobs <- length(private$x[private$w!=0])
                                   return(dplyr::bind_cols(private$..moments))
                                 },
                                 #' @field summary Tibble with estimated parameters, std.errors and statistics
                                 summary = function(){
                                   df <- dplyr::bind_rows(
                                     dplyr::bind_cols(term = names(self$means), estimate = unlist(self$means), std.error = unlist(self$std.errors$means)),
                                     dplyr::bind_cols(term = names(self$sd), estimate = unlist(self$sd), std.error = unlist(self$std.errors$sd)),
                                     dplyr::bind_cols(term = names(self$p), estimate = unlist(self$p), std.error = unlist(self$std.errors$p))
                                   )
                                   dplyr::mutate(df, statistic = estimate / std.error, p.value = (1-pnorm(abs(statistic))) / 2)
                                 }
                               ))

#' Moments of a gaussian mixture
#'
#' Compute the first fourth moments and statistics for a Gaussian Mixture with K components.
#'
#' @examples
#' GM_moments(c(-0.3, 0.8), c(0.4,1), c(0.5, 0.5))
#' GM_moments(c(-0.8, 0.8), c(0.4,1), c(0.5, 0.5))
#'
#' @rdname GM_moments
#' @name GM_moments
#' @export
GM_moments <- function(means, sd, alpha){
  # Mixture first moments (Expectation)
  e_x = sum(means * alpha)
  # Central moments
  delta_k <- means - e_x
  # Mixture second moment
  m2 <- sum((delta_k^2 + sd^2) * alpha)
  # Mixture third moment
  m3 <- sum((3 * delta_k * sd^2 + delta_k^3) * alpha)
  # Fourth moment
  m4 <- sum((3 * sd^4 + 6 * delta_k^2 * sd^2 + delta_k^4) * alpha)
  # Fifth moment
  m5 <- sum((delta_k^5 + 10 * delta_k^3 * sd^2 + 15 * delta_k * sd^4) * alpha)
  # Sixth moment
  m6 <- sum((delta_k^6 + 15 * delta_k^4 * sd^2 + 45 * delta_k^2 * sd^4 + 15 * sd^6)*alpha)
  # Mixture variance
  v_x <- m2 - e_x^2
  # Mixture skewness
  # https://stats.stackexchange.com/questions/54733/skewness-of-a-mixture-density
  sk_x = m3 / (m2^(3/2))
  # Mixture excess kurtosis
  # https://en.wikipedia.org/wiki/Kurtosis#Sample_kurtosis
  # https://en.wikipedia.org/wiki/Normal_distribution#Moments
  kt_x = m4 / (m2^2) - 3

  dplyr::tibble(
    m1 = e_x,
    m2 = m2,
    m3 = m3,
    m4 = m4,
    m5 = m5,
    m6 = m6,
    mean = e_x,
    variance = v_x,
    skewness = sk_x,
    kurtosis = kt_x
  )
}

#' Match the first three moments of a Gaussian Mixture
#'
#' @param d Numeric, distance between the two means.
#' @param m1 Numeric, first target moment.
#' @param m2 Numeric, second target moment.
#' @param m3 Numeric, third target moment.
#' @param p Numeric, probability.
#' @export
GM_moments_match <- function(d, m1 = 0, m2 = 1, m3 = 0, p = 0.5){
  means <- c(mu1 = 0, mu2 = 0)
  sigma <- c(sd1 = 0, sd2 = 0)
  probs <- c(p1 = p, p2 = 1-p)
  # First moment
  means[1] <- m1 + (1-p) * d
  means[2] <- m1 - p * d
  # Second moment
  delta <- (m3 - p * (1 - p) * ((1 - p)^2 - p^2) * d^3) / (3 * p * (1-p) * d)
  sigma[2] <- m2 - p * (1 - p) * d^2 - p * delta
  sigma[1] <- sigma[2] + delta
  sigma <- sqrt(sigma)

  structure(
    list(
      means = means,
      sigma = sigma,
      probs = probs
    )
  )
}

#' Compute the log-likelihood of a Gaussian Mixture
#'
#' @param means description
#' @param sd description
#' @param alpha description
#' @examples
#' GM_loglik(c(-0.8, 0.8), c(0.4,1), c(0.5, 0.5), rnorm(100))
#' @export
GM_loglik <- function(means, sd, alpha, x){
  # Log-likelihood
  if (!missing(x)) {
    lik <- dmixnorm(x, means, sd, alpha)
    loss <- sum(log(lik))
  } else {
    loss <- 0
  }
  return(loss)
}


#' Fit GM
#
#' @export
GM_fit <- function(x, method = c("mclust", "mixtools"),components = 2, maxit = 30000){
  clust <- NULL
  method <- match.arg(method, choices = c("mclust", "mixtools"))
  # Fitted parameters with mclust
  if (method == "mclust") {
    clust <- mclust::Mclust(x, G = components, modelNames = c("V"), verbose = FALSE)
    # Estimated parameters
    means <- clust$parameters$mean
    sd <- sqrt(clust$parameters$variance$scale)
    p <- clust$parameters$pro
  }
  # Fitted parameters with mixtools
  if (method == "mixtools" | is.null(clust)) {
    quiet_EM <- purrr::quietly(mixtools::normalmixEM)
    clust <- quiet_EM(x, maxit = maxit, k = components)$result
    # Estimated parameters
    means <- clust$mu
    sd <- clust$sigma
    p <- clust$lambda
  }
  idx <- order(means)
  means = means[idx]
  sd = sd[idx]
  p = p[idx]
  # Assign a name to the parameters
  names(means) <- paste0("mu", 1:components)
  names(sd) <- paste0("sd", 1:components)
  names(p) <- paste0("p", 1:components)

  structure(
    list(
      means = means,
      sd = sd,
      p = p
    )
  )
}
#' Fit GM robust
#
#' @export
GM_fit_rob <- function(x, B = length(x), method = c("mclust", "mixtools"), components = 2, maxit = 30000){
  n <- length(x)
  idx <- sample(min(c(B, n), n))
  means <- list()
  sd <- list()
  probs <- list()
  for(i in 1:length(idx)){
    x_pert <- x[-idx[i]]
    # Fitted parameters
    fit <- GM_fit(x_pert, method, components, maxit)
    means[[i]] <- dplyr::bind_rows(fit$means)
    sd[[i]] <- dplyr::bind_rows(fit$sd)
    probs[[i]] <- dplyr::bind_rows(fit$p)
  }
  means <- unlist(dplyr::summarise_all(dplyr::bind_rows(means), mean))
  sd <- unlist(dplyr::summarise_all(dplyr::bind_rows(sd), mean))
  probs <- unlist(dplyr::summarise_all(dplyr::bind_rows(probs), mean))
  # Assign a name to the parameters
  names(means) <- paste0("mu", 1:length(means))
  names(sd) <- paste0("sd", 1:length(means))
  names(probs) <- paste0("p", 1:length(means))

  structure(
    list(
      means = means,
      sd = sd,
      p = probs
    )
  )
}


#' Multivariate gaussian mixture
#'
#' @rdname mvgaussianMixture
#' @name mvgaussianMixture
#' @export
mvgaussianMixture <- function(x, means, sd, p, components = 2, maxit = 100, abstol = 10e-15, na.rm = FALSE){

  # Ensure that there are not NAs or NaN observations
  idx_NA <- is.na(x)
  if (any(idx_NA)) {
    x <- na.omit(x)
    wrn <- paste0("Removed ", sum(idx_NA), " NA observations!")
    warning(wrn)
  }
  # Number of observations
  n_w <- nrow(x)
  # Number of variables
  j_w <- ncol(x)
  # Empirical moments
  e_x_hat <- colMeans(x, na.rm = na.rm)
  v_x_hat <- apply(x, 2, var)

  # Default starting means
  if (missing(means) || any(is.na(means))){
    # Initialize a matrix for the components
    means <- matrix(0, nrow = components, ncol = j_w, dimnames = list(1:components, dimnames(x)[[2]]))
    probs <- seq(0.8, 0.2, length.out = components)
    for(k in 1:components){
      means[k,] <- apply(x, 2, quantile, probs = probs[k])
    }
  }

  # Default std. deviations
  if (missing(sd) || any(is.na(sd))){
    sd <- list()
    for(k in 1:components){
      sd[[k]] <- diag(v_x_hat)
    }
  }
  # Default probabilities
  if (missing(p) || any(is.na(p))) {
    p <- rep(1/components, components)
  }

  # Routine
  # 0. Initialization
  log_likelihood <- 0
  previous_log_likelihood <- -Inf
  prev_responsibilities <- matrix(0, nrow = n_w, ncol = components)
  previous_params <- list(mean = means, sd = sd, p = p)
  iteration <- 1
  # EM Algorithm
  for (iteration in 1:maxit) {
    # E-step: posterior probabilities
    responsibilities <- prev_responsibilities
    for (i in 1:n_w) {
      for(k in 1:components){
        responsibilities[i, k] <- previous_params$p[k]*mvtnorm::dmvnorm(x[i,], mean = previous_params$mean[k,], sigma = previous_params$sd[[k]])
      }
      # Normalize the posterior probabilities
      responsibilities[i,] <- (responsibilities[i,])/sum(responsibilities[i,], na.rm = TRUE)
      responsibilities[i,][is.na(responsibilities[i,])] <- 0
    }

    # Optimal parameters
    k <- 1
    params <- previous_params
    # M-step: Update the parameters
    for(k in 1:components){
      # Normalizing factor for each component
      n_k <- sum(responsibilities[, k], na.rm = na.rm)
      # Mean parameters k-component
      params$mean[k,] <- apply(responsibilities[, k]*x, 2, sum)/n_k
      # Covariance matrix k-component
      params$sd[[k]] <- diag(apply(x^2*responsibilities[, k], 2, sum)/n_k - params$mean[k,]^2)
      params$sd[[k]][1,2] <- params$sd[[k]][2,1] <- sum(x[,1]*x[,2]*responsibilities[, k])/n_k - params$mean[k,][1]*params$mean[k,][2]
      # Probability k-component
      params$p[k] <- n_k/n_w
    }

    if(any(params$p > 0.9)){
      warning("Probability greater than 0.9 Break!")
      params <- previous_params
      break
    }

    # Calculate the log-likelihood
    log_likelihood <- 0
    for(i in 1:n_w) {
      ll <- 0
      for(k in 1:components){
        ll <- ll + params$p[k]*mvtnorm::dmvnorm(x[i,], mean = params$mean[,k], sigma = params$sd[[k]])
      }
      log_likelihood <- sum(c(log_likelihood, log(ll)), na.rm = TRUE)
    }

    # Check for convergence
    stop_condition <- abs(log_likelihood - previous_log_likelihood) < abstol
    if (stop_condition) {
      break
    } else {
      # Update log-likelihood
      previous_log_likelihood <- log_likelihood
      # Update parameters
      previous_params <- params
    }
    print(log_likelihood)
    if (iteration == maxit) {
      message("Max iteration reached (", iteration, ")")
    }
  }

  # Final classification of each component
  B_hat <- matrix(0, nrow = n_w, ncol = components)
  for(i in 1:n_w) {
    ll <- c()
    for(k in 1:components){
      ll[k] <- sum(responsibilities[i,k])
    }
    B_hat[i, which.max(ll)] <- 1
  }
  colnames(B_hat) <- paste0("B", 1:components)
  B_hat <- dplyr::as_tibble(B_hat)

  # ML-parameters
  params <- previous_params
  # Reorder the components by decreasing means
  colnames(responsibilities) <- paste0("B", 1:components)
  responsibilities <- dplyr::as_tibble(responsibilities)

  # Log-likelihood on fitted parameters
  # Calculate the log-likelihood
  # Calculate the log-likelihood
  log_likelihood <- 0
  for(i in 1:n_w) {
    ll <- 0
    for(k in 1:components){
      ll <- ll + params$p[k]*mvtnorm::dmvnorm(x[i,], mean = params$mean[,k], sigma = params$sd[[k]])
    }
    log_likelihood <- sum(c(log_likelihood, log(ll)), na.rm = TRUE)
  }

  upd_params <- list()
  upd_params$means <- params$mean
  upd_params$sigma2 <- params$mean
  upd_params$rho <- c(0,0)
  upd_params$p <- params$p
  i <- 1
  for(i in 1:length(params$sd)){
    upd_params$sigma2[i,] <- diag(params$sd[[i]])
    upd_params$rho[i] <-  params$sd[[i]][upper.tri(params$sd[[i]])]/prod(sqrt(diag(params$sd[[i]])))
  }
  structure(
    list(
      B_hat = B_hat,
      iteration = iteration,
      params = upd_params,
      responsibilities = responsibilities,
      log_lik = log_likelihood
    ),
    class = c("mvgaussianMixture")
  )
}



