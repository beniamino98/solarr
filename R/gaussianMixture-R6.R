#' Gaussian mixture
#'
#' Fit the parameters of a gaussian mixture with k-components.
#'
#' @examples
#' means = c(0.5,2)
#' sd = rep(1, 2)
#' p = c(0.2,  0.8)
#' # Grid
#' grid <- seq(-4, 4, 0.01)
#' plot(dmixnorm(grid, means, sd, p))
#' # Simulated sample
#' x <- rmixnorm(5000, means, sd, p)
#' # Gaussian mixture model
#' gm <- gaussianMixture$new(components=2)
#' # Fit the model
#' gm$fit(x$X)
#' gm
#' self <- gm$.__enclos_env__$self
#' private <- gm$.__enclos_env__$private
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
#' @keywords gaussianMixture
#' @note Version 1.0.2
#' @export
gaussianMixture <- R6::R6Class("gaussianMixture",
                               public = list(
                                 #' @field control List to contain custom control parameters.
                                 control = list(),
                                 #' @field maxit Integer, maximum number of iterations.
                                 maxit = 5000,
                                 #' @field maxrestarts Integer, maximum number of restarts.
                                 maxrestarts = 500,
                                 #' @field abstol Numeric, absolute level for convergence.
                                 abstol = 1e-08,
                                 #' @field components Integer, number of mixture components.
                                 components = 2,
                                 #' @field method Character, package used for fitting.
                                 method = "mixtools",
                                 #' @description
                                 #' Initialize a `gaussianMixture` object.
                                 #' @param components (`integer(1)`), number of components.
                                 #' @param method Character, package used to fit the parameters. Can be `mclust` or `mixtools`.
                                 #' @param maxit (`integer(1)`) Numeric, maximum number of iterations.
                                 #' @param abstol (`numeric(1)`) Numeric, absolute level for convergence.
                                 #' @param maxrestarts (`integer(1)`) Numeric, maximum number of restarts.
                                 initialize = function(components = 2, method = "mixtools", maxit = 5000, maxrestarts = 500, abstol = 1e-08){
                                   # Maximum iteration
                                   self$maxit <- maxit
                                   # Absolute level for convergence
                                   self$abstol <- abstol
                                   # Maximum number of restarts
                                   self$maxrestarts <- maxrestarts
                                   # Number of components
                                   self$components <- components
                                   # Method for fit
                                   self$method <- match.arg(method, choices = c("mixtools", "mclust"))
                                   # ************************************************
                                   # Iniitialization
                                   ## Means
                                   init_means <- seq(-3, 3, length.out = components)
                                   names(init_means) <- paste0("mu", 1:components)
                                   ## Std. deviations
                                   init_sd <- rep(1, components)
                                   names(init_sd) <- paste0("sd", 1:components)
                                   ## Probabilities
                                   init_p <- rep(1/components, components)
                                   names(init_p) <- paste0("p", 1:components)
                                   # ************************************************
                                   # Initialize coefficients
                                   private[["..means"]] <- init_means
                                   private[["..sd"]] <- init_sd
                                   private[["..p"]] <- init_p/sum(init_p)
                                   # Initialize std. errors
                                   private[["..std.errors.means"]] <- setNames(rep(NA, components), names(init_means))
                                   private[["..std.errors.sd"]] <- setNames(rep(NA, components), names(init_sd))
                                   private[["..std.errors.p"]] <- setNames(rep(NA, components), names(init_p))
                                 },
                                 #' @description
                                 #' Set the time series internally
                                 #' @param x vector
                                 #' @param weights observations weights, if a weight is equal to zero the observation is excluded, otherwise is included with unitary weight.
                                 set_time_series = function(x, weights){
                                   # Add custom weights
                                   if (missing(weights)) {
                                     w <- rep(1, length(x))
                                   } else {
                                     w <- ifelse(weights == 0, 0, 1)
                                   }
                                   # Ensure no NA
                                   idx_NA <- which(is.na(x))
                                   if (!purrr::is_empty(idx_NA)) {
                                     cli::cli_alert_warning(paste0("Some NA found at indexes: ", idx_NA, collapse = " "))
                                     w[is.na(x)] <- 0
                                   }
                                   # Ensure no infinity
                                   idx_infty <- which(is.infinite(x))
                                   if (!purrr::is_empty(idx_infty)) {
                                     cli::cli_alert_warning(paste0("Some infinity values found at indexes: ", idx_infty, collapse = " "))
                                     w[is.na(x)] <- 0
                                   }
                                   # ************************************************
                                   # Add x
                                   private[["x"]] <- x
                                   # Add weights
                                   private[["w"]] <- w
                                 },
                                 #' @description
                                 #' Compute the log-likelihood.
                                 #' @param x vector
                                 #' @param params Optional. Named list with mixture parameters.
                                 #' @param per_obs Logical, when `TRUE` the log-likelihood is returned per observation,
                                 #' otherwise is summed.
                                 logLik = function(x, params, per_obs = FALSE){
                                   if (missing(params)) {
                                     params <- self$coefficients
                                   }
                                   # Calculate the likelihoods
                                   likelihoods <- matrix(NA, nrow = length(x), ncol = self$components)
                                   for(k in 1:self$components) {
                                     likelihoods[,k] <- params$p[k] * dnorm(x, params$means[k], params$sd[k])
                                   }
                                   # Calculate the total log-likelihood
                                   log_likelihood <- log(rowSums(likelihoods))
                                   if (!per_obs){
                                     log_likelihood <- sum(log_likelihood, na.rm = TRUE)
                                   }
                                   return(log_likelihood)
                                 },
                                 #' @description
                                 #' Compute the posterior probabilities (E-step),
                                 #' @param x Time series to fit.
                                 #' @param params A named list with mixture parameters.
                                 E_step = function(x, params){
                                   if (missing(x)) {
                                     x <- private$x
                                   }
                                   if (missing(params)) {
                                     params <- self$coefficients
                                   }
                                   responsabilities <- matrix(0, nrow = length(x), ncol = self$components)
                                   for(k in 1:self$components) {
                                     responsabilities[, k] <- params$p[k] * dnorm(x, params$means[k], params$sd[k])
                                   }
                                   # Normalize the posterior probabilities
                                   responsabilities <- apply(responsabilities, 2, function(x) ifelse(is.na(x)|is.nan(x), 0, x))
                                   responsabilities <- responsabilities / rowSums(responsabilities)
                                   colnames(responsabilities) <- paste0("B", 1:self$components)
                                   return(dplyr::as_tibble(responsabilities))
                                 },
                                 #' @description
                                 #' Classify the time series in the components with highest likelihood.
                                 #' @param x Time series to fit.
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
                                 #' @param mu_target Numeric, target mean of the mixture to match.
                                 #' @param var_target Numeric, target variance of the mixture to match.
                                 #' When `missing` all the available observations will be used.
                                 fit = function(x, weights, mu_target = NA, var_target = NA){
                                   # Number of observations
                                   n <- length(x)
                                   self$set_time_series(x, weights)
                                   # Add time series
                                   x <- private$x
                                   w <- private$w
                                   # Fitted parameters
                                   clust <- GM_fit(x[w!=0], method = self$method, components = self$components,
                                                   maxit = self$maxit, maxrestarts = self$maxrestarts)
                                   # Moments matching
                                   if (self$components == 2) {
                                     if (!is.na(mu_target) | !is.na(var_target)) {
                                       start_params <- unlist(purrr::flatten(clust)[-6])
                                       clust <- GM_fit_moments_match(x[w!=0], start_params, mu_target = mu_target, var_target = var_target,
                                                                     maxit = self$maxit, abstol = self$abstol)
                                     }
                                   }
                                   # Update the parameters
                                   private$..means <- clust$means
                                   private$..sd <- clust$sd
                                   private$..p <- clust$p
                                   # Assign a name to the parameters
                                   names(private$..means) <- paste0("mu", 1:self$components)
                                   names(private$..sd) <- paste0("sd", 1:self$components)
                                   names(private$..p) <- paste0("p", 1:self$components)
                                   # Filter
                                   self$filter()
                                 },
                                 #' @description
                                 #' Fit the parameters (EM-algorithm)
                                 #' @param x vector
                                 #' @param weights observations weights, if a weight is equal to zero the observation is excluded, otherwise is included with unitary weight.
                                 #' When `missing` all the available observations will be used.
                                 EM = function(x, weights){
                                   if (missing(x)) {
                                     x <- private$x
                                   }
                                   if (missing(weights)) {
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
                                   for (iteration in 1:self$maxit) {
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
                                   # Filter
                                   self$filter()
                                 },
                                 #' @description
                                 #' Update the parameters inside the object.
                                 #' @param means Named vector, new means parameters.
                                 #' @param sd Named vector, new std. deviation parameters.
                                 #' @param p Named vector, new probabilities.
                                 update = function(means, sd, p){
                                   # Update means
                                   if (!missing(means)) {
                                     # Extract old coefficients
                                     new_coefs <- self$means
                                     # Extract names
                                     names_old <- names(new_coefs)
                                     names_new <- names(means)
                                     # Update only if they are present
                                     for(i in 1:length(means)){
                                       condition <- names_new[i] %in% names_old
                                       if (condition) {
                                         old_coef <- private[["..means"]][names_new[i]]
                                         if (old_coef != means[i]){
                                           private[["..std.errors.means"]][names_new[i]] <- NA_integer_
                                           private[["..means"]][names_new[i]] <- means[i]
                                         }
                                       }
                                     }
                                   }
                                   # ************************************************
                                   # Update std. deviations
                                   if (!missing(sd)) {
                                     # Extract old coefficients
                                     new_coefs <- self$sd
                                     # Extract names
                                     names_old <- names(new_coefs)
                                     names_new <- names(sd)
                                     # Update only if they are present
                                     for(i in 1:length(sd)){
                                       condition <- names_new[i] %in% names_old
                                       if (condition) {
                                         old_coef <- private[["..sd"]][names_new[i]]
                                         if (old_coef != sd[i]){
                                           private[["..std.errors.sd"]][names_new[i]] <- NA_integer_
                                           private[["..sd"]][names_new[i]] <- sd[i]
                                         }
                                       }
                                     }
                                   }
                                   # ************************************************
                                   # Update probabilities
                                   if (!missing(p)) {
                                     # Extract old coefficients
                                     new_coefs <- self$p
                                     # Extract names
                                     names_old <- names(new_coefs)
                                     names_new <- names(p)
                                     # Update only if they are present
                                     for(i in 1:length(p)){
                                       condition <- names_new[i] %in% names_old
                                       old_coef <- private[["..p"]][names_new[i]]
                                       if (old_coef != p[i]){
                                         private[["..std.errors.p"]][names_new[i]] <- NA_integer_
                                         private[["..p"]][names_new[i]] <- p[i]
                                       }
                                     }
                                   }
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
                                   private[["..loglik"]] <- self$logLik(x[w!=0])
                                 },
                                 #' @description
                                 #' Update the responsibilities, the log-likelihood, classify again the points and recompute empiric parameters.
                                 #' @details Applied after updating the parameters
                                 filter = function(){
                                   # Update posterior probabilities
                                   private[["..responsabilities"]] <- self$E_step()
                                   # Update log-likelihood
                                   self$update_logLik()
                                   # Final classification of each component
                                   private[["..fitted"]] <- self$classify()
                                   # Compute hessian matrix
                                   self$Hessian()
                                 },
                                 #' @description
                                 #' Hessian matrix `gaussianMixture` class.
                                 Hessian = function(){
                                   # Estimated ML parameters
                                   params <- self$coefficients
                                   # Remove last probability
                                   params$p <- params$p[-c(2)]
                                   # Vectorize
                                   params <- unlist(purrr::flatten(params))
                                   # Time series to use for the computation
                                   y <- private$x[private$w != 0]
                                   # Log likelihood
                                   logLik <- function(params, y, per_obs = FALSE){
                                     par <- list()
                                     par$means <- params[stringr::str_detect(names(params), "mu")]
                                     par$sd <- params[stringr::str_detect(names(params), "sd")]
                                     par$p <- params[stringr::str_detect(names(params), "p")]
                                     par$p <- c(par$p, 1-sum(par$p))
                                     self$logLik(x = y, params = par, per_obs = per_obs)
                                   }
                                   # ************************************************
                                   coefs_names <- names(params)
                                   # Numerical Hessian matrix at MLE
                                   H <- numDeriv::hessian(logLik, x = params, y = y, per_obs = FALSE)
                                   colnames(H) <- rownames(H) <- coefs_names
                                   # Numerical Score at MLE
                                   S <- numDeriv::jacobian(logLik, x = params, y = y, per_obs = TRUE)
                                   colnames(S) <- coefs_names
                                   # Cross products of the score
                                   B <- S[1,] %*% t(S[1,])
                                   for(i in 2:nrow(S)){
                                     B <- B + S[i,] %*% t(S[i,])
                                   }
                                   colnames(B) <- rownames(B) <- coefs_names
                                   # Variance covariance matrix
                                   V <-  solve(-H)
                                   # Sandwitch var-cov matrix
                                   V_rob <- V %*% B %*% V
                                   # Jacobian for last probability
                                   d <- nrow(V_rob)
                                   J <- diag(1, d + 1, d)
                                   J[d + 1, d] <- -1
                                   colnames(J) <- coefs_names
                                   rownames(J) <- names(unlist(c(self$coefficients[[1]], self$coefficients[[2]],self$coefficients[[3]])))
                                   # Sandwitch var-cov matrix (original)
                                   V_rob <- J %*% V_rob %*% t(J)
                                   # Store var-cov matrix
                                   private[["..vcov"]] <- V_rob
                                   # Standard errors
                                   std.errors <-  sqrt(diag(V_rob))
                                   # Store the std. errors
                                   private$..std.errors.means <- std.errors[stringr::str_detect(names(std.errors), "mu")]
                                   private$..std.errors.sd <- std.errors[stringr::str_detect(names(std.errors), "sd")]
                                   private$..std.errors.p <- setNames(std.errors[stringr::str_detect(names(std.errors), "p")], names(self$p))
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
                                   msg_5 <- paste0("Log-Likelihood: ", format(self$loglik, digits = 3), "\n")
                                   msg_6 <- paste0("Version: ", private$version, "\n")
                                   line <- paste0(c(rep("-", ifelse(missing(label), 0, -11) + length(strsplit(msg_0, "")[[1]])), "\n"), collapse = "")
                                   cat(paste0(line, msg_0, line, msg_1, msg_2, msg_3, line,
                                              msg_4, msg_5, line, msg_6))
                                 }
                               ),
                               private  = list(
                                 version = "1.0.2",
                                 x = NA,
                                 w = NA,
                                 ..loglik = NA,
                                 ..fitted = NA,
                                 ..means = NA,
                                 ..sd = NA,
                                 ..p = NA,
                                 ..vcov = NA,
                                 ..std.errors.means = NA,
                                 ..std.errors.sd = NA,
                                 ..std.errors.p = NA,
                                 ..responsabilities = NA
                               ),
                               active = list(
                                 #' @field means Numeric vector containing the location parameter for each component.
                                 means = function(){
                                   private$..means
                                 },
                                 #' @field sd Numeric vector containing the scale parameter for each component.
                                 sd = function(){
                                   private$..sd
                                 },
                                 #' @field p Numeric vector containing the probability for each component.
                                 p = function(){
                                   private$..p
                                 },
                                 #' @field coefficients named list with mixture coefficients.
                                 coefficients = function(){
                                   # Output data
                                   list(means = self$means,
                                        sd = self$sd,
                                        p = self$p)
                                 },
                                 #' @field std.errors named list with mixture parameters.
                                 std.errors = function(){
                                   # Output data
                                   list(means = private$..std.errors.means,
                                        sd = private$..std.errors.sd,
                                        p = private$..std.errors.p)
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
                                   if (!is.na(df)) {
                                     limit_uncertanty <- median(df$uncertanty)
                                     df <- dplyr::mutate(df, uncertanty = ifelse(uncertanty < limit_uncertanty, uncertanty, 0.5))
                                   }
                                   return(df)
                                 },
                                 #' @field responsabilities resp
                                 responsabilities = function(){
                                   private$..responsabilities
                                 },
                                 #' @field moments Tibble with the theoric moments and the number of observations used for fit.
                                 moments = function(){
                                   # Compute the mixture moments
                                   moments <- GM_moments(self$means, self$sd, self$p)
                                   # Add number of observations
                                   moments$nobs <- length(private$x[private$w!=0])
                                   return(moments)
                                 },
                                 #' @field tidy Tibble with estimated parameters, std.errors and statistics
                                 tidy = function(){
                                   df <- dplyr::bind_rows(
                                     dplyr::bind_cols(term = names(self$means), estimate = unlist(self$means), std.error = unlist(self$std.errors$means)),
                                     dplyr::bind_cols(term = names(self$sd), estimate = unlist(self$sd), std.error = unlist(self$std.errors$sd)),
                                     dplyr::bind_cols(term = names(self$p), estimate = unlist(self$p), std.error = unlist(self$std.errors$p))
                                   )
                                   dplyr::mutate(df, statistic = estimate / std.error, p.value = (1-pnorm(abs(statistic))) / 2)
                                 }
                               ))
