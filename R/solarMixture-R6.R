#' Monthly Gaussian mixture with two components
#'
#' @examples
#' # Model fit
#' model <- solarModel$new(spec)
#' model$fit()
#' # Inputs
#' x <- model$data$u_tilde
#' w <- model$data$weights
#' date <- model$data$date
#' # Solar Mixture object
#' sm <- solarMixture$new()
#' sm$fit(x, date, w)
#' params <- sm$parameters
#' sm$std.errors
#' # params[1,]$mu1 <- params[1,]$mu1*0.9
#' # sm$update(means = params[,c(2,3)])
#' @rdname solarMixture
#' @name solarMixture
#' @keywords gaussianMixture
#' @note Version 1.0.0
#' @export
solarMixture <-  R6::R6Class("solarMixture",
                             public = list(
                               #' @field control List to contain custom control parameters.
                               control = list(match.expectation = FALSE, match.variance = FALSE, match.empiric = FALSE),
                               #' @field mu1 Function, see \code{\link{monthlyParams}}.
                               mu1 = NA,
                               #' @field mu2 Function, see \code{\link{monthlyParams}}.
                               mu2 = NA,
                               #' @field sd1 Function, see \code{\link{monthlyParams}}.
                               sd1 = NA,
                               #' @field sd2 Function, see \code{\link{monthlyParams}}.
                               sd2 = NA,
                               #' @field prob Function, see \code{\link{monthlyParams}}.
                               prob = NA,
                               #' @description
                               #' Initialize a `solarMixture` object
                               #' @param components (`integer(1)`), number of components.
                               #' @param method Character, package used to fit the parameters. Can be `mixtools` or `mclust`.
                               #' @param maxit (`integer(1)`) Numeric, maximum number of iterations.
                               #' @param abstol (`numeric(1)`) Numeric, absolute level for convergence.
                               #' @param maxrestarts (`integer(1)`) Numeric, maximum number of restarts.
                               initialize = function(components = 2, method = "mixtools", maxit = 10000, maxrestarts = 1000, abstol = 1e-8){
                                 # Initialize Gaussian mixture objects for all the months
                                 private$..model <- vector("list", 12)
                                 private$date_month <- vector("list", 12)
                                 for(m in 1:12){
                                   private$..model[[m]] <- gaussianMixture$new(components = components,
                                                                               method = method,
                                                                               maxit = maxit,
                                                                               maxrestarts = maxrestarts,
                                                                               abstol = abstol)
                                 }
                                 # Assign standard names
                                 names(private$..model) <- lubridate::month(1:12, label = TRUE)
                                 names(private$date_month) <- lubridate::month(1:12, label = TRUE)
                               },
                               #' @description
                               #' Set the time series internally
                               #' @param x vector
                               #' @param dates date vector
                               #' @param weights observations weights, if a weight is equal to zero the observation is excluded, otherwise is included with unitary weight.
                               set_time_series = function(x, dates, weights){
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
                                 # Add date
                                 private[["date"]] <- dates
                                 # ************************************************
                                 data <- self$data
                                 # Gaussian Mixture parameters
                                 for(m in 1:12){
                                   data_months <- dplyr::filter(data, Month == m)
                                   if (nrow(data_months) <= 1){
                                     cli::cli_alert_warning(paste0("No observations for the month: ", m))
                                     next
                                   }
                                   # Store monthly dates
                                   private$date_month[[m]] <- data_months$date
                                   # Fit the GM model for the m-month
                                   private$..model[[m]]$set_time_series(data_months$x, data_months$w)
                                 }
                               },
                               #' @description
                               #' Fit the parameters with mclust package
                               #' @param x vector
                               #' @param date date vector
                               #' @param weights observations weights, if a weight is equal to zero the observation is excluded, otherwise is included with unitary weight.
                               #' When `missing` all the available observations will be used.
                               #' @param method Character, package used to fit the parameters. Can be `mclust` or `mixtools`.
                               #' @param mu_target Numeric vector with length 12, target mean of the mixture to match.
                               #' @param var_target Numeric vector with length 12, target variance of the mixture to match.
                               fit = function(x, date, weights, method = "mixtools", mu_target = rep(NA, 12), var_target = rep(NA, 12)){
                                 # Set time series internally
                                 self$set_time_series(x, date, weights)
                                 data <- self$data
                                 M <- length(private$..model)
                                 # Gaussian Mixture parameters
                                 for(m in 1:M){
                                   data_months <- dplyr::filter(data, Month == m)
                                   if (nrow(data_months) <= 1){
                                     cli::cli_alert_warning(paste0("No observations for the month: ", m))
                                     next
                                   }
                                   # Fit the GM model for the m-month
                                   private$..model[[m]]$fit(data_months$x, data_months$w, mu_target = mu_target[m], var_target = var_target[m])
                                 }
                                 # Initialize monthly function for Mixture parameters
                                 self$mu1 <- monthlyParams$new(self$coefficients$mu1)
                                 self$mu2 <- monthlyParams$new(self$coefficients$mu2)
                                 self$sd1 <- monthlyParams$new(self$coefficients$sd1)
                                 self$sd2 <- monthlyParams$new(self$coefficients$sd2)
                                 self$prob <- monthlyParams$new(self$coefficients$p1)
                               },
                               #' @description
                               #' Update means, sd, p .
                               #' @param means Numeric matrix of means parameters.
                               #' @param sd Numeric matrix of std. deviation parameters.
                               #' @param p Numeric matrix of probability parameters.
                               update = function(means, sd, p){
                                 # Mean parameters
                                 if (missing(means)) {
                                   means <- self$means
                                 }
                                 # Std. deviations parameters
                                 if (missing(sd)) {
                                   sd <- self$sd
                                 }
                                 # Probability parameters
                                 if (missing(p)) {
                                   p <- self$p
                                 }
                                 # Update parameters
                                 for(m in 1:12){
                                   private$..model[[m]]$update(means = means[m,], sd = sd[m,], p = p[m,])
                                 }
                                 # Update monthly function for Mixture parameters
                                 coefficients <- self$coefficients
                                 self$mu1 <- monthlyParams$new(coefficients$mu1)
                                 self$mu2 <- monthlyParams$new(coefficients$mu2)
                                 self$sd1 <- monthlyParams$new(coefficients$sd1)
                                 self$sd2 <- monthlyParams$new(coefficients$sd2)
                                 self$prob <- monthlyParams$new(coefficients$p1)
                               },
                               #' @description
                               #' Apply the `$update_logLik()` method to all the `gaussianMixture` models
                               update_logLik = function(){
                                 M <- length(private$..model)
                                 # Update log-likelihood
                                 for(m in 1:M){
                                   private$..model[[m]]$update_logLik()
                                 }
                                 return(invisible(NULL))
                               },
                               #' @description
                               #' Apply the `$filter()` method to all the `gaussianMixture` models
                               filter = function(){
                                 M <- length(private$..model)
                                 # Update parameters, fitted data and log-likelihood
                                 for(m in 1:M){
                                   private$..model[[m]]$filter()
                                 }
                                 return(invisible(NULL))
                               },
                               #' @description
                               #' Apply the `$Hessian()` method to all the `gaussianMixture` models
                               Hessian = function(){
                                 M <- length(private$..model)
                                 # Update the Hessian and std. errors
                                 for(m in 1:M){
                                   private$..model[[m]]$Hessian()
                                 }
                                 return(invisible(NULL))
                               },
                               #' @description
                               #' Log-likelihood
                               #' @param x vector
                               #' @param date dates
                               logLik = function(x, date){
                                 M <- length(private$..model)
                                 df <- dplyr::tibble(date = date, Month = lubridate::month(date), x = x)
                                 log_lik <- 0
                                 for(m in 1:M){
                                   data_month <- dplyr::filter(df, Month == m)
                                   log_lik <- log_lik + self$model[[m]]$logLik(filter(df, Month == m)$x)
                                 }
                                 return(log_lik)
                               },
                               #' @description
                               #' Compute the grades
                               #' @param x vector
                               #' @param date dates
                               grades = function(x, date){
                                 M <- length(private$..model)
                                 df <- dplyr::tibble(date = date, Month = lubridate::month(date), x = x, u = NA)
                                 for(m in 1:M){
                                   # Index month
                                   idx_month <- which(df$Month == m)
                                   if (purrr::is_empty(idx_month)) next
                                   # Extract the parameters
                                   means <- self$model[[m]]$means
                                   sd <- self$model[[m]]$sd
                                   p <- self$model[[m]]$p
                                   # Monthly cdf
                                   cdf_M <- function(x) pmixnorm(x, means, sd, p)
                                   # Grades
                                   df$u[idx_month] <- cdf_M(x[idx_month])
                                 }
                                 return(df$u)
                               },
                               #' @description
                               #' Compute the VaR with certain confidence levels
                               #' @param x vector
                               #' @param date dates
                               #' @param alpha confidence levels for the VaR
                               VaR = function(date, alpha = 0.05){
                                 M <- length(private$..model)
                                 df <- dplyr::tibble(date = date, Month = lubridate::month(date))
                                 # Compute VaR
                                 VaR_alpha <- matrix(0, nrow = nrow(df), ncol = length(alpha))
                                 # Assign standard names
                                 colnames(VaR_alpha) <- paste0("VaR_", alpha)
                                 alpha <- alpha[order(alpha)]
                                 for(m in 1:M){
                                   # Index month
                                   idx_month <- which(df$Month == m)
                                   if (purrr::is_empty(idx_month)) next
                                   # Extract the parameters
                                   means <- self$model[[m]]$means
                                   sd <- self$model[[m]]$sd
                                   p <- self$model[[m]]$p
                                   # Monthly quantile
                                   Q_M <- function(x) qmixnorm(x, means, sd, p)
                                   # VaR
                                   VaR_alpha_M <- Q_M(alpha)
                                   for(j in 1:length(idx_month)){
                                     VaR_alpha[idx_month[j],] <- VaR_alpha_M
                                   }
                                 }
                                 return(VaR_alpha)
                               },
                               #' @description
                               #' Compute the VaR with certain confidence levels
                               #' @param x vector
                               #' @param date dates
                               #' @param alpha confidence levels for the VaR
                               ES = function(date, alpha = 0.05){
                                 M <- length(private$..model)
                                 df <- dplyr::tibble(date = date, Month = lubridate::month(date))
                                 # Compute VaR
                                 ES_alpha <- matrix(0, nrow = nrow(df), ncol = length(alpha))
                                 # Assign standard names
                                 colnames(ES_alpha) <- paste0("ES_", alpha)
                                 for(m in 1:M){
                                   # Index month
                                   idx_month <- which(df$Month == m)
                                   if (purrr::is_empty(idx_month)) next
                                   # Extract the parameters
                                   means <- self$model[[m]]$means
                                   sd <- self$model[[m]]$sd
                                   p <- self$model[[m]]$p
                                   # Monthly quantile
                                   Q_M <- function(x) qmixnorm(x, means, sd, p)
                                   VaR_alpha <- Q_M(alpha)
                                   # Expected shortfalls
                                   # Density function
                                   f_R <- function(x) dmixnorm(x, means, sd, p)
                                   # Expected Shortfall
                                   ES <- function(VaR, alpha) integrate(function(x) x*f_R(x), lower = -Inf,
                                                                        upper = VaR, subdivisions = 50L, stop.on.error = FALSE)$value/alpha
                                   ES_alpha_M <- purrr::map2_dbl(VaR_alpha, alpha, ~ES(.x, .y))
                                   for(j in 1:length(idx_month)){
                                     ES_alpha[idx_month[j],] <- ES_alpha_M
                                   }
                                 }
                                 return(ES_alpha)
                               },
                               #' @description
                               #' Print method for `solarMixture` class.
                               print = function(){
                                 M <- length(private$..model)
                                 for(m in 1:M){
                                   nmonth <- as.character(lubridate::month(m, label = TRUE, abbr = FALSE))
                                   self$model[[m]]$print(nmonth)
                                 }
                               }
                             ),
                             private = list(
                               version = "1.0.1",
                               x = NA,
                               w = NA,
                               date = NA,
                               date_month = list(),
                               ..model = list(),
                               deep_clone = function(name, value){
                                 if (name == "..model") {
                                   # Clonazione profonda degli oggetti kernelRegression all'interno della lista ..models
                                   return(lapply(value, function(model) {
                                     if (!is.null(model)) {
                                       model$clone(deep = TRUE)  # Clonazione profonda per ogni oggetto kernelRegression
                                     } else {
                                       NULL
                                     }
                                   }))
                                 } else {
                                   # Per altri campi, usa il comportamento di clonazione predefinito
                                   return(value)
                                 }
                               }
                             ),
                             active = list(
                               #' @field data A tibble with the following columns:
                               #' \describe{
                               #'  \item{date}{Time series of dates.}
                               #'  \item{Month}{Vector of Month.}
                               #'  \item{x}{Time series used for fitting.}
                               #'  \item{w}{Time series of weights.}}
                               data = function(){
                                 dplyr::tibble(date = private$date, Month = lubridate::month(date), x = private$x, w = private$w)
                               },
                               #' @field means Matrix of means where a row represents a month and a column a mixture component.
                               means = function(){
                                 M <- length(private$..model)
                                 components <- private$..model[[1]]$components
                                 means <- matrix(0, nrow = M, ncol = components)
                                 for(m in 1:M){
                                   means[m,] <- unlist(self$model[[m]]$means)
                                 }
                                 colnames(means) <- names(self$model[[1]]$means)
                                 rownames(means) <- lubridate::month(1:M, label = TRUE)
                                 return(means)
                               },
                               #' @field sd Matrix of std. deviations where a row represents a month and a column a mixture component.
                               sd = function(){
                                 M <- length(private$..model)
                                 components <- private$..model[[1]]$components
                                 sd <- matrix(0, nrow = M, ncol = components)
                                 for(m in 1:M){
                                   sd[m,] <- unlist(self$model[[m]]$sd)
                                 }
                                 colnames(sd) <- names(self$model[[1]]$sd)
                                 rownames(sd) <- lubridate::month(1:M, label = TRUE)
                                 return(sd)
                               },
                               #' @field p Matrix of probabilities where a row represents a month and a column a mixture component.
                               p = function(){
                                 M <- length(private$..model)
                                 components <- private$..model[[1]]$components
                                 p <- matrix(0, nrow = M, ncol = components)
                                 for(m in 1:M){
                                   p[m,] <- unlist(self$model[[m]]$p)
                                 }
                                 colnames(p) <- names(self$model[[1]]$p)
                                 rownames(p) <- lubridate::month(1:M, label = TRUE)
                                 return(p)
                               },
                               #' @field model Named List with 12 \code{\link{gaussianMixture}} objects.
                               model = function(){
                                 private$..model
                               },
                               #' @field loglik Numeric, total log-likelihood.
                               loglik = function(){
                                 sum(purrr::map_dbl(self$model, ~.x$loglik))
                               },
                               #' @field fitted A `tibble` with the classified series
                               fitted = function(){
                                 df_fitted <- self$model[[1]]$fitted
                                 if (!is.na(df_fitted)) {
                                   df_fitted <- purrr::map2_df(self$model, private$date_month, ~ dplyr::bind_cols(date = .y, dplyr::select(.x$fitted, B = "B1", uncertanty)))
                                   df_fitted <- dplyr::arrange(df_fitted, date)
                                 }
                                 return(df_fitted)
                               },
                               #' @field moments A `tibble` with the theoric moments. It contains:
                               #' \describe{
                               #'  \item{Month}{Month of the year.}
                               #'  \item{mean}{Theoric monthly expected value of the mixture model.}
                               #'  \item{variance}{Theoric monthly variance of the mixture model.}
                               #'  \item{skewness}{Theoric monthly skewness.}
                               #'  \item{kurtosis}{Theoric monthly kurtosis.}
                               #'  \item{nobs}{Number of observations used for fitting.}
                               #'  \item{loglik}{Monthly log-likelihood.}}
                               moments = function(){
                                 M <- length(private$..model)
                                 dplyr::bind_cols(Month = 1:M, purrr::map_df(self$model, ~.x$moments), loglik = purrr::map_dbl(self$model, ~.x$loglik))
                               },
                               #' @field coefficients A `tibble` with the fitted parameters.
                               coefficients = function(){
                                 M <- length(private$..model)
                                 dplyr::bind_cols(Month = 1:M, purrr::map_df(self$model, ~.x$model))
                               },
                               #' @field std.errors A `tibble` with the fitted std.errors
                               std.errors = function(){
                                 M <- length(private$..model)
                                 dplyr::bind_cols(Month = 1:M, purrr::map_df(self$model, ~dplyr::bind_cols(purrr::map(.x$std.errors, ~dplyr::bind_rows(as.list(.x))))))
                               },
                               #' @field tidy A `tibble` with the fitted std.errors
                               tidy = function(){
                                 M <- length(private$..model)
                                 purrr::map2_df(1:M, self$model, ~dplyr::bind_cols(Month = .x, .y$tidy))
                               }
                             )
)
