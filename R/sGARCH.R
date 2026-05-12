#' Implementation of `rugarch` methods for a GARCH(p,q) as R6 class
#'
#' @rdname sGARCH
#' @name sGARCH
#' @keywords GARCH
#' @note Version 1.0.3
#' @export
sGARCH <- R6::R6Class("sGARCH",
                             public = list(
                               #' @field control List to contain custom control parameters.
                               control = list(mode = "unitOmega"),
                               #' @field model_name Character, model name
                               model_name = "sGARCH(0, 0)",
                               #' @description
                               #' Initialize a standard GARCH model
                               #' @param archOrder Integer scalar, ARCH order.
                               #' @param garchOrder Integer scalar, GARCH order.
                               #' @param mode Character, one of \code{"unitOmega"}, \code{"targetSigma2"}, \code{"freeOmega"}.
                               initialize = function(archOrder = 1, garchOrder = 1, mode = "unitOmega"){
                                 # ARCH order
                                 private[["..archOrder"]] <- archOrder
                                 # GARCH order
                                 private[["..garchOrder"]] <- garchOrder
                                 # Unconditional variance constraint
                                 self$control[["mode"]] <- match.arg(mode, choices =  c("unitOmega","targetSigma2","freeOmega"))
                                 # Matrix components
                                 private[["..b"]] <- GARCH_vector_b(archOrder, garchOrder)
                                 # Model name
                                 self$model_name <- paste0("sGARCH", " (", archOrder, ", ", garchOrder, ")")
                                 # Matrix components
                                 # Companion matrix
                                 private[["..A"]] <- GARCH_companion_matrix(alpha = rep(0, archOrder), beta = rep(0, garchOrder))
                                 # Intercept vector
                                 d <- rep(0, length(private[["..b"]]))
                                 d[1] <- self$omega
                                 d <- matrix(d, ncol = 1)
                                 private[["..d"]] <- matrix(d, ncol = 1)
                               },
                               #' @description
                               #' Fit the GARCH model with `rugarch` function.
                               #' @param x Numeric, vector. Time series to be fitted.
                               #' @param weights Numeric, vector. Optional custom weights.
                               fit = function(x, weights){
                                 # Flexible weights (useful if model do not converge)
                                 if (missing(weights)) {
                                   w <- rep(1, length(x))
                                 } else {
                                   w <- ifelse(weights == 0, 0, 1)
                                 }
                                 # GARCH fit
                                 model <- sGARCH_robust_fit(y = x, weights = w, archOrder = self$archOrder, garchOrder = self$garchOrder,
                                                            mode = self$control[["mode"]], method = "rugarch")
                                 # Update coefficients
                                 private[["..omega"]] <- model$omega
                                 private[["..alpha"]] <- model$alpha
                                 private[["..beta"]]  <- model$beta
                                 # Log-likelihods
                                 private[["..log.likelihoods"]] <- model$loglik
                                 # Std. errors
                                 private[["..std.errors"]] <- model$std.errors
                                 # Matrix components
                                 # Companion matrix
                                 private[["..A"]] <- GARCH_companion_matrix(alpha = self$alpha, beta = self$beta)
                                 # Intercept vector
                                 d <- rep(0, length(private[["..b"]]))
                                 d[1] <- self$omega
                                 d <- matrix(d, ncol = 1)
                                 private[["..d"]] <- matrix(d, ncol = 1)
                               },
                               #' @description
                               #' Filter method from `rugarch` package to compute GARCH variance, residuals and log-likelihoods.
                               #' @param x Numeric, vector. Time series to be filtered.
                               #' @param eps0 Optional numeric initial epsilons to prepend (length p+q).
                               #' @param sigma20 Optional numeric initial variances to prepend (length p+q).
                               filter = function(x, eps0 = NULL, sigma20 = NULL){
                                 sGARCH_filter(x, self$omega, self$alpha, self$beta, eps0, sigma20)
                               },
                               #' @description
                               #' Update the coefficients of the model
                               #' @param coefficients Numeric, named vector. Model's coefficients.
                               update = function(coefficients){
                                 if (missing(coefficients)) {
                                   return(invisible(NULL))
                                 }
                                 # Extract old coefficients
                                 new_coefs <- self$coefficients
                                 # Extract names
                                 names_old <- names(new_coefs)
                                 names_new <- names(coefficients)
                                 # Update only if they are present
                                 for(i in 1:length(coefficients)){
                                   if (names_new[i] %in% names_old) {
                                     old_coef <- new_coefs[names_new[i]]
                                     if (old_coef != coefficients[i]){
                                       new_coefs[names_new[i]] <- coefficients[i]
                                       private[["..std.errors"]][names_new[i]] <- NA_integer_
                                     }
                                   }
                                 }
                                 # Update ARCH parameters
                                 if (self$archOrder > 0) {
                                   private[["..alpha"]] <- new_coefs[stringr::str_detect(names_old, "alpha")]
                                 }
                                 # Update GARCH parameters
                                 if (self$garchOrder > 0) {
                                   private[["..beta"]] <- new_coefs[stringr::str_detect(names_old, "beta")]
                                 }
                                 # Update intercept
                                 if (self$control$mode == "unitOmega") {
                                   private[["..omega"]] <- c(omega = (1 - sum(private[["..alpha"]]) - sum(private[["..beta"]]))[[1]])
                                 } else {
                                   private[["..omega"]] <- new_coefs["omega"]
                                 }
                                 # Matrix components
                                 # Companion matrix
                                 private[["..A"]] <- GARCH_companion_matrix(alpha = self$alpha, beta = self$beta)
                                 # Intercept vector
                                 d <- rep(0, length(private[["..b"]]))
                                 d[1] <- self$omega
                                 d <- matrix(d, ncol = 1)
                                 private[["..d"]] <- matrix(d, ncol = 1)
                               },
                               #' @description
                               #' Numerical computation of the std. errors of the parameters.
                               #' @param std.errors Numeric std. errors.
                               update_std.errors = function(std.errors){
                                 if (missing(std.errors) || purrr::is_empty(std.errors)) {
                                   return(invisible(NULL))
                                 }
                                 # Extract old coefficients
                                 new_std.errors <- private[["..std.errors"]]
                                 # Extract names
                                 names_old <- names(new_std.errors)
                                 names_new <- names(std.errors)
                                 # Update only if they are present
                                 for(i in 1:length(std.errors)){
                                   if (names_new[i] %in% names_old) {
                                     new_std.errors[names_new[i]] <- std.errors[i]
                                   }
                                 }
                                 # Update private std. errors
                                 private[["..std.errors"]] <- new_std.errors
                               },
                               #' @description
                               #' Next step GARCH std. deviation forecast
                               #' @param eps0 Numeric initial epsilons to prepend (length p+q).
                               #' @param sigma20 Numeric initial variances to prepend (length p+q).
                               next_step = function(eps0, sigma20){
                                 sGARCH_next_step(self$omega, self$alpha, self$beta, eps0, sigma20)
                               },
                               #' @description
                               #' Print method for `GARCH_modelR6` class.
                               print = function(){
                                 # Format parameters with std. errors
                                 green <- function(x) paste0("\033[1;32m", x, "\033[0m")
                                 red <- function(x) paste0("\033[1;31m", x, "\033[0m")
                                 msg_col <- function(x) ifelse(x, green(x), red(x))
                                 # Intercept
                                 par.intercept <- format(self$omega, digits = 4)
                                 std.error.intercept <- self$std.errors[names(self$std.errors) == names(par.intercept)]
                                 format.intercept <- paste0(par.intercept, " (", format(std.error.intercept, digits = 3), ")")
                                 # ARCH parameters
                                 par.arch <- format(self$alpha, digits = 4)
                                 std.error.arch <- self$std.errors[names(self$std.errors) %in% names(par.arch)]
                                 format.arch <- paste0(paste0(par.arch, " (", format(std.error.arch, digits = 3), ")"), collapse = " ")
                                 # GARCH parameters
                                 par.garch <- format(self$beta, digits = 4)
                                 std.error.garch <- self$std.errors[names(self$std.errors) %in% names(par.garch)]
                                 format.garch <- paste0(paste0(par.garch, " (", format(std.error.garch, digits = 3), ")"), collapse = " ")
                                 # ********************************************************
                                 msg0 <- paste0("--------------------- ", self$model_name, " ---------------------")
                                 cat(msg0,  "\n")
                                 cat(" ARCH: ", msg_col(self$archOrder != 0), "\n",
                                     "GARCH: ", msg_col(self$garchOrder != 0), "\n",
                                     "Intercept: ", format.intercept, " \n",
                                     "ARCH parameters: ", format.arch, " \n",
                                     "GARCH parameters: ", format.garch, " \n",
                                     "Log-lik: ", format(self$loglik, digits = 3), "\n",
                                     "Version: ", private$version, "\n",
                                     paste0(rep("-", length(strsplit(msg0, "")[[1]])-1), collapse = ""), "\n")
                               }
                             ),
                             private = list(
                               version = "1.0.3",
                               ..archOrder = 0,
                               ..garchOrder = 0,
                               ..omega = c(omega = 1),
                               ..alpha = c(alpha1 = 0),
                               ..beta = c(beta1 = 0),
                               ..A = matrix(),
                               ..b = matrix(),
                               ..d = matrix(),
                               ..log.likelihoods = c(NA),
                               ..std.errors = c(omega = NA, alpha1 = NA, beta1 = NA)
                             ),
                             active = list(
                               #' @field archOrder Numeric scalar, ARCH order.
                               archOrder = function(){
                                 private$..archOrder
                               },
                               #' @field garchOrder Numeric scalar, GARCH order.
                               garchOrder = function(){
                                 private$..garchOrder
                               },
                               #' @field order Numeric named vector, orders of the GARCH model. The first element is the ARCH order, while the second the GARCH order.
                               order = function(){
                                 c(ARCH = private$..archOrder, GARCH = private$..garchOrder)
                               },
                               #' @field omega Numeric scalar, intercept parameter.
                               omega = function(){
                                 private$..omega
                               },
                               #' @field alpha Numeric vector, ARCH parameters.
                               alpha = function(){
                                 private$..alpha
                               },
                               #' @field beta Numeric vector, GARCH parameters.
                               beta = function(){
                                 private$..beta
                               },
                               #' @field coefficients Numeric vector, model's coefficients.
                               coefficients = function(){
                                 c(private$..omega, private$..alpha, private$..beta)
                               },
                               #' @field A Numeric matrix, companion matrix.
                               A = function(){
                                 private$..A
                               },
                               #' @field b Numeric vector.
                               b = function(){
                                 private$..b
                               },
                               #' @field d Numeric vector
                               d = function(){
                                 private$..d
                               },
                               #' @field std.errors Numeric vector, standard errors of the model's coefficients.
                               std.errors = function(){
                                 private$..std.errors
                               },
                               #' @field sigma2_inf Numeric scalar, long-term unconditional std. deviation.
                               sigma2_inf = function(){
                                 private$..omega / (1 - sum(private$..alpha) + sum(private$..beta))
                               },
                               #' @field loglik model log-likelihood
                               loglik = function(){
                                 sum(private$..log.likelihoods)
                               },
                               #' @field tidy Tibble with estimated parameters and relative std. errors.
                               tidy = function(){
                                 # Extract the parameters
                                 params <- self$coefficients
                                 dplyr::tibble(
                                   term = names(params),
                                   estimate = params,
                                   std.error = self$std.errors,
                                 )
                               }
                             )
)
