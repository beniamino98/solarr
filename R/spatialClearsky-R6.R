#' Spatial clear-sky model
#'
#' @description R6 implementation of a clear-sky seasonal model for solar radiation,
#' using extraterrestrial radiation and seasonal harmonics.
#'
#' @docType class
#' @rdname spatialClearsky
#' @name spatialClearsky
#' @keywords clearsky
#' @note Version 1.0.1
#' @export
spatialClearsky <- R6::R6Class("spatialClearsky",
                               public = list(
                                 control = control_spatialClearsky(),
                                 #' @description
                                 #' Initialize a spatialClearsky model
                                 #' @param spec
                                 #' @param control description
                                 initialize = function(spec, control = control_spatialClearsky()){
                                   #' Specification
                                   private$..spec <- spec
                                   # Formula for clear-sky model
                                   formula_Ct <- spec$formula_Ct
                                   formula_params <- spec$formula_params
                                   # Initialize the matrix of parameters
                                   private$..coefficients <- matrix(0, nrow = spec$dim$p, ncol = spec$dim$K)
                                   # Add the function to compute extraterrestrial radiation
                                   private$..ssf <- seasonalSolarFunctions$new(spec$method_H0)
                                   # Control for fit
                                   self$control <- control
                                 },
                                 #' @description
                                 #' Initialize a spatialClearsky model
                                 fit = function(data_list, coords){
                                   # Extract specifications
                                   formula_Ct <- self$spec$formula_Ct
                                   formula_params <- self$spec$formula_params
                                   C_tilde_col <- self$spec$C_tilde_col
                                   R_col <- self$spec$R_col
                                   # Fit control
                                   control <- self$control
                                   lambda = control$lambda
                                   gamma = control$gamma
                                   eps = control$eps
                                   # if TRUE, interpret gamma as normalized weight and set lambda accordingly
                                   normalize_loss = control$normalize_loss
                                   # Type of penalty
                                   penalty = control$penalty
                                   # length p: weights for time-coef blocks (0=no penalty, 1=penalize)
                                   mask_time = control$mask_time
                                   # Settings
                                   settings = control$settings
                                   # Fit the parameters
                                   fit <- spatialClearsky_optimizer_CLS(formula_Ct,
                                                                        formula_params,
                                                                        eps = eps,
                                                                        data_list,
                                                                        coords,
                                                                        C_tilde_col = C_tilde_col,
                                                                        R_col = R_col,
                                                                        lambda = lambda,
                                                                        gamma = gamma,
                                                                        normalize_loss = normalize_loss,
                                                                        penalty = penalty,
                                                                        mask_time = mask_time,
                                                                        settings = settings)

                                   # Fitted parameters
                                   private$..coefficients <- fit$Beta_hat
                                 },
                                 #' @description
                                 #' Predict a clear-sky
                                 predict = function(n, lat, lon, alt = 0){
                                   # Project coordinates
                                   X <- project_CRS(lat, lon)
                                   # Create a dataset with latitude and longitude
                                   dd <- data.frame(lat = lat, lon = lon, x = 0)
                                   # Substitute projected coordinates
                                   dd <- mutate(dd, lat = X$lat, lon = X$lon)
                                   # Matrix of regressors (n x p)
                                   newdata <- data.frame(n = n, H0 = private$..ssf$Hon(n, lat, alt))
                                   newdata[[private$..spec$C_tilde_col]] <- 0
                                   # recycle if needed
                                   if (nrow(dd) == 1L && nrow(newdata) > 1L) dd <- dd[rep(1L, nrow(newdata)), , drop = FALSE]
                                   # Matrix of regressors (n x K)
                                   B_new <- model.matrix(private$..spec$formula_params, data = dd)
                                   # Matrix of regressors (n x p)
                                   X_new <- model.matrix(private$..spec$formula_Ct, data = newdata)     # n x p
                                   # compute C_hat rowwise: C = sum_{j=1..K} B_j * (X %*% Beta[,j])
                                   # (n x p) %*% (p x K) = (n x K), then elementwise * B_new, rowSums
                                   XB <- X_new %*% self$coefficients
                                   drop(rowSums(XB * B_new))
                                 },
                                 #' @description
                                 #' Predict the clear-sky parameters
                                 parameters = function(lat, lon){
                                   # Project coordinates
                                   X <- project_CRS(lat, lon)
                                   # Create a dataset with latitude and longitude
                                   dd <- data.frame(lat = X$lat, lon = X$lon, x = 0)
                                   # Matrix of regressors (1 x K)
                                   B_new <- model.matrix(private$..spec$formula_params, data = dd)
                                   # Forecast clear-sky coefficients
                                   c_hat <- B_new %*% t(self$coefficients)
                                   c_hat
                                 },
                                 #' @description
                                 #' Update the parameters
                                 update = function(coefficients){
                                   # Check dimension
                                   stopifnot(ncol(coefficients) == self$spec$dim$K & nrow(coefficients) == self$spec$dim$p)
                                   private$..coefficients <- coefficients
                                 },
                                 #' @description
                                 #'
                                 print = function(){
                                   msg_0 <- "-------------------------- spatialClearsky --------------------------"
                                   # Extract order
                                   spec <- private$..spec
                                   cat(paste0(msg_0, "\n"))
                                   msg_1 <- paste0("- Order: ", spec$order, "\n",
                                                   "- Order (H0): ", spec$order_H0, "\n",
                                                   "- Period: ", spec$period, "\n",
                                                   "- External regressors: 1 (H0)", "\n",
                                                   "- Linear trend: ", spec$include.trend, "\n",
                                                   "- Parameters: ", as.character(spec$formula_params), "\n",
                                                   "- Method: ", spec$method, "\n",
                                                   "- Version: ", private$version, "\n")
                                   msg_line <- paste0(rep("-", length(strsplit(msg_0, "")[[1]])), collapse = "")
                                   cat(msg_1)
                                   cat(paste0(msg_line, "\n"))
                                   print(self$coefficients)
                                 }
                               ),
                               private = list(
                                 version = "1.0.1",
                                 ..spec = NA,
                                 ..ssf = NA,
                                 ..coefficients = NA
                               ),
                               active = list(
                                 #' @field coefficients description
                                 coefficients = function(){
                                   private$..coefficients
                                 },
                                 #' @field spec description
                                 spec = function(){
                                   private$..spec
                                 }
                               ))
