#' Monthly Gaussian mixture with two components
#'
#' @rdname solarMixture
#' @name solarMixture
#' @export
solarMixture <-  R6::R6Class("solarMixture",
                             public = list(
                               #' @field maxit Integer, maximum number of iterations.
                               maxit = 100,
                               #' @field abstol Numeric, absolute level for convergence.
                               abstol = 10e-5,
                               #' @field components Integer, number of components.
                               components = 2,
                               #' @description
                               #' Initialize a `solarMixture` object
                               #' @param components Integer, number of components.
                               #' @param maxit Integer, maximum number of iterations.
                               #' @param abstol Numeric, absolute level for convergence.
                               initialize = function(components = 2, maxit = 100, abstol = 10e-15){
                                 self$components <- components
                                 self$maxit <- maxit
                                 self$abstol <- abstol
                               },
                               #' @description
                               #' Fit the parameters with mclust package
                               #' @param x vector
                               #' @param date date vector
                               #' @param weights observations weights, if a weight is equal to zero the observation is excluded, otherwise is included with unitary weight.
                               #' When `missing` all the available observations will be used.
                               #' @param EM logical when TRUE will be applied EM algo
                               fit = function(x, date, weights, EM = FALSE){
                                 # Add data
                                 private$x <- x
                                 private$date <- date
                                 # Weights
                                 if (missing(weights)){
                                   w <- rep(1, n)
                                 } else {
                                   w <- ifelse(weights == 0, 0, 1)
                                 }
                                 w[is.na(x)] <- 0
                                 private$w <- w

                                 data <- self$data
                                 # Gaussian Mixture parameters
                                 GM_models <- list()
                                 # Monthly data
                                 data_months <- list()
                                 for(m in 1:12){
                                   data_months[[m]] <- dplyr::filter(data, Month == m)
                                   w <- data_months[[m]]$w
                                   # Monthly data
                                   eps <- data_months[[m]]$x
                                   # Fitted model
                                   GM_models[[m]] <- gaussianMixture$new(maxit = self$maxit, abstol = self$abstol, components = self$components)
                                   GM_models[[m]]$fit(eps, w)
                                   if (EM) GM_models[[m]]$EM(eps, w)
                                   # Add fitted Bernoulli series
                                   data_months[[m]]$B <- GM_models[[m]]$fitted$B1
                                   # Reorder parameters
                                   #if (GM_models[[m]]$means[1] < GM_models[[m]]$means[2]) {
                                   #   data_months[[m]]$B <- GM_models[[m]]$fitted$B2
                                   # }
                                 }
                                 # Add models
                                 private$model1 = GM_models[[1]]$clone(deep = TRUE)
                                 private$model2 = GM_models[[2]]$clone(deep = TRUE)
                                 private$model3 = GM_models[[3]]$clone(deep = TRUE)
                                 private$model4 = GM_models[[4]]$clone(deep = TRUE)
                                 private$model5 = GM_models[[5]]$clone(deep = TRUE)
                                 private$model6 = GM_models[[6]]$clone(deep = TRUE)
                                 private$model7 = GM_models[[7]]$clone(deep = TRUE)
                                 private$model8 = GM_models[[8]]$clone(deep = TRUE)
                                 private$model9 = GM_models[[9]]$clone(deep = TRUE)
                                 private$model10 = GM_models[[10]]$clone(deep = TRUE)
                                 private$model11 = GM_models[[11]]$clone(deep = TRUE)
                                 private$model12 = GM_models[[12]]$clone(deep = TRUE)
                                 # Add fitted values
                                 private$..fitted <- dplyr::select(dplyr::bind_rows(data_months), date, B)
                               }
                             ),
                             private = list(
                               x = NA,
                               w = NA,
                               date = NA,
                               model1 = list(),
                               model2 = list(),
                               model3 = list(),
                               model4 = list(),
                               model5 = list(),
                               model6 = list(),
                               model7 = list(),
                               model8 = list(),
                               model9 = list(),
                               model10 = list(),
                               model11 = list(),
                               model12 = list(),
                               ..fitted = list()
                             ),
                             active = list(
                               #' @field data A tibble with the following columns:
                               #' #' \describe{
                               #'  \item{date}{Time series of dates.}
                               #'  \item{Month}{Vector of Month.}
                               #'  \item{x}{Time series used for fitting.}
                               #'  \item{w}{Time series of weights.}
                               #'}
                               data = function(){
                                 dplyr::tibble(date = private$date, Month = lubridate::month(date), x = private$x, w = private$w)
                               },
                               #' @field model Named List with 12 \code{\link{gaussianMixture}} objects.
                               model = function(){
                                 GM_models <- list(private$model1, private$model2, private$model3, private$model4,
                                                   private$model5, private$model6, private$model7, private$model8,
                                                   private$model9, private$model10, private$model11, private$model12)
                                 names(GM_models) <- lubridate::month(1:12, label = TRUE)
                                 return(GM_models)
                               },
                               #' @field loglik Numeric, total log-likelihood.
                               loglik = function(){
                                 sum(purrr::map_dbl(self$model, ~.x$loglik))
                               },
                               #' @field fitted A `tibble` with the classified series
                               fitted = function(){
                                 private$..fitted
                               },
                               #' @field moments A `tibble` with the theoric moments. It contains:
                               #' #' \describe{
                               #'  \item{Month}{Month of the year.}
                               #'  \item{mean}{Theoric monthly expected value of the mixture model.}
                               #'  \item{variance}{Theoric monthly variance of the mixture model.}
                               #'  \item{skewness}{Theoric monthly skewness.}
                               #'  \item{kurtosis}{Theoric monthly kurtosis.}
                               #'  \item{nobs}{Number of observations used for fitting.}
                               #'  \item{loglik}{Monthly log-likelihood.}
                               #'}
                               moments = function(){
                                 dplyr::bind_cols(Month = 1:12, purrr::map_df(self$model, ~.x$moments), loglik = purrr::map_dbl(self$model, ~.x$loglik))
                               },
                               #' @field parameters A `tibble` with the fitted parameters.
                               parameters = function(){
                                 dplyr::bind_cols(Month = 1:12, purrr::map_df(self$model, ~.x$model))
                               }
                             )
                            )



#' Monthly multivariate Gaussian mixture with two components
#'
#' @param model_Ct arg
#' @param model_GHI arg
#'
#' @rdname solarModel_mvmixture
#' @name solarModel_mvmixture
#' @export
solarModel_mvmixture <- function(model_Ct, model_GHI){

  # Extract a bivariate dataset
  data_GHI <- dplyr::select(model_GHI$data, date, Month, GHI = "u_tilde")
  data_Ct <- dplyr::select(model_Ct$data, date, Month, Ct = "u_tilde")
  data <- dplyr::inner_join(data_GHI, data_Ct, by = c("date", "Month"))
  # Remove outliers
  outliers_date <- c(model_Ct$.__enclos_env__$private$outliers$date, model_GHI$.__enclos_env__$private$outliers$date)
  data <- dplyr::filter(data, !(date %in% outliers_date))

  # Extract Gaussian mixture parameters
  NM_model_GHI <- model_GHI$.__enclos_env__$private$..NM_model
  NM_model_GHI$rho_up <- NM_model_GHI$rho_dw <- 0
  NM_model_Ct <- model_Ct$.__enclos_env__$private$..NM_model
  NM_model_Ct$rho_up <- NM_model_Ct$rho_dw <- 0

  # Gaussian Mixture parameters
  m <- 1
  # Monthly data
  data_months <- list()
  for(m in 1:12){
    data_months[[m]] <- dplyr::filter(data, Month == m)
    # Monthly data
    eps <- data_months[[m]][,c(3,4)]
    # Fitted model
    gm <- mvgaussianMixture(eps, components = 2, na.rm = TRUE)
    # Update Gaussian mixture parameters (GHI)
    NM_model_GHI$mu_up[m] <- gm$params$means[1,2]
    NM_model_GHI$mu_dw[m] <- gm$params$means[2,2]
    NM_model_GHI$sd_up[m] <- sqrt(gm$params$sigma2[1,2])
    NM_model_GHI$sd_dw[m] <- sqrt(gm$params$sigma2[2,2])
    NM_model_GHI$p_up[m] <- gm$params$p[1]
    NM_model_GHI$p_dw[m] <- gm$params$p[2]
    NM_model_GHI$rho_up[m] <- gm$params$rho[1]
    NM_model_GHI$rho_dw[m] <- gm$params$rho[2]
    # Update Gaussian mixture parameters (Ct)
    NM_model_Ct$mu_up[m] <- gm$params$means[1,1]
    NM_model_Ct$mu_dw[m] <- gm$params$means[2,1]
    NM_model_Ct$sd_up[m] <- sqrt(gm$params$sigma2[1,1])
    NM_model_Ct$sd_dw[m] <- sqrt(gm$params$sigma2[2,1])
    NM_model_Ct$p_up[m] <- gm$params$p[1]
    NM_model_Ct$p_dw[m] <- gm$params$p[2]
    NM_model_Ct$rho_up[m] <- gm$params$rho[1]
    NM_model_Ct$rho_dw[m] <- gm$params$rho[2]
    # Add fitted Bernoulli series
    data_months[[m]]$B <- gm$B_hat$B1
  }
  # Fitted series of bernoulli
  data_months <- dplyr::select(dplyr::bind_rows(data_months), date, B)

  # Update data (Ct)
  data_Ct <- model_Ct$.__enclos_env__$private$..data
  data_Ct <- dplyr::left_join(dplyr::select(data_Ct, -B), data_months, by = "date")
  data_Ct$B[is.na(data_Ct$B)] <- 0
  model_Ct$.__enclos_env__$private$..data <- data_Ct
  model_Ct$.__enclos_env__$private$..NM_model <- NM_model_Ct
  # Update data (GHI)
  data_GHI <- model_GHI$.__enclos_env__$private$..data
  data_GHI <- dplyr::left_join(dplyr::select(data_GHI, -B), data_months, by = "date")
  data_GHI$B[is.na(data_GHI$B)] <- 0
  model_GHI$.__enclos_env__$private$..data <- data_GHI
  model_GHI$.__enclos_env__$private$..NM_model <- NM_model_GHI

  structure(
    list(
      model_Ct = model_Ct,
      model_GHI = model_GHI
    ),
    class = c("solarModelMixture", "list")
  )
}


