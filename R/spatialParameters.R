#' `spatialParameters` object
#'
#' @rdname spatialParameters
#' @name spatialParameters
#' @export
spatialParameters <- R6::R6Class("spatialParameters",
                                 public = list(
                                   #' @description
                                   #' Initialize a `spatialParameters` object
                                   #' @param solarModels list of `solarModel` objects.
                                   #' @param models an optional list of models.
                                   #' @param quiet logical
                                   initialize = function(solarModels, models, quiet = FALSE){
                                     private$quiet <- quiet
                                     # Extract a dataset containing all the parameters
                                     private$..data <- purrr::map_df(solarModels, ~solarModel_parameters(.x, as_tibble = TRUE))
                                     # Extract the names of all the parameters
                                     private$params_names <- colnames(private$..data)[-c(1,2,3,4,5)]
                                     if (missing(models)) {
                                       # Initialize a list for all the models
                                       private$..models <- as.list(rep(NA, length(private$params_names)))
                                       # Name the list with the parameter's names
                                       names(private$..models) <- private$params_names
                                     } else {
                                       if (any(names(models) != private$params_names)){
                                         warning("The parameters in `models` do not match the parameters extracted from `solarModel`!")
                                       } else {
                                         # Update the models with the given list
                                         private$..models <- models
                                       }
                                     }
                                     # Initialize the sample dataset to store the predictions
                                     private$sample <- solarModel_parameters(solarModels[[1]], as_tibble = FALSE)
                                   },
                                   #' @description
                                   #' Fit a `kernelRegression` object for a parameter or a group of parameters.
                                   #' @param params list of parameters names to fit. When missing all the parameters will be fitted.
                                   fit = function(params){
                                     if (missing(params)){
                                       params <- private$params_names
                                     } else {
                                       params <- match.arg(params, choices = private$params_names, several.ok = TRUE)
                                     }
                                     if(!private$quiet) message("Fitting ", length(params), " parameters...")
                                     # Fit the kernel regressions for the parameters
                                     i <- 1
                                     for(par in params){
                                       if(!private$quiet) message("Fit parameter: ", par, " ", i, "/", length(params))
                                       private$..models[[par]] <- kernelRegression$new(paste0(par, "~lat+lon"), data = private$..data)
                                       i <- i + 1
                                     }
                                   },
                                   #' @description
                                   #' Predict all the parameters for a specified location.
                                   #' @param lat numeric, latitude in degrees.
                                   #' @param lon numeric, longitude in degrees.
                                   #' @param as_tibble logical, when `TRUE` will be returned a `tibble`.
                                   predict = function(lat, lon, as_tibble = FALSE){
                                     i <- 1
                                     predictions <- list()
                                     for(i in 1:length(lat)){
                                       if(!private$quiet) message("Fitting parameters for latitude: ", lat[i], " and longitude ", lon[i])
                                       # Retrieve the sample dataset
                                       l_params <- private$sample
                                       # Substitute latitude and longitude
                                       l_params$location$place <- ""
                                       l_params$location$lat <- lat[i]
                                       l_params$location$lon <- lon[i]
                                       l_params$location$alt <- 0
                                       # Predict the parameters at the coordinates
                                       for(j in 2:length(l_params)){
                                         for(par in colnames(l_params[[j]])){
                                           l_params[[j]][, par] <- private$..models[[par]]$predict(lat = lat[i], lon = lon[i])
                                         }
                                       }
                                       # Store the predictions
                                       predictions[[i]] <- l_params
                                       # Convert in a dataset
                                       if (as_tibble) {
                                         predictions[[i]] <- dplyr::bind_cols(predictions[[i]])
                                       }
                                     }
                                     # Convert in a dataset
                                     if (as_tibble) {
                                       predictions <- dplyr::bind_rows(predictions)
                                     }
                                     return(predictions)

                                   }
                                 ),
                                 private = list(
                                   quiet = FALSE,
                                   ..data = NA,
                                   params_names = NA,
                                   ..models = NA,
                                   sample = NA
                                 ),
                                 active = list(
                                   models = function(){
                                     private$..models
                                   },
                                   data = function(){
                                     private$..data
                                   }
                                 ))

#' Kernel regression
#'
#' Fit a kernel regression.
#'
#' @param formula formula
#' @param data data
#' @param ... other parameters to be passed to. See \code{np::npreg}.
#'
#' @rdname kernelRegression
#' @name kernelRegression
#' @export
kernelRegression <- R6::R6Class("kernelRegression",
                                public = list(
                                  #' @description
                                  #' Initialize a `kernelRegression` object
                                  #' @param formula formula, an object of class `formula` (or one that can be coerced to that class).
                                  #' @param data 	an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
                                  #' If not found in data, the variables are taken from environment(formula), typically the environment from which `lm` is called.
                                  #' @param ... other parameters to be passed to the function `np::npreg`.
                                  initialize = function(formula, data, ...){
                                    # Model formula
                                    private$formula <- as.formula(formula)
                                    # Fit a kernel regression
                                    private$..model <- np::npreg(private$formula, data = data, ...)
                                    # Store regressors names
                                    private$external_regressors_names <- formula.tools::rhs.vars(private$formula)
                                  },
                                  #' @description
                                  #' Predict method
                                  #' @param ... arguments to fit.
                                  predict = function(...){
                                    l <- list(...)
                                    if (is.null(names(l))) {
                                      names(l) <- private$external_regressors_names
                                    }
                                    np:::predict.npregression(private$..model, newdata = dplyr::bind_rows(l))
                                  }
                                ),
                                private = list(
                                  ..model = NA,
                                  formula = NA,
                                  external_regressors_names = NA
                                ),
                                active = list(
                                  model = function(){
                                    private$..model
                                  }
                                ))

