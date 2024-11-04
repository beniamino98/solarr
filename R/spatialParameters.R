#' `spatialParameters` object
#'
#' @rdname spatialParameters
#' @name spatialParameters
#' @export
spatialParameters <- R6::R6Class("spatialParameters",
                                 public = list(
                                   #' @field quiet Logical
                                   quiet = FALSE,
                                   #' @description
                                   #' Initialize a `spatialParameters` object
                                   #' @param data dataset with spatial parameters and lon, lat.
                                   #' @param params_names Names of the parameters to fit.
                                   #' @param models an optional list of kernelRegression models already fitted.
                                   #' @param sample List of parameter used as sample.
                                   initialize = function(data, params_names, models, sample){
                                     # Dataset containing all the parameters
                                     private$..data <- data
                                     # Names of all the parameters
                                     private$params_names <- params_names
                                     if (missing(models)) {
                                       # Initialize a list for all the models
                                       private$..models <- as.list(rep(NA, length(private$params_names)))
                                       # Name the list with the parameter's names
                                       names(private$..models) <- private$params_names
                                     } else {
                                       if (any(names(models) != private$params_names)){
                                         if (!self$quiet) warning("The parameters in `models` do not match the parameters extracted from `solarModel`!")
                                       } else {
                                         # Update the models with the given list
                                         private$..models <- models
                                       }
                                     }
                                     # Initialize the sample dataset to store the predictions
                                     private$sample <- sample
                                   },
                                   #' @description
                                   #' Fit a `kernelRegression` object for a parameter or a group of parameters.
                                   #' @param params list of parameters names to fit. When missing all the parameters will be fitted.
                                   fit = function(params){
                                     if (missing(params)) {
                                       params <- private$params_names
                                     } else {
                                       params <- match.arg(params, choices = private$params_names, several.ok = TRUE)
                                     }
                                     if (!self$quiet) message("Fitting ", length(params), " parameters...")
                                     # Fit the kernel regressions for the parameters
                                     i <- 1
                                     for(par in params){
                                       if(!self$quiet) message("Fit parameter: ", par, " ", i, "/", length(params))
                                       private$..models[[par]] <- kernelRegression$new()
                                       private$..models[[par]]$fit(paste0(par, "~lat+lon"), data = private$..data[, c(par, "lat","lon")])
                                       i <- i + 1
                                     }
                                   },
                                   #' @description
                                   #' Predict all the parameters for a specified location.
                                   #' @param lat numeric, latitude in degrees.
                                   #' @param lon numeric, longitude in degrees.
                                   #' @param as_tibble logical, when `TRUE` will be returned a `tibble`.
                                   predict = function(lat, lon, as_tibble = FALSE){
                                     predictions <- list()
                                     coords <- dplyr::tibble(lat = lat, lon = lon)
                                     for(i in 1:nrow(coords)){
                                       if(!self$quiet) message("Fitting parameters for latitude: ", coords$lat[i], " and longitude ", coords$lon[i])
                                       # Retrieve the sample dataset
                                       l_params <- private$sample
                                       # Substitute latitude and longitude
                                       l_params$location$place <- ""
                                       l_params$location$lat <- coords$lat[i]
                                       l_params$location$lon <- coords$lon[i]
                                       l_params$location$alt <- 0
                                       # Predict the parameters at the coordinates
                                       for(j in 2:length(l_params)){
                                         for(par in colnames(l_params[[j]])){
                                           newdata = data.frame(lat = coords$lat[i], lon = coords$lon[i])
                                           l_params[[j]][, par] <- private$..models[[par]]$predict(newdata)
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
                                   ..data = NA,
                                   params_names = NA,
                                   ..models = NA,
                                   sample = NA,
                                   deep_clone = function(name, value){
                                     if (name == "..models") {
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
                                   #' @field models list of `kernelRegression` objects
                                   models = function(){
                                     private$..models
                                   },
                                   #' @field data dataset with the parameters used for fitting
                                   data = function(){
                                     private$..data
                                   }
                                 ))
