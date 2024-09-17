#' Spatial kernel regression
#'
#' Fit kernel regression on all the parameters of a list containing
#' `solarModels` at different coordinates.
#'
#' @param solarModels a list containing `solarModels` objects.
#' @param quiet logical
#'
#' @rdname spatialParameters
#' @name spatialParameters
#' @export
spatialParameters <- function(solarModels, quiet = FALSE){

  # Extract a dataset containing all the parameters
  df_params <- purrr::map_df(solarModels, ~solarModel_parameters(.x, as_tibble = TRUE))
  # Extract the names of all the parameters
  params_names <- colnames(df_params)[-c(1,2,3)]

  # Initialize a list for all the models
  params_models <- list()
  # Fit the kernel regressions for all the parameters
  for(i in 1:length(params_names)){
    if(!quiet) message("Fit parameter: ", params_names[i], " ", i, "/", length(params_names))
    params_models[[i]] <- kernelRegression(paste0(params_names[i], "~lat+lon"), data = df_params)
  }
  # Name the list with the parameter's names
  names(params_models) <- params_names

  # Initialize the list version of the dataset to store the functions
  l_params <- solarModel_parameters(solarModels[[1]], as_tibble = FALSE)
  # Store the fit function
  for(i in 2:(length(l_params))){
    for(j in 1:length(l_params[[i]])){
      l_params[[i]][,j][[1]] <- list(f = params_models[[names(l_params[[i]][,j])]]$fit_function)
    }
  }
  # Output
  structure(
    list(
      l_params = l_params,
      models = params_models
    ),
    class = c("spatialParameters", "list")
  )
}

#' Predict method
#'
#' Predict method for the class `spatialParameters`.
#'
#' @param object an object of the class `spatialParameters`. See \code{\link{clearskyModel}}.
#' @param lat numeric latitude of the locations.
#' @param lon numeric longitude of the locations.
#'
#' @rdname spatialParameters_predict
#' @name spatialParameters_predict
#' @export
spatialParameters_predict <- function(object, lat, lon, as_tibble = FALSE, quiet = FALSE){

  i <- 1
  predictions <- list()
  for(i in 1:length(lat)){
    if(!quiet) message("Fitting parameters for latitude: ", lat[i], " and longitude ", lon[i])
    # Substitute latitude and longitude
    l_params <- object$l_params
    l_params$location$place <- ""
    l_params$location$lat <- lat[i]
    l_params$location$lon <- lon[i]
    # Predict the parameters at the coordinates
    for(j in 2:length(l_params)){
      l_params[[j]] <- purrr::map_df(l_params[[j]], ~.x$f(lat = lat[i], lon = lon[i]))
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
kernelRegression <- function(formula, data, ...){

  # Model formula
  formula <- as.formula(formula)
  # Fit a kernel regression
  kernel_model <- np::npreg(formula, data = data, ...)
  # Custom fit function
  fit_function = function(...){
    l <- list(...)
    if (is.null(names(l))) {
      names(l) <- formula.tools::rhs.vars(formula)
    }
    predict(kernel_model, newdata = dplyr::bind_rows(l))
  }

  # Output object
  structure(
    list(
      model = kernel_model,
      fit_function = fit_function
    ),
    class = c("kernelRegression", "list")
  )
}
