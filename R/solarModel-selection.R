#' Fit a grid of solarModels
#'
#' @param spec Specification
#' @param place Reference location
#' @param arOrder Numeric, maximum AR order.
#' @param maOrder Numeric, maximum MA order.
#' @param archOrder Numeric, maximum ARCH order.
#' @param garchOrder Numeric, maximum GARCH order.
#' @examples
#' spec <- solarModel_spec$new()
#' spec$specification("Bologna")
#' models <- solarModel_grid(spec, 1,1,1,1)
#' models
#'
#' @rdname solarModel_grid
#' @name solarModel_grid
#' @keywords solarModel
#' @note Version 1.0.0.
#' @export
solarModel_grid <- function(spec, arOrder = 2, maOrder = 2, archOrder = 1, garchOrder = 1, QMLE = FALSE){

  AR_order <- paste0("ARMA(", 1:arOrder)
  MA_order <- paste0(",", 0:maOrder, ")")
  grid <- expand.grid(x = AR_order, y = MA_order)
  for(i in 1:nrow(grid)){
    grid[i,3] <- paste0(grid[i,1], grid[i,2])
  }
  grid_AR_MA <- grid[,3]

  ARCH_order <- paste0("GARCH(", 0:archOrder)
  GARCH_order <- paste0(",", 0:garchOrder, ")")
  grid <- expand.grid(x = ARCH_order, y = GARCH_order)
  for(i in 1:nrow(grid)){
    if (grid[i,1] == "GARCH(0" & grid[i,2] != ",0)"){
      next
    }
    grid[i,3] <- paste0(grid[i,1], grid[i,2])
  }
  grid_ARCH_GARCH <- na.omit(grid[,3])

  grid <- expand.grid(x = grid_AR_MA, y = grid_ARCH_GARCH)
  for(i in 1:nrow(grid)){
    grid[i,3] <- paste0(grid[i,1], "-", grid[i,2])
  }
  grid_models <- grid[,3]
  # Extract model orders
  order <- stringr::str_extract_all(grid_models, "[0-9],[0-9]")
  order <- purrr::map(order, ~as.numeric(unlist(stringr::str_split(.x, ","))))

  models <- list()
  models_name <- c()
  for(i in 1:length(order)){
    print(paste0("Fitting: ", i, "/", length(order)))
    # Specification
    spec$set_mean.model(arOrder = order[[i]][1], maOrder = order[[i]][2])
    if (order[[i]][3] == 0 & order[[i]][4] == 0){
      spec$set_variance.model(archOrder = order[[i]][3], garchOrder = order[[i]][4])
    } else {
      spec$set_variance.model(archOrder = order[[i]][3], garchOrder = order[[i]][4])
    }
    models_name[i] <- spec$model_name
    print(models_name[i])
    # Initialize the model
    model <- solarModel$new(spec)
    # Model fit
    model$fit()
    print(model$loglik)
    if (QMLE) {
      model <- solarModel_QMLE(model)
      print(model$loglik)
    }
    # Fit result
    models[[i]] <- model$clone(TRUE)
  }
  names(models) <- models_name
  models_name <- names(models)

  dplyr::tibble(
    model.name = models_name,
    link = model$spec$transform$link,
    AR = purrr::map_dbl(order, ~.x[1]),
    MA = purrr::map_dbl(order, ~.x[2]),
    ARCH = purrr::map_dbl(order, ~.x[3]),
    GARCH = purrr::map_dbl(order, ~.x[4]),
    spec = purrr::map(models, ~.x$spec$clone(TRUE)),
    model = models
  )
}

#' Compute the AIC and BIC of a solarModel object
#'
#' @param model solarmodel
#'
#' @rdname solarModel_AIC_BIC
#' @name solarModel_AIC_BIC
#' @keywords solarModel
#' @note Version 1.0.0.
#' @export
solarModel_AIC_BIC <- function(model, target = "GHI", type = "full"){

  type <- match.arg(type, choices = c("train", "test", "full"))

  moments <- model$moments$conditional
  if (type == "train") {
    moments <-  moments[model$data$isTrain,]
  } else if (type == "test"){
    moments <-  moments[!model$data$isTrain,]
  }
  # Log-likelihood
  L <- model$logLik(moments, target = target)
  n <- sum(!is.infinite(L))
  L <- sum(L[!is.infinite(L)])
  # Model's parameters
  params <- unlist(model$coefficients[-c(1,2)])
  # Remove zero params
  params <- params[params!=0]
  # Number of parameters (minus omega)
  k <- length(params) - 1
  # AIC and BIC
  AIC <- -2 * L + 2 * k
  BIC <- log(n) * k - 2 * L

  dplyr::tibble(
    Place = model$place,
    Model = model$spec$model_name,
    p = model$spec$mean.model$arOrder,
    q = model$spec$mean.model$maOrder,
    r = model$spec$variance.model$archOrder,
    s = model$spec$variance.model$garchOrder,
    ARMA = list(ARMA = model$spec$mean.model$order),
    GARCH = list(GARCH = model$spec$variance.model$order),
    Spec = list(model$spec$clone(TRUE)),
    L = L,
    k = k,
    n = n,
    AIC = AIC,
    BIC = BIC
  )
}

