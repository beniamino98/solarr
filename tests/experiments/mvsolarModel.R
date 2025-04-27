mvsolarModel <- R6::R6Class("mvsolarModel",
                          public = list(
                            #' @field place Character, optional name of the location considered.
                            place = NA_character_,
                            #' @field target Character, name of the target variable to model. Can be `GHI` or `clearsky`.
                            target = "GHI",
                            #' @field dates A named list, with the range of dates used in the model.
                            dates = list(),
                            #' @field coords A named list with the coordinates of the location considered. Contains:
                            #' \describe{
                            #'  \item{lat}{Numeric, reference latitude in degrees.}
                            #'  \item{lon}{Numeric, reference longitude in degrees.}
                            #'  \item{alt}{Numeric, reference altitude in metres.}
                            #'}
                            coords = list(lat = NA, lon = NA, alt = NA),
                            #' @description
                            #' Initialize a `solarModel`
                            #' @param spec an object with class `solarModelSpec`. See the function \code{\link{solarModel_spec}} for details.
                            initialize = function(spec){
                              # **************************************************** #
                              # Public components
                              self$place <- spec$place
                              self$dates <- spec$dates
                              self$coords <- spec$coord
                              # **************************************************** #
                              # GHI model
                              spec$target <- "GHI"
                              spec$control$stochastic_clearsky <- TRUE
                              private$..model_GHI <- solarModel$new(spec)
                              # Clearsky model
                              spec$target <- "clearsky"
                              spec$control$stochastic_clearsky <- FALSE
                              private$..model_Ct <- solarModel$new(spec)
                              # Monthly likelihood
                              private$..loglik <- dplyr::tibble(Month = 1:12, loglik = 0)
                            },
                            # ***************************************************************************** #
                            #' @description
                            #' Initialize and fit a \code{\link{solarModel}} object given the specification contained in `$control`.
                            fit = function(){
                              private$..model_GHI$fit()
                              private$..model_Ct$fit()
                              # Compute log-likelihood
                              self$logLik()
                            },
                            #' @description
                            #' Initialize and fit a `solarMixture` object.
                            fit_mixture_model = function(){

                            },
                            # ***************************************************************************** #
                            #' @description
                            #' Update the parameters inside object
                            #' @param params List of parameters. See the slot `$parameters` for a template.
                            update = function(params_GHI, params_Ct){
                              private$..model_GHI$update(params_GHI)
                              private$..model_Ct$update(params_Ct)
                            },
                            #' @description
                            #' Filter the time series when new parameters are supplied in the method `$update(params)`.
                            #' @return Update the slots `$data`, `$seasonal_data`, `$monthly_data`, `$moments$conditional`,
                            #' `$moments$unconditional` and `$loglik`.
                            filter = function(){
                              private$..model_GHI$filter(params_GHI)
                              private$..model_Ct$filter(params_Ct)
                            }
                          ),
                          # ***************************************************************************** #
                          private = list(
                            ..coords = NA,
                            ..model_GHI = NA,
                            ..model_Ct = NA,
                            ..loglik = NA
                          ),
                          # ***************************************************************************** #
                          active = list(
                            #' @field model_GHI A data frame with the fitted data, and the seasonal and monthly parameters.
                            model_GHI = function(){
                              private$..model_GHI
                            },
                            #' @field model_Ct A data frame containing seasonal and monthly parameters.
                            model_Ct = function(){
                              private$..model_Ct
                            },
                            #' @field loglik The log-likelihood computed on train data.
                            loglik = function(){
                              sum(private$..loglik)
                            },
                            #' @field location A data frame with coordinates of the location considered.
                            location = function(){
                              dplyr::bind_cols(place = self$place, target = self$target, dplyr::bind_rows(self$coords))
                            }
                          )
)


mod <- mvsolarModel$new(spec)

mod$fit()

mod$model_GHI

mod$model_Ct
