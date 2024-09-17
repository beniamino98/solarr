#' Solar Model transformation functions
#'
#' @export
solarTransform <- R6::R6Class("solarTransform",
                              public = list(
                                #' @description
                                #' Solar Model transformation functions
                                #' @param alpha bound parameters.
                                #' @param beta bound parameters.
                                #' @export
                                initialize = function(alpha = 0.001, beta = 0.95){
                                  # Control parameters
                                  if (alpha <= 0){
                                    stop("Alpha cannot be lower or equal than zero.")
                                  } else if (alpha + beta >= 1){
                                    stop("`alpha + beta` cannot be greater or equal than one.")
                                  }
                                  # Store tranformation parameters
                                  private$alpha <- alpha
                                  private$beta <- beta
                                },
                                #' @description
                                #' Transformation from Xt to Yt and viceversa
                                #' @param x numeric vector in \eqn{(\alpha, \alpha+\beta)}.
                                #' @param inverse when `TRUE` will compute the inverse transform.
                                #' @details The function computes the transformation:
                                #'
                                #' \deqn{Y_t = \log(\log(\beta) - \log(x - \alpha))}
                                #'
                                #' In case in which `inverse = TRUE` it computes
                                #'
                                #' \deqn{X_t = \alpha + \beta \exp(-\exp(Y_t))}
                                Yt = function(x, inverse = FALSE) {
                                  if (inverse) {
                                    private$alpha + private$beta*exp(-exp(x))
                                  } else {
                                    log(log(private$beta) - log(x - private$alpha))
                                  }

                                },
                                #' @description
                                #' Solar radiation function
                                #' @param x numeric vector in \eqn{(\alpha, \alpha+\beta)}.
                                #' @param Ct clearsky radiation.
                                #' @details The function computes the inverse of `GHI`
                                #'
                                #' \deqn{GHI = C_t(1 - x)}
                                #'
                                #' In case in which `inverse = TRUE` it computes
                                #'
                                #' \deqn{X_t = 1 - \frac{x}{C_t}}
                                GHI = function(x, Ct, inverse = FALSE) {
                                  if (inverse) {
                                    1 - x/Ct
                                  } else {
                                    Ct*(1-x)
                                  }
                                },
                                #' @description
                                #' Compute parameters from a time series
                                #' @param x time series of risk driver
                                #' @param threshold for minimum
                                parameters = function(x, threshold = 0.01){
                                  # Upper and lower bounds
                                  range_Xt <- range(x)
                                  # Approximation parameter
                                  epsilon <- range_Xt[1]*threshold
                                  # Transformation parameters
                                  alpha_ <- range_Xt[1] - epsilon
                                  beta_ <- range_Xt[2] - range_Xt[1] + 2*epsilon
                                  # Store Transform parameters
                                  # Transform parameters
                                  list(alpha = alpha_, beta = beta_, epsilon = epsilon,
                                       Xt_min = range_Xt[1], Xt_max = range_Xt[2])
                                },
                                #' @description
                                #' Update the parameters
                                #' @param alpha bounds parameter.
                                #' @param beta bounds parameter.
                                #' @param threshold for minimum
                                update = function(alpha, beta) {
                                  # Old parameters
                                  if (missing(alpha)){
                                    alpha <- private$alpha
                                  }
                                  if (missing(beta)){
                                    beta <- private$beta
                                  }
                                  # Control consistency
                                  if (alpha < 0){
                                    warning("Error: alpha is lower than 0")
                                    return(invisible(NULL))
                                  }
                                  if (alpha + beta > 1){
                                    warning("Error: alpha + beta is greater than 1")
                                    return(invisible(NULL))
                                  }
                                  # Update parameters
                                  private$alpha <- alpha
                                  private$beta <- beta
                                }
                              ),
                              private = list(
                                alpha = NA,
                                beta = NA
                              )
)



