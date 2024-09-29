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
                                initialize = function(alpha = 0, beta = 1){
                                  # Control parameters
                                  if (alpha < 0){
                                    stop("Alpha cannot be lower than zero.")
                                  } else if (alpha + beta > 1){
                                    stop("`alpha + beta` cannot be greater than one.")
                                  }
                                  # Store tranformation parameters
                                  private$..alpha <- alpha
                                  private$..beta <- beta
                                },
                                #' @description
                                #' Solar radiation function
                                #' @param x numeric vector in \eqn{(\alpha, \alpha+\beta)}.
                                #' @param Ct clear sky radiation.
                                #' @details The function computes:
                                #' \deqn{GHI(x) = C_t(1 - x)}
                                GHI = function(x, Ct) {
                                  Ct*(1-x)
                                },
                                #' @description
                                #' Solar radiation function in terms of y
                                #' @param y numeric vector in \eqn{(-\infty, \infty)}.
                                #' @param Ct clear sky radiation.
                                #' @details The function computes:
                                #' \deqn{GHI(y) = C_t(1 - \alpha-\beta \exp(-\exp(x)))}
                                GHI_y = function(y, Ct) {
                                  Ct*(1-self$iY(y))
                                },
                                #' @description
                                #' Compute the risk driver process for solar radiation
                                #' @param x numeric vector in \eqn{C_t (\alpha, \alpha+\beta)}.
                                #' @param Ct clear sky radiation.
                                #' @details The function computes the inverse of the `GHI`funcion
                                #' \deqn{iGHI(x) = 1 - \frac{x}{C_t}}
                                iGHI = function(x, Ct) {
                                  1 - x/Ct
                                },
                                #' @description
                                #' Transformation function from X to Y
                                #' @param x numeric vector in \eqn{(\alpha, \alpha+\beta)}.
                                #' @param inverse when `TRUE` will compute the inverse transform.
                                #' @details The function computes the transformation:
                                #' \deqn{Y(x) = \log(\log(\beta) - \log(x - \alpha))}
                                Y = function(x) {
                                  log(log(private$..beta) - log(x - private$..alpha))
                                },
                                #' @description
                                #' Inverse transformation from Y to X.
                                #' @param y numeric vector in \eqn{(-\infty, \infty)}.
                                #' @details The function computes the transformation:
                                #' \deqn{iY(y) = \alpha + \beta \exp(-\exp(y))}
                                iY = function(y) {
                                  private$..alpha + private$..beta*exp(-exp(y))
                                },
                                #' @description
                                #' Fit the best parameters from a time series
                                #' @param x time series of solar risk drivers in \eqn{(0, 1)}.
                                #' @param threshold for minimum
                                parameters = function(x, threshold = 0.01){
                                  # Upper and lower bounds
                                  range_Xt <- range(x)
                                  # Approximation parameter
                                  epsilon <- range_Xt[1]*threshold
                                  # Transformation parameters
                                  alpha_ <- range_Xt[1] - epsilon
                                  beta_ <- (range_Xt[2] - range_Xt[1]) + 2*epsilon
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
                                    alpha <- private$..alpha
                                  }
                                  if (missing(beta)){
                                    beta <- private$..beta
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
                                  private$..alpha <- alpha
                                  private$..beta <- beta
                                }
                              ),
                              private = list(
                                ..alpha = NA,
                                ..beta = NA
                              ),
                              active = list(
                                #' @description
                                #' Return the first transformation parameters
                                alpha = function(){
                                  private$..alpha
                                },
                                #' @description
                                #' Return the second transformation parameters
                                beta = function(){
                                  private$..beta
                                }
                              )
)


