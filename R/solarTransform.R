#' @title solarTransform
#' Solar functions
#' @description
#' Solar Model transformation functions
#'
#' @examples
#' st <- solarTransform$new()
#' st$GHI(0.4, 3)
#' st$GHI(st$iGHI(0.4, 3), 3)
#' @export
solarTransform <- R6::R6Class("solarTransform",
                              public = list(
                                #' @description
                                #' Initialize a solarTransform object.
                                #' @param alpha Numeric, transformation parameter.
                                #' @param beta Numeric, transformation parameter.
                                initialize = function(alpha = 0, beta = 1){
                                  # Control parameters
                                  if (alpha < 0) {
                                    stop("Alpha cannot be lower than zero.")
                                  } else if (alpha + beta > 1) {
                                    stop("`alpha + beta` cannot be greater than one.")
                                  }
                                  # Update the parameters
                                  private$..alpha <- alpha
                                  private$..beta <- beta
                                },
                                #' @method GHI solarTransform
                                #' @description
                                #' Solar radiation function
                                #' @param x Numeric values in \eqn{(\alpha, \alpha+\beta)}.
                                #' @param Ct Numeric, clear sky radiation.
                                #' @details The function computes:
                                #' \deqn{GHI(x) = C_t(1 - x)}
                                #' @return Numeric values in \eqn{C_t(1-\alpha-\beta, 1-\alpha)}.
                                GHI = function(x, Ct){
                                  Ct*(1 - x)
                                },
                                #' @method GHI_y solarTransform
                                #' @description
                                #' Solar radiation function in terms of y
                                #' @param y Numeric values in \eqn{(-\infty, \infty)}.
                                #' @param Ct Numeric, clear sky radiation.
                                #' @details The function computes:
                                #' \deqn{GHI(y) = C_t(1 - \alpha-\beta \exp(-\exp(x)))}
                                #' @return Numeric values in \eqn{[C_t(1-\alpha-\beta), C_t(1-\alpha)]}.
                                GHI_y = function(y, Ct){
                                  Ct*(1 - self$iY(y))
                                },
                                #' @method iGHI solarTransform
                                #' @description
                                #' Compute the risk driver process
                                #' @param x Numeric values in \eqn{[C_t(1-\alpha-\beta), C_t(1-\alpha)]}.
                                #' @param Ct Numeric, clear sky radiation.
                                #' @details The function computes the inverse of the `GHI`funcion
                                #' \deqn{iGHI(x) = 1 - \frac{x}{C_t}}
                                #' @return Numeric values in \eqn{[\alpha,\alpha+\beta]}.
                                iGHI = function(x, Ct){
                                  1 - x/Ct
                                },
                                #' @method Y solarTransform
                                #' @description
                                #' Transformation function from X to Y
                                #' @param x numeric vector in \eqn{[\alpha, \alpha+\beta]}.
                                #' @details The function computes:
                                #' \deqn{Y(x) = \log(\log(\beta) - \log(x - \alpha))}
                                #' @return Numeric values in \eqn{[-\infty, \infty]}.
                                Y = function(x){
                                  log(log(self$beta) - log(x - self$alpha))
                                },
                                #' @method iY solarTransform
                                #' @description
                                #' Inverse transformation from Y to X.
                                #' @param y numeric vector in \eqn{[-\infty, \infty]}.
                                #' @details The function computes:
                                #' \deqn{iY(y) = \alpha + \beta \exp(-\exp(y))}
                                #' @return Numeric values in \eqn{[\alpha, \alpha + \beta]}.
                                iY = function(y){
                                  self$alpha + self$beta*exp(-exp(y))
                                },
                                #' @method fit solarTransform
                                #' @description
                                #' Fit the best parameters from a time series
                                #' @param x time series of solar risk drivers in \eqn{(0, 1)}.
                                #' @param threshold for minimum
                                #' @details Return a list that contains:
                                #' \describe{
                                #'  \item{alpha}{Numeric, first transformation parameter.}
                                #'  \item{beta}{Numeric, second transformation parameter.}
                                #'  \item{epsilon}{Numeric, threshold used for fitting.}
                                #'  \item{Xt_min}{Numeric, minimum value of the supplied time series.}
                                #'  \item{Xt_min}{Numeric, maximum value of the supplied time series.}
                                #' }
                                #' @return A named list.
                                fit = function(x, threshold = 0.01){
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
                                #' Compute the bounds for the transformed variables.
                                #' @param target target variable. Available choices are:
                                #' \describe{
                                #'  \item{Xt}{Solar risk driver, the bounds returned are \eqn{[\alpha, \alpha + \beta]}.}
                                #'  \item{Kt}{Clearness index, the bounds returned are \eqn{[1-\alpha-\beta, 1-\alpha]}.}
                                #'  \item{Yt}{Solar transform, the bounds returned are \eqn{[-\infty, \infty]}.}
                                #'}
                                #' @return A numeric vector where the first element is the lower bound and the second the upper bound.
                                bounds = function(target = "Xt"){
                                  target = match.arg(target, choices = c("Xt", "Yt", "Kt"))
                                  lower <- c(Xt_min = self$alpha)
                                  upper <- c(Xt_max = self$alpha + self$beta)
                                  if (target == "Yt") {
                                    lower <-  c(Yt_min = -Inf)
                                    upper <-  c(Yt_max = Inf)
                                  } else if (target == "Kt") {
                                    lower <-  c(Kt_min = 1-self$alpha-self$beta)
                                    upper <-  c(Kt_max = 1-self$alpha)
                                  }
                                  return(c(lower, upper))
                                },
                                #' @method update solarTransform
                                #' @description
                                #' Update the transformation parameters
                                #' @param alpha Numeric, transformation parameter.
                                #' @param beta Numeric, transformation parameter.
                                #' @return Update the slots `$alpha` and `$beta`.
                                update = function(alpha, beta) {
                                  # Old parameters
                                  if (missing(alpha)) {
                                    alpha <- private$..alpha
                                  }
                                  if (missing(beta)) {
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
                                #' @field alpha Numeric, the first transformation parameter.
                                alpha = function(){
                                  private$..alpha
                                },
                                #' @field beta Numeric, the second transformation parameter.
                                beta = function(){
                                  private$..beta
                                }
                              )
                            )
