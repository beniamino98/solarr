#' Create a function of time for monthly parameters
#'
#' @param params Vector of length 12 with the monthly parameters.
#'
#' @examples
#' set.seed(1)
#' params <- runif(12)
#' mp <- monthlyParams$new(params)
#' t_now <- as.Date("2022-01-01")
#' t_hor <- as.Date("2024-12-31")
#' dates <- seq.Date(t_now, t_hor, by = "1 day")
#' plot(mp$predict(dates), type = "l")
#' @export
monthlyParams <- R6::R6Class("monthlyParams",
                             public = list(
                               initialize = function(params){
                                 if (length(params) != 12) {
                                   stop("The length of the vector of parameters must be 12!")
                                 }
                                 private$..parameters <- params
                               },
                               predict = function(x){
                                 nmonths <- lubridate::month(x)
                                 par <- c()
                                 for(i in 1:length(nmonths)){
                                   par[i] <- self$parameters[nmonths[i]]
                                 }
                                 return(par)
                               }
                             ),
                             private = list(
                               ..parameters = NA
                             ),
                             active = list(
                               parameters = function(){
                                 private$..parameters
                               }
                             ))
