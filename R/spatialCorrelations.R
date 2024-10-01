#' spatialCorrelation object
#'
#' @export
spatialCorrelation <- R6::R6Class("spatialCorrelation",
                                  public = list(
                                    initialize = function(binprobs, mixture_cr){
                                      # Store active elements
                                      private$..places = colnames(binprobs$commonprob[[1]])
                                      private$..sigma_B = binprobs$sigma
                                      private$..margprob = purrr::map(binprobs$commonprob, diag)
                                      private$..cr_X = mixture_cr
                                      # Store other variables only in private
                                      private$cr_B = binprobs$cr_B
                                      private$commonprob = binprobs$commonprob
                                      # Extract indexes to subset cr_X
                                      var_names_X <- rownames(mixture_cr)
                                      idx_x1 <- which(stringr::str_detect(var_names_X , "x1_"))
                                      idx_x2 <- which(stringr::str_detect(var_names_X , "x2_"))
                                      names(idx_x1) <- paste0("ID_", stringr::str_remove_all(var_names_X[idx_x1], "x1_"))
                                      names(idx_x2) <- paste0("ID_", stringr::str_remove_all(var_names_X[idx_x2], "x2_"))
                                      # Store the indexes
                                      private$idx = list(x1 = idx_x1, x2 = idx_x2)
                                    },
                                    get_sigma_B = function(places, nmonth = 1){
                                      # Retrieve active elements
                                      sigma_B = self$sigma_B[[nmonth]]
                                      # If places is not missing subset the matrices for such place
                                      if (!missing(places)) {
                                        # Check that places names are correct
                                        places = private$check_places(places)
                                        # Subset the matrix
                                        sigma_B = sigma_B[places, places]
                                      }
                                      return(sigma_B)
                                    },
                                    get_margprob = function(places, nmonth = 1){
                                      # Retrieve active elements
                                      margprob = self$margprob[[nmonth]]
                                      # If places is not missing subset the matrices for such place
                                      if (!missing(places)) {
                                        # Check that places names are correct
                                        places = private$check_places(places)
                                        # Subset the matrix
                                        margprob = margprob[places]
                                      }
                                      return(margprob)
                                    },
                                    get_cr_X = function(places, nmonth = 1){
                                      # Retrieve active elements
                                      cr_X = self$cr_X
                                      # If places is not missing subset the matrices for such place
                                      if (!missing(places)) {
                                        # Check that places names are correct
                                        places = private$check_places(places)
                                        # Subset the matrix
                                        cr_X = cr_X[private$index_X(places), private$index_X(places)]
                                      }
                                      return(cr_X)
                                    },
                                    get = function(places, nmonth = 1, date){
                                      # Compute reference month from date
                                      if (!missing(date)) {
                                        nmonth <- lubridate::month(as.Date(date))
                                      }
                                      # Retrieve the elements
                                      sigma_B = self$get_sigma_B(places, nmonth)
                                      margprob = self$get_margprob(places, nmonth)
                                      cr_X = self$get_cr_X(places, nmonth)
                                      # Output list
                                      list(
                                        sigma_B = sigma_B,
                                        margprob = margprob,
                                        cr_X = cr_X
                                      )
                                    }
                                  ),
                                  private = list(
                                    ..places = NA,
                                    ..sigma_B = NA,
                                    cr_B = NA,
                                    commonprob = NA,
                                    ..margprob = NA,
                                    ..cr_X = NA,
                                    idx = NA,
                                    index_X = function(places) {
                                      idx_x <- c(private$idx$x1[names(private$idx$x1) %in% places],
                                                 private$idx$x2[names(private$idx$x2) %in% places])
                                      return(idx_x[order(idx_x)])
                                    },
                                    check_places = function(places){
                                      # Retrieve active elements
                                      available_places = self$places
                                      # Check that places names are correct
                                      places = match.arg(places, choices = available_places, several.ok = TRUE)
                                      return(places)
                                    }
                                  ),
                                  active = list(
                                    places = function(){
                                      private$..places
                                    },
                                    sigma_B = function(){
                                      private$..sigma_B
                                    },
                                    cr_X = function(){
                                      private$..cr_X
                                    },
                                    margprob = function(){
                                      private$..margprob
                                    }
                                  ))


