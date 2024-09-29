#' spatialCorrelation object
#'
#' @noRd
#' @export
spatialCorrelation_mixture <- function(models, nmonths){

  if (missing(nmonths)){
    filter_data <- function(data) dplyr::filter(data, isTrain & w != 0)
  } else {
    filter_data <- function(data) dplyr::filter(data, Month %in% nmonths & isTrain & w != 0)
  }

  # Number of models
  n_models <- length(models)
  # Extract a dataset containing all U and B variables
  data <- dplyr::bind_cols(models[[1]]$data, w = models[[1]]$outliers$weights)
  data <- filter_data(data)
  # Rename the variables
  data <- dplyr::left_join(data, models[[1]]$NM_model, by = "Month")
  df_u <- dplyr::tibble(date = data$date, u1 = data$u, B1 = data$B, x1_1 = ((u1-data$mu_up)/data$sd_up)*B1, x2_1 = ((u1-data$mu_dw)/data$sd_dw)*(1-B1))
  df_u <- dplyr::select(df_u, -B1)
  # Create a unique dataset
  for(i in 2:n_models){
    message(i, "/", n_models, "\r", appendLF = FALSE)
    data <- dplyr::bind_cols(models[[i]]$data, w = models[[i]]$outliers$weights)
    data <- filter_data(data)
    data <- dplyr::left_join(data, models[[i]]$NM_model, by = "Month")
    data <- dplyr::tibble(date = data$date, u1 = data$u, B1 = data$B, x1 = ((u1-data$mu_up)/data$sd_up)*B1, x2 = ((u1-data$mu_dw)/data$sd_dw)*(1-B1))
    data <- dplyr::select(data, -B1)
    colnames(data) <- c("date", paste0(c("u", "x1_", "x2_"), i))
    df_u <- dplyr::left_join(df_u, data, by = "date")
  }
  # Compute correlations
  message("Computing correlations... ")
  if(missing(nmonths)){
    df_u <- na.omit(df_u)[,-1]
    # Split the datasets
    data_12 <- dplyr::select(df_u, dplyr::contains("x"))
    message("Computing mixture correlations: X1-X1 | X2-X2 | X1-X2 | X2-X1 |")
    # Store the results
    cr_list <- cor(data_12)
  } else {
    m <- 1
    cr_list <- list()
    df_u$Month <- lubridate::month(df_u$date)
    df_u <- na.omit(df_u)[,-1]
    for(m in nmonths){
      message("Computing correlation for Month: ", m, "...")
      # Split the datasets
      data_12 <- dplyr::select(dplyr::filter(df_u, Month == m), dplyr::contains("x"))
      message("Computing mixture correlations: X1-X1 | X2-X2 | X1-X2 | X2-X1 |")
      # Store the results
      cr_list[[m]] <- cor(data_12)
    }
    names(cr_list) <- lubridate::month(nmonths, label = TRUE)
  }
  # Unique names as models names
  return(cr_list)
}

#' spatialCorrelation object
#'
#' @noRd
#' @export
spatialCorrelation_binprobs <- function(models, nmonths = 1:12){
  # Number of models
  n_models <- length(models)
  # Extract a dataset containing all U and B variables
  data <- dplyr::bind_cols(models[[1]]$data, w = models[[1]]$outliers$weights)
  data <- dplyr::filter(data, Month %in% nmonths & isTrain & w != 0)
  # Rename the variables
  df_B <- dplyr::tibble(date = data$date, u1 = data$u, B1 = data$B)
  # Create a unique dataset
  for(i in 2:n_models){
    message(i, "/", n_models, "\r", appendLF = FALSE)
    data <- dplyr::bind_cols(models[[i]]$data, w = models[[i]]$outliers$weights)
    data <- dplyr::filter(data, Month %in% nmonths & isTrain & w != 0)
    data <- dplyr::tibble(date = data$date, u = data$u, B = data$B)
    colnames(data) <- c("date", paste0(c("u", "B"), i))
    df_B <- dplyr::left_join(df_B, data, by = "date")
  }
  df_B <- na.omit(df_B)
  # Add month variable
  df_B$Month <- lubridate::month(df_B$date)
  # Compute correlations
  cr_list <- list()
  cr_list$sigma <- cr_list$commonprob <- cr_list$cr_B <- list()
  m <- 1
  for(m in nmonths){
    message("Computing correlation for Month: ", m, "...")
    # Split the datasets
    data_B <- dplyr::select(dplyr::filter(df_B, Month == m), dplyr::contains("B"))
    message("Computing mixture correlations...", appendLF = FALSE)
    cr_list$cr_B[[m]] <- cor(data_B)
    message("Done!", appendLF = TRUE)
    message("Computing marginal probabilities B...", appendLF = FALSE)
    # Marginal probabilities
    cr_list$commonprob[[m]] <- diag(apply(data_B, 2, mean, na.rm = TRUE))
    message("Done!", appendLF = TRUE)
    message("Computing joint probabilities of event [Bi = 1 & Bj = 1]...")
    # Common probabilities
    for(i in 1:(n_models-1)){
      for(j in (i+1):n_models){
        pairwise_probs_ij <-  mean(data_B[,i] == 1 & data_B[,j] == 1, na.rm = TRUE)
        cr_list$commonprob[[m]][i,j] <- cr_list$commonprob[[m]][j,i] <- pairwise_probs_ij
      }
      if (i == ncol(data_B)) {
        message("Done!")
      } else {
        message(i, "/", ncol(data_B), "\r", appendLF = FALSE)
      }
    }
    message("Adjusting commonprobs...")
    # Convert commonprob into a matrix sigma to improve simulations
    cr_list$commonprob[[m]] <- interpolate_commonprob(cr_list$commonprob[[m]], check.commonprob = FALSE)
    message("Done!", appendLF = TRUE)
    message("Computing implied sigma from commonprob...", appendLF = FALSE)
    # Adjust implied matrix sigma
    cr_list$sigma[[m]] <- bindata::commonprob2sigma(cr_list$commonprob[[m]])
    cr_list$sigma[[m]] <- makeSemiPositive(cr_list$sigma[[m]], neg_values = 1e-3)
    message("Done!", appendLF = TRUE)
    # Unique names as models names
    colnames(cr_list$cr_B[[m]]) <- rownames(cr_list$cr_B[[m]]) <- names(models)
    colnames(cr_list$commonprob[[m]]) <- rownames(cr_list$commonprob[[m]]) <- names(models)
    colnames(cr_list$sigma[[m]]) <- rownames(cr_list$sigma[[m]]) <- names(models)
  }
  names(cr_list$commonprob) <- lubridate::month(nmonths, label = TRUE)
  names(cr_list$cr_B) <- lubridate::month(nmonths, label = TRUE)
  names(cr_list$sigma) <- lubridate::month(nmonths, label = TRUE)
  return(cr_list)
}


#' spatialCorrelation object
#'
#' @export
spatialCorrelation <- R6::R6Class("spatialCorrelation",
                                  public = list(
                                    initialize = function(binprobs, mixture_cr){
                                      private$sigma_B = binprobs$sigma
                                      private$cr_B = binprobs$cr_B
                                      private$commonprob = binprobs$commonprob
                                      private$margprob = purrr::map(binprobs$commonprob, diag)
                                      private$cr_X = mixture_cr
                                      private$places = colnames(binprobs$commonprob[[1]])
                                      # Extract index for subsetting cr_X
                                      var_names_X <- rownames(mixture_cr)
                                      idx_x1 <- which(stringr::str_detect(var_names_X , "x1_"))
                                      idx_x2 <- which(stringr::str_detect(var_names_X , "x2_"))
                                      names(idx_x1) <- paste0("ID_", stringr::str_remove_all(var_names_X[idx_x1], "x1_"))
                                      names(idx_x2) <- paste0("ID_", stringr::str_remove_all(var_names_X[idx_x2], "x2_"))
                                      private$idx = list(x1 = idx_x1, x2 = idx_x2)
                                    },
                                    get = function(places, nmonth = 1, date){

                                      if (!missing(date)){
                                        nmonth <- lubridate::month(as.Date(date))
                                      }

                                      sigma_B = private$sigma_B[[nmonth]]
                                      margprob = private$margprob[[nmonth]]
                                      cr_X = private$cr_X

                                      if (!missing(places)){
                                        places = match.arg(places, choices = private$places, several.ok = TRUE)
                                        sigma_B = sigma_B[places, places]
                                        margprob = margprob[places]
                                        idx_x <- c(private$idx$x1[names(private$idx$x1) %in% places],
                                                   private$idx$x2[names(private$idx$x2) %in% places])
                                        idx_x <- idx_x[order(idx_x)]
                                        cr_X = cr_X[idx_x, idx_x]
                                      }
                                      list(
                                        sigma_B = sigma_B,
                                        margprob = margprob,
                                        cr_X = cr_X
                                      )

                                    }
                                  ),
                                  private = list(
                                    places = NA,
                                    sigma_B = NA,
                                    cr_B = NA,
                                    commonprob = NA,
                                    margprob = NA,
                                    cr_X = NA,
                                    idx = NA
                                  ))














