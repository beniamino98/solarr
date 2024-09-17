#' Kolmogorov distribution and quantile functions
#'
#' @param x vector of quantiles.
#' @param p vector of probabilities.
#' @param k finite value for approximation of infinite sum.
#' @return A probability, a numeric vector in 0, 1.
#' @rdname pks
#' @aliases qks
#' @keywords internal
pks <- function(x, k = 100){
  # Infinite sum function
  # k: approximation for infinity
  inf_sum <- function(x, k = 100){
    isum <- 0
    for(i in 1:k){
      isum <- isum + 2*(-1)^(i-1)*exp(-2*(x^2)*(i^2))
    }
    isum
  }
  # Compute probability
  p <- c()
  for(i in 1:length(x)){
    p[i] <- 1-inf_sum(x[i], k = k)
  }
  p
}

#' @rdname pks
#' @keywords internal
qks <- function(p, k = 100){
  loss <- function(x, p){
    (pks(x, k = k) - p)^2
  }
  x <- c()
  for(i in 1:length(p)){
    x[i] <- optim(par = 2, method = "Brent", lower = 0, upper = 10, fn = loss, p = p[i])$par
  }
  x
}

#' Kolmogorov Smirnov test for a distribution
#'
#' Test against a specific distribution with `ks_test` and
#' perform a two sample invariance test for a time series with `ks_ts_test`
#'
#' @param x a vector.
#' @param pdf a function. The theoric density to use for comparison.
#' @param ci p.value for rejection.
#' @param min_quantile minimum quantile for the grid of values.
#' @param max_quantile maximum quantile for the grid of values.
#' @param k finite value for approximation of infinite sum.
#' @param plot when `TRUE` a plot is returned, otherwise a `tibble`.
#' @param seed random seed for two sample test.
#'
#' @name ks_test
#' @rdname ks_test
#'
#' @return when `plot = TRUE` a plot is returned, otherwise a `tibble`.
#' @aliases ks_test
#' @aliases ks_ts_test
#' @export
ks_test <- function(x, cdf, ci = 0.05, min_quantile = 0.015, max_quantile = 0.985, k = 1000, plot = FALSE){

  # number of observations
  n <- length(x)
  # Set the interval values for upper and lower band
  grid <- seq(quantile(x, min_quantile), quantile(x, max_quantile), 0.01)
  # Empirical cdf
  cdf_x <- ecdf(x)
  # Compute the KS-statistic
  ks_stat <- sqrt(n)*max(abs(cdf_x(grid) - cdf(grid)))
  # Compute the rejection level
  rejection_lev <- qks(1-ci, k = k)

  # ========================== Plot ==========================
  if (plot) {
    y_breaks <- seq(0, 1, 0.2)
    y_labels <- paste0(format(y_breaks*100, digits = 2), "%")
    grid_max <- grid[which.max(abs(cdf_x(grid) - cdf(grid)))]
    plt <- ggplot()+
      geom_ribbon(aes(grid, ymax = cdf_x(grid), ymin = cdf(grid)),
                  alpha = 0.5, fill = "green") +
      geom_line(aes(grid, cdf_x(grid)))+
      geom_line(aes(grid, cdf(grid)), color = "red")+
      geom_segment(aes(x = grid_max, xend = grid_max,
                       y = cdf_x(grid_max), yend = cdf(grid_max)),
                   linetype = "solid", color = "magenta")+
      geom_point(aes(grid_max, cdf_x(grid_max)), color = "magenta")+
      geom_point(aes(grid_max, cdf(grid_max)), color = "magenta")+
      scale_y_continuous(breaks = y_breaks, labels = y_labels)+
      labs(x = "x", y = "cdf")+
      theme_bw()
    return(plt)
  } else {
    kab <- dplyr::tibble(
      n = n,
      alpha = paste0(format(ci*100, digits = 3), "%"),
      KS = ks_stat,
      rejection_lev = rejection_lev,
      H0 = ifelse(KS > rejection_lev, "Rejected", "Non-Rejected"))  %>%
      dplyr::mutate_if(is.numeric, format, digits = 4, scientific = FALSE)
    return(kab)
  }
}


# Two sample Kolmogorov Smirnov test for a time series
#' @rdname ks_test
#' @export
ks_ts_test <- function(x, ci = 0.05, min_quantile = 0.015, max_quantile = 0.985, seed = 1, plot = FALSE){
  # Random seed
  set.seed(seed)
  # number of observations
  n <- length(x)
  # Random split of the time series
  idx_split <- sample(n, 1)
  x1 <- x[1:idx_split]
  x2 <- x[(idx_split+1):n]
  # Number of elements for each sub sample
  n1 <- length(x1)
  n2 <- length(x2)
  # Grid of values for KS-statistic
  grid <- seq(quantile(x, min_quantile), quantile(x, max_quantile), 0.01)
  # Empiric cdfs
  cdf_1 <- ecdf(x1)
  cdf_2 <- ecdf(x2)
  # KS-statistic
  ks_stat <- max(abs(cdf_1(grid) - cdf_2(grid)))
  # Rejection level
  rejection_lev <- sqrt(-0.5*log(ci/2))*sqrt((n1+n2)/(n1*n2))

  # ========================== Plot ==========================
  if (plot) {
    y_breaks <- seq(0, 1, 0.2)
    y_labels <- paste0(format(y_breaks*100, digits = 2), "%")
    grid_max <- grid[which.max(abs(cdf_1(grid) - cdf_2(grid)))]
    plt <- ggplot()+
      geom_ribbon(aes(grid, ymax = cdf_1(grid), ymin = cdf_2(grid)),
                  alpha = 0.5, fill = "green") +
      geom_line(aes(grid, cdf_1(grid)))+
      geom_line(aes(grid, cdf_2(grid)), color = "red")+
      geom_segment(aes(x = grid_max, xend = grid_max,
                       y = cdf_1(grid_max), yend = cdf_2(grid_max)),
                   linetype = "solid", color = "magenta")+
      geom_point(aes(grid_max, cdf_1(grid_max)), color = "magenta")+
      geom_point(aes(grid_max, cdf_2(grid_max)), color = "magenta")+
      scale_y_continuous(breaks = y_breaks, labels = y_labels)+
      labs(x = "x", y = "cdf")+
      theme_bw()
    return(plt)
  } else {
    kab <- dplyr::tibble(
      idx_split = idx_split,
      ci = paste0(ci*100, "%"),
      n1 = n1,
      n2 = n2,
      KS = ks_stat,
      rejection_lev = rejection_lev,
      H0 = ifelse(KS > rejection_lev, "Rejected", "Non-Rejected")
    )
    return(kab)
  }
}
