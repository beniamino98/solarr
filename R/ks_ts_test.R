#' Two sample Kolmogorov Smirnov test for a time series
#'
#' @param x a vector.
#' @param ci p.value for rejection.
#' @param min_quantile minimum quantile for the grid of values.
#' @param max_quantile maximum quantile for the grid of values.
#' @param seed random seed.
#' @param plot when `TRUE` a plot is returned, otherwise a `tibble`.
#'
#' @name ks_ts_test
#' @rdname ks_ts_test
#' @return when `plot = TRUE` a plot is returned, otherwise a `tibble`.
#' @export
ks_ts_test <- function(x, ci = 0.01, min_quantile = 0.015, max_quantile = 0.985, seed = 1, plot = FALSE){
  set.seed(seed) # random seed

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
      labs(x = TeX("x"), y = "cdf")+
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
