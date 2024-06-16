#' Plot solar model
#'
#' @export

plot.solarModel <- function(object, nplot = 1, plot_year = 2019){

  plot_ut <- function(object, M = 1){

    df <- dplyr::filter(object$data, Month == M)
    pdf_emp <- density(df$ut)
    pdf_mod <- dnorm_mix(params = unlist(object$NM_model[M,2:6]))

    ggplot()+
      geom_line(aes(pdf_emp$x, pdf_emp$y, color = "emp"))+
      geom_line(aes(pdf_emp$x, pdf_mod(pdf_emp$x), color = "nm"), linetype = "dashed")+
      theme_bw()+
      scale_color_manual(values = c(emp = "black", nm = "red"),
                         labels = c(emp = "Empiric", nm = "Normal-Mixture"))+
      theme(legend.position = "none")+
      labs(x = NULL, y = NULL, color = NULL, subtitle = lubridate::month(M, label = T, abbr = F))
  }

  plot_1 <- plot_ut(object, 1)
  plot_2 <- plot_ut(object, 2)
  plot_3 <- plot_ut(object, 3)
  plot_4 <- plot_ut(object, 4)
  plot_5 <- plot_ut(object, 5)
  plot_6 <- plot_ut(object, 6)
  plot_7 <- plot_ut(object, 7)
  plot_8 <- plot_ut(object, 8)
  plot_9 <- plot_ut(object, 9)
  plot_10 <- plot_ut(object, 10)
  plot_11 <- plot_ut(object, 11)
  plot_12 <- plot_ut(object, 12)

  plot_ghi <- dplyr::filter(object$data, Year == plot_year) %>%
    ggplot()+
    geom_line(aes(date, GHI, color = "emp"), alpha = 0.8)+
    geom_line(aes(date, GHI_hat, color = "fit"))+
    geom_line(aes(date, GHI_bar, color = "seasonal"))+
    geom_line(aes(date, clearsky, color = "clearsky"))+
    theme_bw()+
    scale_color_manual(values = c(emp = "black", fit = "red",
                                  seasonal = "blue", clearsky = "orange"),
                       labels = c(emp = "Empiric", fit = "Fitted",
                                  seasonal = "Seasonal", clearsky = "Clearsky"))+
    theme(legend.position = "top")+
    labs(x = NULL, y = "GHI", color = NULL,
         title = object$data$place[1])

  if (nplot == 1){
    gridExtra::grid.arrange(plot_1, plot_2, plot_3, plot_4,
                            plot_5, plot_6, plot_7, plot_8,
                            plot_9, plot_10, plot_11, plot_12)
  } else if (nplot == 2){
    plot_ghi
  }
}

#' Plot solar model simulations
#'
#' @export

plot.solarModelSimulation <- function(data, object, nplot = 1, empiric = TRUE){

  df_sim <- dplyr::bind_rows(data$sim)
  df_sim <- dplyr::group_by(df_sim, Year, n) %>% dplyr::mutate(e_x = mean(GHI))
  df_emp <- data$emp
  df_tot <- dplyr::filter(object$data, Month %in% unique(df_emp$Month))
  df_tot <- dplyr::group_by(df_tot, n) %>% dplyr::mutate(e_x = mean(GHI))

  plot_out <- ggplot()
  if (nplot == 1) {
    plot_out <- plot_out +
      geom_line(data = df_sim, aes(date, GHI, group = seed, color = "sim"), alpha = 0.1)+
      geom_line(data = df_sim, aes(date, GHI_bar, group = seed, color = "seasonal"), alpha = 0.1)+
      geom_line(data = df_sim, aes(date, e_x, group = seed, color = "e_x"), linewidth = 0.2)+
      geom_line(data = df_sim, aes(date, Ct, group = seed, color = "skymax"),  linewidth = 0.1)
    if (empiric) {
      plot_out <- plot_out +
        geom_line(data = df_emp, aes(date, GHI, color = "emp"))+
        geom_line(data = filter(df_tot, date %in% df_emp$date), aes(date, e_x, color = "e_x_r"), linewidth = 0.3)
    }

  } else if (nplot == 2){
    plot_out <- plot_out +
      geom_density(data = df_sim, aes(GHI, color = "sim"))

    if (empiric){
      plot_out <- plot_out +
        geom_density(data = df_tot, aes(GHI, color = "emp"))
    }
  }
  plot_out +
    scale_color_manual(values = c(emp = "black", sim = "red",
                                  e_x = "green", e_x_r = "blue", skymax = "orange", seasonal = "purple"),
                       labels = c(emp = "Empiric",
                                  sim = paste(length(data$sim), "simulations")))+
    theme_bw()+
    theme(legend.position = "top")+
    labs(color = NULL)
}

