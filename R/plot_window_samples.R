#### Plot Output ####
# plot sample DLM outputs broken down by window type
plot_window_samples <- function(simulation_path, # path to directory where simulated M values were saved
                                output_path, # path to directory where run_EWDLM output
                                coef_max, # specifies which coef_max should be plotted
                                n = 10) # specifies number of samples
  { 
  true_fits <- read_parquet(paste0(simulation_path, "m_sims.parquet"),
                            col_select = c("site", "window_type", "window_peak", "lag_coefs"))
  
  selected_sites <- true_fits |> 
    dplyr::filter(window_peak == coef_max | window_peak == 0) |> 
    group_by(window_type) |> 
    slice_sample(n = n) |> 
    ungroup()
  
  mat_fits <- open_dataset(output_path, partitioning = c("batch")) |> 
    select(site, mat_coef, mat_ci_lower, mat_ci_upper) |> 
    dplyr::filter(site %in% selected_sites$site) |> 
    collect()
  
  num_lags <- length(mat_fits$mat_coef[[1]])
  
  df_to_plot <- selected_sites |> 
    left_join(select(mat_fits, site, mat_coef), 
              by = "site") |> 
    mutate(lag = list(0:(num_lags - 1))) |> 
    unnest(cols = c("lag_coefs", "mat_coef", "lag"))
  
  est_max <- max(map_dbl(df_to_plot$mat_coef, ~max(as.numeric(.x))))
  y_max <- ceiling(est_max * 1000) / 1000
  
  df_updated <- df_to_plot |> 
    select(site, window_type, lag, mat_coef, lag_coefs) |> 
    mutate(week = lag / 7,
           effect_present = lag_coefs != 0)
  
  shade_key <- tibble(present = c(F, T),
                      color = c(NA, "#D4D4D4"))
  shade_key <- dplyr::filter(shade_key, present %in% df_updated$effect_present)
  
  labels <- c(
    "0" = "Null",
    "1" = "Early",
    "2" = "Middle",
    "3" = "Late",
    "4" = "Full",
    "5" = "Dual")
  
  plot <- df_updated |> 
    ggplot(aes(x = week, y = mat_coef, group = site)) +
    geom_rect(aes(xmin = week - (0.5/7), xmax = week + (0.5/7), fill = effect_present), ymin = -y_max, ymax = y_max, alpha = 0.1) + 
    geom_line(alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_x_continuous(breaks = seq(0, 35, 5)) +
    scale_y_continuous(limits = c(-y_max, y_max), breaks = seq(-y_max, y_max, by = 0.001)) +
    facet_wrap(~ window_type, labeller = labeller(window_type = labels), scales = "free", nrow = 2) +
    theme_classic() +
    scale_fill_manual(values = shade_key$color, na.value = NA, guide = "none") +
    labs(title = "DLM Examples by Window Type", 
         x = "Weeks Since Conception", 
         y = expression(paste("Estimated Change in M Value for a 10 μg/", m^3, " Increase in ", PM[2.5]))) +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 14,
                                    margin = margin(t = 5,
                                                    b = 10)),
          axis.title.x = element_text(size = 12,
                                      margin = margin(t = 10,
                                                      b = 5)),
          axis.title.y = element_text(size = 12))
  
  return(plot)
}