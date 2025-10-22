CIplot_landmark <- function(dataset1, dataset2, group1_name, group2_name, study_name, time_unit = "years", max_y = 60) {
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  
  # Read in the CSV files
  library(readr)
  ctl <- read_csv(dataset1, col_names = c("time", "ci"),
                            show_col_types = FALSE) %>% 
    mutate(across(everything(), as.numeric)) %>% 
    na.omit()
  
  sx <- read_csv(dataset2, col_names = c("time", "ci"),
                 show_col_types = FALSE) %>% 
    mutate(across(everything(), as.numeric)) %>% 
    na.omit()
  
  
  # Set group names
  ctl$group <- group1_name
  sx$group  <- group2_name
  
  # Combine datasets
  df_ci <- bind_rows(ctl, sx)
  
  # Filter out data before 1 year if time unit is years, adjust accordingly if months
  start_time <- if (time_unit == "months") 12 else 1
  df_ci <- df_ci %>% filter(time >= start_time)
  
  # Determine the maximum time for dynamic x-axis scaling
  max_time <- max(df_ci$time)
  
  # Adjust x-axis limits and breaks based on the time unit
  if (time_unit == "months") {
    x_breaks <- seq(start_time, max_time, by = 12)  # yearly intervals in months
    x_label <- "Months since Randomization"
  } else {  # default to years
    x_breaks <- seq(start_time, max_time, by = 1)
    x_label <- "Years since Randomization"
  }
  
  # Define custom colors
  custom_colors <- c("blue", "red")
  names(custom_colors) <- c(group1_name, group2_name)
  
  # Create the cumulative incidence plot
  p_ci <- ggplot(df_ci, aes(x = time, y = ci, color = group)) +
    geom_step(linewidth = 1) +
    scale_x_continuous(limits = c(start_time, max_time), breaks = x_breaks) +
    scale_y_continuous(limits = c(0, max_y), breaks = seq(0, max_y, max_y / 4)) +
    scale_color_manual(values = custom_colors) +
    labs(
      # title = paste("Cumulative Incidence -", study_name),
      title = paste(study_name,  ": 1 Yr Landmark analysis"), #study_name
      x = x_label,
      y = "Cumulative Incidence (%)",
      color = "Treatment"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.title     = element_blank(),
      legend.position  = c(0.7, 0.2),
      panel.grid.minor = element_blank(),
      plot.margin      = margin(1, 1, 1, 1, "cm")
    ) +
    coord_cartesian(clip = "off") 
  
  # Save and return the plot
  ggsave(paste0("Output/Cumulative_Incidence_Plot_", gsub(" ", "_", study_name), ".png"), p_ci, width = 11, height = 8, dpi = 300)
  
  return(p_ci)
}
