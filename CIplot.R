CIplot <- function(dataset1, dataset2,
                   group1_name, group2_name, study_name,
                   time_unit = "years", max_y = 40,
                   base_size = 12) {
  
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(stringr)
  
  ctl <- read.csv(dataset1, header = FALSE)
  sx  <- read.csv(dataset2, header = FALSE)
  
  names(ctl) <- c("time", "ci")
  names(sx)  <- c("time", "ci")
  ctl$group <- group1_name
  sx$group  <- group2_name
  df_ci <- bind_rows(ctl, sx)
  
  max_time <- max(df_ci$time, na.rm = TRUE)
  
  if (time_unit == "months") {
    # show at most ~6 ticks across the whole range
    x_breaks <- extended_breaks(n = 6)(c(0, max_time))
    x_label  <- "Months since Randomization"
  } else {
    # yearly ticks but keep it readable
    raw_breaks <- seq(0, ceiling(max_time), by = 1)
    # if too many, thin them
    x_breaks <- if (length(raw_breaks) > 8) raw_breaks[raw_breaks %% 2 == 0] else raw_breaks
    x_label  <- "Years since Randomization"
  }
  
  # gentle padding so steps/labels don’t clip at edges
  x_expand <- expansion(mult = c(0.00, 0.02))
  y_expand <- expansion(mult = c(0.00, 0.05))
  
  # wrap long study titles so they don’t collide with the plot area
  title_wrapped <- str_wrap(study_name, width = 30)
  
  # use default palette unless you really need specific colors
  p_ci <- ggplot(df_ci, aes(x = time, y = ci, color = factor(group))) +
    geom_step(linewidth = 0.9, lineend = "round") +
    scale_x_continuous(limits = c(0, max_time), breaks = x_breaks, expand = x_expand) +
    scale_y_continuous(limits = c(0, max_y),
                       breaks = pretty_breaks(n = 5),
                       expand = y_expand) +
    labs(
      title = title_wrapped,
      x = x_label,
      y = "Cumulative Incidence (%)",
      color = NULL
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", margin = margin(b = 6)),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.box.margin = margin(t = 2, b = 2),
      legend.key.height = unit(10, "pt"),
      legend.key.width  = unit(16, "pt"),
      panel.grid.minor = element_blank(),
      axis.title.y = element_text(margin = margin(r = 6)),
      axis.title.x = element_text(margin = margin(t = 6)),
      plot.margin = margin(10, 14, 10, 14)  # a bit more right/left padding
    ) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
    coord_cartesian(clip = "off")
  
  # if you still want a file per-plot, save smaller to avoid Word overflow
  try({
    ggsave(
      filename = file.path("Output", paste0("Cumulative_Incidence_Plot_", gsub(" ", "_", study_name), ".png")),
      plot = p_ci, width = 7, height = 4.5, dpi = 300, units = "in", limitsize = FALSE
    )
  }, silent = TRUE)
  
  p_ci
}
