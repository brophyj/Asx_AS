# --- Packages ---
library(brms)
library(tidyverse)
library(tidybayes)
library(ggdist)
library(glue)
library(forcats)

# --- Data (IPD results you posted) ---
dat <- tibble(
  study = c("RECOVERY", "EVOLVED", "AVATAR", "EARLY_TAVR"),
  y = c(0.06, 0.52, 0.54, 1.03),
  lower = c(0.01, 0.27, 0.24, 0.75),
  upper = c(0.46, 1.00, 1.26, 1.41)
)  %>% 
  mutate(
    y_log = log(y),
    sd = (log(upper) - log(lower)) / (2 * 1.96)
  )

# --- Bayesian random-effects meta-analysis (log-HR) ---
bayes_ma <- brm(
  y_log | se(sd) ~ 1 + (1 | study),
  data   = dat,
  family = gaussian(),
  prior  = c(
    prior(normal(0, 1.5), class = "Intercept"),
    prior(normal(0, 0.5), class = "sd")
  ),
  iter = 2000, warmup = 1000, chains = 4, cores = 4,
  refresh = 0, seed = 1234
)

# --- Posterior draws (log scale) ---
out_r_log <- spread_draws(bayes_ma, r_study[study, term], b_Intercept) |>
  filter(term == "Intercept") |>
  transmute(study, effect_log = b_Intercept + r_study)

out_avg_log <- as_draws_df(bayes_ma) |>
  transmute(study = "Average", effect_log = b_Intercept)

out_pred_log <- as_draws_df(bayes_ma) |>
  transmute(
    study      = "Prediction",
    effect_log = rnorm(n(), mean = b_Intercept, sd = sd_study__Intercept)
  )

# --- Combine and move to HR scale ---
study_levels <- c("RECOVERY","EVOLVED","AVATAR","EARLY_TAVR","Average","Prediction")

out_all_hr <- bind_rows(out_r_log, out_avg_log, out_pred_log) |>
  filter(is.finite(effect_log)) |>
  mutate(
    study  = factor(study, levels = study_levels),
    effect = exp(effect_log)
  )

# --- Raw (published) points/intervals (ensure same factor levels) ---
raw_dat <- dat |>
  mutate(study = factor(study, levels = study_levels))

# --- Safe helper to pull HDI bounds regardless of ggdist version ---
get_hdi_bounds <- function(x, width = 0.95) {
  x <- x[is.finite(x)]
  iv <- ggdist::hdi(x, .width = width)
  if (is.numeric(iv)) {
    c(lower = as.numeric(iv[1]), upper = as.numeric(iv[2]))
  } else if (is.data.frame(iv)) {
    # columns may be named ".lower"/".upper" or "lower"/"upper"
    lo <- dplyr::coalesce(iv$.lower %||% iv$lower)
    up <- dplyr::coalesce(iv$.upper %||% iv$upper)
    c(lower = as.numeric(lo[1]), upper = as.numeric(up[1]))
  } else {
    c(lower = NA_real_, upper = NA_real_)
  }
}

# --- Intervals table (single source of truth) ---
hdi_tbl <- out_all_hr |>
  group_by(study) |>
  summarise(
    median = median(effect, na.rm = TRUE),
    # 80% HDI
    h80_l = get_hdi_bounds(effect, 0.80)["lower"],
    h80_u = get_hdi_bounds(effect, 0.80)["upper"],
    # 95% HDI
    h95_l = get_hdi_bounds(effect, 0.95)["lower"],
    h95_u = get_hdi_bounds(effect, 0.95)["upper"],
    .groups = "drop"
  )

# (Optional) print the numbers for the manuscript/table
hdi_tbl |>
  select(study, median, hdi80_lower = h80_l, hdi80_upper = h80_u,
         hdi95_lower = h95_l, hdi95_upper = h95_u) |>
  arrange(match(study, study_levels)) |>
  print(n = Inf)

# Labels for figure (median [95% HDI])
label_tbl <- hdi_tbl |>
  mutate(label = sprintf("%.2f [%.2f, %.2f]", median, h95_l, h95_u))

# Axis limits + label position
x_min  <- quantile(out_all_hr$effect, 0.001, na.rm = TRUE)
x_max  <- quantile(out_all_hr$effect, 0.999, na.rm = TRUE)
x_text <- x_max * 0.55

# --- Figure: densities + our intervals + medians + labels + raw data (nudged down) ---
p <- ggplot() +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3) +
  
  # posterior densities
  stat_slab(
    data = out_all_hr,
    aes(x = effect, y = fct_rev(study)),
    fill = "dodgerblue", color = NA
  ) +
  
  # 95% HDI (thin bar)
  geom_segment(
    data = hdi_tbl,
    aes(y = fct_rev(study), yend = fct_rev(study),
        x = h95_l, xend = h95_u),
    linewidth = 0.8, color = "black"
  ) +
  
  # 80% HDI (thick bar)
  geom_segment(
    data = hdi_tbl,
    aes(y = fct_rev(study), yend = fct_rev(study),
        x = h80_l, xend = h80_u),
    linewidth = 2.3, color = "black"
  ) +
  
  # posterior median (solid dot)
  geom_point(
    data = hdi_tbl,
    aes(x = median, y = fct_rev(study)),
    size = 2.2, color = "black"
  ) +
  
  # raw trial CIs (horizontal) — nudged slightly downward
  geom_errorbarh(
    data = raw_dat,
    aes(y = fct_rev(study), xmin = lower, xmax = upper),
    height = 0.08, color = "black",
    position = position_nudge(y = -0.15)
  ) +
  # raw trial point (open circle) — nudged to match the bar
  geom_point(
    data = raw_dat,
    aes(x = y, y = fct_rev(study)),
    shape = 1, size = 2.8, color = "black", stroke = 0.9,
    position = position_nudge(y = -0.15)
  ) +
  
  # numeric labels (median [95% HDI])
  geom_text(
    data = label_tbl,
    aes(y = fct_rev(study), x = x_text, label = label),
    hjust = 0, size = 3.5
  ) +
  
  scale_x_continuous(
    trans = "log",
    breaks = c(0.05, 0.10, 0.25, 0.50, 1, 2)
  ) +
  coord_cartesian(xlim = c(x_min, x_max)) +
  labs(
    title = "Figure 3- Bayesian (Posterior) Estimates of Primary Outcome from \nOne Year Landmark Analysis of IPD",
    x = "Hazard Ratio (log scale)", y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

print(p)
ggsave("Output/Bayesian_MA_primary_posteriorLandmark.png", width = 10, height = 6, dpi = 300)

