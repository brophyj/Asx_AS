library(brms)
library(tidyverse)
library(tidybayes)
library(ggdist)
library(glue)

# Create dataframe from Published HR column
dat <- data.frame(
  study = c("RECOVERY", "EVOLVED", "AVATAR", "EARLY_TAVR"),
  y = c(0.09, 0.79, 0.46, 0.50),
  lower = c(0.01, 0.44, 0.23, 0.40),
  upper = c(0.67, 1.43, 0.90, 0.63)
)

# Compute standard deviation from 95% CI (log scale)
dat <- dat %>%
  mutate(
    sd = (log(upper) - log(lower)) / (2 * 1.96),
    y_log = log(y)
  )

# Bayesian random-effects meta-analysis
bayes_ma <- brm(
  y_log | se(sd) ~ 1 + (1 | study),
  data = dat,
  family = gaussian(),
  iter = 2000, warmup = 1000, chains = 4, cores = 4, refresh = 0
)

# Extract posterior draws: study-specific
out_r <- spread_draws(bayes_ma, r_study[study, term], b_Intercept) %>%
  mutate(effect = b_Intercept + r_study)

# Average effect
out_f <- spread_draws(bayes_ma, b_Intercept) %>%
  mutate(study = "Average", effect = b_Intercept)

# Combine and control display order
out_all <- bind_rows(out_r, out_f) %>%
  ungroup() %>%
  mutate(study = fct_relevel(study, setdiff(unique(study), "Average"), after = Inf))

# Summary on HR scale (exp)
out_all_sum <- out_all %>%
  group_by(study) %>%
  mean_qi(effect) %>%
  mutate(
    label = glue("{round(exp(effect), 2)} [{round(exp(.lower), 2)}, {round(exp(.upper), 2)}]")
  )

# Raw data (for open points)
dat_clean <- dat %>%
  mutate(study = fct_relevel(study, setdiff(study, "Average"), after = Inf))

# Plotting
ggplot(out_all, aes(x = exp(effect), y = study)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3) +
  stat_halfeye(.width = c(0.80, 0.95), fill = "dodgerblue", color = "black") +
  geom_text(
    data = out_all_sum,
    aes(label = label, x = 0.05),
    hjust = 0,
    size = 3.5
  ) +
  geom_point(
    data = dat_clean,
    aes(x = y, y = study),
    shape = 1,
    position = position_nudge(y = -0.2)
  ) +
  scale_x_continuous(
    transform = "log",
    breaks = c(0.05, 0.1, 0.25, 0.5, 1, 2),
    limits = c(0.04, 2)
  ) +
  labs(
    title = "Posterior Estimates with Observed Hazard Ratios",
    x = "Hazard Ratio (log scale)",
    y = NULL
  ) +
  theme_minimal(base_size = 13)
