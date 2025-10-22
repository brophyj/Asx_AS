# This R script performs a Bayesian hierarchical analysis of crossover rates using the `brms` package.
# It includes data preparation, model fitting, posterior draws extraction, and visualization of results.
# It gives identical results to crossover.R which uses Stan directly
library(brms)
library(dplyr)
library(ggplot2)
library(tidyr)
library(glue)
library(patchwork)

study_data <- data.frame(
  study = c("RECOVERY", "EVOLVED", "AVATAR", "EARLY_TAVR"),
  N = c(73, 113, 79, 446),
  cross = c(14, 31, 29, 208)
)

study_data <- study_data %>%
  mutate(study = factor(study, levels = study))  # preserve order

# Use logit link with binomial family
brms_model <- brm(
  formula = cross | trials(N) ~ 1 + (1 | study),
  data = study_data,
  family = binomial(link = "logit"),
  prior = c(
    prior(normal(0, 1.5), class = "Intercept"),
    prior(cauchy(0, 1), class = "sd")  # tau
  ),
  seed = 1234,
  iter = 2000, warmup = 1000, chains = 4, cores = 4,
  control = list(adapt_delta = 0.95),
  refresh = 0
)

draws <- as_draws_df(brms_model)

# Extract posterior estimates per study (log-odds + group-level effect)
study_levels <- levels(study_data$study)
b0 <- draws$b_Intercept
ranef_mat <- as.matrix(draws[, paste0("r_study[", study_levels, ",Intercept]")])
colnames(ranef_mat) <- study_levels

# Get posterior draws of each study’s probability
draws_long <- as.data.frame(ranef_mat) %>%
  mutate(global = b0) %>%
  pivot_longer(-global, names_to = "study", values_to = "ranef") %>%
  mutate(p = plogis(global + ranef))

subset_studies <- c("AVATAR", "EVOLVED", "RECOVERY")
df_avg <- draws_long %>%
  filter(study %in% subset_studies) %>%
  group_by(.draw = row_number() %% (nrow(.) / length(subset_studies))) %>%
  summarise(p = mean(p)) %>%
  mutate(study = "Pooled average\n(AVATAR, EVOLVED,\nRECOVERY)")

df_all <- bind_rows(
  draws_long %>% select(study, p),
  df_avg
) %>%
  mutate(
    study = factor(study, levels = c("RECOVERY", "EVOLVED", "AVATAR", 
                                     "Pooled average\n(AVATAR, EVOLVED,\nRECOVERY)", 
                                     "EARLY_TAVR"))
  )

g_xover <- ggplot(df_all, aes(x = p, color = study, fill = study)) +
  geom_density(alpha = 0.4) +
  labs(
    title = "Individual study estimates with pooled average from Bayesian\nhierarchical model for AVATAR, EVOLVED, and RECOVERY\n(surgical intervention arms) compared to TAVR intervention trial",
    subtitle = "Light Blue Area compared to Light Purple Area",
    x = "Crossover Probability",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 12, hjust = 0)
  )

# Add plot annotation for the figure label
g_xover +
  plot_annotation(title = 'Distributions of 1 Year Crossover Rates') &
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

ggsave("Output/xover.png", g_xover, width = 8, height = 6)

# Get EARLY_TAVR and pooled difference
p_early <- draws_long %>% filter(study == "EARLY_TAVR") %>% pull(p)
p_avg <- df_avg$p
difference <- p_early - p_avg

thresholds <- seq(0, 0.2, by = 0.025)
prob_values <- sapply(thresholds, function(d) mean(difference > d))

df_plot <- data.frame(
  diff_threshold = thresholds,
  probability = prob_values
)

ggplot(df_plot, aes(x = diff_threshold, y = probability)) +
  geom_line() +
  geom_point(size = 2) +
  labs(
    title = "Probabilities of Differences in Crossover Rates",
    x = "Crossover Rate Difference (EARLY_TAVR − Pooled Surgical Trials)",
    y = "Posterior Probability"
  ) +
  theme_minimal()

summary_table <- df_all %>%
  group_by(study) %>%
  summarise(
    mean_p   = mean(p),
    median_p = median(p),
    sd_p     = sd(p),
    lower95  = quantile(p, 0.025),
    upper95  = quantile(p, 0.975),
    .groups = "drop"
  )

print(summary_table)
