# ====================================================
# Complete R Script for Bayesian Hierarchical Modeling of Crossover Rates
# and Plotting Posterior Densities with a Pooled Average for a Subset of Studies
# (Legend ordered as: RECOVERY, EVOLVED, AVATAR, Pooled average (AVATAR, EVOLVED, RECOVERY), EARLY_TAVR)
# ====================================================

# ----- Step 1: Load Required Libraries -----
library(cmdstanr)
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)

# ----- Step 2: Define Study Data -----
# Data for 4 studies: Study, total N, and number of crossovers.
# AVATAR difficult to get 1 year croovoer 25 total crossovers in early publication 
# with median time 400 days so reasonable to assume 1 year crossover of 12/79
# EARLY_TAVR 1 year crossover from Supplemental Figure 5 47% * 231 at risk

study_data <- data.frame(
  study = c("RECOVERY", "EVOLVED", "AVATAR", "EARLY_TAVR"),
  N = c(73, 113, 79, 231),
  cross = c(9, 31, 12, 109)
)
# We want to calculate the pooled average for the subset {AVATAR, EVOLVED, RECOVERY}.

# ----- Step 3: Prepare Data for Stan -----
stan_data <- list(
  S = nrow(study_data),       # number of studies (4)
  N = study_data$N,           # total subjects per study
  cross = study_data$cross    # number of crossovers per study
)

# ----- Step 4: Write the Stan Model Code (using new array syntax) -----
# Hierarchical model with a logit-link: 
# logit(p[s]) = mu + tau * alpha[s], where alpha[s] ~ normal(0,1)
stan_model_code <- "
data {
  int<lower=1> S;                      // number of studies
  array[S] int<lower=0> cross;         // number of crossovers per study
  array[S] int<lower=1> N;             // total subjects per study
}
parameters {
  real mu;                             // global average on logit scale
  real<lower=0> tau;                   // between-study heterogeneity
  vector[S] alpha;                     // random effects on logit scale
}
transformed parameters {
  vector[S] p;
  for (s in 1:S)
    p[s] = inv_logit(mu + tau * alpha[s]); // p[s] is the crossover probability for study s
}
model {
  // Priors
  mu ~ normal(0, 1.5);
  tau ~ cauchy(0, 1);
  alpha ~ normal(0, 1);
  
  // Likelihood
  for (s in 1:S)
    cross[s] ~ binomial(N[s], p[s]);
}
"
stan_file <- "crossover_model.stan"
writeLines(stan_model_code, con = stan_file)

# ----- Step 5: Compile and Fit the Stan Model -----
mod <- cmdstan_model(stan_file)
fit <- mod$sample(
  data = stan_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)
print(fit$summary())

# ----- Step 6: Extract Posterior Draws for p -----
# 'p' is the vector of crossover probabilities for each study.
draws <- fit$draws("p")  # dimensions: [iterations, chains, S]
df_p <- as_draws_df(draws)

# ----- Step 7: Create a Long-Format Data Frame for Density Plotting -----
df_long <- df_p %>%
  pivot_longer(
    cols = starts_with("p["),
    names_to = "study_index",
    values_to = "p"
  ) %>%
  mutate(
    study_index = as.numeric(gsub("p\\[|\\]", "", study_index)),
    study = study_data$study[study_index]
  )

# ----- Step 8: Compute Pooled Average for a Subset of Studies -----
# Subset: AVATAR, EVOLVED, RECOVERY.
subset_studies <- c("AVATAR", "EVOLVED", "RECOVERY")
subset_indices <- which(study_data$study %in% subset_studies)
if(length(subset_indices) > 0) {
  subset_cols <- paste0("p[", subset_indices, "]")
  p_avg <- rowMeans(df_p[, subset_cols, drop = FALSE])
} else {
  p_avg <- NA
}
df_avg <- data.frame(
  study = "Pooled average\n(AVATAR, EVOLVED,\nRECOVERY)",
  p = p_avg
)

# ----- Step 9: Combine with the Long-Format Data Frame -----
df_long2 <- bind_rows(df_long, df_avg)

# ----- Step 10: Reorder Legend (Factor Levels) -----
# We want the legend order to be: RECOVERY, EVOLVED, AVATAR, Pooled average (AVATAR, EVOLVED, RECOVERY), EARLY_TAVR.
desired_levels <- c("RECOVERY", "EVOLVED", "AVATAR", "Pooled average\n(AVATAR, EVOLVED,\nRECOVERY)", "EARLY_TAVR")
df_long2$study <- factor(df_long2$study, levels = desired_levels)

# ----- Step 11: Plot the Posterior Densities -----
g_xover <- ggplot(df_long2, aes(x = p, color = study, fill = study)) +
  geom_density(alpha = 0.4) +
  labs(
    title = 'Figure 1A - Distributions of 1 Year AVR Crossover Rates in CS Arms',
    x = "Crossover Probability",
    y = "Density"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")


ggsave("Output/xover.png")

# calculate probabilities of differences 

# 1) Identify EARLY_TAVR draws (assuming it's p[4])
p_early <- df_p[["p[4]"]]

# 2) Compute difference: EARLY_TAVR - Pooled
difference <- p_early - p_avg

# 3) Compute probabilities
prob_diff_gt_0  <- mean(difference > 0)
prob_diff_gt_01 <- mean(difference > 0.1)
prob_diff_gt_015 <- mean(difference > 0.15)
prob_diff_gt_02 <- mean(difference > 0.2)
prob_diff_gt_025 <- mean(difference > 0.25)
# 4) Print results
cat("P(EARLY_TAVR - Pooled > 0)    =", prob_diff_gt_0,  "\n")
cat("P(EARLY_TAVR - Pooled > 0.1)  =", prob_diff_gt_01, "\n")
cat("P(EARLY_TAVR - Pooled > 0.15)  =", prob_diff_gt_01, "\n")
cat("P(EARLY_TAVR - Pooled > 0.2)  =", prob_diff_gt_02, "\n")
cat("P(EARLY_TAVR - Pooled > 0.25)  =", prob_diff_gt_025, "\n")

# --------------------------------------------------
# 1) Define thresholds and compute probabilities
# --------------------------------------------------
thresholds <- seq(0, 0.3, by = 0.025)

# For each threshold d, compute P(difference > d)
prob_values <- sapply(thresholds, function(d) {
  mean(difference > d)
})

# Create a data frame for plotting
df_plot <- data.frame(
  diff_threshold = thresholds,
  probability = prob_values
)

# --------------------------------------------------
# 2) Plot the Probability vs. Difference
# --------------------------------------------------

g_xover1 <- ggplot(df_plot, aes(x = diff_threshold, y = probability)) +
  geom_line() +
  geom_point(size = 2) +
  labs(
    title = "Figure 1B - Posterior Probabilities for Crossover Rate Differences\nBetween EARLY_TAVR and Pooled Surgical Trials ",
    x = "Crossover Rate Difference (EARLY TAVR − Pooled Surgical Trials)",
    y = "Posterior Probability"
  ) +
  theme_minimal()

# g_xover1 +
#  patchwork::plot_annotation(title = "Figure 1B - Posterior Probabilities of Differences in Crossover Rates") 

ggsave("Output/xover1.png")


# Table output

summary_by_study <- df_long2 %>%
  group_by(study) %>%
  summarise(
    mean_p   = mean(p),
    median_p = median(p),
    sd_p     = sd(p),
    lower95  = quantile(p, 0.025),
    upper95  = quantile(p, 0.975),
    n        = n()
  )

print(summary_by_study)

# Side by side plots
library(ggplot2)
library(cowplot)     # to wrap images as ggplots
library(patchwork)

p1 <- ggdraw() + draw_image("Output/xover.png")
p2 <- ggdraw() + draw_image("Output/xover1.png")

(p1 | p2) + plot_annotation(tag_levels = "A")  # side by side
ggsave("Output/xover_side_by_side.png", width = 12, height = 5, dpi = 300)

