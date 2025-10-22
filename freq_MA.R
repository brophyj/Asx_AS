# frequentist MA

library(meta)
library(dplyr)
library(ggplot2)
library(ggpubr)

# Data for frequentist meta-analysis
# Cakcuated IPD hazard ratios (HR) and their confidence intervals (CI)
meta_data <- data.frame(
   study = c("RECOVERY", "EVOLVED", "AVATAR", "EARLY_TAVR"),
  HR = c(0.05, 0.71, 0.51, 0.56),
  CI_lower = c(0.01, 0.43, 0.27, 0.45),
  CI_upper = c(0.40, 1.18, 0.94, 0.69)
 )

# Log HR and SE
meta_data$logHR <- log(meta_data$HR)
meta_data$selogHR <- (log(meta_data$CI_upper) - log(meta_data$CI_lower)) / (2 * 1.96)

# Meta-analysis
res <- metagen(
  TE = meta_data$logHR,
  seTE = meta_data$selogHR,
  studlab = meta_data$study,
  sm = "HR",
  method.tau = "REML",
  hakn = TRUE,
  prediction = FALSE,
)

# Forest plot with x-axis limited to 2
forest(res,
       prediction = TRUE,
       comb.random = TRUE,
       print.tau2 = TRUE,
       backtransf = TRUE,
       xlab = "Hazard Ratio",
       leftcols = c("studlab"),
       rightcols = c("effect", "ci"),
       print.subgroup.labels = FALSE,
       xlim = c(0.01, 2))  # Set x-axis limits here

########################


# Step 1: Published HR data
published <- data.frame(
  study = c("RECOVERY", "EVOLVED", "AVATAR", "EARLY_TAVR"),
  HR = c(0.09, 0.79, 0.46, 0.50),
  CI_lower = c(0.01, 0.44, 0.23, 0.40),
  CI_upper = c(0.67, 1.43, 0.90, 0.63)
)

# Step 2: Compute log HR and SE
published <- published %>%
  mutate(
    logHR = log(HR),
    selogHR = (log(CI_upper) - log(CI_lower)) / (2 * 1.96)
  )

# Step 3: Run meta-analysis with prediction interval
res_pub <- metagen(
  TE = logHR,
  seTE = selogHR,
  studlab = study,
  sm = "HR",
  method.tau = "REML",
  method.random.ci = "HK",
  prediction = TRUE,
  data = published
)

# Step 3b: Forest plot with meta
forest(res_pub,
       prediction = FALSE,
       comb.random = TRUE,
       backtransf = TRUE,
       xlab = "Hazard Ratio",
       leftcols = c("studlab"),
       rightcols = c("effect", "ci"),
       print.subgroup.labels = FALSE,
       xlim = c(0.01, 2),
       main = "Meta-analysis of Published Hazard Ratios (Random Effects Model)")

# Step 4: Extract pooled and prediction data
pooled <- tibble(
  study = "Pooled",
  HR = exp(res_pub$TE.random),
  ci_lower = exp(res_pub$lower.random),
  ci_upper = exp(res_pub$upper.random)
)

pi <- tibble(
  study = "Prediction",
  HR = exp(res_pub$TE.random),
  ci_lower = exp(res_pub$lower.predict),
  ci_upper = exp(res_pub$upper.predict)
)

# Step 5: Build full plot data frame
published_plot <- published %>%
  mutate(ci_lower = CI_lower, ci_upper = CI_upper) %>%
  select(study, HR, ci_lower, ci_upper)

plot_data <- bind_rows(published_plot, pooled, pi) %>%
  mutate(
    study = factor(study, levels = rev(c(published$study, "Pooled", "Prediction"))),
    label = sprintf("%.2f [%.2f, %.2f]", HR, ci_lower, ci_upper)
  )

library(meta)
library(dplyr)
library(ggplot2)

# Published data
published <- data.frame(
  study = c("RECOVERY", "EVOLVED", "AVATAR", "EARLY_TAVR"),
  HR = c(0.09, 0.79, 0.46, 0.50),
  CI_lower = c(0.01, 0.44, 0.23, 0.40),
  CI_upper = c(0.67, 1.43, 0.90, 0.63)
)

# Log HR and SE
published <- published %>%
  mutate(
    logHR = log(HR),
    selogHR = (log(CI_upper) - log(CI_lower)) / (2 * 1.96)
  )

# Meta-analysis
res_pub <- metagen(
  TE = logHR,
  seTE = selogHR,
  studlab = study,
  sm = "HR",
  method.tau = "REML",
  method.random.ci = "HK",
  prediction = TRUE,
  data = published
)

# Extract pooled and prediction
pooled <- tibble(
  study = "Pooled",
  HR = exp(res_pub$TE.random),
  ci_lower = exp(res_pub$lower.random),
  ci_upper = exp(res_pub$upper.random)
)

pi <- tibble(
  study = "Prediction",
  HR = exp(res_pub$TE.random),
  ci_lower = exp(res_pub$lower.predict),
  ci_upper = exp(res_pub$upper.predict)
)

# Assemble data
published_plot <- published %>%
  mutate(ci_lower = CI_lower, ci_upper = CI_upper) %>%
  select(study, HR, ci_lower, ci_upper)

plot_data <- bind_rows(published_plot, pooled, pi) %>%
  mutate(
    study = factor(study, levels = rev(c(published$study, "Pooled", "Prediction"))),
    label = sprintf("%.2f [%.2f, %.2f]", HR, ci_lower, ci_upper),
    is_summary = study %in% c("Pooled", "Prediction")
  )

# Final ggplot
ggplot(plot_data, aes(x = study, y = HR)) +
  geom_point(aes(shape = is_summary), size = 3) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.1) +
  geom_text(
    aes(label = label),
    hjust = 0,
    nudge_x = 0.3,
    size = 4.5,
    family = "mono"
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
  scale_y_log10(
    limits = c(0.01, 2.5),
    breaks = c(0.01, 0.05, 0.1, 0.2, 0.5, 1, 2)
  ) +
  coord_flip(clip = "off") +
  labs(
    title = "Published Hazard Ratios with Pooled Estimate and Prediction Interval",
    x = NULL,
    y = "Hazard Ratio (log scale)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.margin = margin(10, 180, 10, 10),
    axis.text.y = element_text(hjust = 1)
  )
