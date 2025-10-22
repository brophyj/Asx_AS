# AVATAR crossover from CS to AVR simulation
# 53 observed events and censoring at day 2930 (with a few earlier censoring times from Supplement), 
# Gamma (scale, shape) distribution best choice for right skwed data 
# Gamma parameters are shape and scale
# The mean is mean = shape × scale, and the variance is var = shape × scale²
# Given median = 476 and IQR = 1098 - 226, we can estimate the shape and scale parameters
# Assmume moderate skewness, we can use a shape parameter of around 2.5
# Can then solve for scale using the median (700) and the quantile function of the gamma distribution

# libraries required for the simulation and plotting
library(dplyr)
library(ggplot2)
library(survival)

# Step 1: Define gamma parameters based on median and IQR
median_time <- 476
iqr <- 1098 - 226
shape_gamma <- 2.5  # Chosen shape to reflect skew
scale_gamma <- median_time / qgamma(0.5, shape = shape_gamma)

# Step 2: Simulate event times for 53 patients
set.seed(123)
event_times <- rgamma(35, shape = shape_gamma, scale = scale_gamma)

# Step 3: Add 19 censored patients
# Specific censoring times you provided:
specific_censor_times <- c(199, 787, 731, 446, 504, 799, 1341,1419,1174)

# Remaining censored at day 2930
remaining_censored <- rep(2400, 9)  # 53-35-9=9

# All censoring times
censor_times <- c(specific_censor_times, remaining_censored)

# Step 4: Combine into single dataset
dat <- tibble(
  time = c(event_times, censor_times),
  status = c(rep(1, length(event_times)), rep(0, length(censor_times)))  # 1 = event, 0 = censored
)

# Step 5: Count how many had events before 180 days
num_events_before_180 <- dat %>% filter(status == 1, time <= 180) %>% nrow()
cat("Number of events before day 180:", num_events_before_180, "\n")

# Step 6: Histogram of event times only
dat_events_only <- dat %>% filter(status == 1)

p1 <- ggplot(dat_events_only, aes(x = time)) +
  geom_histogram(bins = 20, fill = "skyblue", color = "black", boundary = 0) +
  geom_vline(xintercept = 180, linetype = "dashed", color = "red") +
  labs(
    title = "Histogram of Event Times (Excludes Censored)",
    subtitle = "Red dashed line = 180 days",
    x = "Time (days)",
    y = "Number of Events"
  ) +
  theme_minimal()

# Step 7: Histogram with censoring shown as tick marks
p2 <- ggplot(dat_events_only, aes(x = time)) +
  geom_histogram(bins = 20, fill = "skyblue", color = "black", boundary = 0) +
  geom_vline(xintercept = 180, linetype = "dashed", color = "red") +
  geom_rug(data = dat %>% filter(status == 0), aes(x = time), sides = "b", color = "darkgray") +
  labs(
    title = "Event Times with Censoring Shown",
    subtitle = "Censoring shown as tick marks; Red dashed line = 180 days",
    x = "Time (days)",
    y = "Number of Events"
  ) +
  theme_minimal()

# Step 8: Print plots
print(p1)
print(p2)
