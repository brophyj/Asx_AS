# ====================================================
# Stan program but with any plot
# ====================================================

# ----- Step 1: Write the Random-Slope Stan Model -----
stan_model_code <- "
data {
  int<lower=1> N;
  vector[N] time;
  array[N] int<lower=0,upper=1> event;
  int<lower=1> S;
  array[N] int<lower=1, upper=S> study;
  array[N] int<lower=0,upper=1> group;
}
parameters {
  vector[S] alpha;     // study-specific log baseline hazards
  vector[S] gamma;     // study-specific log hazard ratios
  real mu_alpha;
  real<lower=0> sigma_alpha;
  real mu_gamma;       // global mean log hazard ratio
  real<lower=0> sigma_gamma; // between-study heterogeneity
}
model {
  // Priors for alpha
  mu_alpha ~ normal(0, 5);
  sigma_alpha ~ cauchy(0, 2.5);
  alpha ~ normal(mu_alpha, sigma_alpha);
  
  // Priors for gamma
  mu_gamma ~ normal(0, 5);
  sigma_gamma ~ cauchy(0, 0.25); // sigma_gamma ~ cauchy(0, 2.5); too wide for PI
  gamma ~ normal(mu_gamma, sigma_gamma);
  
  // Exponential survival likelihood
  for (i in 1:N) {
    real lambda = exp(alpha[study[i]] + gamma[study[i]] * group[i]);
    if (event[i] == 1)
      target += log(lambda) - lambda * time[i];
    else
      target += - lambda * time[i];
  }
}
"
stan_file <- "random_slope_meta.stan"
writeLines(stan_model_code, con = stan_file)

# ----- Step 2: Load Required Libraries -----
library(cmdstanr)
library(posterior)
library(dplyr)
library(tidyr)
library(survival)
library(broom)
library(ggplot2)

# ----- Step 3: Load and Prepare Data -----
# Data file "Output/dat_primary_landmark.rds" created with datasets.R 
dat <- readRDS("Output/dat_all_km.rds")  # assumed columns: time, event, group, study

dat <- dat %>%
  mutate(
    group_bin = if_else(group == "Int", 1L, 0L),
    study_id  = as.integer(factor(study))
  )

stan_data <- list(
  N     = nrow(dat),
  time  = dat$time,
  event = dat$event,
  S     = length(unique(dat$study_id)),
  study = dat$study_id,
  group = dat$group_bin
)

# ----- Step 4: Compile and Fit the Stan Model -----
mod <- cmdstan_model(stan_file)
fit <- mod$sample(
  data = stan_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)
print(fit$summary())

# ----- Step 5: Extract Posterior Draws for gamma, mu_gamma, sigma_gamma -----
draws_gamma       <- fit$draws("gamma")        # dimensions: [iterations, chains, S]
draws_mu_gamma    <- fit$draws("mu_gamma")     # dimensions: [iterations, chains]
draws_sigma_gamma <- fit$draws("sigma_gamma")  # dimensions: [iterations, chains]

# Flatten iteration and chain dimensions:
gamma_mat       <- as_draws_matrix(draws_gamma)       # now [iterations*chains, S]
mu_gamma_vec    <- as_draws_matrix(draws_mu_gamma)[, 1]
sigma_gamma_vec <- as_draws_matrix(draws_sigma_gamma)[, 1]

study_labels <- sort(unique(dat$study))
S <- length(study_labels)
