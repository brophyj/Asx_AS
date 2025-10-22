
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

