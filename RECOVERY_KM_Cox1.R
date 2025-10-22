# Same results as with RECOVERY_censored_landmark1.R
# Set seed for reproducibility

set.seed(456)

library(survival)
library(dplyr)
library(survminer)
library(readr)

#------------------------------------------------------------
# 1. Read in the cumulative incidence data
#------------------------------------------------------------
data_int <- read_csv("Input/RECOVERY_Sx.csv", col_names = c("time", "incidence"),
                     show_col_types = FALSE) # Adds column names time and incidence (%)
data_clinical <- read_csv("Input/RECOVERY_Ctl.csv", col_names = c("time", "incidence"),
                          show_col_types = FALSE) # Adds column names time and incidence (%)


#------------------------------------------------------------
# 1A. Force time=0, incidence=0 row if missing
#------------------------------------------------------------
add_zero_incidence <- function(df) {
  if (!any(df$time == 0)) {
    df <- rbind(data.frame(time = 0, incidence = 0), df)
  }
  df[order(df$time), ]
}

data_int <- add_zero_incidence(data_int)
data_clinical <- add_zero_incidence(data_clinical)

#------------------------------------------------------------
# 2. At-risk numbers & time points
#------------------------------------------------------------
time_points <- c(0, 2, 4, 6, 8)

at_risk_int  <- c(73, 73, 70, 38, 13)      # Intervention
at_risk_clin <- c(72, 68, 65, 36, 12)      # Clinical
N_int  <- at_risk_int[1]      # 78
N_clin <- at_risk_clin[1]     # 79

#------------------------------------------------------------
# 3. IPD reconstruction function
#------------------------------------------------------------
simulate_ipd_with_risk <- function(df, at_risk, time_points, N, group_label) {
  df <- df %>% arrange(time)
  
  ipd_list <- list()
  
  for (i in seq_len(length(time_points) - 1)) {
    t_start <- time_points[i]
    t_end   <- time_points[i + 1]
    
    # Interpolate
    ci_start <- approx(df$time, df$incidence, xout = t_start, rule = 2)$y
    ci_end   <- approx(df$time, df$incidence, xout = t_end,   rule = 2)$y
    
    # Expected cumulative events
    cum_events_start <- N * ci_start / 100
    cum_events_end   <- N * ci_end   / 100
    raw_events <- cum_events_end - cum_events_start
    
    risk_drop <- at_risk[i] - at_risk[i + 1]
    events_interval <- min(round(raw_events), risk_drop)
    censor_interval <- risk_drop - events_interval
    
    # Offset if t_start == 0
    start_for_draws <- if (t_start == 0) 0.1 else t_start
    
    # Simulate events
    if (events_interval > 0) {
      event_times <- runif(events_interval, min = start_for_draws, max = t_end)
      event_df <- data.frame(time = event_times, event = 1)
    } else {
      event_df <- data.frame(time = numeric(0), event = numeric(0))
    }
    
    # Simulate censoring
    if (censor_interval > 0) {
      censor_times <- runif(censor_interval, min = start_for_draws, max = t_end)
      censor_df <- data.frame(time = censor_times, event = 0)
    } else {
      censor_df <- data.frame(time = numeric(0), event = numeric(0))
    }
    
    ipd_list[[i]] <- bind_rows(event_df, censor_df)
  }
  
  # Final time point
  t_last <- time_points[length(time_points)]
  remaining <- at_risk[length(time_points)]
  if (remaining > 0) {
    last_df <- data.frame(time = rep(t_last, remaining), event = 0)
  } else {
    last_df <- data.frame(time = numeric(0), event = numeric(0))
  }
  
  ipd <- bind_rows(ipd_list, last_df)
  ipd$group <- group_label
  ipd
}

#------------------------------------------------------------
# 4. Reconstruct IPD for both groups
#------------------------------------------------------------
ipd_int <- simulate_ipd_with_risk(data_int, at_risk_int, time_points, N_int, "int")
ipd_clin <- simulate_ipd_with_risk(data_clinical, at_risk_clin, time_points, N_clin, "CS")

ipd_all <- bind_rows(ipd_int, ipd_clin) %>% 
  mutate(study = "RECOVERY")
ipd_all$group <- factor(ipd_all$group, levels = c("CS", "int"))

saveRDS(ipd_all, "Output/RECOVERY_km.rds")
#------------------------------------------------------------
# 5. Cox + KM for full follow-up
#------------------------------------------------------------
cox_model <- coxph(Surv(time, event) ~ group, data = ipd_all)
summary(cox_model)

# Plot with improved x-axis formatting
ggsurvplot(
  fit_overall,
  data        = ipd_all,
  risk.table  = TRUE,
  conf.int    = TRUE,
  xlim        = c(0, 8),          # Show 0 to 5 on the x-axis
  break.time.by = 1,              # Tick marks every 1 year
  xlab        = "Time (years)",   # Custom x-axis label
  ylab        = "Survival Probability",
  ggtheme     = theme_minimal(),
  legend.title = "Group"
)
#------------------------------------------------------------
# 6. Restrict analysis to  1 year by censoring after 1 year
# new dataset (ipd_censor1) where survival data are administratively censored at time = 1 (e.g., 1 year)
#------------------------------------------------------------
ipd_censor1 <- ipd_all %>%
  mutate(
    # If time > 1, censor 
    time_1 = pmin(time, 1),
    event_1 = if_else(time <= 1, event, 0)
  )

cox_1 <- coxph(Surv(time_1, event_1) ~ group, data = ipd_censor1)
summary(cox_1)

# Print HR + CI
hr_1 <- exp(coef(cox_1))
ci_1 <- exp(confint(cox_1))
cat("HR at <= 1 year =", hr_1, "\n95% CI =", ci_1[1], "to", ci_1[1], "\n")
# data insufficient to model

# Plot KM up to 1year
fit_1 <- survfit(Surv(time_1, event_1) ~ group, data = ipd_censor1)
ggsurvplot(
  fit_1, data = ipd_censor1,
  xlim = c(0, 1),
  risk.table = TRUE, conf.int = TRUE,
  ggtheme = theme_minimal(), legend.title = "Group"
)

#------------------------------------------------------------
# 7. Landmark analysis starting from 1 year
#------------------------------------------------------------
ipd_landmark1 <- ipd_all %>%
  # Step 1: Keep those who are still at risk at 1 year
  filter(time > 1) %>%
  # Step 2: Shift time so that 1 year is now "time zero"
  mutate(
    time_post1 = time - 1,
    event_post1 = event  # only include post-1y events
  )
cox_landmark1 <- coxph(Surv(time_post1, event_post1) ~ group, data = ipd_landmark1)
summary(cox_landmark1)
# Print HR + CI
hr_landmark1 <- exp(coef(cox_landmark1))
ci_landmark1 <- exp(confint(cox_landmark1))
cat("HR at > 1 year =", hr_landmark1, "\n95% CI =", ci_landmark1[1], "to", ci_landmark1[2], "\n")


# Save the IPD data for further analysis
saveRDS(ipd_landmark, "Output/RECOVERY_km_landmark.rds")
# Save the IPD data for full follow-up analysis
saveRDS(ipd_all, "Output/RECOVERY_km_full.rds")
# Save the IPD data for 1 year analysis
saveRDS(ipd_censor1, "Output/RECOVERY_km_1year.rds")
# Save the IPD data for 1 year censored analysis
saveRDS(ipd_censor, "Output/RECOVERY_km_censored.rds")
# Save the Cox model results
saveRDS(cox_model, "Output/RECOVERY_cox_full.rds")
# Save the Cox model results for 1 year analysis
saveRDS(cox_1, "Output/RECOVERY_cox_1year.rds")
# Save the Cox model results for censored analysis
saveRDS(cox_censor, "Output/RECOVERY_cox_censored.rds")
# Save the Cox model results for landmark analysis
saveRDS(cox_landmark, "Output/RECOVERY_cox_landmark.rds")
#------------------------------------------------------------



