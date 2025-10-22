# to create consolidated data set of reconstructed data
# file /Output/*_km.rds comes from R object ipd in separate *_KM_Cox1.R programs line 120
temp1 <- readRDS("Output/AVATAR_km.rds") %>% 
  mutate(study = "AVATAR") 
temp2 <- readRDS("Output/EARLY_km.rds") %>% 
  mutate(study = "EARLY")
temp3 <- readRDS("Output/EVOLVED_km.rds") %>% 
  mutate(study = "EVOLVED")
temp4 <- readRDS("Output/RECOVERY_km.rds") %>% 
  mutate(study = "RECOVERY")

dat_all_km <- rbind(temp1, temp2, temp3, temp4)
saveRDS(dat_all_km, "Output/dat_all_km.rds")

# now RECOVERY doesn't include hospitalizations which the other trials do,
# admittedly with slightly different definitions

dat_primary <- dat_all_km %>% 
  filter(study != "RECOVERY")
saveRDS(dat_primary, "Output/dat_primary.rds")

dat_primary_landmark <- dat_primary %>% #landmark dataset primary outcome without RECOVERY
  filter(time > 1)
saveRDS(dat_primary_landmark, "Output/dat_primary_landmark.rds")
