##############################
### Vaccine vs reinfection ###
#### Descriptives Tables #####
##############################

# Table 1 #

# Sample size
nrow(data2[data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate),]) # All residents
nrow(data2[data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate),]) 
nrow(data2[data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate),]) 

nrow(data2[data2$n_tests_delta1 == 1 & (data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate)),]) # Negative test
nrow(data2[data2$n_tests_delta2 == 1 & (data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate)),]) 
nrow(data2[data2$n_tests_omni == 1 & (data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate)),]) 

nrow(data2[data2$flu_vax_delta1 == 1 & (data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate)),]) # Flu vaccine
nrow(data2[data2$flu_vax_delta2 == 1 & (data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate)),]) 
nrow(data2[data2$flu_vax_omni == 1 & (data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate)),]) 

# Outcome variables
table(data2$outcome_delta1[data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate)], exclude = NULL) # All residents
table(data2$outcome_delta2[data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate)], exclude = NULL) 
table(data2$outcome_omni[data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate)], exclude = NULL) 

table(data2$outcome_delta1[data2$n_tests_delta1 == 1 & (data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate)),], exclude = NULL) # Negative test
table(data2$outcome_delta2[[data2$n_tests_delta2 == 1 & (data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate)),], exclude = NULL) 
table(data2$outcome_omni[data2$n_tests_omni == 1 & (data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate)),], exclude = NULL) 

# Vaccination status
table(data2$vax_delta1[data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate)], exclude = NULL) # All residents
table(data2$vax_delta2[data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate)], exclude = NULL)
table(data2$vax_omni[data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate)], exclude = NULL)

# Previous Infection
table(data2$prev_inf_delta1[data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate)], exclude = NULL) # All residents
table(data2$prev_inf_delta2[data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate)], exclude = NULL)
table(data2$prev_inf_omni[data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate)], exclude = NULL)

# Age
mean(data2$Patient.Age[data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate)], na.rm = TRUE) # Mean 
sd(data2$Patient.Age[data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate)]) # Standard deviation
nrow(data2[is.na(data2$Patient.Age) & (data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate))]) # Missing data

mean(data2$Patient.Age[data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate)], na.rm = TRUE) # Mean 
sd(data2$Patient.Age[data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate)]) # Standard deviation
nrow(data2[is.na(data2$Patient.Age) & (data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate))]) # Missing data

mean(data2$Patient.Age[data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate)], na.rm = TRUE) # Mean 
sd(data2$Patient.Age[data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate)]) # Standard deviation
nrow(data2[is.na(data2$Patient.Age) & (data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate))]) # Missing data

# Sex
table(data2$Sex[data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate)], exclude = NULL) # Frequency count
table(data2$Sex[data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate)], exclude = NULL) 
table(data2$Sex[data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate)], exclude = NULL) 

# Ethnicity
table(data2$EthnicMainGroup[data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate)], exclude = NULL) # Frequency count
table(data2$EthnicMainGroup[data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate)], exclude = NULL) 
table(data2$EthnicMainGroup[data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate)], exclude = NULL) 

# Deprivation
# med_vax <- data2[,list(age = median(Patient.Age, na.rm=TRUE), imd = median(imd_score, na.rm=T)), by = "vax_v1"] # If want to stratify
mean(data2$imd_score[data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate)], na.rm = TRUE) # Mean 
sd(data2$imd_score[data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate)], na.rm = TRUE) # Standard deviation
nrow(data2[is.na(data2$imd_score) & (data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate))]) # Missing data

mean(data2$imd_score[data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate)], na.rm = TRUE) # Mean 
sd(data2$imd_score[data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate)], na.rm = TRUE) # Standard deviation
nrow(data2[is.na(data2$imd_score) & (data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate))]) # Missing data

mean(data2$imd_score[data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate)], na.rm = TRUE) # Mean 
sd(data2$imd_score[data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate)], na.rm = TRUE) # Standard deviation
nrow(data2[is.na(data2$imd_score) & (data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate))]) # Missing data

# Health status
table(data2$health_issue[data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate)], exclude = NULL) # Frequency count
table(data2$health_issue[data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate)], exclude = NULL) 
table(data2$health_issue[data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate)], exclude = NULL) 

# Tests - baseline
table(data2$n_tests_delta1[data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate)], exclude = NULL) # Frequency count
table(data2$n_tests_delta2[data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate)], exclude = NULL) 
table(data2$n_tests_omni[data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate)], exclude = NULL) 

# Tests - post baseline
table(data2$n_tests_post_delta1[data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate)], exclude = NULL) # Frequency count
table(data2$n_tests_post_delta2[data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate)], exclude = NULL) 
table(data2$n_tests_post_omni[data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate)], exclude = NULL) 

# Summary stats for negative tests per group
data2[, list(n_tests_post_delta1 = sum(n_tests_post_delta1, na.rm = T), pop = .N), by = c("vax_delta1")] # Vaccination status
data2[, list(n_tests_post_delta2 = sum(n_tests_post_delta2, na.rm = T), pop = .N), by = c("vax_delta2")] 
data2[, list(n_tests_post_omni = sum(n_tests_post_omni, na.rm = T), pop = .N), by = c("vax_omni")] 

data2[, list(n_tests_post_delta1 = sum(n_tests_post_delta1, na.rm = T), pop = .N), by = c("prev_inf_delta1")] # Previous infection
data2[, list(n_tests_post_delta2 = sum(n_tests_post_delta2, na.rm = T), pop = .N), by = c("prev_inf_delta2")] 
data2[, list(n_tests_post_omni = sum(n_tests_post_omni, na.rm = T), pop = .N), by = c("prev_inf_omni")] 

data2[, list(n_tests_post_delta1 = sum(n_tests_post_delta1, na.rm = T), pop = .N), by = c("imd_decile")] # IMD
data2[, list(n_tests_post_delta2 = sum(n_tests_post_delta2, na.rm = T), pop = .N), by = c("imd_decile")] 
data2[, list(n_tests_post_omni = sum(n_tests_post_omni, na.rm = T), pop = .N), by = c("imd_decile")] 

data2[, list(n_tests_post_delta1 = sum(n_tests_post_delta1, na.rm = T), pop = .N), by = c("outcome_delta1")] # Positive test
data2[, list(n_tests_post_delta2 = sum(n_tests_post_delta2, na.rm = T), pop = .N), by = c("outcome_delta2")] 
data2[, list(n_tests_post_omni = sum(n_tests_post_omni, na.rm = T), pop = .N), by = c("outcome_omni")]

# Table 1a #

# Vaccination status
table(data2$vax_delta1[data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate)], data2$outcome_delta1[data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate)], exclude = NULL)
table(data2$vax_delta2[data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate)], data2$outcome_delta2[data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate)], exclude = NULL)
table(data2$vax_omni[data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate)], data2$outcome_omni[data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate)], exclude = NULL)

# Previous infection
table(data2$prev_inf_delta1[data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate)], data2$outcome_delta1[data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate)], exclude = NULL)
table(data2$prev_inf_delta2[data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate)], data2$outcome_delta2[data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate)], exclude = NULL)
table(data2$prev_inf_omni[data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate)], data2$outcome_omni[data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate)], exclude = NULL)

# IMD decile
table(data2$imd_decile[data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate)], data2$outcome_delta1[data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate)], exclude = NULL)
table(data2$imd_decile[data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate)], data2$outcome_delta2[data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate)], exclude = NULL)
table(data2$imd_decile[data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate)], data2$outcome_omni[data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate)], exclude = NULL)
