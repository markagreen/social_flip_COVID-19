##############################
### Vaccine vs reinfection ###
#### Descriptives Tables #####
##############################

# libaries
library(table1)

### Simple way using Rmarkdown (copy and paste into spreadsheet) ###

# Table 1

# Negative test (all residents are included here)
table1(~ factor(outcome_delta1) + factor(vax_delta1) + factor(prev_inf_delta1) + Patient.Age + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) | n_tests_delta1, data=data2[data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate),]) # Delta (3rd June to 1st September period)
table1(~ factor(outcome_delta2) + factor(vax_delta2) + factor(prev_inf_delta2) + Patient.Age + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) | n_tests_delta2, data=data2[data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate),]) # Delta (1st September to 27th November period)
table1(~ factor(outcome_omni) + factor(vax_omni) + factor(prev_inf_omni) + Patient.Age + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) | n_tests_omni, data=data2[data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate),]) # Omicron (13th December to 28th Feb period)

# Flu vaccination (all residents are included here)
table1(~ factor(outcome_delta1) + factor(vax_delta1) + factor(prev_inf_delta1) + Patient.Age + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) | flu_vax_delta1, data=data2[data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate),]) # Delta (3rd June to 1st September period)
table1(~ factor(outcome_delta2) + factor(vax_delta2) + factor(prev_inf_delta2) + Patient.Age + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) | flu_vax_delta2, data=data2[data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate),]) # Delta (1st September to 27th November period)
table1(~ factor(outcome_omni) + factor(vax_omni) + factor(prev_inf_omni) + Patient.Age + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) | flu_vax_omni, data=data2[data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate),]) # Omicron (13th December to 28th Feb period)


# Table 1a #

# Delta 1 #

# Negative test (will give all residents)
table1(~ factor(outcome_delta1) | n_tests_delta1*vax_delta1, data=data2[(data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate))]) 
table1(~ factor(outcome_delta1) | n_tests_delta1*prev_inf_delta1, data=data2[(data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate))]) 
table1(~ factor(outcome_delta1) | n_tests_delta1*imd_decile, data=data2[(data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate))]) 

# Flu vax (will give all residents)
table1(~ factor(outcome_delta1) | flu_vax_delta1*vax_delta1, data=data2[(data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate))]) 
table1(~ factor(outcome_delta1) | flu_vax_delta1*prev_inf_delta1, data=data2[(data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate))]) 
table1(~ factor(outcome_delta1) | flu_vax_delta1*imd_decile, data=data2[(data2$DeathDate >= "2021-09-01" | is.na(data2$DeathDate))]) 

# Delta 2 #

# Negative test (will give all residents)
table1(~ factor(outcome_delta2) | n_tests_delta2*vax_delta2, data=data2[(data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate))]) 
table1(~ factor(outcome_delta2) | n_tests_delta2*prev_inf_delta2, data=data2[(data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate))]) 
table1(~ factor(outcome_delta2) | n_tests_delta2*imd_decile, data=data2[(data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate))]) 

# Flu vax (will give all residents)
table1(~ factor(outcome_delta2) | flu_vax_delta2*vax_delta2, data=data2[(data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate))]) 
table1(~ factor(outcome_delta2) | flu_vax_delta2*prev_inf_delta2, data=data2[(data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate))]) 
table1(~ factor(outcome_delta2) | flu_vax_delta2*imd_decile, data=data2[(data2$DeathDate >= "2021-11-27" | is.na(data2$DeathDate))]) 

# Omicron #

# Negative test (will give all residents)
table1(~ factor(outcome_omni) | n_tests_omni*vax_omni, data=data2[(data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate))]) 
table1(~ factor(outcome_omni) | n_tests_omni*prev_inf_omni, data=data2[(data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate))]) 
table1(~ factor(outcome_omni) | n_tests_omni*imd_decile, data=data2[(data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate))]) 

# Flu vax (will give all residents)
table1(~ factor(outcome_omni) | flu_vax_omni*vax_omni, data=data2[(data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate))]) 
table1(~ factor(outcome_omni) | flu_vax_omni*prev_inf_omni, data=data2[(data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate))]) 
table1(~ factor(outcome_omni) | flu_vax_omni*imd_decile, data=data2[(data2$DeathDate >= "2022-03-01" | is.na(data2$DeathDate))]) 


