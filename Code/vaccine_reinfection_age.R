##############################
### Vaccine vs reinfection ###
### Age-stratified models ####
##############################

# Note: We only run this for the negative test data since the flu vaccine data are concentrated in older ages only


# i. Delta 1

age_delta1_1 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + prev_inf_delta1 + vax_delta1:tstart + prev_inf_delta1:tstart + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_delta1 + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "0-10")
age_table <- tidy(age_delta1_1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
age_table$model <- "0-10" 
age_table$adjusted <- "adjusted" 

age_delta1_1a <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + vax_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "0-10")
test <- tidy(age_delta1_1a, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "0-10" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta1_1b <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_delta1 + prev_inf_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "0-10")
test <- tidy(age_delta1_1b, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "0-10" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta1_1c <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "0-10")
test <- tidy(age_delta1_1c, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "0-10" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 



age_delta1_2 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + prev_inf_delta1 + vax_delta1:tstart + prev_inf_delta1:tstart + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_delta1 + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "11-17")
test <- tidy(age_delta1_2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "11-17" 
test$adjusted <- "adjusted" 
age_table <- rbind(age_table, test) 

age_delta1_2a <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + vax_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "11-17")
test <- tidy(age_delta1_2a, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "11-17" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta1_2b <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_delta1 + prev_inf_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "11-17")
test <- tidy(age_delta1_2b, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "11-17" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta1_2c <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "11-17")
test <- tidy(age_delta1_2c, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "11-17" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 




age_delta1_3 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + prev_inf_delta1 + vax_delta1:tstart + prev_inf_delta1:tstart + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_delta1 + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "18-29")
test <- tidy(age_delta1_3, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "18-29" 
test$adjusted <- "adjusted" 
age_table <- rbind(age_table, test) 

age_delta1_3a <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + vax_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "18-29")
test <- tidy(age_delta1_3a, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "18-29" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta1_3b <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_delta1 + prev_inf_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "18-29")
test <- tidy(age_delta1_3b, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "18-29" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta1_3c <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "18-29")
test <- tidy(age_delta1_3c, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "18-29" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 



age_delta1_4 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + prev_inf_delta1 + vax_delta1:tstart + prev_inf_delta1:tstart + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_delta1 + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "30-39")
test <- tidy(age_delta1_4, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "30-39" 
test$adjusted <- "adjusted" 
age_table <- rbind(age_table, test) 

age_delta1_4a <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + vax_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "30-39")
test <- tidy(age_delta1_4a, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "30-39" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta1_4b <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_delta1 + prev_inf_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "30-39")
test <- tidy(age_delta1_4b, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "30-39" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta1_4c <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "30-39")
test <- tidy(age_delta1_4c, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "30-39" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 



age_delta1_5 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + prev_inf_delta1 + vax_delta1:tstart + prev_inf_delta1:tstart + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_delta1 + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "40-49")
test <- tidy(age_delta1_5, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "40-49" 
test$adjusted <- "adjusted" 
age_table <- rbind(age_table, test) 

age_delta1_5a <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + vax_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "40-49")
test <- tidy(age_delta1_5a, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "40-49" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta1_5b <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_delta1 + prev_inf_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "40-49")
test <- tidy(age_delta1_5b, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "40-49" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta1_5c <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "40-49")
test <- tidy(age_delta1_5c, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "40-49" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 



age_delta1_6 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + prev_inf_delta1 + vax_delta1:tstart + prev_inf_delta1:tstart + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_delta1 + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "50-59")
test <- tidy(age_delta1_6, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "50-59" 
test$adjusted <- "adjusted" 
age_table <- rbind(age_table, test) 

age_delta1_6a <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + vax_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "50-59")
test <- tidy(age_delta1_6a, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "50-59" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta1_6b <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_delta1 + prev_inf_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "50-59")
test <- tidy(age_delta1_6b, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "50-59" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta1_6c <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "50-59")
test <- tidy(age_delta1_6c, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "50-59" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 


age_delta1_7 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + prev_inf_delta1 + vax_delta1:tstart + prev_inf_delta1:tstart + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_delta1 + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "60-69")
test <- tidy(age_delta1_7, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "60-69" 
test$adjusted <- "adjusted" 
age_table <- rbind(age_table, test) 

age_delta1_7a <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + vax_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "60-69")
test <- tidy(age_delta1_7a, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "60-69" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta1_7b <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_delta1 + prev_inf_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "60-69")
test <- tidy(age_delta1_7b, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "60-69" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta1_7c <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "60-69")
test <- tidy(age_delta1_7c, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "60-69" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 



age_delta1_8 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + prev_inf_delta1 + vax_delta1:tstart + prev_inf_delta1:tstart + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_delta1 + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "70+")
test <- tidy(age_delta1_8, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "70+" 
test$adjusted <- "adjusted" 
age_table <- rbind(age_table, test) 

age_delta1_8a <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + vax_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "70+")
test <- tidy(age_delta1_8a, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "70+" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta1_8b <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_delta1 + prev_inf_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "70+")
test <- tidy(age_delta1_8b, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "70+" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta1_8c <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),], subset = age_band == "70+")
test <- tidy(age_delta1_8c, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "70+" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 


write.csv(age_table, "./Vaccine v reinfection paper/coxph_age_delta1.csv") # Save
rm(age_table)



# ii. Delta 2

age_delta2_1 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + prev_inf_delta2 + vax_delta2:tstart + prev_inf_delta2:tstart + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_delta2 + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "0-10")
age_table <- tidy(age_delta2_1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
age_table$model <- "0-10" 
age_table$adjusted <- "adjusted" 

age_delta2_1a <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + vax_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "0-10")
test <- tidy(age_delta2_1a, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "0-10" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta2_1b <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_delta2 + prev_inf_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "0-10")
test <- tidy(age_delta2_1b, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "0-10" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta2_1c <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "0-10")
test <- tidy(age_delta2_1c, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "0-10" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 



age_delta2_2 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + prev_inf_delta2 + vax_delta2:tstart + prev_inf_delta2:tstart + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_delta2 + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "11-17")
test <- tidy(age_delta2_2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "11-17" 
test$adjusted <- "adjusted" 
age_table <- rbind(age_table, test) 

age_delta2_2a <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + vax_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "11-17")
test <- tidy(age_delta2_2a, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "11-17" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta2_2b <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_delta2 + prev_inf_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "11-17")
test <- tidy(age_delta2_2b, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "11-17" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta2_2c <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "11-17")
test <- tidy(age_delta2_2c, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "11-17" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 




age_delta2_3 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + prev_inf_delta2 + vax_delta2:tstart + prev_inf_delta2:tstart + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_delta2 + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "18-29")
test <- tidy(age_delta2_3, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "18-29" 
test$adjusted <- "adjusted" 
age_table <- rbind(age_table, test) 

age_delta2_3a <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + vax_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "18-29")
test <- tidy(age_delta2_3a, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "18-29" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta2_3b <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_delta2 + prev_inf_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "18-29")
test <- tidy(age_delta2_3b, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "18-29" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta2_3c <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "18-29")
test <- tidy(age_delta2_3c, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "18-29" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 



age_delta2_4 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + prev_inf_delta2 + vax_delta2:tstart + prev_inf_delta2:tstart + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_delta2 + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "30-39")
test <- tidy(age_delta2_4, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "30-39" 
test$adjusted <- "adjusted" 
age_table <- rbind(age_table, test) 

age_delta2_4a <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + vax_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "30-39")
test <- tidy(age_delta2_4a, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "30-39" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta2_4b <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_delta2 + prev_inf_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "30-39")
test <- tidy(age_delta2_4b, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "30-39" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta2_4c <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "30-39")
test <- tidy(age_delta2_4c, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "30-39" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 



age_delta2_5 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + prev_inf_delta2 + vax_delta2:tstart + prev_inf_delta2:tstart + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_delta2 + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "40-49")
test <- tidy(age_delta2_5, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "40-49" 
test$adjusted <- "adjusted" 
age_table <- rbind(age_table, test) 

age_delta2_5a <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + vax_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "40-49")
test <- tidy(age_delta2_5a, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "40-49" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta2_5b <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_delta2 + prev_inf_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "40-49")
test <- tidy(age_delta2_5b, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "40-49" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta2_5c <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "40-49")
test <- tidy(age_delta2_5c, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "40-49" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 



age_delta2_6 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + prev_inf_delta2 + vax_delta2:tstart + prev_inf_delta2:tstart + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_delta2 + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "50-59")
test <- tidy(age_delta2_6, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "50-59" 
test$adjusted <- "adjusted" 
age_table <- rbind(age_table, test) 

age_delta2_6a <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + vax_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "50-59")
test <- tidy(age_delta2_6a, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "50-59" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta2_6b <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_delta2 + prev_inf_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "50-59")
test <- tidy(age_delta2_6b, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "50-59" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta2_6c <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "50-59")
test <- tidy(age_delta2_6c, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "50-59" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 


age_delta2_7 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + prev_inf_delta2 + vax_delta2:tstart + prev_inf_delta2:tstart + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_delta2 + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "60-69")
test <- tidy(age_delta2_7, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "60-69" 
test$adjusted <- "adjusted" 
age_table <- rbind(age_table, test) 

age_delta2_7a <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + vax_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "60-69")
test <- tidy(age_delta2_7a, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "60-69" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta2_7b <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_delta2 + prev_inf_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "60-69")
test <- tidy(age_delta2_7b, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "60-69" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta2_7c <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "60-69")
test <- tidy(age_delta2_7c, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "60-69" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 



age_delta2_8 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + prev_inf_delta2 + vax_delta2:tstart + prev_inf_delta2:tstart + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_delta2 + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "70+")
test <- tidy(age_delta2_8, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "70+" 
test$adjusted <- "adjusted" 
age_table <- rbind(age_table, test) 

age_delta2_8a <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + vax_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "70+")
test <- tidy(age_delta2_8a, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "70+" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta2_8b <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_delta2 + prev_inf_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "70+")
test <- tidy(age_delta2_8b, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "70+" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_delta2_8c <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),], subset = age_band == "70+")
test <- tidy(age_delta2_8c, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "70+" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

write.csv(age_table, "./Vaccine v reinfection paper/coxph_age_delta2.csv") # Save
rm(age_table)



# iii. Omicron

age_omni_1 <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + prev_inf_omni + vax_omni:tstart + prev_inf_omni:tstart + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_omni + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "0-10")
age_table <- tidy(age_omni_1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
age_table$model <- "0-10" 
age_table$adjusted <- "adjusted" 

age_omni_1a <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + vax_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "0-10")
test <- tidy(age_omni_1a, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "0-10" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_omni_1b <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_omni + prev_inf_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "0-10")
test <- tidy(age_omni_1b, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "0-10" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_omni_1c <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "0-10")
test <- tidy(age_omni_1c, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "0-10" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 



age_omni_2 <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + prev_inf_omni + vax_omni:tstart + prev_inf_omni:tstart + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_omni + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "11-17")
test <- tidy(age_omni_2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "11-17" 
test$adjusted <- "adjusted" 
age_table <- rbind(age_table, test) 

age_omni_2a <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + vax_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "11-17")
test <- tidy(age_omni_2a, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "11-17" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_omni_2b <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_omni + prev_inf_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "11-17")
test <- tidy(age_omni_2b, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "11-17" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_omni_2c <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "11-17")
test <- tidy(age_omni_2c, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "11-17" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 




age_omni_3 <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + prev_inf_omni + vax_omni:tstart + prev_inf_omni:tstart + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_omni + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "18-29")
test <- tidy(age_omni_3, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "18-29" 
test$adjusted <- "adjusted" 
age_table <- rbind(age_table, test) 

age_omni_3a <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + vax_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "18-29")
test <- tidy(age_omni_3a, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "18-29" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_omni_3b <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_omni + prev_inf_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "18-29")
test <- tidy(age_omni_3b, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "18-29" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_omni_3c <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "18-29")
test <- tidy(age_omni_3c, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "18-29" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 



age_omni_4 <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + prev_inf_omni + vax_omni:tstart + prev_inf_omni:tstart + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_omni + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "30-39")
test <- tidy(age_omni_4, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "30-39" 
test$adjusted <- "adjusted" 
age_table <- rbind(age_table, test) 

age_omni_4a <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + vax_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "30-39")
test <- tidy(age_omni_4a, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "30-39" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_omni_4b <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_omni + prev_inf_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "30-39")
test <- tidy(age_omni_4b, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "30-39" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_omni_4c <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "30-39")
test <- tidy(age_omni_4c, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "30-39" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 



age_omni_5 <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + prev_inf_omni + vax_omni:tstart + prev_inf_omni:tstart + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_omni + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "40-49")
test <- tidy(age_omni_5, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "40-49" 
test$adjusted <- "adjusted" 
age_table <- rbind(age_table, test) 

age_omni_5a <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + vax_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "40-49")
test <- tidy(age_omni_5a, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "40-49" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_omni_5b <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_omni + prev_inf_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "40-49")
test <- tidy(age_omni_5b, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "40-49" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_omni_5c <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "40-49")
test <- tidy(age_omni_5c, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "40-49" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 



age_omni_6 <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + prev_inf_omni + vax_omni:tstart + prev_inf_omni:tstart + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_omni + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "50-59")
test <- tidy(age_omni_6, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "50-59" 
test$adjusted <- "adjusted" 
age_table <- rbind(age_table, test) 

age_omni_6a <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + vax_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "50-59")
test <- tidy(age_omni_6a, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "50-59" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_omni_6b <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_omni + prev_inf_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "50-59")
test <- tidy(age_omni_6b, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "50-59" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_omni_6c <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "50-59")
test <- tidy(age_omni_6c, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "50-59" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 


age_omni_7 <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + prev_inf_omni + vax_omni:tstart + prev_inf_omni:tstart + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_omni + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "60-69")
test <- tidy(age_omni_7, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "60-69" 
test$adjusted <- "adjusted" 
age_table <- rbind(age_table, test) 

age_omni_7a <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + vax_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "60-69")
test <- tidy(age_omni_7a, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "60-69" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_omni_7b <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_omni + prev_inf_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "60-69")
test <- tidy(age_omni_7b, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "60-69" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_omni_7c <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "60-69")
test <- tidy(age_omni_7c, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "60-69" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 



age_omni_8 <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + prev_inf_omni + vax_omni:tstart + prev_inf_omni:tstart + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_omni + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "70+")
test <- tidy(age_omni_8, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "70+" 
test$adjusted <- "adjusted" 
age_table <- rbind(age_table, test) 

age_omni_8a <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + vax_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "70+")
test <- tidy(age_omni_8a, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "70+" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_omni_8b <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_omni + prev_inf_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "70+")
test <- tidy(age_omni_8b, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "70+" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

age_omni_8c <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),], subset = age_band == "70+")
test <- tidy(age_omni_8c, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "70+" 
test$adjusted <- "unadjusted" 
age_table <- rbind(age_table, test) 

write.csv(age_table, "./Vaccine v reinfection paper/coxph_age_omni.csv") # Save
rm(age_table)



### Age plots ###

# Load data and create required measures
delta1_table <- read.csv("./Vaccine v reinfection paper/coxph_age_delta1.csv") # Load
delta1_table$time <- "Delta (6th June - 1st Sept)" # Define time period of analysis
delta1_table$term[delta1_table$term == "prev_inf_delta1"] <- "prev_inf" # Rename term for consistency across models
delta1_table$term[delta1_table$term == "vax_delta11"] <- "1 dose"
delta1_table$term[delta1_table$term == "vax_delta12"] <- "2 doses"
delta1_table$term[delta1_table$term == "vax_delta13"] <- "3 doses"

delta2_table <- read.csv("./Vaccine v reinfection paper/coxph_age_delta2.csv") # Repeat
delta2_table$time <- "Delta (1st Sept - 27th Nov)"
delta2_table$term[delta2_table$term == "prev_inf_delta2"] <- "prev_inf"
delta2_table$term[delta2_table$term == "vax_delta21"] <- "1 dose"
delta2_table$term[delta2_table$term == "vax_delta22"] <- "2 doses"
delta2_table$term[delta2_table$term == "vax_delta23"] <- "3 doses"

omni_table <- read.csv("./Vaccine v reinfection paper/coxph_age_omni.csv") #  Repeat
omni_table$time <- "Omicron (13th Dec - 1st Mar)"
omni_table$term[omni_table$term == "prev_inf_omni"] <- "prev_inf"
omni_table$term[omni_table$term == "vax_omni1"] <- "1 dose"
omni_table$term[omni_table$term == "vax_omni2"] <- "2 doses"
omni_table$term[omni_table$term == "vax_omni3"] <- "3 doses"


## IMD ##

# a. Unadjusted model #

# Subset required data
table_imd <- delta1_table[delta1_table$term == "imd_score" & delta1_table$adjusted == "unadjusted",]
table_imd <- rbind(table_imd, delta2_table[delta2_table$term == "imd_score" & delta2_table$adjusted == "unadjusted",])
table_imd <- rbind(table_imd, omni_table[omni_table$term == "imd_score" & omni_table$adjusted == "unadjusted",])

# Plot
plot_imd_unadj <- ggplot(table_imd, aes(x = model, y = estimate, ymin = conf.low, ymax = conf.high, group = time, color = time)) + 
  geom_point(position=position_dodge(0.4), size = 2) + 
  geom_errorbar(position=position_dodge(0.4), width = 0, size = 0.5) +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_color_viridis_d(begin = 0.1, end = 0.9) +
  labs(
    x = "Age group",
    y = "Hazards Ratio (95% Confidence Intervals)",
    colour = "Period",
  ) +
  coord_cartesian(ylim=c(0.985,1.015)) +
  theme(legend.position="bottom")
plot_imd_unadj

# b. Adjusted model #

# Subset required data
table_imd <- delta1_table[delta1_table$term == "imd_score" & delta1_table$adjusted == "adjusted",]
table_imd <- rbind(table_imd, delta2_table[delta2_table$term == "imd_score" & delta2_table$adjusted == "adjusted",])
table_imd <- rbind(table_imd, omni_table[omni_table$term == "imd_score" & omni_table$adjusted == "adjusted",])

# Plot
plot_imd_adj <- ggplot(table_imd, aes(x = model, y = estimate, ymin = conf.low, ymax = conf.high, group = time, color = time)) + 
  geom_point(position=position_dodge(0.4), size = 2) + 
  geom_errorbar(position=position_dodge(0.4), width = 0, size = 0.5) +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_color_viridis_d(begin = 0.1, end = 0.9) +
  labs(
    x = "Age group",
    y = "Hazards Ratio (95% Confidence Intervals)",
    colour = "Period",
  ) +
  coord_cartesian(ylim=c(0.985,1.015)) +
  theme(legend.position="bottom")
plot_imd_adj



## Previous infection ##


# a. Unadjusted model #

# Subset required data
table_prev <- delta1_table[delta1_table$term == "prev_inf" & delta1_table$adjusted == "unadjusted",]
table_prev <- rbind(table_prev, delta2_table[delta2_table$term == "prev_inf" & delta2_table$adjusted == "unadjusted",])
table_prev <- rbind(table_prev, omni_table[omni_table$term == "prev_inf" & omni_table$adjusted == "unadjusted",])

# Plot
plot_prev_unadj <- ggplot(table_prev, aes(x = model, y = estimate, ymin = conf.low, ymax = conf.high, group = time, color = time)) + 
  geom_point(position=position_dodge(0.4), size = 2) + 
  geom_errorbar(position=position_dodge(0.4), width = 0, size = 0.5) +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_color_viridis_d(begin = 0.1, end = 0.9) +
  labs(
    x = "Age group",
    y = "Hazards Ratio (95% Confidence Intervals)",
    colour = "Period"
  ) +
  coord_cartesian(ylim=c(0,1.01)) +
  theme(legend.position="bottom")
plot_prev_unadj


# b. Adjusted model #

# Subset required data
table_prev <- delta1_table[delta1_table$term == "prev_inf" & delta1_table$adjusted == "adjusted",]
table_prev <- rbind(table_prev, delta2_table[delta2_table$term == "prev_inf" & delta2_table$adjusted == "adjusted",])
table_prev <- rbind(table_prev, omni_table[omni_table$term == "prev_inf" & omni_table$adjusted == "adjusted",])

# Plot
plot_prev_adj <- ggplot(table_prev, aes(x = model, y = estimate, ymin = conf.low, ymax = conf.high, group = time, color = time)) + 
  geom_point(position=position_dodge(0.4), size = 2) + 
  geom_errorbar(position=position_dodge(0.4), width = 0, size = 0.5) +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_color_viridis_d(begin = 0.1, end = 0.9) +
  labs(
    x = "Age group",
    y = "Hazards Ratio (95% Confidence Intervals)",
    colour = "Period"
  ) +
  coord_cartesian(ylim=c(0,1.01)) +
  theme(legend.position="bottom")
plot_prev_adj


## Vaccine status ##

# a. Unadjusted

# Subset required data
table_vax <- delta1_table[(delta1_table$term == "1 dose" | delta1_table$term == "2 doses" | delta1_table$term == "3 doses") & delta1_table$adjusted == "unadjusted",]
table_vax <- rbind(table_vax, delta2_table[(delta2_table$term == "1 dose" | delta2_table$term == "2 doses" | delta2_table$term == "3 doses") & delta2_table$adjusted == "unadjusted",])
table_vax <- rbind(table_vax, omni_table[(omni_table$term == "1 dose" | omni_table$term == "2 doses" | omni_table$term == "3 doses") & omni_table$adjusted == "unadjusted",])

# Remove ages 0-10 as poorly fitting model
table_vax$model <- as.factor(table_vax$model)
table_vax <- table_vax[table_vax$model != "0-10",]

# Plot
plot_vax_unadj <- ggplot(table_vax[table_vax$model != "0-10",], aes(x = model, y = estimate, ymin = conf.low, ymax = conf.high, group = interaction(term,time), colour = time, shape = term)) + 
  geom_point(position=position_dodge(0.6), size = 2) + 
  geom_errorbar(position=position_dodge(0.6), width = 0, size = 0.5) +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_color_viridis_d(begin = 0.1, end = 0.9) +
  labs(
    x = "Age group",
    y = "Hazards Ratio (95% Confidence Intervals)",
    colour = "Period",
    shape = "Vaccination status"
  ) +
  scale_x_discrete(drop=FALSE) + # Plot 0-10 as missing
  coord_cartesian(ylim=c(0.1,3)) +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin()) # Split legend across two rows
plot_vax_unadj


# b. Adjusted

# Subset required data
table_vax <- delta1_table[(delta1_table$term == "1 dose" | delta1_table$term == "2 doses" | delta1_table$term == "3 doses") & delta1_table$adjusted == "adjusted",]
table_vax <- rbind(table_vax, delta2_table[(delta2_table$term == "1 dose" | delta2_table$term == "2 doses" | delta2_table$term == "3 doses") & delta2_table$adjusted == "adjusted",])
table_vax <- rbind(table_vax, omni_table[(omni_table$term == "1 dose" | omni_table$term == "2 doses" | omni_table$term == "3 doses") & omni_table$adjusted == "adjusted",])

# Remove ages 0-10 as poorly fitting model
table_vax$model <- as.factor(table_vax$model)
table_vax <- table_vax[table_vax$model != "0-10",]

# Plot
plot_vax_adj <- ggplot(table_vax[table_vax$model != "0-10",], aes(x = model, y = estimate, ymin = conf.low, ymax = conf.high, group = interaction(term,time), colour = time, shape = term)) + 
  geom_point(position=position_dodge(0.6), size = 2) + 
  geom_errorbar(position=position_dodge(0.6), width = 0, size = 0.5) +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_color_viridis_d(begin = 0.1, end = 0.9) +
  labs(
    x = "Age group",
    y = "Hazards Ratio (95% Confidence Intervals)",
    colour = "Period",
    shape = "Vaccination status"
  ) +
  scale_x_discrete(drop=FALSE) + # Plot 0-10 as missing
  coord_cartesian(ylim=c(0.1,3)) +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin()) # Split legend across two rows
plot_vax_adj


## Combine into a single plots ##

# Libraries
library(patchwork)

# Create blank plot with just y-axis label
p4 <- ggplot(data.frame(l = plot_imd_adj$labels$y, x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), angle = 90) + 
  theme_void() +
  coord_cartesian(clip = "off")

# Set all axis labels as blank
plot_vax_adj$labels$y <- plot_prev_adj$labels$y <- plot_imd_adj$labels$y <- ""
plot_prev_adj$labels$x <- plot_imd_adj$labels$x <- "" # Do same for two plots, but leave bottom one for label
plot_vax_unadj$labels$y <- plot_prev_unadj$labels$y <- plot_imd_unadj$labels$y <- ""
plot_prev_unadj$labels$x <- plot_imd_unadj$labels$x <- "" # Do same for two plots, but leave bottom one for label

# Combine into single plot - undjusted models
plot_age_unadj <- p4 + (plot_imd_unadj / 
  plot_prev_unadj /
  plot_vax_unadj) +
  plot_layout(widths = c(1, 25), guides = 'collect') & theme(legend.position="bottom", legend.box="vertical", legend.margin=margin()) 
plot_age_unadj

ggsave(plot = plot_age_unadj, filename = "./Vaccine v reinfection paper/age_stratified_plots_unadj.jpeg") # Save

# Combine into single plot - adjusted models
plot_age_adj <- p4 + (plot_imd_adj / 
                       plot_prev_adj /
                       plot_vax_adj) +
  plot_layout(widths = c(1, 25), guides = 'collect') & theme(legend.position="bottom", legend.box="vertical", legend.margin=margin()) 
plot_age_adj

ggsave(plot = plot_age_adj, filename = "./Vaccine v reinfection paper/age_stratified_plots_adj.jpeg") # Save





