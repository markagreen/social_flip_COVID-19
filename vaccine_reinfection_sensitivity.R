##############################
### Vaccine vs reinfection ###
#### Sensivitiy analyses #####
##############################


# Note: run vaccine_reinfection_survival.R first


# 1. Run analyses for all residents - vaccination status #

# i. Delta period 1

# Unadjusted
model_i_a1 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + vax_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate),]) # Only for people who survived to end of study period
table <- tidy(model_i_a1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
table$model <- "model_i_a1" # Describe model type

model_i_a2 <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_delta1 + prev_inf_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate),]) # Repeat
test <- tidy(model_i_a2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model_i_a2" # Describe model type
table <- rbind(table, test) # Join onto main table

model_i_a3 <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate),]) # Repeat
test <- tidy(model_i_a3, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model_i_a3" 
table <- rbind(table, test) 

model_i_a4 <- coxph(Surv(tstart, tstop, infection) ~ factor(imd_decile) + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate),]) # Repeat
test <- tidy(model_i_a4, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model_i_a4" 
table <- rbind(table, test) 


# Fully adjusted
model_i_b1 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + prev_inf_delta1 + vax_delta1:tstart + prev_inf_delta1:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) +  n_tests_delta1 + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate),])
test <- tidy(model_i_b1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model_i_b1" 
table <- rbind(table, test) 

model_i_b2 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + prev_inf_delta1 + vax_delta1:tstart + prev_inf_delta1:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + factor(imd_decile) + factor(health_issue) + n_tests_delta1 + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate),])
test <- tidy(model_i_b2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model_i_b2" 
table <- rbind(table, test) 

write.csv(table, "./Vaccine v reinfection paper/coxph_delta1.csv") # Save
rm(table)

# ii. Delta period 2

# Unadjusted
model_ii_a1 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + vax_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate),]) # Only for people who survived to end of study period
table2 <- tidy(model_ii_a1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
table2$model <- "model_ii_a1" 

model_ii_a2 <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_delta2 + prev_inf_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate),]) 
test <- tidy(model_ii_a2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model_ii_a2" 
table2 <- rbind(table2, test) 

model_ii_a3 <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate),]) 
test <- tidy(model_ii_a3, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model_ii_a3" 
table2 <- rbind(table2, test) 

model_ii_a4 <- coxph(Surv(tstart, tstop, infection) ~ factor(imd_decile) + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate),]) 
test <- tidy(model_ii_a4, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model_ii_a4" 
table2 <- rbind(table2, test) 

# Fully adjusted
model_ii_b1 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + prev_inf_delta2 + vax_delta2:tstart + prev_inf_delta2:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_delta2 + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate),])
test <- tidy(model_ii_b1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model_ii_b1" 
table2 <- rbind(table2, test) 

model_ii_b2 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + prev_inf_delta2 + vax_delta2:tstart + prev_inf_delta2:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + factor(imd_decile) + factor(health_issue) + n_tests_delta2 + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate),])
test <- tidy(model_ii_b2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model_ii_b2" 
table2 <- rbind(table2, test) 

write.csv(table2, "./Vaccine v reinfection paper/coxph_delta2.csv") # Save
rm(table2)

# iii. Omnicron

# Unadjusted
model_iii_a1 <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + vax_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate),]) # Only for people who survived to end of study period
table3 <- tidy(model_iii_a1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
table3$model <- "model_iii_a1" 

model_iii_a2 <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_omni + prev_inf_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate),])
test <- tidy(model_iii_a2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model_iii_a2" 
table3 <- rbind(table3, test)

model_iii_a3 <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate),])
test <- tidy(model_iii_a3, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model_iii_a3" 
table3 <- rbind(table3, test)

model_iii_a4 <- coxph(Surv(tstart, tstop, infection) ~ factor(imd_decile) + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate),])
test <- tidy(model_iii_a4, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model_iii_a4" 
table3 <- rbind(table3, test)

# Fully adjusted
model_iii_b1 <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + prev_inf_omni + vax_omni:tstart + prev_inf_omni:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_omni + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate),])
test <- tidy(model_iii_b1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model_iii_b1" 
table3 <- rbind(table3, test)

model_iii_b2 <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + prev_inf_omni + vax_omni:tstart + prev_inf_omni:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + factor(imd_decile) + factor(health_issue) + n_tests_omni + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate),])
test <- tidy(model_iii_b2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model_iii_b2" 
table3 <- rbind(table3, test)

write.csv(table3, "./Vaccine v reinfection paper/coxph_omni.csv") # Save
rm(table3)
gc()


# 2. Checking for immortal time bias in vaccination status model #

# i. Delta period 1

# Unadjusted
sa_model_i_u1 <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_delta1 + prev_inf_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & survival_delta1$immortal_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),]) # Only for people who survived to end of study period
sa_table <- tidy(sa_model_i_u1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
sa_table$model <- "sa_model_i_u_negtest" # Describe model type

sa_model_i_u2 <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_delta1 + prev_inf_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$flu_vax_delta1 == 1 & survival_delta1$immortal_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),]) # Repeat
test <- tidy(sa_model_i_u2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "sa_model_i_u_flu" # Describe model type
sa_table <- rbind(sa_table, test) # Join onto main table


# Fully adjusted
sa_model_i_a1 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + prev_inf_delta1 + vax_delta1:tstart + prev_inf_delta1:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) +  n_tests_delta1 + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & survival_delta1$immortal_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),])
test <- tidy(sa_model_i_a1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "sa_model_i_a_negtest" 
sa_table <- rbind(sa_table, test) 

sa_model_i_a2 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + prev_inf_delta1 + vax_delta1:tstart + prev_inf_delta1:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + factor(imd_decile) + factor(health_issue) + n_tests_delta1 + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$flu_vax_delta1 == 1 & survival_delta1$immortal_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),])
test <- tidy(sa_model_i_a2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model sa_model_i_a2
test$model <- "sa_model_i_a_flu" 
sa_table <- rbind(sa_table, test) 

write.csv(sa_table, "./Vaccine v reinfection paper/coxph_immortal_bias_delta1.csv") # Save
rm(sa_table)


# ii. Delta period 2

# Unadjusted
sa_model_ii_u1 <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_delta2 + prev_inf_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & survival_delta2$immortal_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),]) # Only for people who survived to end of study period
sa_table2 <- tidy(sa_model_ii_u1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
sa_table2$model <- "sa_model_ii_u_negtest" # Describe model type

sa_model_ii_u2 <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_delta2 + prev_inf_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$flu_vax_delta2 == 1 & survival_delta2$immortal_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),]) # Repeat
test <- tidy(sa_model_ii_u2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "sa_model_ii_u_flu" # Describe model type
sa_table2 <- rbind(sa_table2, test) # Join onto main table


# Fully adjusted
sa_model_ii_a1 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + prev_inf_delta2 + vax_delta2:tstart + prev_inf_delta2:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) +  n_tests_delta2 + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & survival_delta2$immortal_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),])
test <- tidy(sa_model_ii_a1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "sa_model_ii_a_negtest" 
sa_table2 <- rbind(sa_table2, test) 

sa_model_ii_a2 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + prev_inf_delta2 + vax_delta2:tstart + prev_inf_delta2:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + factor(imd_decile) + factor(health_issue) + n_tests_delta2 + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$flu_vax_delta2 == 1 & survival_delta2$immortal_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),])
test <- tidy(sa_model_ii_a2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model 
test$model <- "sa_model_ii_a_flu" 
sa_table2 <- rbind(sa_table2, test) 

write.csv(sa_table2, "./Vaccine v reinfection paper/coxph_immortal_bias_delta2.csv") # Save
rm(sa_table2)


# iii. Omnicron

# Unadjusted
sa_model_iii_u1 <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_omni + prev_inf_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & survival_omni$immortal_omni == 1 & (survival_omni$DeathDate >= "2021-11-27" | is.na(survival_omni$DeathDate)),]) # Only for people who survived to end of study period
sa_table3 <- tidy(sa_model_iii_u1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
sa_table3$model <- "sa_model_iii_u_negtest" # Describe model type

sa_model_iii_u2 <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_omni + prev_inf_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$flu_vax_omni == 1 & survival_omni$immortal_omni == 1 & (survival_omni$DeathDate >= "2021-11-27" | is.na(survival_omni$DeathDate)),]) # Repeat
test <- tidy(sa_model_iii_u2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "sa_model_iii_u_flu" # Describe model type
sa_table3 <- rbind(sa_table3, test) # Join onto main table


# Fully adjusted
sa_model_iii_a1 <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + prev_inf_omni + vax_omni:tstart + prev_inf_omni:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) +  n_tests_omni + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & survival_omni$immortal_omni == 1 & (survival_omni$DeathDate >= "2021-11-27" | is.na(survival_omni$DeathDate)),])
test <- tidy(sa_model_iii_a1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "sa_model_iii_a_negtest" 
sa_table3 <- rbind(sa_table3, test) 

sa_model_iii_a2 <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + prev_inf_omni + vax_omni:tstart + prev_inf_omni:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + factor(imd_decile) + factor(health_issue) + n_tests_omni + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$flu_vax_omni == 1 & survival_omni$immortal_omni == 1 & (survival_omni$DeathDate >= "2021-11-27" | is.na(survival_omni$DeathDate)),])
test <- tidy(sa_model_iii_a2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model 
test$model <- "sa_model_iii_a_flu" 
sa_table3 <- rbind(sa_table3, test) 

write.csv(sa_table3, "./Vaccine v reinfection paper/coxph_immortal_bias_omni.csv") # Save
rm(sa_table3)


## 3. Accounting for vaccination type - all residents ##


# i. Delta period 1

# Unadjusted
sa2_model_i_u <- coxph(Surv(tstart, tstop, infection) ~ vax_type_delta1 + vax_type_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate),]) # Only for people who survived to end of study period
sa_table42 <- tidy(sa2_model_i_u, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
sa_table42$model <- "sa2_model_i_u" # Describe model type

# Fully adjusted
sa2_model_i_a <- coxph(Surv(tstart, tstop, infection) ~ vax_type_delta1 + prev_inf_delta1 + vax_type_delta1:tstart + prev_inf_delta1:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) +  n_tests_delta1 + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate),])
test <- tidy(sa2_model_i_a, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "sa2_model_i_a" 
sa_table42 <- rbind(sa_table42, test) 

# ii. Delta period 2

# Unadjusted
sa2_model_ii_u <- coxph(Surv(tstart, tstop, infection) ~ vax_type_delta2 + vax_type_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate),]) # Only for people who survived to end of study period
test <- tidy(sa2_model_ii_u, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "sa2_model_ii_u"
sa_table42 <- rbind(sa_table42, test) 

# Fully adjusted
sa2_model_ii_a <- coxph(Surv(tstart, tstop, infection) ~ vax_type_delta2 + prev_inf_delta2 + vax_type_delta2:tstart + prev_inf_delta2:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_delta2 + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate),])
test <- tidy(sa2_model_ii_a, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "sa2_model_ii_a" 
sa_table42 <- rbind(sa_table42, test)


# iii. Omnicron

# Unadjusted
sa2_model_iii_u <- coxph(Surv(tstart, tstop, infection) ~ vax_type_omni + vax_type_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$DeathDate >= "2021-03-01" | is.na(survival_omni$DeathDate),]) # Only for people who survived to end of study period
test <- tidy(sa2_model_iii_u, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "sa2_model_iii_u" 
sa_table42 <- rbind(sa_table42, test)

# Fully adjusted
sa2_model_iii_a <- coxph(Surv(tstart, tstop, infection) ~ vax_type_omni + prev_inf_omni + vax_type_omni:tstart + prev_inf_omni:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_omni + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$DeathDate >= "2021-03-01" | is.na(survival_omni$DeathDate),])
test <- tidy(sa2_model_iii_a, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "sa2_model_iii_a" 
sa_table42 <- rbind(sa_table42, test)

write.csv(sa_table42, "./Vaccine v reinfection paper/coxph_vax_type.csv") # Save
rm(sa_table42)
gc()




## 4. Matching methods ##

# # Model to predict vax status
# library(nnet)
# 
# 
# 
# # Delta 1
# load("./Data/survival_delta1.RData") # Load data
# for_analysis <- survival_delta1[(survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),] # Subset data
# model_predict <- multinom(vax_delta1 ~ factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + factor(health_issue) + factor(n_tests_delta1) + prev_inf_delta1 + imd_score,
#                           data = for_analysis) # Model
# ps_matrix <- predict(model_predict, type = 'probs') # generate propensity score
# 
# # Get IPW
# w <- rep(0, nrow(for_analysis)) #inititalize weights
# group <- for_analysis$vax_delta1 # Save values for individuals
# for (i in levels(group)) { # Generate IPW values (ATE)
#   w[group == i] <- 1 / ps_matrix[, i]
# }
# for_analysis <- cbind(for_analysis, w) # Join onto main dataset
# 
# # Model
# model4_i <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + prev_inf_delta1 + vax_delta1:tstart + prev_inf_delta1:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + cluster(FK_Patient_Link_ID), data = for_analysis, weights = w)
# sa_table4 <- tidy(model4_i, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
# sa_table4$model <- "delta1" 
# rm(survival_delta1, model_predict, for_analysis) # Tidy
# 
# 
# 
# # Delta 2
# load("./Data/survival_delta2.RData") # Load data
# for_analysis <- survival_delta2[(survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),] # Subset data
# model_predict <- multinom(vax_delta2 ~ factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + factor(health_issue) + factor(n_tests_delta2) + prev_inf_delta2 + imd_score,
#                           data = for_analysis) # Model
# ps_matrix <- predict(model_predict, type = 'probs') # generate propensity score
# 
# # Get IPW
# w <- rep(0, nrow(for_analysis)) #inititalize weights
# group <- for_analysis$vax_delta2 # Save values for individuals
# for (i in levels(group)) { # Generate IPW values
#   w[group == i] <- 1/ps_matrix[,i]
# }
# for_analysis <- cbind(for_analysis, w) # Join onto main dataset
# 
# # Model
# model4_ii <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + prev_inf_delta2 + vax_delta2:tstart + prev_inf_delta2:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + cluster(FK_Patient_Link_ID), data = for_analysis, weights = w)
# test <- tidy(model4_ii, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
# test$model <- "delta2" 
# sa_table4 <- rbind(sa_table4, test)
# rm(survival_delta2, model_predict, for_analysis) # Tidy
# 
# 
# 
# # Omicron
# load("./Data/survival_omni.RData") # Load data
# for_analysis <- survival_omni[(survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),] # Subset data
# model_predict <- multinom(vax_omni ~ factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + factor(health_issue) + factor(n_tests_omni) + prev_inf_omni + imd_score,
#                          data = for_analysis) # Model
# ps_matrix <- predict(model_predict, type = 'probs') # generate propensity score
# 
# # Get IPW
# w <- rep(0, nrow(for_analysis)) #inititalize weights
# group <- for_analysis$vax_omni # Save values for individuals
# for (i in levels(group)) { # Generate IPW values
#   w[group == i] <- 1/ps_matrix[,i]
# }
# for_analysis <- cbind(for_analysis, w) # Join onto main dataset
# 
# # Model
# model4_iii <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + prev_inf_omni + vax_omni:tstart + prev_inf_omni:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + cluster(FK_Patient_Link_ID), data = for_analysis, weights = w)
# test <- tidy(model4_iii, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
# test$model <- "omicron" 
# sa_table4 <- rbind(sa_table4, test)
# rm(survival_omni, model_predict, for_analysis) # Tidy
# 
# write.csv(sa_table4, "./Vaccine v reinfection paper/coxph_matching.csv") # Save
# rm(sa_table4)
# gc()



# 5. Restrict to people who also received the flu vaccine and 65+ #

# i. Delta period 1

# Unadjusted
sa_model_i_u1 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + vax_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$flu_vax_delta1 == 1 & survival_delta1$Patient.Age >= 65 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),]) # Only for people who survived to end of study period
sa_table5 <- tidy(sa_model_i_u1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
sa_table5$model <- "sa_model_i_u1" # Describe model type

sa_model_i_u2 <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_delta1 + prev_inf_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$flu_vax_delta1 == 1 & survival_delta1$Patient.Age >= 65 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),]) # Repeat
test <- tidy(sa_model_i_u2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "sa_model_i_u2" # Describe model type
sa_table5 <- rbind(sa_table5, test) # Join onto main table

sa_model_i_u3 <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$flu_vax_delta1 == 1 & survival_delta1$Patient.Age >= 65 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),]) # Repeat
test <- tidy(sa_model_i_u3, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "sa_model_i_u3" 
sa_table5 <- rbind(sa_table5, test) 

# Fully adjusted
sa_model_i_a1 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + prev_inf_delta1 + vax_delta1:tstart + prev_inf_delta1:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) +  n_tests_delta1 + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$flu_vax_delta1 == 1 & survival_delta1$Patient.Age >= 65 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),])
test <- tidy(sa_model_i_a1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "sa_model_i_a1" 
sa_table5 <- rbind(sa_table5, test) 

write.csv(sa_table5, "./Vaccine v reinfection paper/coxph_flu_vax_65_delta1.csv") # Save
rm(sa_table5)

# ii. Delta period 2

# Unadjusted
model_ii_a1 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + vax_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$flu_vax_delta2 == 1 & survival_delta2$Patient.Age >= 65 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),]) # Only for people who survived to end of study period
sa_table5 <- tidy(model_ii_a1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
sa_table5$model <- "model_ii_a1" 

model_ii_a2 <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_delta2 + prev_inf_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$flu_vax_delta2 == 1 & survival_delta2$Patient.Age >= 65 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),]) 
test <- tidy(model_ii_a2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model_ii_a2" 
sa_table5 <- rbind(sa_table5, test) 

model_ii_a3 <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$flu_vax_delta2 == 1 & survival_delta2$Patient.Age >= 65 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),]) 
test <- tidy(model_ii_a3, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model_ii_a3" 
sa_table5 <- rbind(sa_table5, test) 

# Fully adjusted
model_ii_b1 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + prev_inf_delta2 + vax_delta2:tstart + prev_inf_delta2:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_delta2 + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$flu_vax_delta2 == 1 & survival_delta2$Patient.Age >= 65 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),])
test <- tidy(model_ii_b1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model_ii_b1" 
sa_table5 <- rbind(sa_table5, test) 

write.csv(sa_table5, "./Vaccine v reinfection paper/coxph_flu_vax_65_delta2.csv") # Save
rm(sa_table5)

# iii. Omnicron

# Unadjusted
model_iii_a1 <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + vax_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$flu_vax_omni == 1 & survival_delta2$Patient.Age >= 65 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),]) # Only for people who survived to end of study period
sa_table5 <- tidy(model_iii_a1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
sa_table5$model <- "model_iii_a1" 

model_iii_a2 <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_omni + prev_inf_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$flu_vax_omni == 1 & survival_delta2$Patient.Age >= 65 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),])
test <- tidy(model_iii_a2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model_iii_a2" 
sa_table5 <- rbind(sa_table5, test)

model_iii_a3 <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$flu_vax_omni == 1 & survival_delta2$Patient.Age >= 65 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),])
test <- tidy(model_iii_a3, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model_iii_a3" 
sa_table5 <- rbind(sa_table5, test)


# Fully adjusted
model_iii_b1 <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + prev_inf_omni + vax_omni:tstart + prev_inf_omni:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_omni + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$flu_vax_omni == 1 & survival_delta2$Patient.Age >= 65 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),])
test <- tidy(model_iii_b1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model_iii_b1" 
sa_table5 <- rbind(sa_table5, test)

write.csv(sa_table5, "./Vaccine v reinfection paper/coxph_flu_vax_omni.csv") # Save
rm(sa_table5)



## 6. Accounting for vaccination type - stratified by negative test or flu vaccine ##


# i. Delta period 1

# Unadjusted
sa2_model_i_u1 <- coxph(Surv(tstart, tstop, infection) ~ vax_type_delta1 + vax_type_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),]) # Model for negative tests
sa_table6 <- tidy(sa2_model_i_u1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
sa_table6$model <- "delta1_unajd_negtest" # Describe model type

sa2_model_i_u2 <- coxph(Surv(tstart, tstop, infection) ~ vax_type_delta1 + vax_type_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_omni$flu_vax_omni == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),]) # Model for flu vaccine
test <- tidy(sa2_model_i_u2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "delta1_unajd_flu" # Describe model type
sa_table6 <- rbind(sa_table6, test) 

# Fully adjusted
sa2_model_i_a1 <- coxph(Surv(tstart, tstop, infection) ~ vax_type_delta1 + prev_inf_delta1 + vax_type_delta1:tstart + prev_inf_delta1:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) +  n_tests_delta1 + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),])
test <- tidy(sa2_model_i_a1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "delta1_ajd_negtest" 
sa_table6 <- rbind(sa_table6, test) 

sa2_model_i_a2 <- coxph(Surv(tstart, tstop, infection) ~ vax_type_delta1 + prev_inf_delta1 + vax_type_delta1:tstart + prev_inf_delta1:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) +  n_tests_delta1 + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_omni$flu_vax_omni == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),])
test <- tidy(sa2_model_i_a2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "delta1_ajd_flu" 
sa_table6 <- rbind(sa_table6, test) 

# ii. Delta period 2

# Unadjusted
sa2_model_i_u1 <- coxph(Surv(tstart, tstop, infection) ~ vax_type_delta2 + vax_type_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),]) # Model for negative tests
test <- tidy(sa2_model_i_u1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "delta2_unajd_negtest" # Describe model type
sa_table6 <- rbind(sa_table6, test) 

sa2_model_i_u2 <- coxph(Surv(tstart, tstop, infection) ~ vax_type_delta2 + vax_type_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_omni$flu_vax_omni == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),]) # Model for flu vaccine
test <- tidy(sa2_model_i_u2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "delta2_unajd_flu" # Describe model type
sa_table6 <- rbind(sa_table6, test) 

# Fully adjusted
sa2_model_i_a1 <- coxph(Surv(tstart, tstop, infection) ~ vax_type_delta2 + prev_inf_delta2 + vax_type_delta2:tstart + prev_inf_delta2:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) +  n_tests_delta2 + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),])
test <- tidy(sa2_model_i_a1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "delta2_ajd_negtest" 
sa_table6 <- rbind(sa_table6, test) 

sa2_model_i_a2 <- coxph(Surv(tstart, tstop, infection) ~ vax_type_delta2 + prev_inf_delta2 + vax_type_delta2:tstart + prev_inf_delta2:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) +  n_tests_delta2 + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_omni$flu_vax_omni == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),])
test <- tidy(sa2_model_i_a2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "delta2_ajd_flu" 
sa_table6 <- rbind(sa_table6, test) 


# iii. Omnicron

# Unadjusted
sa2_model_i_u1 <- coxph(Surv(tstart, tstop, infection) ~ vax_type_omni + vax_type_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),]) # Model for negative tests
test <- tidy(sa2_model_i_u1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "omni_unajd_negtest" # Describe model type
sa_table6 <- rbind(sa_table6, test) 

sa2_model_i_u2 <- coxph(Surv(tstart, tstop, infection) ~ vax_type_omni + vax_type_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$flu_vax_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),]) # Model for flu vaccine
test <- tidy(sa2_model_i_u2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "omni_unajd_flu" # Describe model type
sa_table6 <- rbind(sa_table6, test) 

# Fully adjusted
sa2_model_i_a1 <- coxph(Surv(tstart, tstop, infection) ~ vax_type_omni + prev_inf_omni + vax_type_omni:tstart + prev_inf_omni:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) +  n_tests_omni + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),])
test <- tidy(sa2_model_i_a1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "omni_ajd_negtest" 
sa_table6 <- rbind(sa_table6, test) 

sa2_model_i_a2 <- coxph(Surv(tstart, tstop, infection) ~ vax_type_omni + prev_inf_omni + vax_type_omni:tstart + prev_inf_omni:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) +  n_tests_omni + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$flu_vax_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),])
test <- tidy(sa2_model_i_a2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "omni_ajd_flu" 
sa_table6 <- rbind(sa_table6, test) 

write.csv(sa_table6, "./Vaccine v reinfection paper/coxph_vax_type_sensitivity.csv") # Save
rm(sa_table6)
gc()

