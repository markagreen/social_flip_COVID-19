##############################
### Vaccine vs reinfection ###
####### Paper analyses #######
##############################

# Purpose: To produce a series of analyses that are included in the paper.

# In case internet drops, redo this when back
# setwd("/Volumes/MAST/Smart_LCR_report")

# Libraries
library(data.table) # For data wrangling (plus those below)
library(lubridate) 
library(tidyr)
library(plyr)
library(ggplot2) # For visualisations
library(viridis)
library(zoo) # For rolling average
library(ggeffects) # For marginal effects
library(sf) # For spatial modelling
library(spdep)

### 1. Tidy data ###

#  Load data
#load("/Volumes/Sharing Folder/2022-01-31/VaccineAnalysisDatasets.RData") # Mac
load("Q:/2022-01-31/VaccineAnalysisDatasets.RData") # Windows

# 1a. Vaccination data #

# Subset only COVID-19 vaccines (flu jabs are in here for example) and people who completed the vaccine
data1 <- data.table(data1) # Convert data type
data1 <- data1[data1$Status == "completed" & (data1$VaccinationProcedure == "1324681000000101" | data1$VaccinationProcedure == "1324691000000104" | data1$VaccinationProcedure == "1362591000000103")] # Subset

# Tidy vaccination data and link to main dataset
data1$t_date <- ymd_hms(data1$DateTimeAdministered) # Convert to time-date
data1 <- data1[data1$t_date < "2021-11-27"] # Drop vaccinations after time period
data1[order(FK_Patient_Link_ID, t_date), num_vax:=1:.N, by=.(FK_Patient_Link_ID)] # Count number of vaccine doses
vax <- data1[num_vax <= 3] # Store upto three vaccine doses only
vaccinated <- dcast(vax, FK_Patient_Link_ID ~ num_vax, value.var = "t_date") # Reshape to wide format to list when got each vaccine
names(vaccinated) <- c("FK_Patient_Link_ID", "vax1_date", "vax2_date", "vax3_date") # Rename variables to help R
data2 <- merge(data2, vaccinated, by = "FK_Patient_Link_ID", all.x = TRUE) # Join onto main dataset
rm(data1, vax, vaccinated) # Save space

# 1b. Test and Trace data #

# Tidy data
data3 <- data.table(data3) # Convert format
data3$t_date <- ymd_hms(data3$SpecimenDate) # Convert to time-date
data3 <- data3[data3$t_date < "2021-11-27"] # Drop vaccinations after time period
data3$day <- as_date(data3$t_date) # Define day
data3 <- data3[order(FK_Patient_Link_ID, day), num_inf:=1:.N, by=.(FK_Patient_Link_ID)] # Count number of positive tests
data3 <- data3[order(FK_Patient_Link_ID, day), tdiff := difftime(day, shift(day, fill=day[1L]), units="days"), by=FK_Patient_Link_ID] # Calculate time inbetween positive tests

# Define first infection
data3$first_inf <- 0
data3$first_inf[data3$num_inf == 1] <- 1 # num_inf = 1 - always gives this tdiff = 0 so is first infection

# Define reinfection(s)
data3$reinf <- 0 
data3$reinf[data3$num_inf > 1 & data3$tdiff > 90] <- 1 # If num_inf > 1 and tdiff is 90+ apart, then define as reinfection 

# Calculate number of infections per person
hold <- data3[data3$first_inf == 1 | data3$reinf == 1] # Select only first infection or any reinfections only
hold <- hold[order(FK_Patient_Link_ID, day), num_inf:=1:.N, by=.(FK_Patient_Link_ID)] # Count actual number of infections
infections <- dcast(hold, FK_Patient_Link_ID ~ num_inf, value.var = "day") # Reshape to wide
names(infections) <- c("FK_Patient_Link_ID", "inf1_date", "inf2_date", "inf3_date") # Rename variables to help R
data2 <- merge(data2, infections, by = "FK_Patient_Link_ID", all.x = TRUE) # Join onto main dataset
rm(hold, data3, infections) # Drop objects not required

# 1c. Final tidying #

# Drop data from analyses
data2 <- data.table(data2) # Convert format
data2 <- data2[data2$FK_Patient_Link_ID != -1,] # Drop missing person ID (n=21)
data2 <- data2[!is.na(data2$FK_Patient_Link_ID),] # Drop missing person ID (n=2)
data2 <- data2[data2$Deceased == "N",] # Exclude deceased people (n=50698)

# Load and join on IMD decile
imd <- fread("/Volumes/MAST/SMART_LCR_report/Data/imd2019.csv") # Mac
#imd <- fread("R:/SMART_LCR_report/Data/imd2019.csv") # Windows
data2 <- merge(data2, imd, by.x = "LSOA_Code", by.y = "lsoa11", all.x = TRUE) # join together
rm(imd)



## 2. Generate variables for analysis ##

# 2a. Outcome variable #

# Tested positive from 1st September
data2$outcome <- 0 # Create blank variable
data2$outcome[data2$inf1_date >= "2021-09-01 00:00:00.0000000" | data2$inf2_date >= "2021-09-01 00:00:00.0000000" | data2$inf3_date >= "2021-09-01 00:00:00.0000000"] <- 1 # Tested positive over this time period

# Date of positive test
data2$outcome_date <- NA # Create blank variable
data2$outcome_date <- as.Date(data2$outcome_date) # Define as date variable
data2$outcome_date[data2$inf1_date >= "2021-09-01 00:00:00.0000000" & !is.na(data2$inf1_date)] <- data2$inf1_date[data2$inf1_date >= "2021-09-01 00:00:00.0000000" & !is.na(data2$inf1_date)] # If have tested positive, then store date
data2$outcome_date[data2$inf2_date >= "2021-09-01 00:00:00.0000000" & !is.na(data2$inf2_date)] <- data2$inf2_date[data2$inf2_date >= "2021-09-01 00:00:00.0000000" & !is.na(data2$inf2_date)]  # Repeat (note people can't test positive more than once in this study period)
data2$outcome_date[data2$inf3_date >= "2021-09-01 00:00:00.0000000" & !is.na(data2$inf3_date)] <- data2$inf3_date[data2$inf3_date >= "2021-09-01 00:00:00.0000000" & !is.na(data2$inf3_date)]  # Repeat

# 2b. Vaccination status #

# Fully vaccinated -6 months to -3 weeks of 1st September
data2$vax_v1 <- 0 # Create blank variable
data2$vax_v1[data2$vax2_date < "2021-08-11 00:00:00.0000000"] <- 1 # Vaccinated as per definition

# Fully vaccinated -6 months to -3 weeks of 1st September with booster
data2$vax_v1b <- 0 # Create blank variable
data2$vax_v1b[data2$vax2_date < "2021-08-11 00:00:00.0000000"] <- 1 # Vaccinated as per
data2$vax_v1b[data2$vax3_date < "2021-08-11 00:00:00.0000000"] <- 2 # Define if boosted
data2$vax_v1b <- factor(data2$vax_v1b) # Store as factor

# Above but add in first dose
data2$vax_v1c <- 0 # Create blank variable
data2$vax_v1c[data2$vax1_date < "2021-08-11 00:00:00.0000000"] <- 1 # First dose
data2$vax_v1c[data2$vax2_date < "2021-08-11 00:00:00.0000000"] <- 2 # Second dose
data2$vax_v1c[data2$vax3_date < "2021-08-11 00:00:00.0000000"] <- 3 # Define if boosted
data2$vax_v1c <- factor(data2$vax_v1c) # Store as factor

# Fully vaccinated -6 months to -3 months (2) / -3 months to -3 weeks of 1st September (1)
data2$vax_v2 <- 0 # Create blank variable
data2$vax_v2[data2$vax2_date >= "2021-06-01 00:00:00.0000000" & data2$vax2_date < "2021-08-11 00:00:00.0000000"] <- 1 # Vaccinated as per definition
data2$vax_v2[data2$vax2_date < "2021-06-01 00:00:00.0000000"] <- 2 # Vaccinated as per definition
data2$vax_v2 <- factor(data2$vax_v2) # Store as factor

# Fully vaccinated -6 months to -3 months (2) / -3 months to -3 weeks of 1st September (1), with booster (3)
data2$vax_v2b <- 0 # Create blank variable
data2$vax_v2b[data2$vax2_date >= "2021-06-01 00:00:00.0000000" & data2$vax2_date < "2021-08-11 00:00:00.0000000"] <- 1 # Vaccinated as per definition
data2$vax_v2b[data2$vax2_date < "2021-06-01 00:00:00.0000000"] <- 2 # Vaccinated as per definition
data2$vax_v2b[data2$vax3_date < "2021-08-11 00:00:00.0000000"] <- 3 # Define if boosted
data2$vax_v2b <- factor(data2$vax_v2b) # Store as factor

# Received first or second dose -6 months to -3 weeks of 1st September
data2$vax_v3 <- 0 # Create blank variable
data2$vax_v3[data2$vax1_date < "2021-08-11 00:00:00.0000000"] <- 1 # Vaccinated as per definition

# Received third dose by -3 weeks of 1st September
data2$boosted <- 0 # Create blank variable
data2$boosted[data2$vax3_date < "2021-08-11 00:00:00.0000000"] <- 1 # Vaccinated as per definition

# Received third dose any point
data2$boosted_any <- 0 # Create blank variable
data2$boosted_any[!is.na(data2$vax3_date)] <- 1 # Vaccinated as per definition

# Received any vaccine after inclusion criteria
data2$exclude <- 0 # Create blank variable
data2$exclude[data2$vax3_date >= "2021-08-11 00:00:00.0000000"] <- 1 # Boosted
data2$exclude[data2$vax2_date >= "2021-08-11 00:00:00.0000000"] <- 2 # Received second dose
data2$exclude[data2$vax1_date >= "2021-08-11 00:00:00.0000000"] <- 3 # Received first dose

# 2c. Previous infection #

# Previous infection -6 months to -2 weeks
data2$prev_inf_v1 <- 0 # Create blank variable
data2$prev_inf_v1[data2$inf1_date < "2021-08-18 00:00:00.0000000" | data2$inf2_date < "2021-08-18 00:00:00.0000000" | data2$inf3_date < "2021-08-18 00:00:00.0000000"] <- 1 # Tested positive in this time period

# Previous infection -6 months to -3 months (2), -3 months to -2 weeks (1)
data2$prev_inf_v2 <- 0 # Create blank variable
data2$prev_inf_v2[(data2$inf1_date >= "2021-06-01 00:00:00.0000000" & data2$inf1_date < "2021-08-18 00:00:00.0000000") | (data2$inf3_date >= "2021-06-01 00:00:00.0000000" & data2$inf2_date < "2021-08-18 00:00:00.0000000") | (data2$inf3_date >= "2021-06-01 00:00:00.0000000" & data2$inf3_date < "2021-08-18 00:00:00.0000000")] <- 1 # Tested positive in first period
data2$prev_inf_v2[data2$inf1_date < "2021-06-01 00:00:00.0000000" | data2$inf2_date < "2021-06-01 00:00:00.0000000" | data2$inf3_date < "2021-06-01 00:00:00.0000000"] <- 2 # Vaccinated as per definition
data2$prev_inf_v2 <- factor(data2$prev_inf_v2) 

# More than one reinfection (max is 3)
data2$inf_3times <- 0
data2$inf_3times[!is.na(data2$inf3_date)] <- 0

# 2d. Control variables #

# Patient age squared
data2$Patient.Age2 <- data2$Patient.Age ^ 2

# Age bands
data2[, age_group:=cut(Patient.Age, breaks = c(0,10,17,29,39,49,59,69,120), include.lowest=T)] # Split into 10-year age bands
data2[age_group=="[0,10]", age_band:="0-10"] # Add labels
data2[age_group=="(10,17]", age_band:="11-17"]
data2[age_group=="(17,29]", age_band:="18-29"]
data2[age_group=="(29,39]", age_band:="30-39"]
data2[age_group=="(39,49]", age_band:="40-49"]
data2[age_group=="(49,59]", age_band:="50-59"]
data2[age_group=="(59,69]", age_band:="60-69"]
data2[age_group=="(69,120]", age_band:="70+"]
data2[, age_group:=NULL] # Drop
data2$age_band <- as.factor(data2$age_band) # Set as factor
data2$age_band <- relevel(data2$age_band, ref = "18-29") # Define White as reference

# IMD
data2$imd_score2 <- data2$imd_score ^ 2 # IMD squared

# Sex
data2$Sex[data2$Sex == "U"] <- NA # Drop missing sex (n=117))

# Ethnicity
# What to do with NULL ethnicity?
data2$EthnicMainGroup <- as.factor(data2$EthnicMainGroup) # Set as factor
data2$EthnicMainGroup <- relevel(data2$EthnicMainGroup, ref = "White") # Define White as reference

# Previous number of tests
load("./Data/number_tests.RData") # Load data - Mac
#load("R:/SMART_LCR_report/Data/number_tests.RData") # Windows
data2 <- merge(data2, tests_6mo, by = "FK_Patient_Link_ID", all.x = TRUE) # Join onto main dataset
data2$lft_6mo[is.na(data2$lft_6mo)] <- 0 # Assume if missing then have had no tests so fiil in with 0s (and repeat below)
data2$pcr_6mo[is.na(data2$pcr_6mo)] <- 0 
data2$lft_3mo[is.na(data2$lft_3mo)] <- 0 
data2$pcr_3mo[is.na(data2$pcr_3mo)] <- 0 
data2$lft_1mo[is.na(data2$lft_1mo)] <- 0 
data2$pcr_1mo[is.na(data2$pcr_1mo)] <- 0 
rm(tests_6mo) # Save space

# Create explanatory var
data2$test_1mo <- 0
data2$test_1mo[data2$test_1mo > 0 | data2$pcr_1mo > 0] <- 1


### 3. Descriptive tables and plots ###

# To stop reporting numbers by e+10 etc
options(scipen=999) 

# Table 1 #
table(data2$outcome)
table(data2$vax_v1)
table(data2$prev_inf_v1)
mean(data2$Patient.Age)
sd(data2$Patient.Age)
table(data2$Sex)
table(data2$EthnicMainGroup)
# (table(data2$EthnicMainGroup)/2811086)*100 # Percentages
mean(data2$imd_score, na.rm=T)
sd(data2$imd_score, na.rm=T)
# med_vax <- data2[,list(age = median(Patient.Age, na.rm=TRUE), imd = median(imd_score, na.rm=T)), by = "vax_v1"] # If want to stratify

# Summary stats by group
total_pos <- data2[, list(positive = sum(outcome, na.rm = T), pop = .N), by = c("vax_v1", "prev_inf_v1")] # Aggregate to get total numbers for positive tests 

# Figure 1a - vax vs reinfection#

# Aggregate to number of tests per day
trends_day <- data2[, list(positive = sum(outcome, na.rm = T)), by = c("outcome_date", "vax_v1", "prev_inf_v1")] # Aggregate to get total numbers for positive tests 
pop <- data2[, list(pop = .N), by = c("vax_v1", "prev_inf_v1")] # Get population sizes
trends_day <- merge(trends_day, pop, by = c("vax_v1", "prev_inf_v1"), all.x = TRUE) # Join together

# Drop missing dates
trends_day <- trends_day[!is.na(trends_day$outcome_date)]

# Percentage of tests per day that were positive
trends_day$positive_rate <- (trends_day$positive / (trends_day$pop)) * 10000

# Create 7 day rolling averages for plotting purposes
trends_day <- trends_day %>%
  dplyr::arrange(outcome_date) %>% # Order by day (ascending order)
  dplyr::group_by(vax_v1, prev_inf_v1) %>% # Group by characteristics
  dplyr::mutate(positive_rate_7day = zoo::rollmean(positive_rate, k = 7, fill = NA)) # Generate 7 day average

# Plot labels
labels <- c("Unvaccinated", "Vaccinated") # Define facet labels
names(labels) <- c("0", "1")

# Plot
fig1a <- ggplot(trends_day) +
  geom_line(aes(x = outcome_date, y = positive_rate_7day, group = factor(prev_inf_v1), color = factor(prev_inf_v1))) + # 7 day rolling average
  facet_wrap(~ vax_v1, labeller = labeller(vax_v1 = labels)) +
  scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9, labels = c("No", "Yes")) + # Make colour blind friendly (ribbon)
  scale_x_date(date_labels = "%b", date_breaks = "1 month") + # Add monthly labels
  xlab("Date") +
  ylab("Incidence per 100,000") +
  labs(color = "Previous positive") +
  theme(text = element_text(size = 12)) + # Change text size
  ylim(0,NA) +
  theme(legend.position="bottom")
fig1a
ggsave(plot = fig1a, file = "./Vaccine v reinfection paper/figure1a.jpeg")


# Figure 1b - age band #

# Aggregate to number of tests per day
trends_day_age <- data2[, list(positive = sum(outcome, na.rm = T)), by = c("outcome_date", "vax_v1", "age_band")] # Aggregate to get total numbers for positive tests 
pop_age <- data2[, list(pop = .N), by = c("vax_v1", "age_band")] # Get population sizes
trends_day_age <- merge(trends_day_age, pop_age, by = c("vax_v1", "age_band"), all.x = TRUE) # Join together

# Drop missing dates
trends_day_age <- trends_day_age[!is.na(trends_day_age$outcome_date)]

# Percentage of tests per day that were positive
trends_day_age$positive_rate <- (trends_day_age$positive / (trends_day_age$pop)) * 10000

# Create 7 day rolling averages for plotting purposes
trends_day_age <- trends_day_age %>%
  dplyr::arrange(outcome_date) %>% # Order by day (ascending order)
  dplyr::group_by(vax_v1, age_band) %>% # Group by characteristics
  dplyr::mutate(positive_rate_7day = zoo::rollmean(positive_rate, k = 7, fill = NA)) # Generate 7 day average

# Plot labels
labels <- c("Unvaccinated", "Vaccinated") # Define facet labels
names(labels) <- c("0", "1")

# Set order of age band
trends_day_age$age_band <- ordered(trends_day_age$age_band, levels = c("0-10", "11-17", "18-29", "30-39", "40-49", "50-59", "60-69", "70+"))

# Plot
fig1b <- ggplot(trends_day_age) +
  geom_line(aes(x = outcome_date, y = positive_rate_7day, group = factor(age_band), color = factor(age_band))) + # 7 day rolling average
  facet_wrap(~ vax_v1, labeller = labeller(vax_v1 = labels)) +
  scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9) + # Make colour blind friendly (ribbon)
  scale_x_date(date_labels = "%b", date_breaks = "1 month") + # Add monthly labels
  xlab("Date") +
  ylab("Incidence per 100,000") +
  labs(color = "Age band") +
  theme(text = element_text(size = 12)) + # Change text size
  ylim(0,NA) +
  theme(legend.position="bottom")
fig1b
ggsave(plot = fig1b, file = "./Vaccine v reinfection paper/figure1b.jpeg")


# Figure 1c - ethnicity #

# Aggregate to number of tests per day
trends_day_eth <- data2[, list(positive = sum(outcome, na.rm = T)), by = c("outcome_date", "vax_v1", "EthnicMainGroup")] # Aggregate to get total numbers for positive tests 
pop_eth <- data2[, list(pop = .N), by = c("vax_v1", "EthnicMainGroup")] # Get population sizes
trends_day_eth <- merge(trends_day_eth, pop_eth, by = c("vax_v1", "EthnicMainGroup"), all.x = TRUE) # Join together

# Drop missing data
trends_day_eth <- trends_day_eth[!is.na(trends_day_eth$outcome_date)]

# Percentage of tests per day that were positive
trends_day_eth$positive_rate <- (trends_day_eth$positive / (trends_day_eth$pop)) * 10000

# Create 7 day rolling averages for plotting purposes
trends_day_eth <- trends_day_eth %>%
  dplyr::arrange(outcome_date) %>% # Order by day (ascending order)
  dplyr::group_by(vax_v1, EthnicMainGroup) %>% # Group by characteristics
  dplyr::mutate(positive_rate_7day = zoo::rollmean(positive_rate, k = 7, fill = NA)) # Generate 7 day average

# Plot labels
labels <- c("Unvaccinated", "Vaccinated") # Define facet labels
names(labels) <- c("0", "1")

# Set order of ethnicity
trends_day_eth$EthnicMainGroup <- ordered(trends_day_eth$EthnicMainGroup, levels = c("Asian or Asian British", "Black or Black British", "Mixed", "NULL", "Other Ethnic Groups", "White"))
levels(trends_day_eth$EthnicMainGroup)[levels(trends_day_eth$EthnicMainGroup)=="Other Ethnic Groups"] <- "Other" # Rename for plotting

# Plot
fig1c <- ggplot(trends_day_eth) +
  geom_line(aes(x = outcome_date, y = positive_rate_7day, group = factor(EthnicMainGroup), color = factor(EthnicMainGroup))) + # 7 day rolling average
  facet_wrap(~ vax_v1, labeller = labeller(vax_v1 = labels)) +
  scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9) + # Make colour blind friendly (ribbon)
  scale_x_date(date_labels = "%b", date_breaks = "1 month") + # Add monthly labels
  xlab("Date") +
  ylab("Incidence per 100,000") +
  labs(color = "Ethnicity") +
  theme(text = element_text(size = 12)) + # Change text size
  ylim(0,NA) +
  theme(legend.position="bottom")
fig1c
ggsave(plot = fig1c, file = "./Vaccine v reinfection paper/figure1c.jpeg")


# Figure 1d - deprivation #

# Aggregate to number of tests per day
trends_day_imd <- data2[, list(positive = sum(outcome, na.rm = T)), by = c("outcome_date", "vax_v1", "imd_decile")] # Aggregate to get total numbers for positive tests 
pop_imd <- data2[, list(pop = .N), by = c("vax_v1", "imd_decile")] # Get population sizes
trends_day_imd <- merge(trends_day_imd, pop_imd, by = c("vax_v1", "imd_decile"), all.x = TRUE) # Join together

# Drop missing data
trends_day_imd <- trends_day_imd[!is.na(trends_day_imd$outcome_date)]
trends_day_imd <- trends_day_imd[!is.na(trends_day_imd$imd_decile),]

# Percentage of tests per day that were positive
trends_day_imd$positive_rate <- (trends_day_imd$positive / (trends_day_imd$pop)) * 10000

# Create 7 day rolling averages for plotting purposes
trends_day_imd <- trends_day_imd %>%
  dplyr::arrange(outcome_date) %>% # Order by day (ascending order)
  dplyr::group_by(vax_v1, imd_decile) %>% # Group by characteristics
  dplyr::mutate(positive_rate_7day = zoo::rollmean(positive_rate, k = 7, fill = NA)) # Generate 7 day average

# Plot labels
labels <- c("Unvaccinated", "Vaccinated") # Define facet labels
names(labels) <- c("0", "1")

# Plot
fig1d <- ggplot(trends_day_imd) +
  geom_line(aes(x = outcome_date, y = positive_rate_7day, group = factor(imd_decile), color = factor(imd_decile))) + # 7 day rolling average
  facet_wrap(~ vax_v1, labeller = labeller(vax_v1 = labels)) +
  scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9) + # Make colour blind friendly (ribbon)
  scale_x_date(date_labels = "%b", date_breaks = "1 month") + # Add monthly labels
  xlab("Date") +
  ylab("Incidence per 100,000") +
  labs(color = "IMD decile") +
  theme(text = element_text(size = 12)) + # Change text size
  ylim(0,NA) +
  theme(legend.position="bottom")
fig1d
ggsave(plot = fig1d, file = "./Vaccine v reinfection paper/figure1d.jpeg")


# Put all plots together into one #

library(patchwork)
figure1 <- (fig1a + fig1b) / 
  (fig1c + fig1d) +
  plot_annotation(tag_levels = 'A')
figure1
ggsave(plot = figure1, file = "./Vaccine v reinfection paper/figure1.jpeg")

# Tidy
rm(figure1, fig1a, fig1b, fig1c, fig1d)


### 4. Regression analyses ###


# 4a. Main explanatory variables #

# Unadjusted - vaccination status
model1a <- glm(outcome ~ vax_v1, data = data2, family = "binomial") # Model
table <- as.data.table(cbind(exp(model1a$coefficients), exp(confint(model1a)), summary(model1a)$coefficients[,"Pr(>|z|)"]), keep.rownames = TRUE) # Save key results (estimates and CIs as odds ratios)
names(table)[names(table) == "V1"] <- "estimate" # Rename columns
names(table)[names(table) == "V4"] <- "p-value"
names(table)[names(table) == "rn"] <- "variable"
names(table)[names(table) == "2.5 %"] <- "lower"
names(table)[names(table) == "97.5 %"] <- "upper"
table$model <- "a" # Make note of which model

# Unadjusted - previous infection status
model1b <- glm(outcome ~ prev_inf_v1, data = data2, family = "binomial")
hold <- as.data.table(cbind(exp(model1b$coefficients), exp(confint(model1b)), summary(model1b)$coefficients[,"Pr(>|z|)"]), keep.rownames = TRUE) # Save key results (estimates and CIs as odds ratios)
hold$model <- "b" # Make note of which model
names(hold)[names(hold) == "V1"] <- "estimate" # Rename columns
names(hold)[names(hold) == "V4"] <- "p-value"
names(hold)[names(hold) == "rn"] <- "variable"
names(hold)[names(hold) == "2.5 %"] <- "lower"
names(hold)[names(hold) == "97.5 %"] <- "upper"
table <- rbind(table, hold) # join onto main table of results

# Unadjusted - vaccination and previous infection status
model1c <- glm(outcome ~ prev_inf_v1 + vax_v1, data = data2, family = "binomial")
hold <- as.data.table(cbind(exp(model1c$coefficients), exp(confint(model1c)), summary(model1c)$coefficients[,"Pr(>|z|)"]), keep.rownames = TRUE) # Save key results (estimates and CIs as odds ratios)
hold$model <- "c" # Make note of which model
names(hold)[names(hold) == "V1"] <- "estimate" # Rename columns
names(hold)[names(hold) == "V4"] <- "p-value"
names(hold)[names(hold) == "rn"] <- "variable"
names(hold)[names(hold) == "2.5 %"] <- "lower"
names(hold)[names(hold) == "97.5 %"] <- "upper"
table <- rbind(table, hold) # join onto main table of results


# Fully adjusted model
model1d <- glm(outcome ~ vax_v1 + prev_inf_v1 + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + imd_score + test_1mo, data = data2, family = "binomial")
hold <- as.data.table(cbind(exp(model1e$coefficients), exp(confint(model1e)), summary(model1e)$coefficients[,"Pr(>|z|)"]), keep.rownames = TRUE) # Save key results (estimates and CIs as odds ratios)
hold$model <- "d" # Make note of which model
names(hold)[names(hold) == "V1"] <- "estimate" # Rename columns
names(hold)[names(hold) == "V4"] <- "p-value"
names(hold)[names(hold) == "rn"] <- "variable"
names(hold)[names(hold) == "2.5 %"] <- "lower"
names(hold)[names(hold) == "97.5 %"] <- "upper"
table <- rbind(table, hold) # join onto main table of results
write.csv(table, "./Vaccine v reinfection paper/model1_results.csv")
rm(hold) # Save space

# Marginal effects at the mean
pred <- ggpredict(model1d, terms = c("prev_inf_v1", "vax_v1b")) # Get margins
# Plot
plot1 <- ggplot(pred, aes(x = factor(x), y = predicted, ymin = conf.low, ymax = conf.high, color = factor(group))) + 
  geom_point(position=position_dodge(0.1)) + 
  geom_errorbar(position=position_dodge(0.1), width = 0.2) +
  scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9, labels = c("No", "Yes", "Yes + booster")) + # Make colour blind friendly (line)
  scale_x_discrete(labels = c("No", "Yes")) +
  labs(
    x = "Previously tested positive for COVID-19",
    y = "Predicted probability of COVID-19",
    colour = "Fully vaccinated"
  ) +
  ylim(0,NA) +
  theme(legend.position="bottom")
plot1
ggsave(plot = plot1, file = "./Vaccine v reinfection paper/model1.jpeg")

# Save space
rm(model1a, model1b, model1c, model1d, table, pred, plot1)


# 4b. Add time since vax/positive into measure #

# Unadjusted - vaccination status
model2a <- glm(outcome ~ vax_v2, data = data2, family = "binomial") # Model
table2 <- as.data.table(cbind(exp(model2a$coefficients), exp(confint(model2a)), summary(model2a)$coefficients[,"Pr(>|z|)"]), keep.rownames = TRUE) # Save key results (estimates and CIs as odds ratios)
names(table2)[names(table2) == "V1"] <- "estimate" # Rename columns
names(table2)[names(table2) == "V4"] <- "p-value"
names(table2)[names(table2) == "rn"] <- "variable"
names(table2)[names(table2) == "2.5 %"] <- "lower"
names(table2)[names(table2) == "97.5 %"] <- "upper"
table2$model <- "a" # Make note of which model

# Unadjusted - vaccination and previous infection status
model2b <- glm(outcome ~ prev_inf_v1 + vax_v2, data = data2, family = "binomial")
hold <- as.data.table(cbind(exp(model2c$coefficients), exp(confint(model2c)), summary(model2c)$coefficients[,"Pr(>|z|)"]), keep.rownames = TRUE) # Save key results (estimates and CIs as odds ratios)
hold$model <- "b" # Make note of which model
names(hold)[names(hold) == "V1"] <- "estimate" # Rename columns
names(hold)[names(hold) == "V4"] <- "p-value"
names(hold)[names(hold) == "rn"] <- "variable"
names(hold)[names(hold) == "2.5 %"] <- "lower"
names(hold)[names(hold) == "97.5 %"] <- "upper"
table2 <- rbind(table2, hold) # join onto main table of results


# Fully adjusted model
model2c <- glm(outcome ~ vax_v2 + prev_inf_v1 + age_band + factor(Sex) + factor(EthnicMainGroup) + imd_score + test_1mo, data = data2, family = "binomial")
hold <- as.data.table(cbind(exp(model2e$coefficients), exp(confint(model2e)), summary(model2e)$coefficients[,"Pr(>|z|)"]), keep.rownames = TRUE) # Save key results (estimates and CIs as odds ratios)
hold$model <- "c" # Make note of which model
names(hold)[names(hold) == "V1"] <- "estimate" # Rename columns
names(hold)[names(hold) == "V4"] <- "p-value"
names(hold)[names(hold) == "rn"] <- "variable"
names(hold)[names(hold) == "2.5 %"] <- "lower"
names(hold)[names(hold) == "97.5 %"] <- "upper"
table2 <- rbind(table2, hold) # join onto main table of results
write.csv(table2, "./Vaccine v reinfection paper/model2_results.csv")
rm(hold) # Save space

# Marginal effects at the mean
pred2 <- ggpredict(model2c, terms = c("prev_inf_v1", "vax_v2b")) # Get margins
# Plot
plot2 <- ggplot(pred2, aes(x = factor(x), y = predicted, ymin = conf.low, ymax = conf.high, color = factor(group))) + 
  geom_point(position=position_dodge(0.1)) + 
  geom_errorbar(position=position_dodge(0.1), width = 0.2) +
  scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9, labels = c("No", "-3 weeks to 3 months", "3+ months", "Boosted")) + # Make colour blind friendly (line)
  scale_x_discrete(labels = c("No", "Yes")) +
  labs(
    x = "Previously tested positive for COVID-19",
    y = "Predicted probability of COVID-19",
    colour = "Fully vaccinated"
  ) +
  ylim(0,NA) +
  theme(legend.position="bottom")
plot2
ggsave(plot = plot2, file = "./Vaccine v reinfection paper/model2.jpeg")

# Save space
rm(model2a, model2b, model2c, model2d, model2e, table2, pred2, plot2)


# 4c. Stratify model by age group #


# Stratify by age group
model3a <- glm(outcome ~ vax_v1b + prev_inf_v1 + factor(Sex) + factor(EthnicMainGroup) + imd_score + test_1mo, data = data2[data2$Patient.Age >=0 & data2$Patient.Age < 11], family = "binomial")
table3 <- as.data.table(cbind(exp(model3a$coefficients), exp(confint.default(model3a)), summary(model3a)$coefficients[,"Pr(>|z|)"]), keep.rownames = TRUE) # Save key results (estimates and CIs as odds ratios)
table3$model <- "0-10" # Make note of which model
names(table3)[names(table3) == "V1"] <- "estimate" # Rename columns
names(table3)[names(table3) == "V4"] <- "p-value"
names(table3)[names(table3) == "rn"] <- "variable"
names(table3)[names(table3) == "2.5 %"] <- "lower"
names(table3)[names(table3) == "97.5 %"] <- "upper"

model3b <- glm(outcome ~ vax_v1b + prev_inf_v1 + factor(Sex) + factor(EthnicMainGroup) + imd_score + test_1mo, data = data2[data2$Patient.Age >=11 & data2$Patient.Age < 18], family = "binomial")
hold <- as.data.table(cbind(exp(model3b$coefficients), exp(confint.default(model3b)), summary(model3b)$coefficients[,"Pr(>|z|)"]), keep.rownames = TRUE) # Save key results (estimates and CIs as odds ratios)
hold$model <- "11-17" # Make note of which model
names(hold)[names(hold) == "V1"] <- "estimate" # Rename columns
names(hold)[names(hold) == "V4"] <- "p-value"
names(hold)[names(hold) == "rn"] <- "variable"
names(hold)[names(hold) == "2.5 %"] <- "lower"
names(hold)[names(hold) == "97.5 %"] <- "upper"
table3 <- rbind(table3, hold) # join onto main table of results

model3c <- glm(outcome ~ vax_v1b + prev_inf_v1 + factor(Sex) + factor(EthnicMainGroup) + imd_score + test_1mo, data = data2[data2$Patient.Age >=18 & data2$Patient.Age < 30], family = "binomial")
hold <- as.data.table(cbind(exp(model3c$coefficients), exp(confint.default(model3c)), summary(model3c)$coefficients[,"Pr(>|z|)"]), keep.rownames = TRUE) # Save key results (estimates and CIs as odds ratios)
hold$model <- "18-29" # Make note of which model
names(hold)[names(hold) == "V1"] <- "estimate" # Rename columns
names(hold)[names(hold) == "V4"] <- "p-value"
names(hold)[names(hold) == "rn"] <- "variable"
names(hold)[names(hold) == "2.5 %"] <- "lower"
names(hold)[names(hold) == "97.5 %"] <- "upper"
table3 <- rbind(table3, hold) # join onto main table of results

model3d <- glm(outcome ~ vax_v1b + prev_inf_v1 + factor(Sex) + factor(EthnicMainGroup) + imd_score + test_1mo, data = data2[data2$Patient.Age >=30 & data2$Patient.Age < 40], family = "binomial")
hold <- as.data.table(cbind(exp(model3d$coefficients), exp(confint.default(model3d)), summary(model3d)$coefficients[,"Pr(>|z|)"]), keep.rownames = TRUE) # Save key results (estimates and CIs as odds ratios)
hold$model <- "30-39" # Make note of which model
names(hold)[names(hold) == "V1"] <- "estimate" # Rename columns
names(hold)[names(hold) == "V4"] <- "p-value"
names(hold)[names(hold) == "rn"] <- "variable"
names(hold)[names(hold) == "2.5 %"] <- "lower"
names(hold)[names(hold) == "97.5 %"] <- "upper"
table3 <- rbind(table3, hold) # join onto main table of results

model3e <- glm(outcome ~ vax_v1b + prev_inf_v1 + factor(Sex) + factor(EthnicMainGroup) + imd_score + test_1mo, data = data2[data2$Patient.Age >=40 & data2$Patient.Age < 50], family = "binomial")
hold <- as.data.table(cbind(exp(model3e$coefficients), exp(confint.default(model3e)), summary(model3e)$coefficients[,"Pr(>|z|)"]), keep.rownames = TRUE) # Save key results (estimates and CIs as odds ratios)
hold$model <- "40-49" # Make note of which model
names(hold)[names(hold) == "V1"] <- "estimate" # Rename columns
names(hold)[names(hold) == "V4"] <- "p-value"
names(hold)[names(hold) == "rn"] <- "variable"
names(hold)[names(hold) == "2.5 %"] <- "lower"
names(hold)[names(hold) == "97.5 %"] <- "upper"
table3 <- rbind(table3, hold) # join onto main table of results

model3f <- glm(outcome ~ vax_v1b + prev_inf_v1 + factor(Sex) + factor(EthnicMainGroup) + imd_score + test_1mo, data = data2[data2$Patient.Age >=50 & data2$Patient.Age < 60], family = "binomial")
hold <- as.data.table(cbind(exp(model3f$coefficients), exp(confint.default(model3f)), summary(model3f)$coefficients[,"Pr(>|z|)"]), keep.rownames = TRUE) # Save key results (estimates and CIs as odds ratios)
hold$model <- "50-59" # Make note of which model
names(hold)[names(hold) == "V1"] <- "estimate" # Rename columns
names(hold)[names(hold) == "V4"] <- "p-value"
names(hold)[names(hold) == "rn"] <- "variable"
names(hold)[names(hold) == "2.5 %"] <- "lower"
names(hold)[names(hold) == "97.5 %"] <- "upper"
table3 <- rbind(table3, hold) # join onto main table of results

model3g <- glm(outcome ~ vax_v1b + prev_inf_v1 + factor(Sex) + factor(EthnicMainGroup) + imd_score + test_1mo, data = data2[data2$Patient.Age >=60 & data2$Patient.Age < 70], family = "binomial")
hold <- as.data.table(cbind(exp(model3g$coefficients), exp(confint.default(model3g)), summary(model3g)$coefficients[,"Pr(>|z|)"]), keep.rownames = TRUE) # Save key results (estimates and CIs as odds ratios)
hold$model <- "60-69" # Make note of which model
names(hold)[names(hold) == "V1"] <- "estimate" # Rename columns
names(hold)[names(hold) == "V4"] <- "p-value"
names(hold)[names(hold) == "rn"] <- "variable"
names(hold)[names(hold) == "2.5 %"] <- "lower"
names(hold)[names(hold) == "97.5 %"] <- "upper"
table3 <- rbind(table3, hold) # join onto main table of results

model3h <- glm(outcome ~ vax_v1b + prev_inf_v1 + factor(Sex) + factor(EthnicMainGroup) + imd_score + test_1mo, data = data2[data2$Patient.Age >=70], family = "binomial")
hold <- as.data.table(cbind(exp(model3h$coefficients), exp(confint.default(model3h)), summary(model3h)$coefficients[,"Pr(>|z|)"]), keep.rownames = TRUE) # Save key results (estimates and CIs as odds ratios)
hold$model <- "70+" # Make note of which model
names(hold)[names(hold) == "V1"] <- "estimate" # Rename columns
names(hold)[names(hold) == "V4"] <- "p-value"
names(hold)[names(hold) == "rn"] <- "variable"
names(hold)[names(hold) == "2.5 %"] <- "lower"
names(hold)[names(hold) == "97.5 %"] <- "upper"
table3 <- rbind(table3, hold) # join onto main table of results
write.csv(table3, "./Vaccine v reinfection paper/age_group_analysis.csv")

# Plot
table3_to_plot <- table3[table3$variable == "vax_v1b1" | table3$variable == "vax_v1b2" | table3$variable == "prev_inf_v1"]
plot3 <- ggplot(table3_to_plot, aes(x = model, y = estimate, ymin = lower, ymax = upper, color = factor(variable))) + 
  geom_point(position=position_dodge(0.1)) + 
  geom_errorbar(position=position_dodge(0.1), width = 0.2) +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9, labels = c("Previous infection", "Fully vaccinated", "Fully vaccinated + booster")) + # Make colour blind friendly (line)
  labs(
    x = "Age group",
    y = "Odds Ratio",
    colour = "Coefficient"
  ) +
  coord_cartesian(ylim=c(0,2.5)) +
  theme(legend.position="bottom")
plot3
ggsave(plot = plot3, file = "./Vaccine v reinfection paper/model3.jpeg")

# Plot for IMD
table3_to_plot2 <- table3[table3$variable == "imd_score"]
plot3a <- ggplot(table3_to_plot2, aes(x = model, y = estimate, ymin = lower, ymax = upper)) + 
  geom_point(position=position_dodge(0.1)) + 
  geom_errorbar(position=position_dodge(0.1), width = 0.2) +
  geom_hline(yintercept = 1, linetype = 2) +
  labs(
    x = "Age group",
    y = "Odds Ratio"
  ) +
  #coord_cartesian(ylim=c(0,2.5)) +
  theme(legend.position="bottom")
plot3a
ggsave(plot = plot3, file = "./Vaccine v reinfection paper/model3_ses.jpeg")

# Save space
rm(model3a, model3b, model3c, model3d, model3e, model3f, model3g, model3h, model3i)


# 4d. Sensitivity analyses #

# Had one dose of vaccine
model4a <- glm(outcome ~ vax_v3 + prev_inf_v1 + age_band + factor(Sex) + factor(EthnicMainGroup) + imd_score + test_1mo, data = data2, family = "binomial")
table4 <- as.data.table(cbind(exp(model4a$coefficients), exp(confint(model4a)), summary(model4a)$coefficients[,"Pr(>|z|)"]), keep.rownames = TRUE) # Save key results (estimates and CIs as odds ratios)
names(table4)[names(table4) == "V1"] <- "estimate" # Rename columns
names(table4)[names(table4) == "V4"] <- "p-value"
names(table4)[names(table4) == "rn"] <- "variable"
names(table4)[names(table4) == "2.5 %"] <- "lower"
names(table4)[names(table4) == "97.5 %"] <- "upper"
table4$model <- "One dose"

# IMD decile
model4b <- glm(outcome ~ vax_v1b + prev_inf_v1 + age_band + factor(Sex) + factor(EthnicMainGroup) + factor(imd_decile) + test_1mo, data = data2, family = "binomial")
hold <- as.data.table(cbind(exp(model4b$coefficients), exp(confint(model4b)), summary(model4b)$coefficients[,"Pr(>|z|)"]), keep.rownames = TRUE) # Save key results (estimates and CIs as odds ratios)
names(hold)[names(hold) == "V1"] <- "estimate" # Rename columns
names(hold)[names(hold) == "V4"] <- "p-value"
names(hold)[names(hold) == "rn"] <- "variable"
names(hold)[names(hold) == "2.5 %"] <- "lower"
names(hold)[names(hold) == "97.5 %"] <- "upper"
hold$model <- "IMD Decile"
table4 <- rbind(table4, hold) # join onto main table of results

# Remove boosted any point in time
model4c <- glm(outcome ~ vax_v1 + prev_inf_v1 + age_band + factor(Sex) + factor(EthnicMainGroup) + imd_score + test_1mo, data = data2[data2$boosted_any == 0], family = "binomial")
hold <- as.data.table(cbind(exp(model4c$coefficients), exp(confint(model4c)), summary(model4c)$coefficients[,"Pr(>|z|)"]), keep.rownames = TRUE) # Save key results (estimates and CIs as odds ratios)
names(hold)[names(hold) == "V1"] <- "estimate" # Rename columns
names(hold)[names(hold) == "V4"] <- "p-value"
names(hold)[names(hold) == "rn"] <- "variable"
names(hold)[names(hold) == "2.5 %"] <- "lower"
names(hold)[names(hold) == "97.5 %"] <- "upper"
hold$model <- "remove boosted all - model 1"
table4 <- rbind(table4, hold) # join onto main table of results

# Remove boosted any point in time
model4d <- glm(outcome ~ vax_v2 + prev_inf_v1 + age_band + factor(Sex) + factor(EthnicMainGroup) + imd_score + test_1mo, data = data2[data2$boosted_any == 0], family = "binomial")
hold <- as.data.table(cbind(exp(model4d$coefficients), exp(confint(model4d)), summary(model4d)$coefficients[,"Pr(>|z|)"]), keep.rownames = TRUE) # Save key results (estimates and CIs as odds ratios)
names(hold)[names(hold) == "V1"] <- "estimate" # Rename columns
names(hold)[names(hold) == "V4"] <- "p-value"
names(hold)[names(hold) == "rn"] <- "variable"
names(hold)[names(hold) == "2.5 %"] <- "lower"
names(hold)[names(hold) == "97.5 %"] <- "upper"
hold$model <- "remove boosted all - model 2"
table4 <- rbind(table4, hold) # join onto main table of results

# Received any vaccinated since September
data2$vac_sept <- 0
data2$vac_sept[data2$vax1_date >= "2021-09-01"] <- 1
data2$vac_sept[data2$vax2_date >= "2021-09-01"] <- 1
data2$vac_sept[data2$vax3_date >= "2021-09-01"] <- 1
model4e <- glm(outcome ~ vax_v1 + prev_inf_v1 + age_band + factor(Sex) + factor(EthnicMainGroup) + imd_score + test_1mo, data = data2[data2$vac_sept == 0], family = "binomial")
hold <- as.data.table(cbind(exp(model4e$coefficients), exp(confint(model4e)), summary(model4e)$coefficients[,"Pr(>|z|)"]), keep.rownames = TRUE) # Save key results (estimates and CIs as odds ratios)
names(hold)[names(hold) == "V1"] <- "estimate" # Rename columns
names(hold)[names(hold) == "V4"] <- "p-value"
names(hold)[names(hold) == "rn"] <- "variable"
names(hold)[names(hold) == "2.5 %"] <- "lower"
names(hold)[names(hold) == "97.5 %"] <- "upper"
hold$model <- "vax since sept - model 1"
table4 <- rbind(table4, hold) # join onto main table of results

model4f <- glm(outcome ~ vax_v2 + prev_inf_v1 + age_band + factor(Sex) + factor(EthnicMainGroup) + imd_score + test_1mo, data = data2[data2$vac_sept == 0], family = "binomial")
hold <- as.data.table(cbind(exp(model4f$coefficients), exp(confint(model4f)), summary(model4f)$coefficients[,"Pr(>|z|)"]), keep.rownames = TRUE) # Save key results (estimates and CIs as odds ratios)
names(hold)[names(hold) == "V1"] <- "estimate" # Rename columns
names(hold)[names(hold) == "V4"] <- "p-value"
names(hold)[names(hold) == "rn"] <- "variable"
names(hold)[names(hold) == "2.5 %"] <- "lower"
names(hold)[names(hold) == "97.5 %"] <- "upper"
hold$model <- "vax since sept - model2"
table4 <- rbind(table4, hold) # join onto main table of results
write.csv(table4, "./Vaccine v reinfection paper/sensitivity_analyses.csv")

# Save space
rm(model4a, model4b, model4c, model4d, model4e)


# 4e. Why the social gradient #

# Add in interactions for vaccine and previous infection to IMD
model5 <- glm(outcome ~ vax_v1b + prev_inf_v1 + vax_v1b*prev_inf_v1*imd_score + age_band + factor(Sex) + factor(EthnicMainGroup) + imd_score + test_1mo, data = data2, family = "binomial")

# Marginal effects at the mean
pred3 <- ggpredict(model5, terms = c("imd_score", "vax_v1b", "prev_inf_v1")) # Get margins
# Plot labels
labels <- c("New infection", "Previously tested positive") # Define facet labels
names(labels) <- c("0", "1")
# Plot
plot5 <- ggplot(pred3, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, color = group, fill = group)) + 
  geom_line() + 
  geom_ribbon(alpha = 0.2, colour = NA) +
  facet_wrap(~facet, labeller = labeller(facet = labels)) +
  scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9, labels = c("No", "Yes", "Yes + boosted")) + # Make colour blind friendly (line)
  scale_fill_viridis_d(option = "turbo", begin = 0.1, end = 0.9, labels = c("No", "Yes", "Yes + boosted")) + # Make colour blind friendly (ribbon)
  labs(
    x = "2019 IMD Score (0 = Low Deprivation, 100 = High Deprivation)",
    y = "Predicted probability of COVID-19",
    colour = "Fully vaccinated",
    fill = "Fully vaccinated"
  ) +
  ylim(0,NA) +
  theme(legend.position="bottom")
plot5
ggsave(plot = plot5, file = "./Vaccine v reinfection paper/model_ses.jpeg")


### 5. Spatial analyses ###

# # Tidy data
# lsoa_pos <- data2[, list(positive = sum(sept_pos, na.rm = TRUE), pop = .N), by = c("LSOA_Code", "vax_v1", "prev_inf_v1")] # Aggregate to LSOA
# lsoa_pos$rate <- (lsoa_pos$positive / lsoa_pos$pop)
# lsoa_pos_wd <- dcast(lsoa_pos, LSOA_Code ~ vax_v1 + prev_inf_v1, value.var = "rate", fill = 0)
# 
# # Join to spatial
# lcr_lsoas <- read_sf(dsn = "./Shapefiles/E47000004.shp") # Load LSOA shapefiles - windows
# # lcr_lsoas <- read_sf(dsn = "R:/SMART_LCR_report/Shapefiles/E47000004.shp") # Load LSOA shapefiles - windows
# lcr_lsoas <- merge(lcr_lsoas, lsoa_pos_wd, by.x = "lsoa11cd", by.y = "LSOA_Code", all.x = TRUE)
# lcr_lsoas <- st_transform(lcr_lsoas, 4326) # Set CRS
# lcr_lsoas.nb <- poly2nb(lcr_lsoas, snap=0.0002) # Identify neighbours of each LSOA
# 
# # Define adjacency using a row-standardised matrix
# lw <- nb2listw(lcr_lsoas.nb)
# W <- as(as_dgRMatrix_listw(lw), "CsparseMatrix")
# 
# # Repeat for Getis-Ord statistic
# localGi <- localG(lcr_lsoas$`Vaccinated_New infection`, listw = lw) # Calculate Gi statistics 
# lcr_lsoas <- cbind(lcr_lsoas, data.matrix(localGi)) # Join results onto shapefile
# names(lcr_lsoas)[names(lcr_lsoas) == "data.matrix.localGi."] <- "gstat_vax_new" # Rename column
# 
# localGi <- localG(lcr_lsoas$`Unvaccinated_New.infection`, listw = lw) # Calculate Gi statistics 
# lcr_lsoas <- cbind(lcr_lsoas, data.matrix(localGi)) # Join results onto shapefile
# names(lcr_lsoas)[names(lcr_lsoas) == "data.matrix.localGi."] <- "gstat_unvax_new" # Rename column
# 
# localGi <- localG(lcr_lsoas$`Vaccinated_Reinfection`, listw = lw) # Calculate Gi statistics 
# lcr_lsoas <- cbind(lcr_lsoas, data.matrix(localGi)) # Join results onto shapefile
# names(lcr_lsoas)[names(lcr_lsoas) == "data.matrix.localGi."] <- "gstat_vax_reinf" # Rename column
# 
# localGi <- localG(lcr_lsoas$`Unvaccinated_Reinfection`, listw = lw) # Calculate Gi statistics 
# lcr_lsoas <- cbind(lcr_lsoas, data.matrix(localGi)) # Join results onto shapefile
# names(lcr_lsoas)[names(lcr_lsoas) == "data.matrix.localGi."] <- "gstat_unvax_reinf" # Rename column
# 
# # Make long
# lcr_lsoas_long <- melt(lcr_lsoas, id.vars = c("lsoa11cd", "geometry"), measure = 6:9)
# lcr_lsoas_long = st_as_sf(lcr_lsoas_long, sf_column_name = 'geometry') # Convert back to sf object
# 
# # Time period labels
# labels <- c("Unvaccinated & new infection", "Vaccinated & new infection", "Unvaccinated & reinfection", "Vaccinated & reinfection")
# names(labels) <- c("gstat_unvax_new", "gstat_vax_new", "gstat_unvax_reinf", "gstat_vax_reinf")
# 
# # Map
# map <- ggplot() +
#   geom_sf(data = lcr_lsoas_long, aes(fill = value), lwd = 0) +
#   scale_fill_viridis(option = "turbo") + # Make colour blind friendly
#   facet_wrap(~variable, labeller = labeller(variable = labels)) +
#   xlab("Longitude") + # Add x-axis label
#   ylab("Latitude") + # Add y-axis label
#   labs(title = "COVID-19 incidence since 1st September 2021", # Edit plot title
#        fill = "Gi statistic") # Edit legend title (note must match fill as that is what we are plotting)
# map
# ggsave(plot = map, filename = "/Volumes/MAST/SMART_LCR_report/update outputs/clusters_cases_since_sept.jpeg")
# 
