##############################
### Vaccine vs reinfection ###
##### Survival analyses ######
##############################

# Purpose: To produce a series of survival analyses that are included in the paper.

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
library(survival) # For survival analyses
library(broom)

### 1. Tidy data ###

#  Load data
# load("/Volumes/Sharing Folder/2022-03-02/VaccineAnalysisDatasets.RData") # Mac
# load("/Volumes/Sharing Folder/2022-03-02/VaccineAnalysisDatasets_Data4a.RData") 
# load("/Volumes/Sharing Folder/2022-03-02/VaccineAnalysisDatasets_Data4b.RData") 
load("Q:/2022-03-02/VaccineAnalysisDatasets.RData") # Windows
load("Q:/2022-03-02/VaccineAnalysisDatasets_Data4a.RData")
load("Q:/2022-03-02/VaccineAnalysisDatasets_Data4b.RData")
data4 <- rbind(data4a, data4b) # Combined parts 1 and 2 of data4 (PHE positives)
rm(cohort1, data4a, data4b)
gc()

# 1a. Vaccination data #

# Subset only COVID-19 vaccines (flu jabs are in here for example) and people who completed the vaccine
data1 <- data.table(data1) # Convert data type
data1$t_date <- ymd_hms(data1$DateTimeAdministered) # Convert to time-date
covid_vax <- data1[data1$Status == "completed" & (data1$VaccinationProcedure == "1324681000000101" | data1$VaccinationProcedure == "1324691000000104" | data1$VaccinationProcedure == "1362591000000103")] # Subset
# Vaccine SNOMED codes are in "VaccinationProcedure". 1st dose=1324681000000101, 2nd dose=1324691000000104, Booster=1362591000000103. Completed = people who attended and were received their vaccine (not all receive them after attending)

# Flu codes 822851000000102 (Seasonal influenza vaccination), 884861000000100 (Administration of first intranasal seasonal influenza vaccination), 1037351000000105 (first inactivated seasonal influenza vaccination given by pharmacist) and 1037371000000101 (second inactivated seasonal influenza vaccination given by pharmacist)

# Tidy vaccination data and link to main dataset
covid_vax <- unique(covid_vax, by = c("FK_Patient_Link_ID", "t_date")) # Remove repeated vaccines (just over 60k cases are repeated, so remove one of each)
covid_vax[order(FK_Patient_Link_ID, t_date), num_vax:=1:.N, by=.(FK_Patient_Link_ID)] # Count number of vaccine doses
vax <- covid_vax[num_vax <= 3] # Store upto three vaccine doses only
vaccinated <- dcast(vax, FK_Patient_Link_ID ~ num_vax, value.var = c("t_date", "Manufacturer")) # Reshape to wide format to list when got each vaccine
names(vaccinated) <- c("FK_Patient_Link_ID", "vax1_date", "vax2_date", "vax3_date", "vax1_type", "vax2_type", "vax3_type") # Rename variables to help R
data2 <- merge(data2, vaccinated, by = "FK_Patient_Link_ID", all.x = TRUE) # Join onto main dataset
rm(vax, vaccinated, covid_vax) # Save space

# Recode vax type into mRNA or not
data2$vax1_mRNA <- NA
data2$vax1_mRNA[data2$vax1_type == "AstraZeneca" | data2$vax1_type == "AstraZeneca UK Ltd" | data2$vax1_type == "Janssen-Cilag" | data2$vax1_type == "Johnson & Johnson" | data2$vax1_type == "Oxford Astra Zeneca"] <- 0
data2$vax1_mRNA[data2$vax1_type == "Moderna" | data2$vax1_type == "Moderna, Inc" | data2$vax1_type == "Pfizer" | data2$vax1_type == "Pfizer-BioNTech" | data2$vax1_type == "Pfizer Ltd" | data2$vax1_type == "Pfizer/BioNTech" | data2$vax1_type == "Pfizer-BioNTech COVID-19 mRNA Vaccine BNT162b2" | data2$vax1_type == "Baxter Oncology GmbH"] <- 1

data2$vax2_mRNA <- NA
data2$vax2_mRNA[data2$vax2_type == "AstraZeneca" | data2$vax2_type == "AstraZeneca UK Ltd" | data2$vax2_type == "Janssen-Cilag" | data2$vax2_type == "Johnson & Johnson" | data2$vax2_type == "Oxford Astra Zeneca"] <- 0
data2$vax2_mRNA[data2$vax2_type == "Moderna" | data2$vax2_type == "Moderna, Inc" | data2$vax2_type == "Pfizer" | data2$vax2_type == "Pfizer-BioNTech" | data2$vax2_type == "Pfizer Ltd" | data2$vax2_type == "Pfizer/BioNTech" | data2$vax2_type == "Pfizer-BioNTech COVID-19 mRNA Vaccine BNT162b2" | data2$vax2_type == "Baxter Oncology GmbH"] <- 1

data2$vax3_mRNA <- NA
data2$vax3_mRNA[data2$vax3_type == "AstraZeneca" | data2$vax3_type == "AstraZeneca UK Ltd" | data2$vax3_type == "Janssen-Cilag" | data2$vax3_type == "Johnson & Johnson" | data2$vax3_type == "Oxford Astra Zeneca"] <- 0
data2$vax3_mRNA[data2$vax3_type == "Moderna" | data2$vax3_type == "Moderna, Inc" | data2$vax3_type == "Pfizer" | data2$vax3_type == "Pfizer-BioNTech" | data2$vax3_type == "Pfizer Ltd" | data2$vax3_type == "Pfizer/BioNTech" | data2$vax3_type == "Pfizer-BioNTech COVID-19 mRNA Vaccine BNT162b2" | data2$vax3_type == "Baxter Oncology GmbH"] <- 1

# Redo but for flu vaccinations
data1 <- data1[data1$Status == "completed" & (data1$VaccinationProcedure == "822851000000102" | data1$VaccinationProcedure == "884861000000100" | data1$VaccinationProcedure == "1037351000000105" | data1$VaccinationProcedure == "1037371000000101")] # Subset
# Flu snomed codes 822851000000102 (Seasonal influenza vaccination), 884861000000100 (Administration of first intranasal seasonal influenza vaccination), 1037351000000105 (first inactivated seasonal influenza vaccination given by pharmacist) and 1037371000000101 (second inactivated seasonal influenza vaccination given by pharmacist)
data1 <- unique(data1, by = c("FK_Patient_Link_ID", "t_date")) # Remove repeated vaccines (just over 54k cases are repeated, so remove one of each)
data1[order(FK_Patient_Link_ID, t_date), num_vax:=1:.N, by=.(FK_Patient_Link_ID)] # Count number of vaccine doses
vax <- data1[num_vax <= 4] # Store upto four vaccine doses only (gives us 99.97%)
vaccinated <- dcast(vax, FK_Patient_Link_ID ~ num_vax, value.var = c("t_date")) # Reshape to wide format to list when got each vaccine
names(vaccinated) <- c("FK_Patient_Link_ID", "flu1_date", "flu2_date", "flu3_date", "flu4_date") # Rename variables to help R
data2 <- merge(data2, vaccinated, by = "FK_Patient_Link_ID", all.x = TRUE) # Join onto main dataset
rm(vax, vaccinated, data1) # Save space
gc()


# 1b. Test and Trace data #

# Tidy data
data3 <- data.table(data3) # Convert format
data3$t_date <- ymd_hms(data3$SpecimenDate) # Convert to time-date
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
names(infections) <- c("FK_Patient_Link_ID", "inf1_date", "inf2_date", "inf3_date", "inf4_date") # Rename variables to help R
data2 <- merge(data2, infections, by = "FK_Patient_Link_ID", all.x = TRUE) # Join onto main dataset
rm(hold, data3, infections) # Drop objects not required
gc()

# 1c. Final tidying #

# Drop data from analyses
data2 <- data.table(data2) # Convert format
data2 <- data2[data2$FK_Patient_Link_ID != -1,] # Drop missing person ID (n=21)
data2 <- data2[!is.na(data2$FK_Patient_Link_ID),] # Drop missing person ID (n=2)

# Load and join on IMD decile
# imd <- fread("/Volumes/MAST/SMART_LCR_report/Data/imd2019.csv") # Mac
imd <- fread("R:/SMART_LCR_report/Data/imd2019.csv") # Windows
data2 <- merge(data2, imd, by.x = "LSOA_Code", by.y = "lsoa11", all.x = TRUE) # join together
rm(imd)

# Subset only people registered as living in Cheshire and Merseyside
# cipha <- fread("/Volumes/MAST/SMART_LCR_report/Data/cipha_areas.csv") # Mac
cipha <- fread("R:/SMART_LCR_report/Data/cipha_areas.csv") # Windows
data2 <- merge(data2, cipha, by.x = "LSOA_Code", by.y = "LSOA11CD", all.x = TRUE) # join together
data2 <- data2[data2$cipha_region == 1] # Drop people outside of CIPHA
rm(cipha)

## 2. Generate variables for analysis ##

# 2a. Outcome variables #

# i. Delta period 1 - 3rd June 2021 to 1st September 2021 # [3rd June is when PHE said was 99% dominant] - or 13th May when most cases ~50% according to PHE

# Tested positive (yes/no)
data2$outcome_delta1 <- 0 # Create blank variable
data2$outcome_delta1[(data2$inf1_date >= "2021-06-03 00:00:00.0000000" & data2$inf1_date < "2021-09-01 00:00:00.0000000") | (data2$inf2_date >= "2021-06-03 00:00:00.0000000" & data2$inf2_date < "2021-09-01 00:00:00.0000000") | (data2$inf3_date >= "2021-06-03 00:00:00.0000000" & data2$inf3_date < "2021-09-01 00:00:00.0000000") | (data2$inf4_date >= "2021-06-03 00:00:00.0000000" & data2$inf4_date < "2021-09-01 00:00:00.0000000")] <- 1 # Tested positive over this time period

# Date of positive test
data2$outcome_date_delta1 <- NA # Create blank variable
data2$outcome_date_delta1 <- as.Date(data2$outcome_date_delta1) # Define as date variable
data2$outcome_date_delta1[data2$inf1_date >= "2021-06-03 00:00:00.0000000" & data2$inf1_date < "2021-09-01 00:00:00.0000000" & !is.na(data2$inf1_date)] <- data2$inf1_date[data2$inf1_date >= "2021-06-03 00:00:00.0000000" & data2$inf1_date < "2021-09-01 00:00:00.0000000" & !is.na(data2$inf1_date)] # If have tested positive, then store date
data2$outcome_date_delta1[data2$inf2_date >= "2021-06-03 00:00:00.0000000" & data2$inf2_date < "2021-09-01 00:00:00.0000000" & !is.na(data2$inf2_date)] <- data2$inf2_date[data2$inf2_date >= "2021-06-03 00:00:00.0000000" & data2$inf2_date < "2021-09-01 00:00:00.0000000" & !is.na(data2$inf2_date)]  # Repeat (note people can't test positive more than once in this study period)
data2$outcome_date_delta1[data2$inf3_date >= "2021-06-03 00:00:00.0000000" & data2$inf3_date < "2021-09-01 00:00:00.0000000" & !is.na(data2$inf3_date)] <- data2$inf3_date[data2$inf3_date >= "2021-06-03 00:00:00.0000000" & data2$inf3_date < "2021-09-01 00:00:00.0000000" & !is.na(data2$inf3_date)]  # Repeat
data2$outcome_date_delta1[data2$inf4_date >= "2021-06-03 00:00:00.0000000" & data2$inf4_date < "2021-09-01 00:00:00.0000000" & !is.na(data2$inf4_date)] <- data2$inf4_date[data2$inf4_date >= "2021-06-03 00:00:00.0000000" & data2$inf4_date < "2021-09-01 00:00:00.0000000" & !is.na(data2$inf4_date)]  # Repeat

# ii. Delta period 2 - 1st September 2021 to 27th November 2021 #

# Tested positive
data2$outcome_delta2 <- 0 # Create blank variable
data2$outcome_delta2[(data2$inf1_date >= "2021-09-01 00:00:00.0000000" & data2$inf1_date < "2021-11-27 00:00:00.0000000") | (data2$inf2_date >= "2021-09-01 00:00:00.0000000" & data2$inf2_date < "2021-11-27 00:00:00.0000000") | (data2$inf3_date >= "2021-09-01 00:00:00.0000000" & data2$inf3_date < "2021-11-27 00:00:00.0000000") | (data2$inf4_date >= "2021-09-01 00:00:00.0000000" & data2$inf4_date < "2021-11-27 00:00:00.0000000")] <- 1 # Tested positive over this time period

# Date of positive test
data2$outcome_date_delta2 <- NA # Create variable
data2$outcome_date_delta2 <- as.Date(data2$outcome_date_delta2) # Define as date variable
data2$outcome_date_delta2[data2$inf1_date >= "2021-09-01 00:00:00.0000000" & data2$inf1_date < "2021-11-27 00:00:00.0000000" & !is.na(data2$inf1_date)] <- data2$inf1_date[data2$inf1_date >= "2021-09-01 00:00:00.0000000" & data2$inf1_date < "2021-11-27 00:00:00.0000000" & !is.na(data2$inf1_date)] # If have tested positive, then store date
data2$outcome_date_delta2[data2$inf2_date >= "2021-09-01 00:00:00.0000000" & data2$inf2_date < "2021-11-27 00:00:00.0000000" & !is.na(data2$inf2_date)] <- data2$inf2_date[data2$inf2_date >= "2021-09-01 00:00:00.0000000" & data2$inf2_date < "2021-11-27 00:00:00.0000000" & !is.na(data2$inf2_date)]  # Repeat (note people can't test positive more than once in this study period)
data2$outcome_date_delta2[data2$inf3_date >= "2021-09-01 00:00:00.0000000" & data2$inf3_date < "2021-11-27 00:00:00.0000000" & !is.na(data2$inf3_date)] <- data2$inf3_date[data2$inf3_date >= "2021-09-01 00:00:00.0000000" & data2$inf3_date < "2021-11-27 00:00:00.0000000" & !is.na(data2$inf3_date)]  # Repeat
data2$outcome_date_delta2[data2$inf4_date >= "2021-09-01 00:00:00.0000000" & data2$inf4_date < "2021-11-27 00:00:00.0000000" & !is.na(data2$inf4_date)] <- data2$inf4_date[data2$inf4_date >= "2021-09-01 00:00:00.0000000" & data2$inf4_date < "2021-11-27 00:00:00.0000000" & !is.na(data2$inf4_date)]  # Repeat

# Alpha - 8th December 2020 - 22 February 2021 (first person to be vaccinated to when delta first detected in UK - 22 February 2021, although delta not majority until 13th May]

# iii. Omicron - 13th December 2021 to 31st January 2022 

# Tested positive 
data2$outcome_omni <- 0 # Create blank variable
data2$outcome_omni[(data2$inf1_date >= "2021-12-13 00:00:00.0000000" & data2$inf1_date < "2022-03-01 00:00:00.0000000") | (data2$inf2_date >= "2021-12-13 00:00:00.0000000" & data2$inf2_date < "2022-03-01 00:00:00.0000000") | (data2$inf3_date >= "2021-12-13 00:00:00.0000000" & data2$inf3_date < "2022-03-01 00:00:00.0000000") | (data2$inf4_date >= "2021-12-13 00:00:00.0000000" & data2$inf4_date < "2022-03-01 00:00:00.0000000")] <- 1 # Tested positive over this time period

# Date of positive test
data2$outcome_date_omni <- NA # Create variable
data2$outcome_date_omni <- as.Date(data2$outcome_date_omni) # Define as date variable
data2$outcome_date_omni[data2$inf1_date >= "2021-12-13 00:00:00.0000000" & data2$inf1_date < "2022-03-01 00:00:00.0000000" & !is.na(data2$inf1_date)] <- data2$inf1_date[data2$inf1_date >= "2021-12-13 00:00:00.0000000" & data2$inf1_date < "2022-03-01 00:00:00.0000000" & !is.na(data2$inf1_date)] # If have tested positive, then store date
data2$outcome_date_omni[data2$inf2_date >= "2021-12-13 00:00:00.0000000" & data2$inf2_date < "2022-03-01 00:00:00.0000000" & !is.na(data2$inf2_date)] <- data2$inf2_date[data2$inf2_date >= "2021-12-13 00:00:00.0000000" & data2$inf2_date < "2022-03-01 00:00:00.0000000" & !is.na(data2$inf2_date)]  # Repeat (note people can't test positive more than once in this study period)
data2$outcome_date_omni[data2$inf3_date >= "2021-12-13 00:00:00.0000000" & data2$inf3_date < "2022-03-01 00:00:00.0000000" & !is.na(data2$inf3_date)] <- data2$inf3_date[data2$inf3_date >= "2021-12-13 00:00:00.0000000" & data2$inf3_date < "2022-03-01 00:00:00.0000000" & !is.na(data2$inf3_date)]  # Repeat
data2$outcome_date_omni[data2$inf4_date >= "2021-12-13 00:00:00.0000000" & data2$inf4_date < "2022-03-01 00:00:00.0000000" & !is.na(data2$inf4_date)] <- data2$inf4_date[data2$inf4_date >= "2021-12-13 00:00:00.0000000" & data2$inf4_date < "2022-03-01 00:00:00.0000000" & !is.na(data2$inf4_date)]  # Repeat


# 2b. Vaccination status #

# i. Delta period 1

# Fully vaccinated before baseline (3rd June 2021) - 2 weeks for 1 dose, 1 week for 2/3 doses
data2$vax_delta1 <- 0 # Create blank variable
data2$vax_delta1[data2$vax1_date < "2021-05-20 00:00:00.0000000"] <- 1 # First dose
data2$vax_delta1[data2$vax2_date < "2021-05-27 00:00:00.0000000"] <- 2 # Second dose
data2$vax_delta1[data2$vax3_date < "2021-05-27 00:00:00.0000000"] <- 3 # Define if boosted
data2$vax_delta1 <- factor(data2$vax_delta1) # Store as factor

# Vax type
data2$vax_type_delta1 <- 0 # Unvaccinated
data2$vax_type_delta1[data2$vax_delta1 == 1 & data2$vax1_mRNA == 0] <- 1 # 1 dose not mRNA
data2$vax_type_delta1[data2$vax_delta1 == 1 & data2$vax1_mRNA == 1] <- 2 # 1 dose mRNA
data2$vax_type_delta1[data2$vax_delta1 == 2 & (data2$vax1_mRNA == 0 & data2$vax2_mRNA == 0)] <- 3 # 2 doses both not mRNA
data2$vax_type_delta1[data2$vax_delta1 == 2 & (data2$vax1_mRNA == 1 & data2$vax2_mRNA == 1)] <- 4 # 2 doses both mRNA
data2$vax_type_delta1[data2$vax_delta1 == 2 & ((data2$vax1_mRNA == 0 & data2$vax2_mRNA == 1) | (data2$vax1_mRNA == 1 & data2$vax2_mRNA == 0))] <- 5 # 2 doses mixture
data2$vax_type_delta1[data2$vax_delta1 == 3 & data2$vax1_mRNA == 0 & data2$vax2_mRNA == 0 & data2$vax3_mRNA == 0] <- 6 # 3 doses all not mRNA
data2$vax_type_delta1[data2$vax_delta1 == 3 & data2$vax1_mRNA == 1 & data2$vax2_mRNA == 1 & data2$vax3_mRNA == 1] <- 7 # 3 does all mNA
data2$vax_type_delta1[data2$vax_delta1 == 3 & data2$vax1_mRNA == 0 & data2$vax2_mRNA == 0 & data2$vax3_mRNA == 1] <- 8 # 3 doses where two does not mRNA and third dose mRNA
data2$vax_type_delta1[data2$vax_delta1 == 3 & data2$vax_type_delta1 != 6 & data2$vax_type_delta1 != 7 & data2$vax_type_delta1 != 8] <- 9 # 3 doses mixture (other)
data2$vax_type_delta1[data2$vax_delta1 == 1 & is.na(data2$vax1_mRNA)] <- NA # Seen some data errors that need fixing
data2$vax_type_delta1 <- factor(data2$vax_type_delta1) # Store as factor

# Received flu vax 1 year before start of period 
data2$flu_vax_delta1 <- 0 # Create blank variable
data2$flu_vax_delta1[(data2$flu1_date > "2020-06-03 00:00:00.0000000" & data2$flu1_date < "2021-06-03 00:00:00.0000000") | (data2$flu2_date > "2020-06-03 00:00:00.0000000" & data2$flu2_date < "2021-06-03 00:00:00.0000000") | (data2$flu3_date > "2020-06-03 00:00:00.0000000" & data2$flu3_date < "2021-06-03 00:00:00.0000000") | (data2$flu4_date > "2020-06-03 00:00:00.0000000" & data2$flu4_date < "2021-06-03 00:00:00.0000000")] <- 1

# ii. Delta period 2

# Fully vaccinated 3 weeks before baseline (1st September 2021) or 2 weeks for booster
data2$vax_delta2 <- 0 # Create blank variable
data2$vax_delta2[data2$vax1_date < "2021-08-18 00:00:00.0000000"] <- 1 # First dose
data2$vax_delta2[data2$vax2_date < "2021-08-25 00:00:00.0000000"] <- 2 # Second dose
data2$vax_delta2[data2$vax3_date < "2021-08-25 00:00:00.0000000"] <- 3 # Define if boosted
data2$vax_delta2 <- factor(data2$vax_delta2) # Store as factor

# Vax type
data2$vax_type_delta2 <- 0 # Unvaccinated
data2$vax_type_delta2[data2$vax_delta2 == 1 & data2$vax1_mRNA == 0] <- 1 # 1 dose not mRNA
data2$vax_type_delta2[data2$vax_delta2 == 1 & data2$vax1_mRNA == 1] <- 2 # 1 dose mRNA
data2$vax_type_delta2[data2$vax_delta2 == 2 & (data2$vax1_mRNA == 0 & data2$vax2_mRNA == 0)] <- 3 # 2 doses both not mRNA
data2$vax_type_delta2[data2$vax_delta2 == 2 & (data2$vax1_mRNA == 1 & data2$vax2_mRNA == 1)] <- 4 # 2 doses both mRNA
data2$vax_type_delta2[data2$vax_delta2 == 2 & ((data2$vax1_mRNA == 0 & data2$vax2_mRNA == 1) | (data2$vax1_mRNA == 1 & data2$vax2_mRNA == 0))] <- 5 # 2 doses mixture
data2$vax_type_delta2[data2$vax_delta2 == 3 & data2$vax1_mRNA == 0 & data2$vax2_mRNA == 0 & data2$vax3_mRNA == 0] <- 6 # 3 doses all not mRNA
data2$vax_type_delta2[data2$vax_delta2 == 3 & data2$vax1_mRNA == 1 & data2$vax2_mRNA == 1 & data2$vax3_mRNA == 1] <- 7 # 3 does all mNA
data2$vax_type_delta2[data2$vax_delta2 == 3 & data2$vax1_mRNA == 0 & data2$vax2_mRNA == 0 & data2$vax3_mRNA == 1] <- 8 # 3 doses where two does not mRNA and third dose mRNA
data2$vax_type_delta2[data2$vax_delta2 == 3 & data2$vax_type_delta2 != 6 & data2$vax_type_delta2 != 7 & data2$vax_type_delta2 != 8] <- 9 # 3 doses mixture (other)
data2$vax_type_delta2[data2$vax_delta2 == 1 & is.na(data2$vax1_mRNA)] <- NA # Seen some data errors that need fixing
data2$vax_type_delta2 <- factor(data2$vax_type_delta2) # Store as factor

# Received flu vax 1 year before start of period 
data2$flu_vax_delta2 <- 0 # Create blank variable
data2$flu_vax_delta2[(data2$flu1_date > "2020-09-01 00:00:00.0000000" & data2$flu1_date < "2021-09-01 00:00:00.0000000") | (data2$flu2_date > "2020-09-01 00:00:00.0000000" & data2$flu2_date < "2021-09-01 00:00:00.0000000") | (data2$flu3_date > "2020-09-01 00:00:00.0000000" & data2$flu3_date < "2021-09-01 00:00:00.0000000") | (data2$flu4_date > "2020-09-01 00:00:00.0000000" & data2$flu4_date < "2021-09-01 00:00:00.0000000")] <- 1

# iii. Omnicron

# Fully vaccinated 3 weeks before baseline (13th December 2021) or 2 weeks for booster
data2$vax_omni <- 0 # Create blank variable
data2$vax_omni[data2$vax1_date < "2021-11-29 00:00:00.0000000"] <- 1 # First dose
data2$vax_omni[data2$vax2_date < "2021-12-06 00:00:00.0000000"] <- 2 # Second dose
data2$vax_omni[data2$vax3_date < "2021-12-06 00:00:00.0000000"] <- 3 # Define if boosted
data2$vax_omni <- factor(data2$vax_omni) # Store as factor

# Vax type
data2$vax_type_omni <- 0 # Unvaccinated
data2$vax_type_omni[data2$vax_omni == 1 & data2$vax1_mRNA == 0] <- 1 # 1 dose not mRNA
data2$vax_type_omni[data2$vax_omni == 1 & data2$vax1_mRNA == 1] <- 2 # 1 dose mRNA
data2$vax_type_omni[data2$vax_omni == 2 & (data2$vax1_mRNA == 0 & data2$vax2_mRNA == 0)] <- 3 # 2 doses both not mRNA
data2$vax_type_omni[data2$vax_omni == 2 & (data2$vax1_mRNA == 1 & data2$vax2_mRNA == 1)] <- 4 # 2 doses both mRNA
data2$vax_type_omni[data2$vax_omni == 2 & ((data2$vax1_mRNA == 0 & data2$vax2_mRNA == 1) | (data2$vax1_mRNA == 1 & data2$vax2_mRNA == 0))] <- 5 # 2 doses mixture
data2$vax_type_omni[data2$vax_omni == 3 & data2$vax1_mRNA == 0 & data2$vax2_mRNA == 0 & data2$vax3_mRNA == 0] <- 6 # 3 doses all not mRNA
data2$vax_type_omni[data2$vax_omni == 3 & data2$vax1_mRNA == 1 & data2$vax2_mRNA == 1 & data2$vax3_mRNA == 1] <- 7 # 3 does all mNA
data2$vax_type_omni[data2$vax_omni == 3 & data2$vax1_mRNA == 0 & data2$vax2_mRNA == 0 & data2$vax3_mRNA == 1] <- 8 # 3 doses where two does not mRNA and third dose mRNA
data2$vax_type_omni[data2$vax_omni == 3 & data2$vax_type_omni != 6 & data2$vax_type_omni != 7 & data2$vax_type_omni != 8] <- 9 # 3 doses mixture (other)
data2$vax_type_omni[data2$vax_omni == 1 & is.na(data2$vax1_mRNA)] <- NA # Seen some data errors that need fixing
data2$vax_type_omni <- factor(data2$vax_type_omni) # Store as factor

# Received flu vax 1 year before start of period 
data2$flu_vax_omni <- 0 # Create blank variable
data2$flu_vax_omni[(data2$flu1_date > "2020-12-13 00:00:00.0000000" & data2$flu1_date < "2021-12-13 00:00:00.0000000") | (data2$flu2_date > "2020-12-13 00:00:00.0000000" & data2$flu2_date < "2021-12-13 00:00:00.0000000") | (data2$flu3_date > "2020-12-13 00:00:00.0000000" & data2$flu3_date < "2021-12-13 00:00:00.0000000") | (data2$flu4_date > "2020-12-13 00:00:00.0000000" & data2$flu4_date < "2021-12-13 00:00:00.0000000")] <- 1

# 2c. Previous infection #

# i. Delta period 1

# Previous infection 2 weeks before baseline
data2$prev_inf_delta1 <- 0 # Create blank variable
data2$prev_inf_delta1[data2$inf1_date < "2021-05-20 00:00:00.0000000" | data2$inf2_date < "2021-05-20 00:00:00.0000000" | data2$inf3_date < "2021-05-20 00:00:00.0000000" | data2$inf4_date < "2021-05-20 00:00:00.0000000"] <- 1 # Tested positive in this time period

# Immortal time bias - sensitivity analysis (i.e., only keep people who could have tested positive in the entire period - so previous positive 90 days before 3rd June 2021)
data2$immortal_delta1 <- 0
data2$immortal_delta1[data2$inf1_date < "2021-03-05 00:00:00.0000000" | data2$inf2_date < "2021-03-05 00:00:00.0000000" | data2$inf3_date < "2021-03-05 00:00:00.0000000" | data2$inf4_date < "2021-03-05 00:00:00.0000000"] <- 1 # Test positive before cut off
data2$immortal_delta1[is.na(data2$inf1_date) | data2$inf1_date >= "2021-06-03 00:00:00.0000000"] <- 1 # Never tested positive or tested positive for first time in outcome period

# ii. Delta period 2

# Previous infection 2 weeks before baseline
data2$prev_inf_delta2 <- 0 # Create blank variable
data2$prev_inf_delta2[data2$inf1_date < "2021-08-18 00:00:00.0000000" | data2$inf2_date < "2021-08-18 00:00:00.0000000" | data2$inf3_date < "2021-08-18 00:00:00.0000000" | data2$inf4_date < "2021-08-18 00:00:00.0000000"] <- 1 # Tested positive in this time period

# Immortal time bias - repeat
data2$immortal_delta2 <- 0
data2$immortal_delta2[data2$inf1_date < "2021-06-03 00:00:00.0000000" | data2$inf2_date < "2021-06-03 00:00:00.0000000" | data2$inf3_date < "2021-06-03 00:00:00.0000000" | data2$inf4_date < "2021-06-03 00:00:00.0000000"] <- 1
data2$immortal_delta2[is.na(data2$inf1_date) | data2$inf1_date >= "2021-09-01 00:00:00.0000000"] <- 1 # Never tested positive or tested positive for first time in outcome period

# ii. Omnicron

# Previous infection 2 weeks before baseline
data2$prev_inf_omni <- 0 # Create blank variable
data2$prev_inf_omni[data2$inf1_date < "2021-11-29 00:00:00.0000000" | data2$inf2_date < "2021-11-29 00:00:00.0000000" | data2$inf3_date < "2021-11-29 00:00:00.0000000" | data2$inf4_date < "2021-11-29 00:00:00.0000000"] <- 1 # Tested positive in this time period

# Immortal time bias - repeat
data2$immortal_omni <- 0
data2$immortal_omni[data2$inf1_date < "2021-09-14 00:00:00.0000000" | data2$inf2_date < "2021-09-14 00:00:00.0000000" | data2$inf3_date < "2021-09-14 00:00:00.0000000" | data2$inf4_date < "2021-09-14 00:00:00.0000000"] <- 1
data2$immortal_omni[is.na(data2$inf1_date) | data2$inf1_date >= "2021-12-14 00:00:00.0000000"] <- 1 # Never tested positive or tested positive for first time in outcome period


# 2d. Control variables #


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

# Sex
data2$Sex[data2$Sex == "U"] <- NA # Drop missing sex (n=117))

# Ethnicity
# What to do with NULL ethnicity?
data2$EthnicMainGroup <- as.factor(data2$EthnicMainGroup) # Set as factor
data2$EthnicMainGroup <- relevel(data2$EthnicMainGroup, ref = "White") # Define White as reference

# Health status
data2$health_issue <- 1 # Create blank variable (people with health issues on EDC)
data2$health_issue[data2$EDC_Codes == "NULL"] <- 0 # Has issue defined
# Expanded Diagnosis Clusters codes are all of the EDC codes assigned to this patient, separated by spaces. The EDC taxonomy identifies patients with specific diseases or symptoms that are treated in ambulatory and inpatient settings.


# Previous number of tests

# Subset tests by time period
delta1_tests <- data.table(data4[(data4$SpecimenDate < "2021-06-03" & data4$SpecimenDate >= "2021-05-03"),]) # 1 month before delta period 1 baseline (3rd June 2021))
delta2_tests <- data.table(data4[(data4$SpecimenDate < "2021-09-01" & data4$SpecimenDate >= "2021-08-01"),]) # 1 month before delta period 2 baseline (1st September 2021)
omni_tests <- data.table(data4[(data4$SpecimenDate < "2021-12-13" & data4$SpecimenDate >= "2021-11-13"),]) # 1 month before omnicron period baseline (13th December 2021)

# Calculate number of tests per ID/person
delta1_tests <- delta1_tests[, list(n_tests_delta1 = .N), by = "FK_Patient_Link_ID"]
delta2_tests <- delta2_tests[, list(n_tests_delta2 = .N), by = "FK_Patient_Link_ID"]
omni_tests <- omni_tests[, list(n_tests_omni = .N), by = "FK_Patient_Link_ID"]

# Join onto main dataset
data2 <- merge(data2, delta1_tests, by = "FK_Patient_Link_ID", all.x = TRUE) # Join tests onto main spline
data2$n_tests_delta1[is.na(data2$n_tests_delta1)] <- 0  # If missing, then had no tests so reclassify
data2$n_tests_delta1[data2$n_tests_delta1 > 1] <- 1 # Recode variable for if tested or not
data2 <- merge(data2, delta2_tests, by = "FK_Patient_Link_ID", all.x = TRUE) # Repeat
data2$n_tests_delta2[is.na(data2$n_tests_delta2)] <- 0
data2$n_tests_delta2[data2$n_tests_delta2 > 1] <- 1
data2 <- merge(data2, omni_tests, by = "FK_Patient_Link_ID", all.x = TRUE) # Repeat
data2$n_tests_omni[is.na(data2$n_tests_omni)] <- 0
data2$n_tests_omni[data2$n_tests_omni > 1] <- 1
rm(delta1_tests, delta2_tests, omni_tests)


# Number of tests in period

# Subset tests by time period
delta1_tests <- data.table(data4[(data4$SpecimenDate >= "2021-06-03" & data4$SpecimenDate < "2021-09-01"),]) # 1 month before delta period 1 baseline (3rd June 2021))
delta2_tests <- data.table(data4[(data4$SpecimenDate >= "2021-09-01" & data4$SpecimenDate < "2021-11-27"),]) # 1 month before delta period 2 baseline (1st September 2021)
omni_tests <- data.table(data4[(data4$SpecimenDate >= "2021-12-13" & data4$SpecimenDate < "2022-03-01"),]) # 1 month before omnicron period baseline (13th December 2021)
#rm(data4)

# Calculate number of tests per ID/person
delta1_tests <- delta1_tests[, list(n_tests_post_delta1 = .N), by = "FK_Patient_Link_ID"]
delta2_tests <- delta2_tests[, list(n_tests_post_delta2 = .N), by = "FK_Patient_Link_ID"]
omni_tests <- omni_tests[, list(n_tests_post_omni = .N), by = "FK_Patient_Link_ID"]

# Join onto main dataset
data2 <- merge(data2, delta1_tests, by = "FK_Patient_Link_ID", all.x = TRUE) # Join tests onto main spline
data2$n_tests_post_delta1[is.na(data2$n_tests_post_delta1)] <- 0  # If missing, then had no tests so reclassify
data2$n_tests_post_delta1[data2$n_tests_post_delta1 > 1] <- 1 # Recode variable for if tested or not
data2 <- merge(data2, delta2_tests, by = "FK_Patient_Link_ID", all.x = TRUE) # Repeat
data2$n_tests_post_delta2[is.na(data2$n_tests_post_delta2)] <- 0
data2$n_tests_post_delta2[data2$n_tests_post_delta2 > 1] <- 1
data2 <- merge(data2, omni_tests, by = "FK_Patient_Link_ID", all.x = TRUE) # Repeat
data2$n_tests_post_omni[is.na(data2$n_tests_post_omni)] <- 0
data2$n_tests_post_omni[data2$n_tests_post_omni > 1] <- 1
rm(delta1_tests, delta2_tests, omni_tests, data4) # Tidy


# 2e. Reshape data into survival analysis format #

# i. Delta period 1

# Calculate total time period for exposure
data2$full_time_delta1 <- interval(ymd_hms("2021-06-03 00:00:00"), ymd_hms("2021-09-01 00:00:00")) / days(1)

# Count number of days to infection (if infected, else NA)
data2$days_to_infection_delta1 <- NA # Create variable and set as missing to represent those not infected
data2$days_to_infection_delta1 <- interval(ymd_hms("2021-06-03 00:00:00"), data2$outcome_date_delta1) / days(1) # Number of days to infection
data2$days_to_infection_delta1[data2$days_to_infection_delta1 > data2$full_time_delta1] <- NA # Vaccinated after period of analysis

# Days to first vaccination (of received a vaccination post baseline period)
data2$days_to_vax1_delta1 <- NA # As above
data2$days_to_vax1_delta1 <- interval(ymd_hms("2021-06-03 00:00:00"), data2$vax1_date) / days(1) # Number of days to vaccination
data2$days_to_vax1_delta1[data2$days_to_vax1_delta1 < 0] <- NA # If <0 then record as missing as captured at baseline 
data2$days_to_vax1_delta1[data2$days_to_vax1_delta1 > data2$full_time_delta1] <- NA # Vaccinated after period of analysis

# Days to second vaccination (of received a vaccination post baseline period)
data2$days_to_vax2_delta1 <- NA # As above
data2$days_to_vax2_delta1 <- interval(ymd_hms("2021-06-03 00:00:00"), data2$vax2_date) / days(1) # Number of days to vaccination
data2$days_to_vax2_delta1[data2$days_to_vax2_delta1 < 0] <- NA # If <0 then record as missing as captured at baseline 
data2$days_to_vax2_delta1[data2$days_to_vax2_delta1 > data2$full_time_delta1] <- NA # Vaccinated after period of analysis

# Days to third vaccination (of received a vaccination post baseline period)
data2$days_to_vax3_delta1 <- NA # As above
data2$days_to_vax3_delta1 <- interval(ymd_hms("2021-06-03 00:00:00"), data2$vax3_date) / days(1) # Number of days to vaccination
data2$days_to_vax3_delta1[data2$days_to_vax3_delta1 < 0] <- NA # If <0 then record as missing as captured at baseline 
data2$days_to_vax3_delta1[data2$days_to_vax3_delta1 > data2$full_time_delta1] <- NA # Vaccinated after period of analysis


# Reshape data
survival_delta1 <- tmerge(data1 = data2[, c("FK_Patient_Link_ID", "Patient.Age", "age_band", "Sex", "EthnicMainGroup", "imd_score", "imd_decile", "n_tests_delta1", "full_time_delta1", "vax_delta1", "vax1_date", "vax2_date", "vax3_date", "vax_type_delta1", "outcome_date_delta1", "prev_inf_delta1", "inf1_date", "vax1_mRNA", "vax2_mRNA", "vax3_mRNA", "Deceased", "DeathDate", "health_issue", "immortal_delta1", "n_tests_post_delta1", "flu_vax_delta1")], data2 = data2, id = FK_Patient_Link_ID, tstop = full_time_delta1) # Calculate each persons total days in study (to capture people not infected)
survival_delta1 <- tmerge(data1 = survival_delta1, data2 = data2, id = FK_Patient_Link_ID, infection = event(days_to_infection_delta1)) # Add in time period pre- and -infection
survival_delta1 <- tmerge(data1 = survival_delta1, data2 = data2, id = FK_Patient_Link_ID, vax1_new = tdc(days_to_vax1_delta1)) # Time varying vaccination status - first dose
survival_delta1 <- tmerge(data1 = survival_delta1, data2 = data2, id = FK_Patient_Link_ID, vax2_new = tdc(days_to_vax2_delta1)) # Time varying vaccination status - second dose
survival_delta1 <- tmerge(data1 = survival_delta1, data2 = data2, id = FK_Patient_Link_ID, vax3_new = tdc(days_to_vax3_delta1)) # Time varying vaccination status - third dose
survival_delta1 <- tmerge(data1 = survival_delta1, data2 = survival_delta1, id = FK_Patient_Link_ID, n_episodes = cumtdc(tstart)) # Create number of period per person

# Recode vaccination data so time-varying
survival_delta1$days_vax1_outcome <- interval(survival_delta1$vax1_date, survival_delta1$outcome_date_delta1) / days(1) # Number of days between first dose and infection - will be positive if infection is after vax, negative if vax after infection 
survival_delta1$vax_delta1[(survival_delta1$vax1_new == 1) & (survival_delta1$infection == 1) & survival_delta1$days_vax1_outcome >= 14] <- 1 # If received their second vaccine dose two weeks before infection, then recode vaccination status - if infected
survival_delta1$vax_delta1[(survival_delta1$vax1_new == 1) & (survival_delta1$infection == 0)] <- 1 # If received their second vaccine and not infected

survival_delta1$days_vax2_outcome <- interval(survival_delta1$vax2_date, survival_delta1$outcome_date_delta1) / days(1) # Number of days between second dose and infection
survival_delta1$vax_delta1[(survival_delta1$vax2_new == 1) & (survival_delta1$infection == 1) & (survival_delta1$days_vax2_outcome >= 7)] <- 2 # If received their second vaccine dose one week before infection, then recode vaccination status - if infected
survival_delta1$vax_delta1[(survival_delta1$vax2_new == 1) & (survival_delta1$infection == 0)] <- 2 # If received their second vaccine and not infected

survival_delta1$days_vax3_outcome <- interval(survival_delta1$vax3_date, survival_delta1$outcome_date_delta1) / days(1) # Number of days between third dose and infection
survival_delta1$vax_delta1[(survival_delta1$vax3_new == 1) & (survival_delta1$infection == 1) & (survival_delta1$days_vax3_outcome >= 7)] <- 3 # If received their second vaccine dose one week before infection, then recode vaccination status - if infected
survival_delta1$vax_delta1[(survival_delta1$vax3_new == 1) & (survival_delta1$infection == 0)] <- 3 # If received their second vaccine and not infected

# Recode vaccination type so time-varying (sorry not sorry I have not commented every line - lazy)
survival_delta1$vax_type_delta1[(survival_delta1$vax1_new == 1 & survival_delta1$infection == 1 & survival_delta1$days_vax1_outcome >= 14 & survival_delta1$vax1_mRNA == 0)] <- 1
survival_delta1$vax_type_delta1[(survival_delta1$vax1_new == 1 & survival_delta1$infection == 0 & survival_delta1$vax1_mRNA == 0)] <- 1

survival_delta1$vax_type_delta1[(survival_delta1$vax1_new == 1 & survival_delta1$infection == 1 & survival_delta1$days_vax1_outcome >= 14 & survival_delta1$vax1_mRNA == 1)] <- 2 
survival_delta1$vax_type_delta1[(survival_delta1$vax1_new == 1 & survival_delta1$infection == 0 & survival_delta1$vax1_mRNA == 1)] <- 2 

survival_delta1$vax_type_delta1[(survival_delta1$vax2_new == 1 & survival_delta1$infection == 1 & survival_delta1$days_vax2_outcome >= 7 & (survival_delta1$vax1_mRNA == 0 & survival_delta1$vax2_mRNA == 0))] <- 3
survival_delta1$vax_type_delta1[(survival_delta1$vax2_new == 1 & survival_delta1$infection == 0 & (survival_delta1$vax1_mRNA == 0 & survival_delta1$vax2_mRNA == 0))] <- 3

survival_delta1$vax_type_delta1[(survival_delta1$vax2_new == 1 & survival_delta1$infection == 1 & survival_delta1$days_vax2_outcome >= 7 & (survival_delta1$vax1_mRNA == 1 & survival_delta1$vax2_mRNA == 1))] <- 4
survival_delta1$vax_type_delta1[(survival_delta1$vax2_new == 1 & survival_delta1$infection == 0 & (survival_delta1$vax1_mRNA == 1 & survival_delta1$vax2_mRNA == 1))] <- 4

survival_delta1$vax_type_delta1[(survival_delta1$vax2_new == 1 & survival_delta1$infection == 1 & survival_delta1$days_vax2_outcome >= 7 & ((survival_delta1$vax1_mRNA == 0 & survival_delta1$vax2_mRNA == 1) | (survival_delta1$vax1_mRNA == 1 & survival_delta1$vax2_mRNA == 0)))] <- 5
survival_delta1$vax_type_delta1[(survival_delta1$vax2_new == 1 & survival_delta1$infection == 0 & ((survival_delta1$vax1_mRNA == 0 & survival_delta1$vax2_mRNA == 1) | (survival_delta1$vax1_mRNA == 1 & survival_delta1$vax2_mRNA == 0)))] <- 5

survival_delta1$vax_type_delta1[(survival_delta1$vax3_new == 1 & survival_delta1$infection == 1 & survival_delta1$days_vax3_outcome >= 7 & (survival_delta1$vax1_mRNA == 0 & survival_delta1$vax2_mRNA == 0 & survival_delta1$vax3_mRNA == 0))] <- 6
survival_delta1$vax_type_delta1[(survival_delta1$vax3_new == 1 & survival_delta1$infection == 0 & (survival_delta1$vax1_mRNA == 0 & survival_delta1$vax2_mRNA == 0 & survival_delta1$vax3_mRNA == 0))] <- 6

survival_delta1$vax_type_delta1[(survival_delta1$vax3_new == 1 & survival_delta1$infection == 1 & survival_delta1$days_vax3_outcome >= 7 & (survival_delta1$vax1_mRNA == 1 & survival_delta1$vax2_mRNA == 1 & survival_delta1$vax3_mRNA == 1))] <- 7
survival_delta1$vax_type_delta1[(survival_delta1$vax3_new == 1 & survival_delta1$infection == 0 & (survival_delta1$vax1_mRNA == 1 & survival_delta1$vax2_mRNA == 1 & survival_delta1$vax3_mRNA == 1))] <- 7

survival_delta1$vax_type_delta1[(survival_delta1$vax3_new == 1 & survival_delta1$infection == 1 & survival_delta1$days_vax3_outcome >= 7 & (survival_delta1$vax1_mRNA == 0 & survival_delta1$vax2_mRNA == 0 & survival_delta1$vax3_mRNA == 1))] <- 8
survival_delta1$vax_type_delta1[(survival_delta1$vax3_new == 1 & survival_delta1$infection == 0 & (survival_delta1$vax1_mRNA == 0 & survival_delta1$vax2_mRNA == 0 & survival_delta1$vax3_mRNA == 1))] <- 8

survival_delta1$vax_type_delta1[(survival_delta1$vax3_new == 1 & survival_delta1$infection == 1 & survival_delta1$days_vax3_outcome >= 7 & (survival_delta1$vax_type_delta1 != 6 & survival_delta1$vax_type_delta1 != 7 & survival_delta1$vax_type_delta1 != 8))] <- 9
survival_delta1$vax_type_delta1[(survival_delta1$vax3_new == 1 & survival_delta1$infection == 0 & (survival_delta1$vax_type_delta1 != 6 & survival_delta1$vax_type_delta1 != 7 & survival_delta1$vax_type_delta1 != 8))] <- 9

# # Recode previous infection so time-varying
# survival_delta1$days_previnf_outcome <- interval(survival_delta1$inf1_date, survival_delta1$outcome_date_delta1) / days(1) # Number of days between first infection and current infection
# survival_delta1$prev_inf_delta1[(survival_delta1$prev_inf_delta1 == 0) & (survival_delta1$days_previnf_outcome >= 90)] <- 1 # Recode status if not previous infection and current infection greater than 90 days



# ii. Delta period 2

# Calculate total time period for exposure
data2$full_time_delta2 <- interval(ymd_hms("2021-09-01 00:00:00"), ymd_hms("2021-11-27 00:00:00")) / days(1)

# Count number of days to infection from 1st September (if infected, else NA)
data2$days_to_infection_delta2 <- NA # Create variable and set as missing to represent those not infected
data2$days_to_infection_delta2 <- interval(ymd_hms("2021-09-01 00:00:00"), data2$outcome_date_delta2) / days(1) # Number of days to infection
data2$days_to_infection_delta2[data2$days_to_infection_delta2 > data2$full_time_omni] <- NA # Vaccinated after period of analysis

# Days to first vaccination (of received a vaccination post baseline period)
data2$days_to_vax1_delta2 <- NA # As above
data2$days_to_vax1_delta2 <- interval(ymd_hms("2021-09-01 00:00:00"), data2$vax1_date) / days(1) # Number of days to vaccination
data2$days_to_vax1_delta2[data2$days_to_vax1_delta2 < 0] <- NA # If <0 then record as missing as captured at baseline 
data2$days_to_vax1_delta2[data2$days_to_vax1_delta2 > data2$full_time_delta2] <- NA # Vaccinated after period of analysis

# Days to second vaccination (of received a vaccination post baseline period)
data2$days_to_vax2_delta2 <- NA # As above
data2$days_to_vax2_delta2 <- interval(ymd_hms("2021-09-01 00:00:00"), data2$vax2_date) / days(1) # Number of days to vaccination
data2$days_to_vax2_delta2[data2$days_to_vax2_delta2 < 0] <- NA # If <0 then record as missing as captured at baseline 
data2$days_to_vax2_delta2[data2$days_to_vax2_delta2 > data2$full_time_delta2] <- NA # Vaccinated after period of analysis

# Days to third vaccination (of received a vaccination post baseline period)
data2$days_to_vax3_delta2 <- NA # As above
data2$days_to_vax3_delta2 <- interval(ymd_hms("2021-09-01 00:00:00"), data2$vax3_date) / days(1) # Number of days to vaccination
data2$days_to_vax3_delta2[data2$days_to_vax3_delta2 < 0] <- NA # If <0 then record as missing as captured at baseline 
data2$days_to_vax3_delta2[data2$days_to_vax3_delta2 > data2$full_time_delta2] <- NA # Vaccinated after period of analysis


# Reshape data
survival_delta2 <- tmerge(data1 = data2[, c("FK_Patient_Link_ID", "Patient.Age", "age_band", "Sex", "EthnicMainGroup", "imd_score", "imd_decile", "n_tests_delta2", "full_time_delta2", "vax_delta2", "vax1_date", "vax2_date", "vax3_date", "vax_type_delta2", "outcome_date_delta2", "prev_inf_delta2", "inf1_date", "vax1_mRNA", "vax2_mRNA", "vax3_mRNA", "Deceased", "DeathDate", "health_issue", "immortal_delta2", "n_tests_post_delta2", "flu_vax_delta2")], data2 = data2, id = FK_Patient_Link_ID, tstop = full_time_delta2) # Calculate each persons total days in study (to capture people not infected)
survival_delta2 <- tmerge(data1 = survival_delta2, data2 = data2, id = FK_Patient_Link_ID, infection = event(days_to_infection_delta2)) # Add in time period pre- and -infection
survival_delta2 <- tmerge(data1 = survival_delta2, data2 = data2, id = FK_Patient_Link_ID, vax1_new = tdc(days_to_vax1_delta2)) # Time varying vaccination status - first dose
survival_delta2 <- tmerge(data1 = survival_delta2, data2 = data2, id = FK_Patient_Link_ID, vax2_new = tdc(days_to_vax2_delta2)) # Time varying vaccination status - second dose
survival_delta2 <- tmerge(data1 = survival_delta2, data2 = data2, id = FK_Patient_Link_ID, vax3_new = tdc(days_to_vax3_delta2)) # Time varying vaccination status - third dose
survival_delta2 <- tmerge(data1 = survival_delta2, data2 = survival_delta2, id = FK_Patient_Link_ID, n_episodes = cumtdc(tstart)) # Create number of period per person

# Recode vaccination data so time-varying
survival_delta2$days_vax1_outcome <- interval(survival_delta2$vax1_date, survival_delta2$outcome_date_delta2) / days(1) # Number of days between first dose and infection - will be positive if infection is after vax, negative if vax after infection 
survival_delta2$vax_delta2[(survival_delta2$vax1_new == 1) & (survival_delta2$infection == 1) & survival_delta2$days_vax1_outcome >= 14] <- 1 # If received their second vaccine dose two weeks before infection, then recode vaccination status - if infected
survival_delta2$vax_delta2[(survival_delta2$vax1_new == 1) & (survival_delta2$infection == 0)] <- 1 # If received their second vaccine and not infected

survival_delta2$days_vax2_outcome <- interval(survival_delta2$vax2_date, survival_delta2$outcome_date_delta2) / days(1) # Number of days between second dose and infection
survival_delta2$vax_delta2[(survival_delta2$vax2_new == 1) & (survival_delta2$infection == 1) & (survival_delta2$days_vax2_outcome >= 7)] <- 2 # If received their second vaccine dose one week before infection, then recode vaccination status - if infected
survival_delta2$vax_delta2[(survival_delta2$vax2_new == 1) & (survival_delta2$infection == 0)] <- 2 # If received their second vaccine and not infected

survival_delta2$days_vax3_outcome <- interval(survival_delta2$vax3_date, survival_delta2$outcome_date_delta2) / days(1) # Number of days between third dose and infection
survival_delta2$vax_delta2[(survival_delta2$vax3_new == 1) & (survival_delta2$infection == 1) & (survival_delta2$days_vax3_outcome >= 7)] <- 3 # If received their second vaccine dose one week before infection, then recode vaccination status - if infected
survival_delta2$vax_delta2[(survival_delta2$vax3_new == 1) & (survival_delta2$infection == 0)] <- 3 # If received their second vaccine and not infected

# Recode vaccination type so time-varying (sorry not sorry I have not commented every line - lazy)
survival_delta2$vax_type_delta2[(survival_delta2$vax1_new == 1 & survival_delta2$infection == 1 & survival_delta2$days_vax1_outcome >= 14 & survival_delta2$vax1_mRNA == 0)] <- 1
survival_delta2$vax_type_delta2[(survival_delta2$vax1_new == 1 & survival_delta2$infection == 0 & survival_delta2$vax1_mRNA == 0)] <- 1

survival_delta2$vax_type_delta2[(survival_delta2$vax1_new == 1 & survival_delta2$infection == 1 & survival_delta2$days_vax1_outcome >= 14 & survival_delta2$vax1_mRNA == 1)] <- 2 
survival_delta2$vax_type_delta2[(survival_delta2$vax1_new == 1 & survival_delta2$infection == 0 & survival_delta2$vax1_mRNA == 1)] <- 2 

survival_delta2$vax_type_delta2[(survival_delta2$vax2_new == 1 & survival_delta2$infection == 1 & survival_delta2$days_vax2_outcome >= 7 & (survival_delta2$vax1_mRNA == 0 & survival_delta2$vax2_mRNA == 0))] <- 3
survival_delta2$vax_type_delta2[(survival_delta2$vax2_new == 1 & survival_delta2$infection == 0 & (survival_delta2$vax1_mRNA == 0 & survival_delta2$vax2_mRNA == 0))] <- 3

survival_delta2$vax_type_delta2[(survival_delta2$vax2_new == 1 & survival_delta2$infection == 1 & survival_delta2$days_vax2_outcome >= 7 & (survival_delta2$vax1_mRNA == 1 & survival_delta2$vax2_mRNA == 1))] <- 4
survival_delta2$vax_type_delta2[(survival_delta2$vax2_new == 1 & survival_delta2$infection == 0 & (survival_delta2$vax1_mRNA == 1 & survival_delta2$vax2_mRNA == 1))] <- 4

survival_delta2$vax_type_delta2[(survival_delta2$vax2_new == 1 & survival_delta2$infection == 1 & survival_delta2$days_vax2_outcome >= 7 & ((survival_delta2$vax1_mRNA == 0 & survival_delta2$vax2_mRNA == 1) | (survival_delta2$vax1_mRNA == 1 & survival_delta2$vax2_mRNA == 0)))] <- 5
survival_delta2$vax_type_delta2[(survival_delta2$vax2_new == 1 & survival_delta2$infection == 0 & ((survival_delta2$vax1_mRNA == 0 & survival_delta2$vax2_mRNA == 1) | (survival_delta2$vax1_mRNA == 1 & survival_delta2$vax2_mRNA == 0)))] <- 5

survival_delta2$vax_type_delta2[(survival_delta2$vax3_new == 1 & survival_delta2$infection == 1 & survival_delta2$days_vax3_outcome >= 7 & (survival_delta2$vax1_mRNA == 0 & survival_delta2$vax2_mRNA == 0 & survival_delta2$vax3_mRNA == 0))] <- 6
survival_delta2$vax_type_delta2[(survival_delta2$vax3_new == 1 & survival_delta2$infection == 0 & (survival_delta2$vax1_mRNA == 0 & survival_delta2$vax2_mRNA == 0 & survival_delta2$vax3_mRNA == 0))] <- 6

survival_delta2$vax_type_delta2[(survival_delta2$vax3_new == 1 & survival_delta2$infection == 1 & survival_delta2$days_vax3_outcome >= 7 & (survival_delta2$vax1_mRNA == 1 & survival_delta2$vax2_mRNA == 1 & survival_delta2$vax3_mRNA == 1))] <- 7
survival_delta2$vax_type_delta2[(survival_delta2$vax3_new == 1 & survival_delta2$infection == 0 & (survival_delta2$vax1_mRNA == 1 & survival_delta2$vax2_mRNA == 1 & survival_delta2$vax3_mRNA == 1))] <- 7

survival_delta2$vax_type_delta2[(survival_delta2$vax3_new == 1 & survival_delta2$infection == 1 & survival_delta2$days_vax3_outcome >= 7 & (survival_delta2$vax1_mRNA == 0 & survival_delta2$vax2_mRNA == 0 & survival_delta2$vax3_mRNA == 1))] <- 8
survival_delta2$vax_type_delta2[(survival_delta2$vax3_new == 1 & survival_delta2$infection == 0 & (survival_delta2$vax1_mRNA == 0 & survival_delta2$vax2_mRNA == 0 & survival_delta2$vax3_mRNA == 1))] <- 8

survival_delta2$vax_type_delta2[(survival_delta2$vax3_new == 1 & survival_delta2$infection == 1 & survival_delta2$days_vax3_outcome >= 7 & (survival_delta2$vax_type_delta2 != 6 & survival_delta2$vax_type_delta2 != 7 & survival_delta2$vax_type_delta2 != 8))] <- 9
survival_delta2$vax_type_delta2[(survival_delta2$vax3_new == 1 & survival_delta2$infection == 0 & (survival_delta2$vax_type_delta2 != 6 & survival_delta2$vax_type_delta2 != 7 & survival_delta2$vax_type_delta2 != 8))] <- 9

# # Recode previous infection so time-varying
# survival_delta2$days_previnf_outcome <- interval(survival_delta2$inf1_date, survival_delta2$outcome_date_delta2) / days(1) # Number of days between first infection and current infection
# survival_delta2$prev_inf_delta2[(survival_delta2$prev_inf_delta2 == 0) & (survival_delta2$days_previnf_outcome >= 90)] <- 1 # Recode status if not previous infection and current infection greater than 90 days


# iii. Omnicron 

# Calculate total time period for exposure
data2$full_time_omni <- interval(ymd_hms("2021-12-13 00:00:00"), ymd_hms("2022-01-10 00:00:00")) / days(1)

# Count number of days to infection from 1st September (if infected, else NA)
data2$days_to_infection_omni <- NA # Create variable and set as missing to represent those not infected
data2$days_to_infection_omni <- interval(ymd_hms("2021-12-13 00:00:00"), data2$outcome_date_omni) / days(1) # Number of days to infection
data2$days_to_infection_omni[data2$days_to_infection_omni > data2$full_time_omni] <- NA # Vaccinated after period of analysis

# Days to first vaccination (of received a vaccination post baseline period)
data2$days_to_vax1_omni <- NA # As above
data2$days_to_vax1_omni <- interval(ymd_hms("2021-12-13 00:00:00"), data2$vax1_date) / days(1) # Number of days to vaccination
data2$days_to_vax1_omni[data2$days_to_vax1_omni < 0] <- NA # If <0 then record as missing as captured at baseline 
data2$days_to_vax3_omni[data2$days_to_vax3_omni > data2$full_time_omni] <- NA # Vaccinated after period of analysis

# Days to second vaccination (of received a vaccination post baseline period)
data2$days_to_vax2_omni <- NA # As above
data2$days_to_vax2_omni <- interval(ymd_hms("2021-12-13 00:00:00"), data2$vax2_date) / days(1) # Number of days to vaccination
data2$days_to_vax2_omni[data2$days_to_vax2_omni < 0] <- NA # If <0 then record as missing as captured at baseline 
data2$days_to_vax3_omni[data2$days_to_vax3_omni > data2$full_time_omni] <- NA # Vaccinated after period of analysis

# Days to third vaccination (of received a vaccination post baseline period)
data2$days_to_vax3_omni <- NA # As above
data2$days_to_vax3_omni <- interval(ymd_hms("2021-12-13 00:00:00"), data2$vax3_date) / days(1) # Number of days to vaccination
data2$days_to_vax3_omni[data2$days_to_vax3_omni < 0] <- NA # If <0 then record as missing as captured at baseline 
data2$days_to_vax3_omni[data2$days_to_vax3_omni > data2$full_time_omni] <- NA # Vaccinated after period of analysis


# Reshape data
survival_omni <- tmerge(data1 = data2[, c("FK_Patient_Link_ID", "Patient.Age", "age_band", "Sex", "EthnicMainGroup", "imd_score", "imd_decile", "n_tests_omni", "full_time_omni", "vax_omni", "vax1_date", "vax2_date", "vax3_date", "vax_type_omni", "outcome_date_omni", "prev_inf_omni", "inf1_date", "vax1_mRNA", "vax2_mRNA", "vax3_mRNA", "Deceased", "DeathDate", "health_issue", "immortal_omni", "n_tests_post_omni", "flu_vax_omni")], data2 = data2, id = FK_Patient_Link_ID, tstop = full_time_omni) # Calculate each persons total days in study (to capture people not infected)
survival_omni <- tmerge(data1 = survival_omni, data2 = data2, id = FK_Patient_Link_ID, infection = event(days_to_infection_omni)) # Add in time period pre- and -infection
survival_omni <- tmerge(data1 = survival_omni, data2 = data2, id = FK_Patient_Link_ID, vax1_new = tdc(days_to_vax1_omni)) # Time varying vaccination status - first dose
survival_omni <- tmerge(data1 = survival_omni, data2 = data2, id = FK_Patient_Link_ID, vax2_new = tdc(days_to_vax2_omni)) # Time varying vaccination status - second dose
survival_omni <- tmerge(data1 = survival_omni, data2 = data2, id = FK_Patient_Link_ID, vax3_new = tdc(days_to_vax3_omni)) # Time varying vaccination status - third dose
survival_omni <- tmerge(data1 = survival_omni, data2 = survival_omni, id = FK_Patient_Link_ID, n_episodes = cumtdc(tstart)) # Create number of period per person

# Recode vaccination data so time-varying
survival_omni$days_vax1_outcome <- interval(survival_omni$vax1_date, survival_omni$outcome_date_omni) / days(1) # Number of days between first dose and infection - will be positive if infection is after vax, negative if vax after infection 
survival_omni$vax_omni[(survival_omni$vax1_new == 1) & (survival_omni$infection == 1) & survival_omni$days_vax1_outcome >= 14] <- 1 # If received their second vaccine dose two weeks before infection, then recode vaccination status - if infected
survival_omni$vax_omni[(survival_omni$vax1_new == 1) & (survival_omni$infection == 0)] <- 1 # If received their second vaccine and not infected

survival_omni$days_vax2_outcome <- interval(survival_omni$vax2_date, survival_omni$outcome_date_omni) / days(1) # Number of days between second dose and infection
survival_omni$vax_omni[(survival_omni$vax2_new == 1) & (survival_omni$infection == 1) & (survival_omni$days_vax2_outcome >= 7)] <- 2 # If received their second vaccine dose one week before infection, then recode vaccination status - if infected
survival_omni$vax_omni[(survival_omni$vax2_new == 1) & (survival_omni$infection == 0)] <- 2 # If received their second vaccine and not infected

survival_omni$days_vax3_outcome <- interval(survival_omni$vax3_date, survival_omni$outcome_date_omni) / days(1) # Number of days between third dose and infection
survival_omni$vax_omni[(survival_omni$vax3_new == 1) & (survival_omni$infection == 1) & (survival_omni$days_vax3_outcome >= 7)] <- 3 # If received their second vaccine dose one week before infection, then recode vaccination status - if infected
survival_omni$vax_omni[(survival_omni$vax3_new == 1) & (survival_omni$infection == 0)] <- 3 # If received their second vaccine and not infected

# Recode vaccination type so time-varying (sorry not sorry I have not commented every line - lazy)
survival_omni$vax_type_omni[(survival_omni$vax1_new == 1 & survival_omni$infection == 1 & survival_omni$days_vax1_outcome >= 14 & survival_omni$vax1_mRNA == 0)] <- 1
survival_omni$vax_type_omni[(survival_omni$vax1_new == 1 & survival_omni$infection == 0 & survival_omni$vax1_mRNA == 0)] <- 1

survival_omni$vax_type_omni[(survival_omni$vax1_new == 1 & survival_omni$infection == 1 & survival_omni$days_vax1_outcome >= 14 & survival_omni$vax1_mRNA == 1)] <- 2 
survival_omni$vax_type_omni[(survival_omni$vax1_new == 1 & survival_omni$infection == 0 & survival_omni$vax1_mRNA == 1)] <- 2 

survival_omni$vax_type_omni[(survival_omni$vax2_new == 1 & survival_omni$infection == 1 & survival_omni$days_vax2_outcome >= 7 & (survival_omni$vax1_mRNA == 0 & survival_omni$vax2_mRNA == 0))] <- 3
survival_omni$vax_type_omni[(survival_omni$vax2_new == 1 & survival_omni$infection == 0 & (survival_omni$vax1_mRNA == 0 & survival_omni$vax2_mRNA == 0))] <- 3

survival_omni$vax_type_omni[(survival_omni$vax2_new == 1 & survival_omni$infection == 1 & survival_omni$days_vax2_outcome >= 7 & (survival_omni$vax1_mRNA == 1 & survival_omni$vax2_mRNA == 1))] <- 4
survival_omni$vax_type_omni[(survival_omni$vax2_new == 1 & survival_omni$infection == 0 & (survival_omni$vax1_mRNA == 1 & survival_omni$vax2_mRNA == 1))] <- 4

survival_omni$vax_type_omni[(survival_omni$vax2_new == 1 & survival_omni$infection == 1 & survival_omni$days_vax2_outcome >= 7 & ((survival_omni$vax1_mRNA == 0 & survival_omni$vax2_mRNA == 1) | (survival_omni$vax1_mRNA == 1 & survival_omni$vax2_mRNA == 0)))] <- 5
survival_omni$vax_type_omni[(survival_omni$vax2_new == 1 & survival_omni$infection == 0 & ((survival_omni$vax1_mRNA == 0 & survival_omni$vax2_mRNA == 1) | (survival_omni$vax1_mRNA == 1 & survival_omni$vax2_mRNA == 0)))] <- 5

survival_omni$vax_type_omni[(survival_omni$vax3_new == 1 & survival_omni$infection == 1 & survival_omni$days_vax3_outcome >= 7 & (survival_omni$vax1_mRNA == 0 & survival_omni$vax2_mRNA == 0 & survival_omni$vax3_mRNA == 0))] <- 6
survival_omni$vax_type_omni[(survival_omni$vax3_new == 1 & survival_omni$infection == 0 & (survival_omni$vax1_mRNA == 0 & survival_omni$vax2_mRNA == 0 & survival_omni$vax3_mRNA == 0))] <- 6

survival_omni$vax_type_omni[(survival_omni$vax3_new == 1 & survival_omni$infection == 1 & survival_omni$days_vax3_outcome >= 7 & (survival_omni$vax1_mRNA == 1 & survival_omni$vax2_mRNA == 1 & survival_omni$vax3_mRNA == 1))] <- 7
survival_omni$vax_type_omni[(survival_omni$vax3_new == 1 & survival_omni$infection == 0 & (survival_omni$vax1_mRNA == 1 & survival_omni$vax2_mRNA == 1 & survival_omni$vax3_mRNA == 1))] <- 7

survival_omni$vax_type_omni[(survival_omni$vax3_new == 1 & survival_omni$infection == 1 & survival_omni$days_vax3_outcome >= 7 & (survival_omni$vax1_mRNA == 0 & survival_omni$vax2_mRNA == 0 & survival_omni$vax3_mRNA == 1))] <- 8
survival_omni$vax_type_omni[(survival_omni$vax3_new == 1 & survival_omni$infection == 0 & (survival_omni$vax1_mRNA == 0 & survival_omni$vax2_mRNA == 0 & survival_omni$vax3_mRNA == 1))] <- 8

survival_omni$vax_type_omni[(survival_omni$vax3_new == 1 & survival_omni$infection == 1 & survival_omni$days_vax3_outcome >= 7 & (survival_omni$vax_type_omni != 6 & survival_omni$vax_type_omni != 7 & survival_omni$vax_type_omni != 8))] <- 9
survival_omni$vax_type_omni[(survival_omni$vax3_new == 1 & survival_omni$infection == 0 & (survival_omni$vax_type_omni != 6 & survival_omni$vax_type_omni != 7 & survival_omni$vax_type_omni != 8))] <- 9

# # Recode previous infection so time-varying
# survival_omni$days_previnf_outcome <- interval(survival_omni$inf1_date, survival_omni$outcome_date_omni) / days(1) # Number of days between first infection and current infection
# survival_omni$prev_inf_omni[(survival_omni$prev_inf_omni == 0) & (survival_omni$days_previnf_outcome >= 90)] <- 1 # Recode status if not previous infection and current infection greater than 90 days

# # # Save survival analysis datasets for faster load times
# save(survival_delta1, file = "./Data/survival_delta1.RData") # Save
# save(survival_delta2, file = "./Data/survival_delta2.RData")
# save(survival_omni, file = "./Data/survival_omni.RData")
# save(data2, file = "./Data/cleaned_data.RData")
load("./Data/survival_delta1.RData") # Load
load("./Data/survival_delta2.RData")
load("./Data/survival_omni.RData")
load("./Data/cleaned_data.RData")

## 3. Descriptive tables and plots ##

# Generate descriptive tables
# options(scipen=999) # To stop reporting numbers by e+10 etc
source(vaccine_reinfection_descriptives.R) # Run code

# To get plots 
# source(vaccine_reinfection_plots.R) # Run separately as loads in data again

## 4. Survival analysis - vaccination status ##

# a. Only people who test #

# i. Delta period 1

# Unadjusted
model1_i_u1 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + vax_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),]) # Only for people who survived to end of study period
negtest_results <- tidy(model1_i_u1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
negtest_results$model <- "model1_i_u1" # Describe model type

model1_i_u2 <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_delta1 + prev_inf_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),]) # Repeat
test <- tidy(model1_i_u2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model1_i_u2" # Describe model type
negtest_results <- rbind(negtest_results, test) # Join onto main table

model1_i_u3 <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),]) # Repeat
test <- tidy(model1_i_u3, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model1_i_u3" 
negtest_results <- rbind(negtest_results, test) 

model1_i_u4 <- coxph(Surv(tstart, tstop, infection) ~ factor(imd_decile) + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),]) # Repeat
test <- tidy(model1_i_u4, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model1_i_u4" 
negtest_results <- rbind(negtest_results, test) 

# Fully adjusted
model1_i_a1 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + prev_inf_delta1 + vax_delta1:tstart + prev_inf_delta1:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue)  + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),])
test <- tidy(model1_i_a1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model1_i_a1" 
negtest_results <- rbind(negtest_results, test) 

model1_i_a2 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + prev_inf_delta1 + vax_delta1:tstart + prev_inf_delta1:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + factor(imd_decile) + factor(health_issue)  + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$n_tests_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),])
test <- tidy(model1_i_a2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model1_i_a2" 
negtest_results <- rbind(negtest_results, test) 

write.csv(negtest_results, "./Vaccine v reinfection paper/coxph_negtest_delta1.csv") # Save
rm(negtest_results)

# ii. Delta period 2

# Unadjusted
model1_ii_a1 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + vax_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),]) # Only for people who survived to end of study period
negtest_results <- tidy(model1_ii_a1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
negtest_results$model <- "model1_ii_a1" 

model1_ii_a2 <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_delta2 + prev_inf_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),]) 
test <- tidy(model1_ii_a2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model1_ii_a2" 
negtest_results <- rbind(negtest_results, test) 

model1_ii_a3 <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),]) 
test <- tidy(model1_ii_a3, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model1_ii_a3" 
negtest_results <- rbind(negtest_results, test) 

model1_ii_a4 <- coxph(Surv(tstart, tstop, infection) ~ factor(imd_decile) + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),]) 
test <- tidy(model1_ii_a4, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model1_ii_a4" 
negtest_results <- rbind(negtest_results, test) 

# Fully adjusted
model1_ii_b1 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + prev_inf_delta2 + vax_delta2:tstart + prev_inf_delta2:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),])
test <- tidy(model1_ii_b1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model1_ii_b1" 
negtest_results <- rbind(negtest_results, test) 

model1_ii_b2 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + prev_inf_delta2 + vax_delta2:tstart + prev_inf_delta2:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + factor(imd_decile) + factor(health_issue) + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$n_tests_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),])
test <- tidy(model1_ii_b2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model1_ii_b2" 
negtest_results <- rbind(negtest_results, test) 

write.csv(negtest_results, "./Vaccine v reinfection paper/coxph_negtest_delta2.csv") # Save
rm(negtest_results)


# iii. Omnicron

# Unadjusted
model1_iii_a1 <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + vax_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),]) # Only for people who survived to end of study period
negtest_results <- tidy(model1_iii_a1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
negtest_results$model <- "model1_iii_a1" 

model1_iii_a2 <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_omni + prev_inf_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),])
test <- tidy(model1_iii_a2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model1_iii_a2" 
negtest_results <- rbind(negtest_results, test)

model1_iii_a3 <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),])
test <- tidy(model1_iii_a3, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model1_iii_a3" 
negtest_results <- rbind(negtest_results, test)

model1_iii_a4 <- coxph(Surv(tstart, tstop, infection) ~ factor(imd_decile) + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),])
test <- tidy(model1_iii_a4, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model1_iii_a4" 
negtest_results <- rbind(negtest_results, test)

# Fully adjusted
model1_iii_b1 <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + prev_inf_omni + vax_omni:tstart + prev_inf_omni:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),])
test <- tidy(model1_iii_b1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model1_iii_b1" 
negtest_results <- rbind(negtest_results, test)

model1_iii_b2 <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + prev_inf_omni + vax_omni:tstart + prev_inf_omni:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + factor(imd_decile) + factor(health_issue) + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$n_tests_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),])
test <- tidy(model1_iii_b2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model1_iii_b2" 
negtest_results <- rbind(negtest_results, test)

write.csv(negtest_results, "./Vaccine v reinfection paper/coxph_negtest_omni.csv") # Save
rm(negtest_results)


# b. Restrict to people who also received the flu vaccine #

# i. Delta period 1

# Unadjusted
model2_i_u1 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + vax_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$flu_vax_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),]) # Only for people who survived to end of study period
table <- tidy(model2_i_u1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
table$model <- "model2_i_u1" # Describe model type

model2_i_u2 <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_delta1 + prev_inf_delta1:tstart + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$flu_vax_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),]) # Repeat
test <- tidy(model2_i_u2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model2_i_u2" # Describe model type
table <- rbind(table, test) # Join onto main table

model2_i_u3 <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$flu_vax_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),]) # Repeat
test <- tidy(model2_i_u3, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model2_i_u3" 
table <- rbind(table, test) 

model2_i_u4 <- coxph(Surv(tstart, tstop, infection) ~ factor(imd_decile) + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$flu_vax_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),]) # Repeat
test <- tidy(model2_i_u4, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model2_i_u4" 
table <- rbind(table, test) 


# Fully adjusted
model2_i_a1 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + prev_inf_delta1 + vax_delta1:tstart + prev_inf_delta1:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) +  n_tests_delta1 + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$flu_vax_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),])
test <- tidy(model2_i_a1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model2_i_a1" 
table <- rbind(table, test) 

model2_i_a2 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta1 + prev_inf_delta1 + vax_delta1:tstart + prev_inf_delta1:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + factor(imd_decile) + factor(health_issue) +  n_tests_delta1 + cluster(FK_Patient_Link_ID), data = survival_delta1[survival_delta1$flu_vax_delta1 == 1 & (survival_delta1$DeathDate >= "2021-09-01" | is.na(survival_delta1$DeathDate)),])
test <- tidy(model2_i_a2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model2_i_a2" 
table <- rbind(table, test) 

write.csv(table, "./Vaccine v reinfection paper/coxph_flu_vax_delta1.csv") # Save
rm(table)

# ii. Delta period 2

# Unadjusted
model2_ii_a1 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + vax_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$flu_vax_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),]) # Only for people who survived to end of study period
table <- tidy(model2_ii_a1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
table$model <- "model_ii_a1" 

model2_ii_a2 <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_delta2 + prev_inf_delta2:tstart + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$flu_vax_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),]) 
test <- tidy(model2_ii_a2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model2_ii_a2" 
table <- rbind(table, test) 

model2_ii_a3 <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$flu_vax_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),]) 
test <- tidy(model2_ii_a3, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model2_ii_a3" 
table <- rbind(table, test) 

model2_ii_a4 <- coxph(Surv(tstart, tstop, infection) ~ factor(imd_decile) + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$flu_vax_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),]) 
test <- tidy(model2_ii_a4, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model2_ii_a4" 
table <- rbind(table, test) 

# Fully adjusted
model2_ii_b1 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + prev_inf_delta2 + vax_delta2:tstart + prev_inf_delta2:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_delta2 + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$flu_vax_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),])
test <- tidy(model2_ii_b1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model2_ii_b1" 
table <- rbind(table, test) 

model2_ii_b2 <- coxph(Surv(tstart, tstop, infection) ~ vax_delta2 + prev_inf_delta2 + vax_delta2:tstart + prev_inf_delta2:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + factor(imd_decile) + factor(health_issue) + n_tests_delta2 + cluster(FK_Patient_Link_ID), data = survival_delta2[survival_delta2$flu_vax_delta2 == 1 & (survival_delta2$DeathDate >= "2021-11-27" | is.na(survival_delta2$DeathDate)),])
test <- tidy(model2_ii_b2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model2_ii_b2" 
table <- rbind(table, test) 

write.csv(table, "./Vaccine v reinfection paper/coxph_flu_vax_delta2.csv") # Save
rm(table)

# iii. Omnicron

# Unadjusted
model2_iii_a1 <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + vax_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$flu_vax_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),]) # Only for people who survived to end of study period
table <- tidy(model2_iii_a1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
table$model <- "model2_iii_a1" 

model2_iii_a2 <- coxph(Surv(tstart, tstop, infection) ~ prev_inf_omni + prev_inf_omni:tstart + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$flu_vax_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),])
test <- tidy(model2_iii_a2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model2_iii_a2" 
table <- rbind(table, test)

model2_iii_a3 <- coxph(Surv(tstart, tstop, infection) ~ imd_score + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$flu_vax_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),])
test <- tidy(model2_iii_a3, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model2_iii_a3" 
table <- rbind(table, test)

model2_iii_a4 <- coxph(Surv(tstart, tstop, infection) ~ factor(imd_decile) + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$flu_vax_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),])
test <- tidy(model2_iii_a4, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model2_iii_a4" 
table <- rbind(table, test)

# Fully adjusted
model2_iii_b1 <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + prev_inf_omni + vax_omni:tstart + prev_inf_omni:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + imd_score + factor(health_issue) + n_tests_omni + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$flu_vax_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),])
test <- tidy(model2_iii_b1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model2_iii_b1" 
table <- rbind(table, test)

model2_iii_b2 <- coxph(Surv(tstart, tstop, infection) ~ vax_omni + prev_inf_omni + vax_omni:tstart + prev_inf_omni:tstart + factor(age_band) + factor(Sex) + factor(EthnicMainGroup) + factor(imd_decile) + factor(health_issue) + n_tests_omni + cluster(FK_Patient_Link_ID), data = survival_omni[survival_omni$flu_vax_omni == 1 & (survival_omni$DeathDate >= "2022-03-01" | is.na(survival_omni$DeathDate)),])
test <- tidy(model2_iii_b2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) # Save model results
test$model <- "model2_iii_b2" 
table <- rbind(table, test)

write.csv(table, "./Vaccine v reinfection paper/coxph_flu_vax_omni.csv") # Save
rm(table)


# 4. Sensitivity analyses #

source("./Scripts/vaccine_reinfection_sensitivity.R")


# 5. Age stratified analyses #

source("./Scripts/vaccine_reinfection_age.R")

