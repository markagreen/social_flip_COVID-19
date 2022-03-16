##############################
### Vaccine vs reinfection ###
### Descriptives and Plots ###
##############################

# Purpose: To produce a series of survival analyses that are included in the paper.

# In case internet drops, redo this when back
# setwd("/Volumes/MAST/Smart_LCR_report")

# Libraries
library(data.table) # For data wrangling (plus those below)
library(lubridate) 
library(tidyr)
library(dplyr)
library(plyr)
library(ggplot2) # For visualisations
library(viridis)
library(zoo) # For rolling average


### 1. Tidy data ###

#  Load data
load("/Volumes/Sharing Folder/2022-03-02/VaccineAnalysisDatasets.RData") # Mac
#load("Q:/2022-03-02/VaccineAnalysisDatasets.RData") # Windows
rm(cohort1)

# 1a. Vaccination data #

# Subset only COVID-19 vaccines (flu jabs are in here for example) and people who completed the vaccine
data1 <- data.table(data1) # Convert data type
data1 <- data1[data1$Status == "completed" & (data1$VaccinationProcedure == "1324681000000101" | data1$VaccinationProcedure == "1324691000000104" | data1$VaccinationProcedure == "1362591000000103")] # Subset
# Vaccine SNOMED codes are in "VaccinationProcedure". 1st dose=1324681000000101, 2nd dose=1324691000000104, Booster=1362591000000103. Completed = people who attended and were received their vaccine (not all receive them after attending)

# Tidy vaccination data and link to main dataset
data1$t_date <- ymd_hms(data1$DateTimeAdministered) # Convert to time-date
data1 <- unique(data1, by = c("FK_Patient_Link_ID", "t_date")) # Remove repeated vaccines (just over 60k cases are repeated, so remove one of each)
data1[order(FK_Patient_Link_ID, t_date), num_vax:=1:.N, by=.(FK_Patient_Link_ID)] # Count number of vaccine doses
vax <- data1[num_vax <= 3] # Store upto three vaccine doses only
vaccinated <- dcast(vax, FK_Patient_Link_ID ~ num_vax, value.var = c("t_date", "Manufacturer")) # Reshape to wide format to list when got each vaccine
names(vaccinated) <- c("FK_Patient_Link_ID", "vax1_date", "vax2_date", "vax3_date", "vax1_type", "vax2_type", "vax3_type") # Rename variables to help R
data2 <- merge(data2, vaccinated, by = "FK_Patient_Link_ID", all.x = TRUE) # Join onto main dataset
rm(data1, vax, vaccinated) # Save space

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

# 1c. Final tidying #

# Drop data from analyses
data2 <- data.table(data2) # Convert format
data2 <- data2[data2$FK_Patient_Link_ID != -1,] # Drop missing person ID (n=21)
data2 <- data2[!is.na(data2$FK_Patient_Link_ID),] # Drop missing person ID (n=2)


# Load and join on IMD decile
imd <- fread("/Volumes/MAST/SMART_LCR_report/Data/imd2019.csv") # Mac
#imd <- fread("R:/SMART_LCR_report/Data/imd2019.csv") # Windows
data2 <- merge(data2, imd, by.x = "LSOA_Code", by.y = "lsoa11", all.x = TRUE) # join together
rm(imd)



# 1d. Quick and dirty descriptive plot #

# Generate data for plotting
new_inf <- data2[, list(New = .N), by = "inf1_date"] # Aggregate counts of infections for new infection
names(new_inf)[names(new_inf) == "inf1_date"] <- "day" # Rename variable
reinf1 <- data2[, list(Reinfection1 = .N), by = "inf2_date"] # Aggregate counts of reinfections (first reinfection)
names(reinf1)[names(reinf1) == "inf2_date"] <- "day"
reinf2 <- data2[, list(Reinfection2 = .N), by = "inf3_date"] # Aggregate counts of reinfections (second)
names(reinf2)[names(reinf2) == "inf3_date"] <- "day"
reinf3 <- data2[, list(Reinfection3 = .N), by = "inf4_date"] # Aggregate counts of reinfections (third)
names(reinf3)[names(reinf3) == "inf4_date"] <- "day"
reinf <- merge(reinf1, reinf2, by = "day", all.x = TRUE) # Join reinfections together
reinf <- merge(reinf, reinf3, by = "day", all.x = TRUE) 
reinf$Reinfection1[is.na(reinf$Reinfection1)] <- 0 # If missing, then should be 0
reinf$Reinfection2[is.na(reinf$Reinfection2)] <- 0 
reinf$Reinfection3[is.na(reinf$Reinfection3)] <- 0 
reinf$Reinfection <- reinf$Reinfection1 + reinf$Reinfection3 + reinf$Reinfection3 # Sum together all reinfections
reinf <- reinf[,c("day", "Reinfection")] # Drop variables not needed
rm(reinf1, reinf2, reinf3)

# Combine all data together and tidy
trend_data <- merge(new_inf, reinf, by = "day", all.x = TRUE) # Join reinfections onto new infections
trend_data$Reinfection[is.na(trend_data$Reinfection)] <- 0 # NAs are actually 0s
trend_data <- trend_data[!is.na(trend_data$day)] # Drop missing
trend_data <- trend_data[trend_data$day != "1990-01-01"] # Drop missing
trend_data$New_7day <- rollmean(trend_data$New, k = 7, fill = NA) # 7 day average
trend_data$Reinfection_7day <- rollmean(trend_data$Reinfection, k = 7, fill = NA)
trend_data_long <- gather(data = trend_data, key = Infection, value = freq, New_7day:Reinfection_7day) # Convert to long format as easier for plotting

# Plot
plot1 <- ggplot() +
  geom_line(data = trend_data_long, aes(x = day, y = freq, group = Infection, color = Infection)) +
  geom_point(data = trend_data_long[trend_data_long$Infection == "New_7day",], aes(x = day, y = New, group = Infection, color = Infection), alpha = 0.1) +
  geom_point(data = trend_data_long[trend_data_long$Infection == "Reinfection_7day",], aes(x = day, y = Reinfection, group = Infection, color = Infection), alpha = 0.1) +
  xlab("Date") +
  ylab("Number of registered positive tests") +
  scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9, labels = c("First positive", "Further positives")) +
  labs(title = "Number of COVID-19 registered positive tests",
       subtitle = "The line represents the 7-day moving average, with points the daily value.") +
  theme(text = element_text(size = 16), plot.title = element_text(size = 20, face = "bold"),)  +
  scale_x_date(breaks = "4 months", minor_breaks = "1 month", date_labels = "%b %Y")
plot1
ggsave(plot = plot1, filename = "./Vaccine v reinfection paper/reinfection_trends.jpeg") # Save
ggsave(plot = plot1, filename = "./Vaccine v reinfection paper/reinfection_trends.svg") 

# Calculate percentage of reinfections for a specific period
(sum(trend_data$Reinfection[trend_data$day >= "2021-12-13"]) / (sum(trend_data$Reinfection[trend_data$day >= "2021-12-13"]) + sum(trend_data$New[trend_data$day >= "2021-12-13"]))) * 100

# Calculate a denominator (add in deaths)
data2$DeathDate <- ymd(data2$DeathDate) # Convert to time-date 
data2$death_noinf <- 0 # Classify if deaths before or after infection
data2$death_noinf[data2$DeathDate < data2$inf1_date] <- 1
data2$death_reinf <- 0
data2$death_reinf[data2$DeathDate >= data2$inf1_date] <- 1
deaths <- data2[, list(death_noinf = sum(data2$death_noinf, na.rm=T), death_reinf = sum(death_reinf, na.rm=T)), by = "DeathDate"] # Aggregate counts of deaths
names(deaths)[names(deaths) == "DeathDate"] <- "day"
trend_data <- merge(trend_data, deaths, by = "day", all.x = TRUE) # Join deaths onto data
trend_data$death_noinf[is.na(trend_data$death_noinf)] <- 0 # NAs are actually 0s
trend_data$death_reinf[is.na(trend_data$death_reinf)] <- 0

# Calculate rates
trend_data$cum_new_inf <- cumsum(trend_data$New) # Create cumulative sum of infections 
trend_data$pop_new <- nrow(data2) - trend_data$cum_new_inf - trend_data$death_noinf # subtract above from total pop and deaths
trend_data$pop_reinf <- trend_data$cum_new_inf - trend_data$death_reinf # Get pop for reinfections (infection people - deaths)
trend_data$new_rate <- (trend_data$New_7day / trend_data$pop_new) * 100000 # Calculate rates
trend_data$reinf_rate <- (trend_data$Reinfection_7day / trend_data$pop_reinf) * 100000
trend_data_long <- gather(data = trend_data, key = Infection, value = rate, new_rate:reinf_rate) # Convert to long format as easier for plotting

# Plot
plot2 <- ggplot() +
  geom_line(data = trend_data_long, aes(x = day, y = rate, group = Infection, color = Infection)) +
  xlab("Date") +
  ylab("Number of registered positive tests per 100,000") +
  scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9, labels = c("First positive", "Further positives")) +
  labs(title = "Number of COVID-19 registered positive tests per capita",
       subtitle = "The line represents the 7-day moving average.") +
  theme(text = element_text(size = 16), plot.title = element_text(size = 20, face = "bold"),)  +
  scale_x_date(breaks = "4 months", minor_breaks = "1 month", date_labels = "%b %Y")
plot2
ggsave(plot = plot2, filename = "./Vaccine v reinfection paper/reinfection_trends_rate.jpeg") # Save
ggsave(plot = plot2, filename = "./Vaccine v reinfection paper/reinfection_trends_rate.svg") 

# Redo by deprivation #

# Calculate counts of infections and reinfections
new_inf <- data2[, list(New = .N), by = c("inf1_date", "imd_decile")] # Aggregate counts of infections for new infection
names(new_inf)[names(new_inf) == "inf1_date"] <- "day" # Rename variable
reinf1 <- data2[, list(Reinfection1 = .N), by = c("inf2_date", "imd_decile")] # Aggregate counts of reinfections (first reinfection)
names(reinf1)[names(reinf1) == "inf2_date"] <- "day"
reinf2 <- data2[, list(Reinfection2 = .N), by = c("inf3_date", "imd_decile")] # Repeat
names(reinf2)[names(reinf2) == "inf3_date"] <- "day"
reinf3 <- data2[, list(Reinfection3 = .N), by = c("inf4_date", "imd_decile")] # Repeat
names(reinf3)[names(reinf3) == "inf4_date"] <- "day"
reinf <- merge(reinf1, reinf2, by = c("day", "imd_decile"), all.x = TRUE) # Join reinfections together
reinf <- merge(reinf, reinf3, by = c("day", "imd_decile"), all.x = TRUE) 
reinf$Reinfection1[is.na(reinf$Reinfection1)] <- 0 # If missing, then should be 0
reinf$Reinfection2[is.na(reinf$Reinfection2)] <- 0 
reinf$Reinfection3[is.na(reinf$Reinfection3)] <- 0 
reinf$Reinfection <- reinf$Reinfection1 + reinf$Reinfection3 + reinf$Reinfection3 # Sum together all reinfections
reinf <- reinf[,c("day", "imd_decile", "Reinfection")] # Drop variables not needed
reinf <- reinf[!is.na(reinf$day)] # Drop missing day
rm(reinf1, reinf2, reinf3)


# Get into single dataset
trend_data2 <- merge(new_inf, reinf, by = c("day", "imd_decile"), all.x = TRUE) # Join reinfections onto new infections
trend_data2$Reinfection[is.na(trend_data2$Reinfection)] <- 0 # If missing, then should be 0
trend_data2 <- trend_data2[!is.na(trend_data2$day)] # Drop missing
trend_data2 <- trend_data2[trend_data2$day != "1990-01-01"] # Drop missing
trend_data2 <- trend_data2[!is.na(trend_data2$imd_decile)] # Drop missing
trend_data2 <- trend_data2[order(trend_data2$imd_decile, trend_data2$day)] # Order by day and IMD
# trend_data2 <- trend_data2 %>% dplyr::group_by(imd_decile) %>% # Calculate 7 day averages by imd
#   mutate(New_7day=rollmean(New, k = 7, fill = NA, align="center")) %>% 
#   mutate(Reinfection_7day=rollmean(Reinfection, k = 7, fill = NA))
trend_data2 <- data.table(trend_data2) # For next step
trend_data2[, New_7day := rollmean(New, k = 7, fill = NA, align="center"), by = .(imd_decile)] # Calculate 7 day moving average by imd for new infections
trend_data2[, Reinfection_7day := rollmean(Reinfection, k = 7, fill = NA, align="center"), by = .(imd_decile)] # Repeat for reinfections

# Calculate a denominator (add in deaths)
deaths_imd <- data2[, list(death_noinf = sum(data2$death_noinf), death_reinf = sum(death_reinf)), by = c("DeathDate", "imd_decile")] # Aggregate counts of deaths
names(deaths_imd)[names(deaths_imd) == "DeathDate"] <- "day"
trend_data2 <- merge(trend_data2, deaths_imd, by = c("day", "imd_decile"), all.x = TRUE) # Join deaths onto data
trend_data2$death_noinf[is.na(trend_data2$death_noinf)] <- 0 # NAs are actually 0s
trend_data2$death_reinf[is.na(trend_data2$death_reinf)] <- 0

# Population (sort out!)
pop <- data2[, list(pop = .N), by = "imd_decile"] # Aggregate population by imd
trend_data2 <- merge(trend_data2, pop, by = "imd_decile", all.x = TRUE) #  Population at start
trend_data2 <- trend_data2[order(trend_data2$imd_decile, trend_data2$day),] # Order by date
trend_data2 <- data.table(trend_data2) # For next step
trend_data2[, cum_new_inf := cumsum(New), by = .(imd_decile)] # Calculate cumulative infections by imd

# Calculate measures to plot
trend_data2$pop_new <- trend_data2$pop - trend_data2$cum_new_inf - trend_data2$death_noinf # subtract above from total pop and deaths
trend_data2$pop_reinf <- trend_data2$cum_new_inf - trend_data2$death_reinf # Get pop for reinfections (infection people - deaths)
trend_data2$new_rate <- (trend_data2$New_7day / trend_data2$pop_new) * 100000 # Calculate rates
trend_data2$reinf_rate <- (trend_data2$Reinfection_7day / trend_data2$pop_reinf) * 100000
trend_data2_long <- gather(data = trend_data2, key = Infection, value = rate, new_rate:reinf_rate) # Convert to long format as easier for plotting
trend_data2$pc_reinf <- (trend_data2$Reinfection_7day / (trend_data2$New_7day + trend_data2$Reinfection_7day)) * 100 # Calculate % reinfection


# Plot 1: Infection rate over time

# Labels
labels <- c("First positive", "Further positives")
names(labels) <- c("new_rate", "reinf_rate")

# Plot
ses_plot1 <- ggplot(data = trend_data2_long[(trend_data2_long$imd_decile == 1 | trend_data2_long$imd_decile == 10),], aes(x = day, y = rate, group = factor(imd_decile), color = factor(imd_decile))) +
  geom_line() +
  facet_wrap(~Infection, scales = "fixed", labeller = labeller(Infection = labels)) +
  scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9, labels = c("1: Most deprived", "10: Least Deprived")) +
  labs(color = "Deprivation decile") +
  xlab("Date") +
  ylab("Registered positive tests per 100,000 people") +
  ylim(0,NA) +
  theme(legend.position="bottom") +
  scale_x_date(breaks = "7 months", minor_breaks = "1 month", date_labels = "%b %Y")
ses_plot1
ggsave(plot = ses_plot1, filename = "./Vaccine v reinfection paper/ses_trends_whole.jpeg") # Save
ggsave(plot = ses_plot1, filename = "./Vaccine v reinfection paper/ses_trends_whole.svg") 

# Plot 2: Fix date for period

# Plot
ses_plot2 <- ggplot(data = trend_data2_long[trend_data2_long$day >= "2021-06-03" & (trend_data2_long$imd_decile == 1 | trend_data2_long$imd_decile == 10),], aes(x = day, y = rate, group = factor(imd_decile), color = factor(imd_decile))) +
  geom_line() +
  facet_wrap(~Infection, scales = "fixed", labeller = labeller(Infection = labels)) +
  scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9, labels = c("1: Most deprived", "10: Least Deprived")) +
  labs(color = "Deprivation decile") +
  xlab("Date") +
  ylab("Registered positive tests per 100,000 people") +
  ylim(0,NA) +
  theme(legend.position="bottom")
ses_plot2
ggsave(plot = ses_plot2, filename = "./Vaccine v reinfection paper/ses_trends_study_period.jpeg") # Save
ggsave(plot = ses_plot2, filename = "./Vaccine v reinfection paper/ses_trends_study_period.svg") 

# Plot 3: percent of reinfections by decile

# Plot
ses_plot3 <- ggplot(data = trend_data2[trend_data2$day >= "2021-06-03",], aes(x = day, y = pc_reinf, group = factor(imd_decile), color = factor(imd_decile))) + # Plot
  geom_line() +
  scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9) +
  labs(color = "Deprivation decile") +
  xlab("Date") +
  ylab("Percent of all positive tests that are further positive tests (%)") +
  ylim(0,NA)
ses_plot3
ggsave(plot = ses_plot3, filename = "./Vaccine v reinfection paper/ses_percent_reinfections.jpeg") # Save
ggsave(plot = ses_plot3, filename = "./Vaccine v reinfection paper/ses_percent_reinfections.svg") 


