##############################
### Vaccine vs reinfection ###
###### Spatial analyses ######
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
library(patchwork)
library(zoo) # For rolling average
library(ggeffects) # For marginal effects
library(sf) # For spatial modelling
library(spdep)
library(spatialreg)

sf::sf_use_s2(FALSE) # Else throws errors later

### 1. Tidy data ###

#  Load data
#load("/Volumes/Sharing Folder/2022-01-31/VaccineAnalysisDatasets.RData") # Mac
load("Q:/2022-03-02/VaccineAnalysisDatasets.RData") # Windows
rm(cohort1, data4)

# 1a. Vaccination data #

# Subset only COVID-19 vaccines (flu jabs are in here for example) and people who completed the vaccine
data1 <- data.table(data1) # Convert data type
data1 <- data1[data1$Status == "completed" & (data1$VaccinationProcedure == "1324681000000101" | data1$VaccinationProcedure == "1324691000000104" | data1$VaccinationProcedure == "1362591000000103")] # Subset

# Tidy vaccination data and link to main dataset
data1$t_date <- ymd_hms(data1$DateTimeAdministered) # Convert to time-date
data1 <- unique(data1, by = c("FK_Patient_Link_ID", "t_date")) # Remove repeated vaccines (just over 60k cases are repeated, so remove one of each)
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

## 2. Generate variables for mapping ##

# 2a. Positive tests #

# i. Delta period 1 - 3rd June 2021 to 1st September 2021 # [3rd June is when PHE said was 99% dominant] - or 13th May when most cases ~50% according to PHE

# Tested positive (yes/no)
data2$outcome_delta1 <- 0 # Create blank variable
data2$outcome_delta1[(data2$inf1_date >= "2021-06-03 00:00:00.0000000" & data2$inf1_date < "2021-09-01 00:00:00.0000000") | (data2$inf2_date >= "2021-06-03 00:00:00.0000000" & data2$inf2_date < "2021-09-01 00:00:00.0000000") | (data2$inf3_date >= "2021-06-03 00:00:00.0000000" & data2$inf3_date < "2021-09-01 00:00:00.0000000") | (data2$inf4_date >= "2021-06-03 00:00:00.0000000" & data2$inf4_date < "2021-09-01 00:00:00.0000000")] <- 1 # Tested positive over this time period

# ii. Delta period 2 - 1st September 2021 to 27th November 2021 #

data2$outcome_delta2 <- 0 # Create blank variable
data2$outcome_delta2[(data2$inf1_date >= "2021-09-01 00:00:00.0000000" & data2$inf1_date < "2021-11-27 00:00:00.0000000") | (data2$inf2_date >= "2021-09-01 00:00:00.0000000" & data2$inf2_date < "2021-11-27 00:00:00.0000000") | (data2$inf3_date >= "2021-09-01 00:00:00.0000000" & data2$inf3_date < "2021-11-27 00:00:00.0000000") | (data2$inf4_date >= "2021-09-01 00:00:00.0000000" & data2$inf4_date < "2021-11-27 00:00:00.0000000")] <- 1 # Tested positive over this time period

# iii. Omicron

data2$outcome_omni <- 0 # Create blank variable
data2$outcome_omni[(data2$inf1_date >= "2021-12-13 00:00:00.0000000" & data2$inf1_date < "2022-03-01 00:00:00.0000000") | (data2$inf2_date >= "2021-12-13 00:00:00.0000000" & data2$inf2_date < "2022-03-01 00:00:00.0000000") | (data2$inf3_date >= "2021-12-13 00:00:00.0000000" & data2$inf3_date < "2022-03-01 00:00:00.0000000") | (data2$inf4_date >= "2021-12-13 00:00:00.0000000" & data2$inf4_date < "2022-03-01 00:00:00.0000000")] <- 1 # Tested positive over this time period


# 2b. Vaccination status #

# i. Delta period 1

# Fully vaccinated 3 weeks before baseline (3rd June 2021) or 2 weeks for booster
data2$vax_delta1 <- 0 # Create blank variable
data2$vax_delta1[data2$vax1_date < "2021-05-13 00:00:00.0000000"] <- 1 # First dose
data2$vax_delta1[data2$vax2_date < "2021-05-13 00:00:00.0000000"] <- 2 # Second dose
data2$vax_delta1[data2$vax3_date < "2021-05-20 00:00:00.0000000"] <- 3 # Define if boosted

# ii. Delta period 2

# Fully vaccinated 3 weeks before baseline (1st September 2021) or 2 weeks for booster
data2$vax_delta2 <- 0 # Create blank variable
data2$vax_delta2[data2$vax1_date < "2021-08-11 00:00:00.0000000"] <- 1 # First dose
data2$vax_delta2[data2$vax2_date < "2021-08-11 00:00:00.0000000"] <- 2 # Second dose
data2$vax_delta2[data2$vax3_date < "2021-08-18 00:00:00.0000000"] <- 3 # Define if boosted

# iii. Omnicron

# Fully vaccinated 3 weeks before baseline (13th December 2021) or 2 weeks for booster
data2$vax_omni <- 0 # Create blank variable
data2$vax_omni[data2$vax1_date < "2021-11-22 00:00:00.0000000"] <- 1 # First dose
data2$vax_omni[data2$vax2_date < "2021-11-22 00:00:00.0000000"] <- 2 # Second dose
data2$vax_omni[data2$vax3_date < "2021-11-29 00:00:00.0000000"] <- 3 # Define if boosted


# 2c. Previous infection #

# i. Delta period 1

# Previous infection 2 weeks before baseline
data2$prev_inf_delta1 <- 0 # Create blank variable
data2$prev_inf_delta1[data2$inf1_date < "2021-05-20 00:00:00.0000000" | data2$inf2_date < "2021-05-20 00:00:00.0000000" | data2$inf3_date < "2021-05-20 00:00:00.0000000" | data2$inf4_date < "2021-05-20 00:00:00.0000000"] <- 1 # Tested positive in this time period

# ii. Delta period 2

# Previous infection 2 weeks before baseline
data2$prev_inf_delta2 <- 0 # Create blank variable
data2$prev_inf_delta2[data2$inf1_date < "2021-08-18 00:00:00.0000000" | data2$inf2_date < "2021-08-18 00:00:00.0000000" | data2$inf3_date < "2021-08-18 00:00:00.0000000" | data2$inf4_date < "2021-08-18 00:00:00.0000000"] <- 1 # Tested positive in this time period

# ii. Omnicron

# Previous infection 2 weeks before baseline
data2$prev_inf_omni <- 0 # Create blank variable
data2$prev_inf_omni[data2$inf1_date < "2021-11-29 00:00:00.0000000" | data2$inf2_date < "2021-11-29 00:00:00.0000000" | data2$inf3_date < "2021-11-29 00:00:00.0000000" | data2$inf4_date < "2021-11-29 00:00:00.0000000"] <- 1 # Tested positive in this time period


# 2d. Deaths (all-causes) #

# i. Delta period 1

data2$death_delta1 <- NA # Create variable
data2$death_delta1[data2$Deceased == "N" | (data2$Deceased == "Y" & data2$DeathDate >= "2021-06-03 00:00:00.0000000")] <- 0 # Alive at baseline (removes people who died after this period and therefore recorded as died)
data2$death_delta1[data2$DeathDate >= "2021-06-03 00:00:00.0000000" & data2$DeathDate < "2021-09-01 00:00:00.0000000"] <- 1 # Died in time period

# ii. Delta period 2

data2$death_delta2 <- NA # Create variable
data2$death_delta2[data2$Deceased == "N" | (data2$Deceased == "Y" & data2$DeathDate >= "2021-09-01 00:00:00.0000000")] <- 0 # Alive at baseline (removes people who died after this period and therefore recorded as died)
data2$death_delta2[data2$DeathDate >= "2021-09-01 00:00:00.0000000" & data2$DeathDate < "2021-11-27 00:00:00.0000000"] <- 1 # Died in time period

# ii. Omnicron

data2$death_omni <- NA # Create variable
data2$death_omni[data2$Deceased == "N" | (data2$Deceased == "Y" & data2$DeathDate >= "2021-12-13 00:00:00.0000000")] <- 0 # Alive at baseline (removes people who died after this period and therefore recorded as died)
data2$death_omni[data2$DeathDate >= "2021-12-13 00:00:00.0000000" & data2$DeathDate < "2022-03-01 00:00:00.0000000"] <- 1 # Died in time period


# # Deaths within 28 days of a positive COVID-19 test
# # # Time between death and any positive test
# # data2$died_covid[data2$Deceased == "N"] <- 0 # Did not die
# # data2$died_covid[data2$Deceased == "Y" & data2$DeathDate >= "2021-09-01 00:00:00.0000000"] <- 1 # Did not die
# data2$t_death <- ymd_hms(data2$DeathDate) # Convert to time-date
# data2$day_death <- as_date(data2$DeathDate) # Define day
# data2$days_between1 <- data2$day_death - data2$inf1_date # Number of days died after tested positive (first positive)
# data2$days_between1[data2$days_between1 < 0] <- NA # If tested positive after then say is missing
# data2$days_between2 <- data2$day_death - data2$inf2_date # Number of days died after tested positive (second positive)
# data2$days_between2[data2$days_between2 < 0] <- NA # If tested positive after then say is missing
# data2$days_between3 <- data2$day_death - data2$inf3_date # Number of days died after tested positive (third positive)
# data2$days_between3[data2$days_between3 < 0] <- NA # If tested positive after then say is missing
# data2$died_covid[data2$Deceased == "N"] <- 0 # Did not die
# data2$died_covid[data2$Deceased == "Y" & data2$DeathDate >= "2021-09-01 00:00:00.0000000" & (data2$days_between1 <= 28 | data2$days_between2 <= 28 | data2$days_between3 <= 28)] <- 1 # Died in time period within 28 days of +ive test
# 
# data2$died_noncovid[data2$Deceased == "N"] <- 0 # Did not die
# data2$died_noncovid[data2$Deceased == "Y" & data2$DeathDate >= "2021-09-01 00:00:00.0000000" & (data2$days_between1 > 28 | data2$days_between2 > 28 | data2$days_between3 > 28)] <- 1 # Died in time period within 28 days of +ive test


# 2e. Other

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


### 3. Spatial analyses ###

# 3a. Positive tests # 

# Define if person is fully vaccinated and/or previous infection - delta 1
data2$delta1_status <- NA # Create variable
data2$delta1_status[data2$vax_delta1 >= 2 & data2$prev_inf_delta1 == 0] <- 1 # Fully vaccinated and no previous infection
data2$delta1_status[data2$vax_delta1 >= 2 & data2$prev_inf_delta1 == 1] <- 2 # Fully vaccinated and previous infection
data2$delta1_status[data2$vax_delta1 < 2 & data2$prev_inf_delta1 == 0] <- 3 # Not fully vaccinated and no previous infection
data2$delta1_status[data2$vax_delta1 < 2 & data2$prev_inf_delta1 == 1] <- 4 # Not fully vaccinated and previous infection

# Define if person is fully vaccinated and/or previous infection - delta 2
data2$delta2_status <- NA # Create variable
data2$delta2_status[data2$vax_delta2 >= 2 & data2$prev_inf_delta2 == 0] <- 1 # Fully vaccinated and no previous infection
data2$delta2_status[data2$vax_delta2 >= 2 & data2$prev_inf_delta2 == 1] <- 2 # Fully vaccinated and previous infection
data2$delta2_status[data2$vax_delta2 < 2 & data2$prev_inf_delta2 == 0] <- 3 # Not fully vaccinated and no previous infection
data2$delta2_status[data2$vax_delta2 < 2 & data2$prev_inf_delta2 == 1] <- 4 # Not fully vaccinated and previous infection

# Define if person is fully vaccinated and/or previous infection - delta 1
data2$omni_status <- NA # Create variable
data2$omni_status[data2$vax_omni >= 2 & data2$prev_inf_omni == 0] <- 1 # Fully vaccinated and no previous infection
data2$omni_status[data2$vax_omni >= 2 & data2$prev_inf_omni == 1] <- 2 # Fully vaccinated and previous infection
data2$omni_status[data2$vax_omni < 2 & data2$prev_inf_omni == 0] <- 3 # Not fully vaccinated and no previous infection
data2$omni_status[data2$vax_omni < 2 & data2$prev_inf_omni == 1] <- 4 # Not fully vaccinated and previous infection

# Calculate observed positive tests per LSOA
lsoa_inf_delta1 <- data2[, list(outcome_delta1 = sum(outcome_delta1, na.rm = TRUE), pop = .N), by = c("LSOA_Code", "delta1_status")] # Aggregate to LSOA
lsoa_inf_delta1$rate_delta1 <- lsoa_inf_delta1$outcome_delta1 / lsoa_inf_delta1$pop #  Calculate rate of infections
lsoa_inf_delta1 <- lsoa_inf_delta1[, c(1,2,5)] # Drop variables not needed
lsoa_inf_delta1 <- spread(data = lsoa_inf_delta1, key = delta1_status, value = rate_delta1, fill = 0) # Reshape to wide format

lsoa_inf_delta2 <- data2[, list(outcome_delta2 = sum(outcome_delta2, na.rm = TRUE), pop = .N), by = c("LSOA_Code", "delta2_status")] # Repeat
lsoa_inf_delta2$rate_delta2 <- lsoa_inf_delta2$outcome_delta2 / lsoa_inf_delta2$pop
lsoa_inf_delta2 <- lsoa_inf_delta2[, c(1,2,5)] 
lsoa_inf_delta2 <- spread(data = lsoa_inf_delta2, key = delta2_status, value = rate_delta2, fill = 0)

lsoa_inf_omni <- data2[, list(outcome_omni = sum(outcome_omni, na.rm = TRUE), pop = .N), by = c("LSOA_Code", "omni_status")] 
lsoa_inf_omni$rate_omni <- (lsoa_inf_omni$outcome_omni / lsoa_inf_omni$pop) * 100000
lsoa_inf_omni <- lsoa_inf_omni[, c(1,2,5)] 
lsoa_inf_omni <- spread(data = lsoa_inf_omni, key = omni_status, value = rate_omni, fill = 0)

# Delta 1 = spatial analysis

# Join to spatial
# lcr_lsoas <- read_sf(dsn = "./Shapefiles/E47000004.shp") # Load LSOA shapefiles - windows
lcr_lsoas <- read_sf(dsn = "R:/SMART_LCR_report/Shapefiles/E47000004.shp") # Load LSOA shapefiles - windows
lcr_lsoas <- merge(lcr_lsoas, lsoa_inf_delta1, by.x = "lsoa11cd", by.y = "LSOA_Code", all.x = TRUE)
lcr_lsoas <- st_transform(lcr_lsoas, 4326) # Set CRS
lcr_lsoas.nb <- poly2nb(lcr_lsoas, snap=0.0002) # Identify neighbours of each LSOA

# Define adjacency using a row-standardised matrix
lw <- nb2listw(lcr_lsoas.nb)
W <- as(as_dgRMatrix_listw(lw), "CsparseMatrix")

# Calculate Getis-Ord statistic for SMR
localGi <- localG(lcr_lsoas$`1`, listw = lw) # Calculate Gi statistics - vax and no previous infection
lcr_lsoas <- cbind(lcr_lsoas, data.matrix(localGi)) # Join results onto shapefile
names(lcr_lsoas)[names(lcr_lsoas) == "data.matrix.localGi."] <- "gstat_delta1_v1" # Rename column

localGi <- localG(lcr_lsoas$X2, listw = lw) # repeat process - vax and previous infection (for some reason R changes variables from 1 to X1)
lcr_lsoas <- cbind(lcr_lsoas, data.matrix(localGi))
names(lcr_lsoas)[names(lcr_lsoas) == "data.matrix.localGi."] <- "gstat_delta1_v2" 

localGi <- localG(lcr_lsoas$X3, listw = lw) # repeat process - no vax and no previous infection
lcr_lsoas <- cbind(lcr_lsoas, data.matrix(localGi))
names(lcr_lsoas)[names(lcr_lsoas) == "data.matrix.localGi."] <- "gstat_delta1_v3" 

localGi <- localG(lcr_lsoas$X4, listw = lw) # repeat process - no vax and previous infection
lcr_lsoas <- cbind(lcr_lsoas, data.matrix(localGi))
names(lcr_lsoas)[names(lcr_lsoas) == "data.matrix.localGi."] <- "gstat_delta1_v4" 

# Make long format
lcr_lsoas_long <- melt(lcr_lsoas, id.vars = c("lsoa11cd", "geometry"), measure = 6:9)
lcr_lsoas_long = st_as_sf(lcr_lsoas_long, sf_column_name = 'geometry') # Convert back to sf object


# Time period labels
labels <- c("Fully vaccinated & no previous positive", "Fully vaccinated & previous positive", "Not fully vaccinated & no previous positive", "Not fully vaccinated & previous positive")
names(labels) <- c("gstat_delta1_v1", "gstat_delta1_v2", "gstat_delta1_v3", "gstat_delta1_v4")

# Map
map_inf_delta1 <- ggplot() +
  geom_sf(data = lcr_lsoas_long, aes(fill = value), color = NA) +
  scale_fill_viridis(option = "turbo", # Make colour blind friendly
                     begin = 0.1, end = 0.9, # Set colour range
                     limits = c(-4.25, 4.25), oob = scales::squish) + # Define legend values to show on plot
  facet_wrap(~variable, labeller = labeller(variable = labels)) +
  xlab("Longitude") + # Add x-axis label
  ylab("Latitude") + # Add y-axis label
  labs(title = "Delta (3rd June - 1st September 2021)", # Edit plot title
       fill = "Gi statistic") + # Edit legend title (note must match fill as that is what we are plotting)
  theme(legend.position="bottom") +
  theme_void()


# Delta 2 = spatial analysis

# Join to spatial
# lcr_lsoas <- read_sf(dsn = "./Shapefiles/E47000004.shp") # Load LSOA shapefiles - windows
lcr_lsoas <- read_sf(dsn = "R:/SMART_LCR_report/Shapefiles/E47000004.shp") # Load LSOA shapefiles - windows
lcr_lsoas <- merge(lcr_lsoas, lsoa_inf_delta2, by.x = "lsoa11cd", by.y = "LSOA_Code", all.x = TRUE)
lcr_lsoas <- st_transform(lcr_lsoas, 4326) # Set CRS
lcr_lsoas.nb <- poly2nb(lcr_lsoas, snap=0.0002) # Identify neighbours of each LSOA

# Define adjacency using a row-standardised matrix
lw <- nb2listw(lcr_lsoas.nb)
W <- as(as_dgRMatrix_listw(lw), "CsparseMatrix")

# Calculate Getis-Ord statistic for SMR
localGi <- localG(lcr_lsoas$`1`, listw = lw) # Calculate Gi statistics - vax and no previous infection
lcr_lsoas <- cbind(lcr_lsoas, data.matrix(localGi)) # Join results onto shapefile
names(lcr_lsoas)[names(lcr_lsoas) == "data.matrix.localGi."] <- "gstat_delta2_v1" # Rename column

localGi <- localG(lcr_lsoas$X2, listw = lw) # repeat process - vax and previous infection (for some reason R changes variables from 1 to X1)
lcr_lsoas <- cbind(lcr_lsoas, data.matrix(localGi))
names(lcr_lsoas)[names(lcr_lsoas) == "data.matrix.localGi."] <- "gstat_delta2_v2" 

localGi <- localG(lcr_lsoas$X3, listw = lw) # repeat process - no vax and no previous infection
lcr_lsoas <- cbind(lcr_lsoas, data.matrix(localGi))
names(lcr_lsoas)[names(lcr_lsoas) == "data.matrix.localGi."] <- "gstat_delta2_v3" 

localGi <- localG(lcr_lsoas$X4, listw = lw) # repeat process - no vax and previous infection
lcr_lsoas <- cbind(lcr_lsoas, data.matrix(localGi))
names(lcr_lsoas)[names(lcr_lsoas) == "data.matrix.localGi."] <- "gstat_delta2_v4" 

# Make long format
lcr_lsoas_long <- melt(lcr_lsoas, id.vars = c("lsoa11cd", "geometry"), measure = 6:9)
lcr_lsoas_long = st_as_sf(lcr_lsoas_long, sf_column_name = 'geometry') # Convert back to sf object


# Time period labels
labels <- c("Fully vaccinated & no previous positive", "Fully vaccinated & previous positive", "Not fully vaccinated & no previous positive", "Not fully vaccinated & previous positive")
names(labels) <- c("gstat_delta2_v1", "gstat_delta2_v2", "gstat_delta2_v3", "gstat_delta2_v4")

# Map
map_inf_delta2 <- ggplot() +
  geom_sf(data = lcr_lsoas_long, aes(fill = value), color = NA) +
  scale_fill_viridis(option = "turbo", # Make colour blind friendly
                     begin = 0.1, end = 0.9, # Set colour range
                     limits = c(-4.25, 4.25), oob = scales::squish) + # Define legend values to show on plot
  facet_wrap(~variable, labeller = labeller(variable = labels)) +
  xlab("Longitude") + # Add x-axis label
  ylab("Latitude") + # Add y-axis label
  labs(title = "Delta (1st September - 27th Novemeber 2021)", # Edit plot title
       fill = "Gi statistic") + # Edit legend title (note must match fill as that is what we are plotting)
  theme(legend.position="bottom") +
  theme_void()


# Omnicron = spatial analysis

# Join to spatial
# lcr_lsoas <- read_sf(dsn = "./Shapefiles/E47000004.shp") # Load LSOA shapefiles - windows
lcr_lsoas <- read_sf(dsn = "R:/SMART_LCR_report/Shapefiles/E47000004.shp") # Load LSOA shapefiles - windows
lcr_lsoas <- merge(lcr_lsoas, lsoa_inf_omni, by.x = "lsoa11cd", by.y = "LSOA_Code", all.x = TRUE)
lcr_lsoas <- st_transform(lcr_lsoas, 4326) # Set CRS
lcr_lsoas.nb <- poly2nb(lcr_lsoas, snap=0.0002) # Identify neighbours of each LSOA

# Define adjacency using a row-standardised matrix
lw <- nb2listw(lcr_lsoas.nb)
W <- as(as_dgRMatrix_listw(lw), "CsparseMatrix")

# Calculate Getis-Ord statistic for SMR
localGi <- localG(lcr_lsoas$`1`, listw = lw) # Calculate Gi statistics - vax and no previous infection
lcr_lsoas <- cbind(lcr_lsoas, data.matrix(localGi)) # Join results onto shapefile
names(lcr_lsoas)[names(lcr_lsoas) == "data.matrix.localGi."] <- "gstat_omni_v1" # Rename column

localGi <- localG(lcr_lsoas$X2, listw = lw) # repeat process - vax and previous infection (for some reason R changes variables from 1 to X1)
lcr_lsoas <- cbind(lcr_lsoas, data.matrix(localGi))
names(lcr_lsoas)[names(lcr_lsoas) == "data.matrix.localGi."] <- "gstat_omni_v2" 

localGi <- localG(lcr_lsoas$X3, listw = lw) # repeat process - no vax and no previous infection
lcr_lsoas <- cbind(lcr_lsoas, data.matrix(localGi))
names(lcr_lsoas)[names(lcr_lsoas) == "data.matrix.localGi."] <- "gstat_omni_v3" 

localGi <- localG(lcr_lsoas$X4, listw = lw) # repeat process - no vax and previous infection
lcr_lsoas <- cbind(lcr_lsoas, data.matrix(localGi))
names(lcr_lsoas)[names(lcr_lsoas) == "data.matrix.localGi."] <- "gstat_omni_v4" 

# Make long format
lcr_lsoas_long <- melt(lcr_lsoas, id.vars = c("lsoa11cd", "geometry"), measure = 6:9)
lcr_lsoas_long = st_as_sf(lcr_lsoas_long, sf_column_name = 'geometry') # Convert back to sf object


# Time period labels
labels <- c("Fully vaccinated & no previous positive", "Fully vaccinated & previous positive", "Not fully vaccinated & no previous positive", "Not fully vaccinated & previous positive")
names(labels) <- c("gstat_omni_v1", "gstat_omni_v2", "gstat_omni_v3", "gstat_omni_v4")

# Map
map_inf_omni <- ggplot() +
  geom_sf(data = lcr_lsoas_long, aes(fill = value), color = NA) +
  scale_fill_viridis(option = "turbo", # Make colour blind friendly
                     begin = 0.1, end = 0.9, # Set colour range
                     limits = c(-4.25, 4.25), oob = scales::squish) + # Define legend values to show on plot
  facet_wrap(~variable, labeller = labeller(variable = labels)) +
  xlab("Longitude") + # Add x-axis label
  ylab("Latitude") + # Add y-axis label
  labs(title = "Omicron (13th December 2021 - 2nd February 2022)", # Edit plot title
       fill = "Gi statistic") + # Edit legend title (note must match fill as that is what we are plotting)
  theme(legend.position="bottom") +
  theme_void()


# Plot
map_inf_delta1
ggsave(plot = map_inf_delta1, filename = "./Vaccine v reinfection paper/clusters_infections_delta1.jpeg") # Save
ggsave(plot = map_inf_delta1, filename = "./Vaccine v reinfection paper/clusters_infections_delta1.svg")
map_inf_delta2
ggsave(plot = map_inf_delta2, filename = "./Vaccine v reinfection paper/clusters_infections_delta2.jpeg") # Save
ggsave(plot = map_inf_delta2, filename = "./Vaccine v reinfection paper/clusters_infections_delta2.svg")
map_inf_omni
ggsave(plot = map_inf_omni, filename = "./Vaccine v reinfection paper/clusters_infections_omni.jpeg") # Save
ggsave(plot = map_inf_omni, filename = "./Vaccine v reinfection paper/clusters_infections_omni.svg")

# ggsave(plot = inf_map, filename = "./Vaccine v reinfection paper/clusters_infections.jpeg") # Save


# 3b. Deaths (under 75) # 

# Calculate observed deaths per LSOA
hold <- data2 # Save in case need
data2 <- data2[data2$Patient.Age < 75] # Remove if want to have all ages
lsoa_deaths <- data2[, list(death_delta1 = sum(death_delta1, na.rm = TRUE), death_delta2 = sum(death_delta2, na.rm = TRUE), death_omni = sum(death_omni, na.rm = TRUE), pop = .N), by = "LSOA_Code"] # Aggregate to LSOA

# Create expected count (age standardised)
lsoa_age <- data2[, list(pop = .N), by = c("LSOA_Code", "age_band")] # Pop by age band for LSOAs
std_pop <- data2[, list(death_delta1 = sum(death_delta1, na.rm = TRUE), death_delta2 = sum(death_delta2, na.rm = TRUE), death_omni = sum(death_omni, na.rm = TRUE), pop = .N), by = c("age_band")] # Standard pop
std_pop$delta1_rate <- (std_pop$death_delta1 / std_pop$pop) # Calculate rates
std_pop$delta2_rate <- (std_pop$death_delta2 / std_pop$pop) 
std_pop$omni_rate <- (std_pop$death_omni / std_pop$pop) 
std_pop$pop <- NULL # Drop as not needed
lsoa_age <- merge(lsoa_age, std_pop, by = "age_band", all.x = TRUE) # Join together
lsoa_age$exp_delta1_deaths <- lsoa_age$pop * lsoa_age$delta1_rate # Expected deaths by age band and LSOA
lsoa_age$exp_delta2_deaths <- lsoa_age$pop * lsoa_age$delta2_rate
lsoa_age$exp_omni_deaths <- lsoa_age$pop * lsoa_age$omni_rate
exp_deaths <- lsoa_age[, list(exp_delta1_deaths = sum(exp_delta1_deaths, na.rm = TRUE), exp_delta2_deaths = sum(exp_delta2_deaths, na.rm = TRUE), exp_omni_deaths = sum(exp_omni_deaths, na.rm = TRUE)), by = c("LSOA_Code")] # Sum expected deaths for LSOAs
lsoa_deaths <- merge(lsoa_deaths, exp_deaths, by = "LSOA_Code", all.x = TRUE) # Join together
lsoa_deaths$delta1_smr <- (lsoa_deaths$death_delta1 / lsoa_deaths$exp_delta1_deaths) * 100 # Create SMRs
lsoa_deaths$delta2_smr <- (lsoa_deaths$death_delta2 / lsoa_deaths$exp_delta2_deaths) * 100
lsoa_deaths$omni_smr <- (lsoa_deaths$death_omni / lsoa_deaths$exp_omni_deaths) * 100
rm(std_pop, lsoa_age, exp_deaths) # Tidy

# Join to spatial
# lcr_lsoas <- read_sf(dsn = "./Shapefiles/E47000004.shp") # Load LSOA shapefiles - windows
lcr_lsoas <- read_sf(dsn = "R:/SMART_LCR_report/Shapefiles/E47000004.shp") # Load LSOA shapefiles - windows
lcr_lsoas <- merge(lcr_lsoas, lsoa_deaths, by.x = "lsoa11cd", by.y = "LSOA_Code", all.x = TRUE)
lcr_lsoas <- st_transform(lcr_lsoas, 4326) # Set CRS
lcr_lsoas.nb <- poly2nb(lcr_lsoas, snap=0.0002) # Identify neighbours of each LSOA

# Define adjacency using a row-standardised matrix
lw <- nb2listw(lcr_lsoas.nb)
W <- as(as_dgRMatrix_listw(lw), "CsparseMatrix")

# Calculate Getis-Ord statistic for SMR
localGi <- localG(lcr_lsoas$delta1_smr, listw = lw) # Calculate Gi statistics - delta 1
lcr_lsoas <- cbind(lcr_lsoas, data.matrix(localGi)) # Join results onto shapefile
names(lcr_lsoas)[names(lcr_lsoas) == "data.matrix.localGi."] <- "gstat_smr_delta1" # Rename column
localGi <- localG(lcr_lsoas$delta2_smr, listw = lw) # repeat process - delta 2
lcr_lsoas <- cbind(lcr_lsoas, data.matrix(localGi))
names(lcr_lsoas)[names(lcr_lsoas) == "data.matrix.localGi."] <- "gstat_smr_delta2" 
localGi <- localG(lcr_lsoas$omni_smr, listw = lw) # repeat process - omnicron
lcr_lsoas <- cbind(lcr_lsoas, data.matrix(localGi))
names(lcr_lsoas)[names(lcr_lsoas) == "data.matrix.localGi."] <- "gstat_smr_omni" 


# Map
map_deaths_delta1 <- ggplot() +
  geom_sf(data = lcr_lsoas, aes(fill = gstat_smr_delta1), color = NA) +
  scale_fill_viridis(option = "turbo", # Make colour blind friendly
                     limits = c(-3, 4.5), oob = scales::squish) + # Define legend values to show on plot
  xlab("Longitude") + # Add x-axis label
  ylab("Latitude") + # Add y-axis label
  labs(title = "Delta (3/6 - 1/9/21)", # Edit plot title
       fill = "Gi statistic") + # Edit legend title (note must match fill as that is what we are plotting)
  theme_void()

map_deaths_delta2 <- ggplot() +
  geom_sf(data = lcr_lsoas, aes(fill = gstat_smr_delta2), color = NA) +
  scale_fill_viridis(option = "turbo", # Make colour blind friendly
                     limits = c(-3, 4.5), oob = scales::squish) + # Define legend values to show on plot
  xlab("Longitude") + # Add x-axis label
  ylab("Latitude") + # Add y-axis label
  labs(title = "Delta (1/9 - 27/11/21)", # Edit plot title
       fill = "Gi statistic") + # Edit legend title (note must match fill as that is what we are plotting)
theme_void()

map_deaths_omni <- ggplot() +
  geom_sf(data = lcr_lsoas, aes(fill = gstat_smr_omni), color = NA) +
  scale_fill_viridis(option = "turbo", # Make colour blind friendly
                     limits = c(-3, 4.5), oob = scales::squish) + # Define legend values to show on plot
  xlab("Longitude") + # Add x-axis label
  ylab("Latitude") + # Add y-axis label
  labs(title = "Omicron (13/12/21 - 2/2/22)", # Edit plot title
       fill = "Gi statistic") + # Edit legend title (note must match fill as that is what we are plotting)
theme_void()

# ggsave(plot = map_deaths_omni, filename = "./Vaccine v reinfection paper/clusters_deaths_omicron.jpeg") # Save

deaths_map <- map_deaths_delta1 + map_deaths_delta2 + map_deaths_omni + # Join maps together
  plot_layout(guides = "collect") + # Use same legends
  plot_annotation(title = "Deaths in Cheshire and Merseyside",
    subtitle = "Spatial clustering of small area standardised mortality ratios for all-cause mortality (<75 years).")
deaths_map

ggsave(plot = deaths_map, filename = "./Vaccine v reinfection paper/clusters_deaths.jpeg") # Save
ggsave(plot = deaths_map, filename = "./Vaccine v reinfection paper/clusters_deaths.svg")

# # 3b. Deaths by covid or not # 
# 
# # Calculate observed deaths per LSOA
# lsoa_deaths_v2 <- hold[, list(covid_deaths = sum(died_covid, na.rm = TRUE), noncovid_deaths = sum(died_noncovid, na.rm = TRUE), pop = .N), by = "LSOA_Code"] # Aggregate to LSOA
# 
# # Create expected count (age standardised)
# lsoa_age <- hold[, list(pop = .N), by = c("LSOA_Code", "age_band")] # Pop by age band for LSOAs
# std_pop <- hold[, list(covid_deaths = sum(died_covid, na.rm = TRUE), noncovid_deaths = sum(died_noncovid, na.rm = TRUE), pop = .N), by = c("age_band")] # Standard pop
# std_pop$rate_covid <- (std_pop$covid_deaths / std_pop$pop) # Calculate rate
# std_pop$rate_noncovid <- (std_pop$noncovid_deaths / std_pop$pop) 
# lsoa_age <- merge(lsoa_age, std_pop, by = "age_band", all.x = TRUE) # Join together
# lsoa_age$exp_covid_deaths <- lsoa_age$pop.x * lsoa_age$rate_covid # Expected covid deaths
# lsoa_age$exp_noncovid_deaths <- lsoa_age$pop.x * lsoa_age$rate_noncovid # Expected non covid deaths
# exp_deaths <- lsoa_age[, list(exp_covid_deaths = sum(exp_covid_deaths, na.rm = TRUE), exp_noncovid_deaths = sum(exp_noncovid_deaths, na.rm = TRUE)), by = c("LSOA_Code")] # Sum expected deaths for LSOAs
# lsoa_deaths_v2 <- merge(lsoa_deaths_v2, exp_deaths, by = "LSOA_Code", all.x = TRUE) # Join together
# lsoa_deaths_v2$smr_covid <- (lsoa_deaths_v2$covid_deaths / lsoa_deaths_v2$exp_covid_deaths) * 100# Create SMR
# lsoa_deaths_v2$smr_noncovid <- (lsoa_deaths_v2$noncovid_deaths / lsoa_deaths_v2$exp_noncovid_deaths) * 100# Create SMR
# rm(hold, std_pop, lsoa_age, exp_deaths) # Tidy
# 
# # Join to spatial
# lcr_lsoas <- read_sf(dsn = "./Shapefiles/E47000004.shp") # Load LSOA shapefiles - windows
# # lcr_lsoas <- read_sf(dsn = "R:/SMART_LCR_report/Shapefiles/E47000004.shp") # Load LSOA shapefiles - windows
# lcr_lsoas <- merge(lcr_lsoas, lsoa_deaths_v2, by.x = "lsoa11cd", by.y = "LSOA_Code", all.x = TRUE)
# lcr_lsoas <- st_transform(lcr_lsoas, 4326) # Set CRS
# lcr_lsoas.nb <- poly2nb(lcr_lsoas, snap=0.0002) # Identify neighbours of each LSOA
# 
# # Define adjacency using a row-standardised matrix
# lw <- nb2listw(lcr_lsoas.nb)
# W <- as(as_dgRMatrix_listw(lw), "CsparseMatrix")
# 
# # Calculate Getis-Ord statistic for SMR - covid
# localGi <- localG(lcr_lsoas$smr_covid, listw = lw) # Calculate Gi statistics
# lcr_lsoas <- cbind(lcr_lsoas, data.matrix(localGi)) # Join results onto shapefile
# names(lcr_lsoas)[names(lcr_lsoas) == "data.matrix.localGi."] <- "gstat_covid" # Rename column
# 
# # Calculate Getis-Ord statistic for SMR - non-covid
# localGi <- localG(lcr_lsoas$smr_noncovid, listw = lw) # Calculate Gi statistics
# lcr_lsoas <- cbind(lcr_lsoas, data.matrix(localGi)) # Join results onto shapefile
# names(lcr_lsoas)[names(lcr_lsoas) == "data.matrix.localGi."] <- "gstat_noncovid" # Rename column
# 
# # Make long
# lcr_lsoas_long <- melt(lcr_lsoas, id.vars = c("lsoa11cd", "geometry"), measure = 9:10)
# lcr_lsoas_long = st_as_sf(lcr_lsoas_long, sf_column_name = 'geometry') # Convert back to sf object
# 
# # Time period labels
# labels <- c("COVID-19 deaths", "Non COVID-19 deaths")
# names(labels) <- c("gstat_covid", "gstat_noncovid")
# 
# # Map
# map2 <- ggplot() +
#   geom_sf(data = lcr_lsoas_long, aes(fill = value), lwd = 0) +
#   scale_fill_viridis(option = "turbo") + # Make colour blind friendly
#   facet_wrap(~variable, labeller = labeller(variable = labels)) +
#   xlab("Longitude") + # Add x-axis label
#   ylab("Latitude") + # Add y-axis label
#   labs(title = "Deaths since 1st September 2021 (Age-standardised)", # Edit plot title
#        fill = "Gi statistic") + # Edit legend title (note must match fill as that is what we are plotting)
#   theme(legend.position="bottom")
# map2
# ggsave(plot = map2, filename = "./Vaccine v reinfection paper/clusters_deaths_since_sept_bytype.jpeg")
# 


