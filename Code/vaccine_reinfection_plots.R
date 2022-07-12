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
library(patchwork)
library(zoo) # For rolling average


### 1. Tidy data ###

#  Load data
# load("./Data/cleaned_data.RData") 


### 2. Plots ####


# Figure 1 #

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

# Define plot labels for study period 
text_plot <- data.frame(text = c("Delta 1", "Delta 2", "Omicron"), start_date = as.Date(c("2021-06-03", "2021-09-01", "2021-12-13")), end_date = as.Date(c("2021-08-31", "2021-11-26", "2022-02-28")), ymin = c(0, 0, 0), ymax = c(max(trend_data_long$freq, na.rm=T), max(trend_data_long$freq, na.rm=T), max(trend_data_long$freq, na.rm=T)), stringsAsFactors = FALSE)

# Plot
plot1 <- ggplot() +
  geom_rect(data = text_plot, mapping = aes(xmin = start_date, xmax = end_date, ymin = ymin, ymax = ymax), color = "grey", alpha = 0.1) + # Add box for period
  geom_text(mapping = aes(x = start_date, y = 500+ymax, label = text, hjust = 0), size = 3.4, data = text_plot) + # Add in text labels for period (I have added some extra distance to the y value so it is spaced out vs the rectangles above)
  geom_line(data = trend_data_long[trend_data_long$day < "2022-03-01",], aes(x = day, y = freq, group = Infection, color = Infection)) +
  # geom_point(data = trend_data_long[trend_data_long$Infection == "New_7day",], aes(x = day, y = New, group = Infection, color = Infection), alpha = 0.1) +
  # geom_point(data = trend_data_long[trend_data_long$Infection == "Reinfection_7day",], aes(x = day, y = Reinfection, group = Infection, color = Infection), alpha = 0.1) +
  xlab("Date") +
  ylab("Frequency") +
  scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9, labels = c("First positive", "Further positives")) +
  #labs(title = "Number of COVID-19 registered positive tests",
  #     subtitle = "The line represents the 7-day moving average, with points the daily value.") +
  theme(text = element_text(size = 16), plot.title = element_text(size = 20, face = "bold"),)  +
  scale_x_date(breaks = "3 months", minor_breaks = "1 month", date_labels = "%b %y") +
  ylim(0, 500+max(trend_data_long$freq, na.rm=T))
plot1
# ggsave(plot = plot1, filename = "./Vaccine v reinfection paper/reinfection_trends.jpeg") # Save
# ggsave(plot = plot1, filename = "./Vaccine v reinfection paper/reinfection_trends.svg") 

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

# Define plot labels for study period 
text_plot2 <- data.frame(text = c("Delta 1", "Delta 2", "Omicron"), start_date = as.Date(c("2021-06-03", "2021-09-01", "2021-12-13")), end_date = as.Date(c("2021-08-31", "2021-11-26", "2022-02-28")), ymin = c(0, 0, 0), ymax = c(max(trend_data_long$rate, na.rm=T), max(trend_data_long$rate, na.rm=T), max(trend_data_long$rate, na.rm=T)), stringsAsFactors = FALSE)

# Plot
plot2 <- ggplot() +
  geom_rect(data = text_plot2, mapping = aes(xmin = start_date, xmax = end_date, ymin = ymin, ymax = ymax), color = "grey", alpha = 0.1) + # Add box for period
  geom_text(mapping = aes(x = start_date, y = 20+ymax, label = text, hjust = 0), size = 3.4, data = text_plot2) + # Add in text labels for period
  geom_line(data = trend_data_long[trend_data_long$day < "2022-03-01",], aes(x = day, y = rate, group = Infection, color = Infection)) +
  xlab("Date") +
  ylab("Rate per 100,000") +
  scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9, labels = c("First positive", "Further positives")) +
  #labs(title = "Number of COVID-19 registered positive tests per capita",
  #     subtitle = "The line represents the 7-day moving average.") +
  theme(text = element_text(size = 16), plot.title = element_text(size = 20, face = "bold"))  +
  scale_x_date(breaks = "3 months", minor_breaks = "1 month", date_labels = "%b %y") +
  ylim(0, 20+max(trend_data_long$rate, na.rm=T))
plot2
# ggsave(plot = plot2, filename = "./Vaccine v reinfection paper/reinfection_trends_rate.jpeg") # Save
# ggsave(plot = plot2, filename = "./Vaccine v reinfection paper/reinfection_trends_rate.svg") 
# hjust = -0.05, vjust = -0.5

# Combine plot
fig1 <- plot1 / 
  plot2 +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") 
fig1
ggsave(plot = fig1, filename = "./Vaccine v reinfection paper/figure1.jpeg") # Save
ggsave(plot = fig1, filename = "./Vaccine v reinfection paper/figure1_hires.tiff", dpi = 300)
ggsave(plot = fig1, filename = "./Vaccine v reinfection paper/figure1.svg")


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
ses_plot1 <- ggplot(data = trend_data2_long[(trend_data2_long$imd_decile == 1 | trend_data2_long$imd_decile == 10) & trend_data2_long$day < "2022-03-01",], aes(x = day, y = rate, group = factor(imd_decile), color = factor(imd_decile))) +
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
ggsave(plot = ses_plot1, filename = "./Vaccine v reinfection paper/ses_trends_whole_hires.tiff", dpi = 300) 
ggsave(plot = ses_plot1, filename = "./Vaccine v reinfection paper/ses_trends_whole.svg") 

# Plot 2: Fix date for period

# Plot
ses_plot2 <- ggplot() +
  geom_rect(data = text_plot2, mapping = aes(xmin = start_date, xmax = end_date, ymin = ymin, ymax = ymax), color = "grey", alpha = 0.1) + # Add box for period
  geom_text(mapping = aes(x = start_date, y = 400, label = text, hjust = 0), size = 5, data = text_plot2) + # Add in text labels for period
  geom_line(data = trend_data2_long[(trend_data2_long$day >= "2021-06-03" & trend_data2_long$day < "2022-03-01") & (trend_data2_long$imd_decile == 1 | trend_data2_long$imd_decile == 10),], aes(x = day, y = rate, group = factor(imd_decile), color = factor(imd_decile)), size = 0.8) +
  facet_wrap(~Infection, scales = "fixed", labeller = labeller(Infection = labels)) +
  scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9, labels = c("1: Most deprived", "10: Least Deprived")) +
  labs(color = "Deprivation decile") +
  xlab("Date") +
  ylab("Registered positive tests per 100,000 people") +
  ylim(0, 410) +
  scale_x_date(breaks = "1 month", date_labels = "%b") + # Define labels for x-axis
  theme(text = element_text(size = 16), plot.title = element_text(size = 20, face = "bold")) + # Make text size bigger
  theme(legend.position="bottom")
ses_plot2
ggsave(plot = ses_plot2, filename = "./Vaccine v reinfection paper/figure3.jpeg") # Save
ggsave(plot = ses_plot2, filename = "./Vaccine v reinfection paper/figure3.tiff", dpi = 300) 
ggsave(plot = ses_plot2, filename = "./Vaccine v reinfection paper/figure3.svg") 

# Plot 3: Percent of reinfections by decile


# Define plot labels for study period 
text_plot3 <- data.frame(text = c("Delta 1", "Delta 2", "Omicron"), start_date = as.Date(c("2021-06-03", "2021-09-01", "2021-12-13")), end_date = as.Date(c("2021-08-31", "2021-11-26", "2022-02-28")), ymin = c(0, 0, 0), ymax = c(20, 20, 20), stringsAsFactors = FALSE)

# Define transparency for lines
trend_data2$alpha <- NA # Create blank variable
trend_data2$alpha[trend_data2$imd_decile == 1 | trend_data2$imd_decile == 10] <- 1 # Make so stand out
trend_data2$alpha[trend_data2$imd_decile > 1 & trend_data2$imd_decile < 10] <- 0.8 # Make see through

# Plot
ses_plot3 <- ggplot() + # Plot
  geom_rect(data = text_plot3, mapping = aes(xmin = start_date, xmax = end_date, ymin = ymin, ymax = 19), color = "grey", alpha = 0.1) + # Add box for period
  geom_text(mapping = aes(x = start_date, y = ymax, label = text, hjust = 0), size = 5, data = text_plot3) + # Add in text labels for period
  geom_line(data = trend_data2[trend_data2$day >= "2021-06-03" & trend_data2$day < "2022-03-01",], aes(x = day, y = pc_reinf, group = factor(imd_decile), color = factor(imd_decile), alpha = alpha), size = 0.8) +
  scale_color_viridis_d(option = "turbo", begin = 0.1, end = 0.9) + # Make colouur blind friendly
  labs(color = "Deprivation decile") + # Define labels for plot
  xlab("Date") +
  ylab("Percent of all positive tests that are further positive tests (%)") +
  ylim(0,20) + # Set limits of y-axis
  scale_x_date(breaks = "1 month", date_labels = "%b") + # Define labels for x-axis
  theme(text = element_text(size = 16), plot.title = element_text(size = 20, face = "bold")) + # Make text size bigger
  scale_alpha_continuous(guide = "none") # Drop alpha legend
ses_plot3
ggsave(plot = ses_plot3, filename = "./Vaccine v reinfection paper/figure2.jpeg") # Save
ggsave(plot = ses_plot3, filename = "./Vaccine v reinfection paper/figure2.tiff", dpi = 300)
ggsave(plot = ses_plot3, filename = "./Vaccine v reinfection paper/figure2.svg") 


