# Joe Chudzik, Reid Holben
# joseph.chudzik@marquette.edu, reid.holben@marquette.edu
# Data Science Capstone Project - GoKidGoWeb
# keywoods: simulation, data wrangling

# Clean the global library.
rm(list=ls())

# NOTE:
#   Set your working directory to the base location of the repository.
setwd('~/Documents/School/DataScienceCapstone/')

set.seed(42)

# Build a simulated dataset of about 400 individuals

# Name
names <- read.csv('Data/DataForSimulation/randomNames.csv')

# Home address
# Original dataset contains over 300,000 entries (taken from Milwaukee's OpenData Master Address Index (MAI)).
# Following code splits the dataset in half to avoid home addresses the same as destination and vice versa.
#   then creates a random sample of 400 for both home and destination addresses from the master file.
# Code saves the new dataset to be later used for analysis.

#addresses <- read.csv('Data/DataForSimulation/addressList_master.csv')
#halfway <- round(nrow(addresses)/2)
#first_half_addresses <- addresses[1:halfway,]
#second_half_addresses <- addresses[halfway+1:nrow(addresses),]
#second_half_addresses <- na.omit(second_half_addresses)

#home_addr <- first_half_addresses[sample(nrow(first_half_addresses), 400),]
#destination_addr <- second_half_addresses[sample(nrow(second_half_addresses), 400),]

# Save the newly partitioned addresses to load quicker.
#write.csv(home_addr, '~/Documents/School/DataScienceCapstone/Data/DataForSimulation/home_addr.csv')
#write.csv(destination_addr, '~/Documents/School/DataScienceCapstone/Data/DataForSimulation/destination_addr.csv')

# Load the simulated home addresses and clean the data to correct address format.
homes <- read.csv('Data/DataForSimulation/home_addr.csv')
home_addr <- data.frame('HSE_NBR_home' = homes$HSE_NBR, 'STREET_home' = homes$STREET, 'STTYPE_home' = homes$STTYPE, 'UNIT_NBR_home' = homes$UNIT_NBR, 'ZIP_CODE_home' = homes$ZIP_CODE)

# Load in the destinations.
destinations <- read.csv('Data/DataForSimulation/destinations.csv')

# Repeat new destinations 20 times to match the 400 people in the names dataset. 
destinations <- destinations[rep(1:nrow(destinations),each=20),]

# Shuffle the row indeces of the destinations dataframe.
rows <- sample(nrow(destinations))

# Use the random vector above to reorder the destinations dataset.
destinations <- destinations[rows,]

# Time to pick up from home. First simulate some time data from 7am-10am then create dataframe with only the time (date may not be important).
# hpu = home pick up

homePickup_start <- as.POSIXct('07:00z', format='%H:%M')
homePickup_end <- as.POSIXct('10:00z', format='%H:%M')
seconds_hpu <- difftime(homePickup_end, homePickup_start, units='secs')
difference_hpu <- sample(1:seconds_hpu, 400, replace=T)
homePickup_pickupTime <- homePickup_start + difference_hpu

homePickup <- data.frame('homePickup' = format(strptime(homePickup_pickupTime, format='%Y-%m-%d %H:%M:%S'), '%H:%M'))

# Time to drop off at destination.
# ddo = destination drop off

destDropOff_start <- as.POSIXct('09:00z', format='%H:%M')
destDropOff_end <- as.POSIXct('11:00z', format='%H:%M')
seconds_ddo <- difftime(destDropOff_end, destDropOff_start, units='secs')
difference_ddo <- sample(1:seconds_ddo, 400, replace=T)
destination_dropoffTime <- destDropOff_start + difference_ddo

destinationDropoff <- data.frame('destinationDropOff' = format(strptime(destination_dropoffTime, format='%Y-%m-%d %H:%M:%S'), '%H:%M'))

# Time to pick up from destination.
# dpu = destination pick up

destPickUp_start <- as.POSIXct('13:00z', format='%H:%M')
destPickUp_end <- as.POSIXct('17:00z', format='%H:%M')
seconds_dpu <- difftime(destPickUp_end, destPickUp_start, units='secs')
difference_dpu <- sample(1:seconds_dpu, 400, replace=T)
destination_pickupTime <- destPickUp_start + difference_dpu

destinationPickup <- data.frame('destinationPickup' = format(strptime(destination_pickupTime, format='%Y-%m-%d %H:%M:%S'), '%H:%M'))

# Time to drop off at home. (not sure if this is needed)
# hdo = destination drop off

homeDropoff_start <- as.POSIXct('13:15z', format='%H:%M')
homeDropoff_end <- as.POSIXct('18:00z', format='%H:%M')
seconds_hdo <- difftime(homeDropoff_end, homeDropoff_start, units='secs')
difference_hdo <- sample(1:seconds_hdo, 400, replace=T)
home_dropoffTime <- homeDropoff_start + difference_hdo

homeDropoff <- data.frame('homeDropOff' = format(strptime(home_dropoffTime, format='%Y-%m-%d %H:%M:%S'), '%H:%M'))


# Putting everything together.
#   names, home_addr, dest_addr, homePickup, destinationDropoff, destinationPickup, homeDropoff

simulatedData <- cbind(names, home_addr, destinations, homePickup, destinationDropoff, destinationPickup, homeDropoff)

# Write the newly simulated data to data folder.
#write.csv(simulatedData, '~/Documents/GitHub/DataScienceCapstone/Data/simulatedData.csv')







