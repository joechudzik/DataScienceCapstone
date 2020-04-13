library(mapsapi)
library(leaflet)
library(Matrix)
library(blockcluster)

# Used to load the R file with the custom hclust method
# source('hierarchicalClust.R')

# Clear the global library.
rm(list=ls())

data <- read.csv('~/Documents/GitHub/DataScienceCapstone/Data/simulatedData.csv')


# Simulating better times for improved algorithm
# 
ept_range <- seq(as.POSIXct('7:00z', format='%H:%M'), as.POSIXct('10:00z', format='%H:%M'), by=(60*15))
edt_range <- seq(as.POSIXct('12:00z', format='%H:%M'), as.POSIXct('17:00z', format='%H:%M'), by=(60*15))

# repeat ept_range & edt_range and sample it to randomize. Limit the vector to 400 elements
ept <- sample(rep(ept_range, 31)); ept <- ept[1:400]
edt <- sample(rep(edt_range, 31)); edt <- edt[1:400]

# initialize latest pickup time and latest dropoff time with corresponding earlier values (just for memory allocation)
lpt <- ept
ldt <- edt

# offset created by selecting a random number between 0:12 (multiple of 5 minutes within an hour), 
#   multiply by amount of second in a minute, and add it to ept & edt to create lpt or ldt
for(x in 1:length(ept)){
  offset_pickup <- sample(0:12, 1)*5*60
  offset_dropoff <- sample(0:12, 1)*5*60
  lpt[x] <- ept[x] + offset_pickup
  ldt[x] <- edt[x] + offset_pickup
}

dataForSpecificTimeAlgorithm <- cbind(data[,1:14], ept, lpt, edt, ldt)
#write.csv(dataForSpecificTimeAlgorithm, '~/Documents/GitHub/DataScienceCapstone/Data/dataForSpecificTimeAlgorithm.csv')


data <- dataForSpecificTimeAlgorithm


#
# AS OF NOW, API DOES NOT ALLOW MORE THAN 10 ENTRIES FOR DISTANCE CALCULATIONS IN MATRIX
#
# Subset original simulated data for only the fierst 10 entries
data <- data[1:10,]

# Create vectors for pickup and dropoff locations
origin_vector <- paste0(data$HSE_NBR_home,' ', data$STREET_home, ' ', data$STTYPE_home, ', Milwaukee, WI ', data$ZIP_CODE_home)
destination_vector <- paste0(data$HSE_NBR,' ', data$STREET, ' ', data$STTYPE, ', Milwaukee, WI ', data$ZIP_CODE)

# Build matrices from Google's distance API
# Naming convention = "rowByColumn"
originByOriginMatrix <- mp_matrix(
  origins = origin_vector,
  destinations = origin_vector,
  key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
)
destinationByDestinationMatrix <- mp_matrix(
  origins = destination_vector,
  destinations = destination_vector,
  key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
)
originByDestinationMatrix <- mp_matrix(
  origins = origin_vector,
  destinations = destination_vector,
  key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
)
destinationByOriginMatrix <- mp_matrix(
  origins = destination_vector,
  destinations = origin_vector,
  key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
)
# Naming convention = "rowByColumn"
originByOriginMatrix <- mp_get_matrix(originByOriginMatrix)
destinationByDestinationMatrix <- mp_get_matrix(destinationByDestinationMatrix)
originByDestinationMatrix <- mp_get_matrix(originByDestinationMatrix)
destinationByOriginMatrix <- mp_get_matrix(destinationByOriginMatrix)


# Combine the matrices above into a 20x20 matrix
# p1 . . p10    d1 .  . . d10
# .
# .
# p10       p10 d1 .  . . d10
# d1        p1  d1 .  . . d10
# .
# .
# d10       p10 d1 .  . . d10

# Get the matrix size from the number of entries being used. This variable is created to change if we figure out how to 
#   calculate more than 10 entries.
matSize <- nrow(data)

# Initialize the correlation distance matrix.
correlationDistanceMatrix <- matrix(0, nrow=matSize*2, ncol=matSize*2)

# Initialize the placement matrix list.
placementMatrix <- list()

# Create a list of 4 lists (one list per distance matrix found above).
for(i in 1:4){
  placementMatrix[[i]] <- list()
}

# Add the respective distance matrices into specific quadrants.
# Naming convention = "rowByColumn"
placementMatrix[[1]][[1]] <- originByOriginMatrix
placementMatrix[[1]][[2]] <- originByDestinationMatrix
placementMatrix[[2]][[1]] <- destinationByOriginMatrix
placementMatrix[[2]][[2]] <- destinationByDestinationMatrix

# Create the master distance matrix with the allocated distance matrices above
for(i in 1:2){
  for(j in 1:2){
    correlationDistanceMatrix[(i-1)*matSize+1:matSize, (j-1)*matSize+1:matSize] <- placementMatrix[[i]][[j]]
  }
}

# Add person and destination labels on the rows and columns for easier viewing
# Refer to data file for actual people/destinations
colnames(correlationDistanceMatrix) <- c(paste0('p', rep(1:matSize), sep=''), paste0('d', rep(1:matSize), sep=''))
rownames(correlationDistanceMatrix) <- c(paste0('p', rep(1:matSize), sep=''), paste0('d', rep(1:matSize), sep=''))

#write.csv(correlationDistanceMatrix, '~Documents/GitHub/DataScienceCapstone/Data/pdist.csv')


###############################
#### BUILDING PTIME MATRIX ####
###############################

# Building ptime matrix
origByOrigTimeMat <- matrix(0, nrow=matSize, ncol=matSize, 
                            dimnames=list(origin_vector, origin_vector))
destByDestTimeMat <- matrix(0, nrow=matSize, ncol=matSize, 
                            dimnames=list(origin_vector, origin_vector))
origByDestTimeMat <- matrix(0, nrow=matSize, ncol=matSize, 
                            dimnames=list(origin_vector, origin_vector))
destByOrigTimeMat <- matrix(0, nrow=matSize, ncol=matSize, 
                            dimnames=list(origin_vector, origin_vector))

for(i in 1:length(origin_vector)){
  for(j in 1:length(origin_vector)){
    originByOriginTime <- mp_directions(
      origin = origin_vector[i],
      destination = origin_vector[j],
      alternatives = FALSE,
      key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
    )
    originByOriginTime <- mp_get_routes(originByOriginTime)
    origByOrigTimeMat[i,j] <- originByOriginTime$duration_s
    
    destinationByDestinationTime <- mp_directions(
      origin = destination_vector[i],
      destination = destination_vector[j],
      alternatives = FALSE,
      key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
    )
    destinationByDestinationTime <- mp_get_routes(destinationByDestinationTime)
    destByDestTimeMat[i,j] <- destinationByDestinationTime$duration_s
    
    originByDestinationTime <- mp_directions(
      origin = origin_vector[i],
      destination = destination_vector[j],
      alternatives = FALSE,
      key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
    )
    originByDestinationTime <- mp_get_routes(originByDestinationTime)
    origByDestTimeMat[i,j] <- originByDestinationTime$duration_s
    
    destinationByOriginTime <- mp_directions(
      origin = destination_vector[i],
      destination = origin_vector[j],
      alternatives = FALSE,
      key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
    )
    destinationByOriginTime <- mp_get_routes(destinationByOriginTime)
    destByOrigTimeMat[i,j] <- destinationByOriginTime$duration_s
  }
}

correlationTimeMatrix <- matrix(0, nrow=matSize*2, ncol=matSize*2)

# Initialize the placement matrix list.
placementTimeMatrix <- list()

# Create a list of 4 lists (one list per distance matrix found above).
for(i in 1:4){
  placementTimeMatrix[[i]] <- list()
}

# Add the respective distance matrices into specific quadrants.
# Naming convention = "rowByColumn"
placementTimeMatrix[[1]][[1]] <- origByOrigTimeMat
placementTimeMatrix[[1]][[2]] <- origByDestTimeMat
placementTimeMatrix[[2]][[1]] <- destByOrigTimeMat
placementTimeMatrix[[2]][[2]] <- destByDestTimeMat

# Create the master distance matrix with the allocated distance matrices above
for(i in 1:2){
  for(j in 1:2){
    correlationTimeMatrix[(i-1)*matSize+1:matSize, (j-1)*matSize+1:matSize] <- placementTimeMatrix[[i]][[j]]
  }
}

# Add person and destination labels on the rows and columns for easier viewing
# Refer to data file for actual people/destinations
colnames(correlationTimeMatrix) <- c(paste0('p', rep(1:matSize), sep=''), paste0('d', rep(1:matSize), sep=''))
rownames(correlationTimeMatrix) <- c(paste0('p', rep(1:matSize), sep=''), paste0('d', rep(1:matSize), sep=''))

#write.csv(correlationTimeMatrix, '~/Documents/GitHub/DataScienceCapstone/Data/ptime.csv')







###################################
###### PERFORM CO CLUSTERING ######
###################################

# Using the above created correlation distance matrix...

out <- coclusterContinuous(correlationDistanceMatrix, nbcocluster=c(2,2))


clusters <- hclust(dist(correlationDistanceMatrix, method='euclidean'), method='average')
plot(clusters)

clustersWCars <- cutree(clusters, k=3)
plot(clustersWCars)

rect.hclust(clusters, k=3, border=2:6)
abline(h=3, col='red')


