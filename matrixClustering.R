library(mapsapi)
library(leaflet)
library(Matrix)
library(blockcluster)
library(NbClust)

# Clear the global library.
rm(list=ls())

data <- read.csv('~/Documents/GitHub/DataScienceCapstone/Data/simulatedData.csv')

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



###################################
###### PERFORM CO CLUSTERING ######
###################################

# Using the above created correlation distance matrix...

out <- coclusterContinuous(correlationDistanceMatrix, nbcocluster=c(2,2))


clusters <- hclust(dist(correlationDistanceMatrix, method='euclidean'), method='average')
plot(clusters)

clustersWCars <- cutree(clusters, k=3)
plot(clustersWCars)

rect.hclust()
