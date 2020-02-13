library(mapsapi)
library(leaflet)
library(Matrix)

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

originByOriginMatrix <- mp_get_matrix(originByOriginMatrix)
destinationByDestinationMatrix <- mp_get_matrix(destinationByDestinationMatrix)
originByDestinationMatrix <- mp_get_matrix(originByDestinationMatrix)
destinationByOriginMatrix <- mp_get_matrix(destinationByOriginMatrix)


# Combine the matrices above into a 100x100 matrix then perform some type of clustering
#    p1 . . p10 d1 .  . . d10
# p1
# .
# .
# p10       p10 d1
# d1        p1  d1 .  . . d10
# .
# .
# d10       p10 d1 .  . . d10


matrix <- matrix(0, nrow=4, ncol=4)
test1 <- matrix(1, nrow=2, ncol=2)
test2 <- matrix(2, nrow=2, ncol=2)
test3 <- matrix(3, nrow=2, ncol=2)
test4 <- matrix(4, nrow=2, ncol=2)

testlist <- list(test1, test2, test3, test4)

# How do i put this into a loop?
  matrix[1,1] <- testlist[[1]][1,1]
  matrix[1,2] <- testlist[[1]][1,2]
  matrix[2,1] <- testlist[[1]][2,1]
  matrix[2,2] <- testlist[[1]][2,2]
  
  matrix[1,3] <- testlist[[2]][1,1]
  matrix[1,4] <- testlist[[2]][1,2]
  matrix[2,3] <- testlist[[2]][2,1]
  matrix[2,4] <- testlist[[2]][2,2]
  
  matrix[3,1] <- testlist[[3]][1,1]
  matrix[3,2] <- testlist[[3]][1,2]
  matrix[3,3] <- testlist[[4]][1,1]
  matrix[3,4] <- testlist[[4]][1,2]
  
  matrix[4,1] <- testlist[[3]][2,1]
  matrix[4,2] <- testlist[[3]][2,2]
  matrix[4,3] <- testlist[[4]][2,1]
  matrix[4,4] <- testlist[[4]][2,2]
  
  
for(i in 1:nrow(matrix)){
  for(j in 1:ncol(matrix)){
    matrix[i,j] <- testlist[[j]][1,1]
  }
}

