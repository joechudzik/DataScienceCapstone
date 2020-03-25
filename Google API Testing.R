library(mapsapi)
library(leaflet)

rm(list=ls())

data <- read.csv('~/Documents/GitHub/DataScienceCapstone/Data/simulatedData.csv')
data <- data[1:10,]

# Vectors for origins and destinations.
origin_vector <- paste0(data$HSE_NBR_home,' ', data$STREET_home, ' ', data$STTYPE_home, ', Milwaukee, WI ', data$ZIP_CODE_home)
destination_vector <- paste0(data$HSE_NBR,' ', data$STREET, ' ', data$STTYPE, ', Milwaukee, WI ', data$ZIP_CODE)

doc=mp_directions(
  origin = origin_vector[1],
  destination = destination_vector[1],
  alternatives=FALSE,
  key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
)

routes <- mp_get_routes(doc)



for(x in 2:length(origin_vector)){
  doc <- mp_directions(
    origin = origin_vector[x],
    destination = destination_vector[x],
    alternatives = FALSE,
    key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
  )
  route <- mp_get_routes(doc)
  
  routes <- rbind(route, routes)
}

pal <- colorFactor(palette='Dark2', domain=rownames(routes))
leaflet() %>%
  addProviderTiles('CartoDB.Positron') %>%
  addPolylines(data=routes[1], opacity=1, weight=7, color=~pal(rownames(routes)))


# Build destination matrix

doc2 = mp_matrix(
  origins = origin_vector,
  destinations = destination_vector,
  key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
)

dest_matrix = mp_get_matrix(doc2, value = "distance_m") #Distance is meters
colnames(dest_matrix) = paste0('d', rep(1:10))
rownames(dest_matrix) = paste0('h', rep(1:10))

dest_matrix

doc3=mp_directions(
  origin="11500 E Lave Avenue, Englewood, CO 80111",
  destination="9659 E Mississippi Ave Aurora, CO  80247",
  alternatives=TRUE,
  key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
)

w=mp_get_routes(doc3)
w

pal = colorFactor(palette = "Dark2", domain = r$alternative_id)
leaflet() %>% 
  addProviderTiles("CartoDB.DarkMatter") %>%
  addPolylines(data = w, opacity = 1, weight = 7, color = ~pal(alternative_id))
