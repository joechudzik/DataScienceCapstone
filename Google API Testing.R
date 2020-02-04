install.packages("mapsapi")
install.packages("leaflet")
library(mapsapi)
library(leaflet)

doc=mp_directions(
  origin="1313 W Wisconsin Ave, Milwaukee, WI 53233",
  destination="1500 W Wells Street, Milwaukee, WI 53233",
  alternatives=TRUE,
  key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
)

r=mp_get_routes(doc)

r

pal = colorFactor(palette = "Dark2", domain = r$alternative_id)
leaflet() %>% 
  addProviderTiles("CartoDB.DarkMatter") %>%
  addPolylines(data = r, opacity = 1, weight = 7, color = ~pal(alternative_id))

locations = c("1313 W Wisconsin Ave, Milwaukee, WI 53233", "1500 W Wells Street, Milwaukee, WI 53233", "1111 Vel R Phillips Avenue Milwaukee, WI 53203")

doc2 = mp_matrix(
  origins = locations,
  destinations = locations,
  key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
)

m = mp_get_matrix(doc2, value = "distance_m") #Distance is meters
colnames(m) = locations
rownames(m) = locations

m

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
