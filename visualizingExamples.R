
# Clean the global library

# Geocode API key = AIzaSyAVPhF7x_vfjIxHAlMry0k6M5tgC3ZjYeI
# Maps API key = AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg

source('~/Documents/GitHub/DataScienceCapstone/hclust ideas/rclust.R')
library(mapsapi)
library(ggmap)
library(ggrepel)

data <- read.csv('~/Documents/GitHub/DataScienceCapstone/Data/simulatedData.csv')
data <- data[1:10,]

pdis <- read.csv('https://raw.githubusercontent.com/joechudzik/DataScienceCapstone/master/Data/pdist.csv')
row.names(pdis) <- pdis[,1]; pdis <- pdis[,-1]
ptim <- read.csv('https://raw.githubusercontent.com/joechudzik/DataScienceCapstone/master/Data/ptime.csv')
row.names(ptim) <- ptim$X; ptim <- ptim[,-1]


# Plotting pickups and dropoffs on a map

origin_vector <- paste0(data$HSE_NBR_home,' ', data$STREET_home, ' ', data$STTYPE_home, ', Milwaukee, WI ', data$ZIP_CODE_home)
destination_vector <- paste0(data$HSE_NBR,' ', data$STREET, ' ', data$STTYPE, ', Milwaukee, WI ', data$ZIP_CODE)

register_google(key='AIzaSyAVPhF7x_vfjIxHAlMry0k6M5tgC3ZjYeI')
origin_locs <- geocode(origin_vector); origin_locs$type <- 'Pickup'; origin_locs$index <- paste0('p', rep(1:10), sep='')
destination_locs <- geocode(destination_vector); destination_locs$type <- 'Dropoff'; destination_locs$index <- c('d1', 'd1', 'd3', 'd4', 'd5', 'd6', 'd7', 'd8', 'd5', 'd9')
locs <- rbind(origin_locs, destination_locs)

mke <- c(lon = -87.9065, lat = 43.0389)
ggmap(get_googlemap(center = mke, zoom=11)) + 
  geom_point(aes(lon, lat, color=type), data=locs) +
  labs(x='Longitude',y='Latitude') + 
  geom_label_repel(data=locs, aes(x=lon, y=lat, label=index), size=3, vjust=1.25, hjust=-1)



# Newly designed hclust algorithm

# low time importance
lowTimeImp <- rhclust(pdis, ptim, 0)
cutree(lowTimeImp,3)
lowTimeImp$merge.route[5:7]
lowTimeImp$height[5:7]
plot(lowTimeImp, main='Low time importance')
rect.hclust(lowTimeImp, k=3, border=2:6)

# medium time importance
medTimeImp <- rhclust(pdis, ptim, 0.5)
cutree(medTimeImp,3)
plot(medTimeImp, main='Medium time importance')
rect.hclust(medTimeImp, k=3, border=2:6)

# high time importance
highTimeImp <- rhclust(pdis, ptim, 1)
cutree(highTimeImp,3)
plot(highTimeImp, main='High time importance')
rect.hclust(highTimeImp, k=3, border=2:6)
