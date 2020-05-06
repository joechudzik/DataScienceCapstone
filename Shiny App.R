library(mapsapi)
library(leaflet)
library(shiny)
library(dendextend)
library(ggmap)
library(ggrepel)
library(viridis)
library(data.table)
library(rsconnect)

cost <- function(path, pdist=pdist, ptime=ptime, weight, data) {
  cost <- NULL; 
  # shiny app will have a slider that controls lambda (how important ptime will be)
  # for (i in 1:(length(path)-1)) cost[i] <- pdist[path[i],path[i+1]]
  #for (i in 1:(length(path)-1)) cost[i] <- (1-weight)*pdist[path[i],path[i+1]] + weight*ptime[path[i],path[i+1]]
  
  # Build dtime vector
  dtime <- rep(0, length(path))
  for(i in 1:(length(path)-1)) dtime[i+1] <- ptime[path[i],path[i+1]]
  
  # Build routeTime table
  routeTime <- data.table(path=path, et=c(0), lt=c(0))
  
  for (i in 1:length(path)){
    # Retrieve the location type (pickup or dropoff) and the corresponding index
    locationType <- substr(path[i], 1, 1)
    locationIndex <- as.integer(substr(path[i], 2, length(path[i])+1))
    
    # If location type is a pickup, allocate corresponding ept and lpt to the correct row in routeTime
    #   Otherwise if an error occurred, just fill with NA
    switch(locationType,
           'p' = {routeTime$et[i] <- data$ept[locationIndex]; routeTime$lt[i] <- data$lpt[locationIndex]},
           'd' = {routeTime$et[i] <- data$edt[locationIndex]; routeTime$lt[i] <- data$ldt[locationIndex]},
           {routeTimes$et[i] <- NA; routeTimes$lt[i] <- NA})
  }
  
  # Find the cumulative distance
  est <- routeTime$et - cumsum(dtime)
  lst <- routeTime$lt - cumsum(dtime)
  
  # Logic of function (testing if path works)
  # if( max(est, na.rm=T)>min(lst, na.rm=T) && weight!=0) {
  cumin.lst <- NA; 
  cumin.lst[!is.na(lst)] <- rev(cummin(na.omit(rev(lst))))
  if (sum(!(est < cumin.lst),na.rm=T)  && weight!=0){
    return(list(cost=Inf))
  }
  else{
    st <- min(lst, max(est), na.rm=T)
    trip_time <- st + dtime
    wait_time <- routeTime$et - trip_time
    for(i in 1:length(wait_time)){ if( is.na(wait_time[i]) ) wait_time[i] <- 0}
    trip_time <- trip_time + wait_time
    new_time <- wait_time + dtime
    for (i in 1:(length(path)-1)) cost[i] <- (1-weight)*pdist[path[i],path[i+1]] + weight*trip_time[i]
    return(list(cost=sum(cost),st=st,wait_time=wait_time,trip_time=trip_time))
  }
}

# Searching among possible routes by variety of combination of "a" and "b" and pick the optimal one.
rcost <- function(a,b,pdist=pdist,ptime=ptime,weight,data) {
  rcosts <- list()
  if (is.numeric(a) && is.numeric(b)) {
    routs <- list(c(paste0("p",a),paste0("p",b),paste0("d",a),paste0("d",b)),
                  c(paste0("p",a),paste0("p",b),paste0("d",b),paste0("d",a)),
                  c(paste0("p",a),paste0("d",a),paste0("p",b),paste0("d",b)),
                  c(paste0("p",b),paste0("d",b),paste0("p",a),paste0("d",a)),
                  c(paste0("p",b),paste0("p",a),paste0("d",a),paste0("d",b)),
                  c(paste0("p",b),paste0("p",a),paste0("d",b),paste0("d",a)))
  } else if (is.numeric(a) || is.numeric(b)) {
    if (is.numeric(a)) {tem <- a; a <- b; b <- tem; rm("tem")}
    b <- paste0(c("p","d"),b)
    pdis <- pdist[b,a]; i <- which.min(pdis[1,]); j <- which.min(pdis[2,])
    if (i < j) {
      routs <- list(c(a[0:(i-1)],b[1],a[i:(j-1)],b[2],a[j:length(a)]),
                    c(a[0:(i-1)],b[1],a[i:j],b[2],a[-(1:j)]),
                    c(a[1:i],b[1],a[intersect((i+1):length(a),1:(j-1))],b[2],a[j:length(a)]),
                    c(a[1:i],b[1],a[(i+1):j],b[2],a[-(1:j)]))
    } else {
      routs <- unique(list(c(a[0:(i-1)],b,a[i:length(a)]), c(a[1:i],b,a[-(1:i)]),
                           c(a[0:(j-1)],b,a[j:length(a)]), c(a[1:j],b,a[-(1:j)])))
    }
  } else {
    routs <- list(c(a,b),c(b,a))
  }
  for (i in 1:length(routs)) rcosts[[i]] <- cost(routs[[i]], pdist, ptime, weight, data)
  # May need to add some logic to remove the very high or NA value from rcosts
  indx <- which.min(unlist(rcosts)[names(unlist(rcosts))=="cost"])
  return(list(rcost=rcosts[[indx]], route=routs[[indx]]))
}

# Calculate all of the pairwise costs needed to start the algorithm
pcosts <- function(pdist=pdist,ptime=ptime,weight, data) {
  N <- dim(pdist)[1]/2
  pcosts <- list(rcosts=matrix(NA,nr=N,nc=N), routes=array(list(),dim=c(N,N)))
  for (i in 1:(N-1)) for (j in (i+1):N) {
    pcost.ij <- rcost(i,j,pdist,ptime,weight,data)
    pcosts$rcosts[i,j] <- pcosts$rcosts[j,i] <- pcost.ij$rcost$cost
    pcosts$routes[i,j] <- pcosts$routes[j,i] <- list(pcost.ij$route)
  }
  colnames(pcosts$rcosts) <- rownames(pcosts$rcosts) <- colnames(pcosts$routes) <- rownames(pcosts$routes) <- paste0("-",1:N)
  return(pcosts)
}

# Route hclust function
rhclust <- function(pdist,ptime,weight,data) {
  temp <- list(pcosts(pdist,ptime,weight,data)); N <- dim(pdist)[1]/2;
  rhclust <- list(merge=matrix(0,nr=N-1,nc=2), merge.route=list())
  for (i in 1:(N-2)) {
    indx <- rownames(which(temp[[i]]$rcosts==min(temp[[i]]$rcosts, na.rm = T), arr.ind = T))[1:2]
    rhclust$merge[i,] <- as.numeric(indx); rhclust$merge.route[i] <- temp[[i]]$routes[indx[1],indx[2]]; 
    rhclust$height[i] <- temp[[i]]$rcosts[indx[1],indx[2]]
    ind <- !colnames(temp[[i]]$rcosts)%in%indx
    temp[[(i+1)]] <- list(rcosts=as.matrix(temp[[i]]$rcosts[ind,ind]), routes=as.matrix(temp[[i]]$routes[ind,ind]))
    inds <- as.list(colnames(temp[[i]]$rcosts)[ind])
    if (sum(ind)==1) colnames(temp[[(i+1)]][[1]]) <- rownames(temp[[(i+1)]][[1]]) <- colnames(temp[[(i+1)]][[2]]) <- rownames(temp[[(i+1)]][[2]]) <- inds[[1]]
    cost.j <- list(rcosts=array(NA,length(inds)+1), routes=array(list(),length(inds)+1))
    for (j in 1:length(inds)) {
      if (sign(as.numeric(inds[[j]]))==-1) inds[[j]] <- -as.numeric(inds[[j]]) else inds[[j]] <- rhclust$merge.route[[as.numeric(inds[[j]])]]
      tmp <- rcost(inds[[j]], rhclust$merge.route[[i]], pdist, ptime, weight, data)
      cost.j$rcosts[j] <- tmp$rcost$cost; cost.j$routes[[j]] <- tmp$route
    }
    temp[[(i+1)]]$rcosts <- cbind(rbind(temp[[(i+1)]]$rcosts,cost.j$rcosts[-length(cost.j$rcosts)]),cost.j$rcosts)
    temp[[(i+1)]]$routes <- cbind(rbind(temp[[(i+1)]]$routes,cost.j$routes[-length(cost.j$routes)]),cost.j$routes)
    tmp <- colnames(temp[[(i+1)]]$rcosts); tmp[length(tmp)] <- i
    colnames(temp[[(i+1)]]$rcosts) <- rownames(temp[[(i+1)]]$rcosts) <- colnames(temp[[(i+1)]]$routes) <- rownames(temp[[(i+1)]]$routes) <- tmp
  }
  rhclust$merge[(i+1),] <- as.numeric(tmp); 
  rhclust$merge.route[(i+1)] <- temp[[(i+1)]]$routes[1,2]; 
  rhclust$height[(i+1)] <- temp[[(i+1)]]$rcosts[1,2]; rhclust$height[rhclust$height==Inf] <- max(rhclust$height[rhclust$height<Inf])
  rhclust$labels <- NULL; rhclust$method <- "pool"; rhclust$dist.method <- "pdist"; rhclust$call <- sys.call()
  rhclust$order <- 1:N; class(rhclust) <- "hclust"
  rhclust$order <- unlist(as.dendrogram(rhclust))
  return(rhclust)
}
data <- read.csv('https://raw.githubusercontent.com/joechudzik/DataScienceCapstone/master/Data/dataForSpecificTimeAlgorithm.csv',stringsAsFactors=FALSE)
data <- data[1:10,]
pdis <- read.csv('https://raw.githubusercontent.com/joechudzik/DataScienceCapstone/master/Data/pdist.csv')
row.names(pdis) <- pdis[,1]; pdis <- pdis[,-1]
ptim <- read.csv('https://raw.githubusercontent.com/joechudzik/DataScienceCapstone/master/Data/ptime.csv')
row.names(ptim) <- ptim$X; ptim <- ptim[,-1]

# Vectors for origins and destinations.
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

r_colors <- rgb(t(col2rgb(colors()) / 255))
names(r_colors) <- colors()



ui <- fluidPage(
  leafletOutput("routemap"),
  p(),
  titlePanel("Actions"),
  
  sidebarLayout(
    fluidRow(
      column(5,sliderInput("integer", "Time Factor:",
                           min = 0, max = 1,
                           value = 0,step=.5),
             
             tableOutput("values")),
      
      column(5,sliderInput("cluster", "Clusters:",
                           min = 2, max = 5,
                           value = 3)),
      
      column(5,radioButtons("radio", h3("Routes (As many as clusters)"),
                            choices = list("1" = 1, "2" = 2,
                                           "3" = 3,"4" = 4,
                                           "5" = 5),selected = 1)),
      
      #column(3,plotOutput("plot")),
      
    ),
    # Main panel for displaying outputs ----
    mainPanel(
      plotOutput("plot"),
      leafletOutput("mymap")
    )
  )
)
server <- function(input, output, session) {
  
  points <- eventReactive(input$recalc, {
    cbind(rnorm(40) * 2 + 13, rnorm(40) + 48)
  }, ignoreNULL = FALSE)
  

  output$mymap <- renderLeaflet({
      leaflet() %>%
        addProviderTiles("CartoDB.Positron") %>%
        addCircleMarkers(data = origin_locs, label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs, label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data = routes[1], opacity = 1, weight = 7, color = 'purple')
    })
  
  output$plot <- renderPlot({
    input$newplot
    dend <- rhclust(pdis, ptim, input$integer, data)
    cutree(dend,input$cluster)
    plot(dend, main='Time Dendogram')
    rect.hclust(dend, k=input$cluster, border=2:6)
  })
  
  output$routemap <- renderLeaflet({
    if(input$radio==1 & input$integer==0 & input$cluster==2){
      temp <- rhclust(pdis, ptim, 0, data)
      call=substr(temp$merge.route[[6]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[6]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR,' ', data[row[2],col[2]:(col[2]+4)]$STREET, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE)
      ori[2] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR_home,' ', data[row[3],col[3]:(col[3]+4)]$STREET_home, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE_home)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = dest[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = ori[2],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      leaflet() %>%
        addProviderTiles("CartoDB.Positron") %>%
        addCircleMarkers(data = origin_locs[4,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = origin_locs[8,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[4,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = destination_locs[8,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') }
    else if(input$radio==2 & input$integer==0 & input$cluster==2){
      temp <- rhclust(pdis, ptim, 0, data)
      call=substr(temp$merge.route[[8]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[8]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      ori[2] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR_home,' ', data[row[2],col[2]:(col[2]+4)]$STREET_home, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR,' ', data[row[3],col[3]:(col[3]+4)]$STREET, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      ori[3] <- paste0(data[row[5],col[5]:(col[5]+4)]$HSE_NBR_home,' ', data[row[5],col[5]:(col[5]+4)]$STREET_home, ' ', data[row[5],col[5]:(col[5]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[5],col[5]:(col[5]+4)]$ZIP_CODE_home)
      ori[4] <- paste0(data[row[6],col[6]:(col[6]+4)]$HSE_NBR_home,' ', data[row[6],col[6]:(col[6]+4)]$STREET_home, ' ', data[row[6],col[6]:(col[6]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[6],col[6]:(col[6]+4)]$ZIP_CODE_home)
      dest[3] <- paste0(data[row[7],col[7]:(col[7]+4)]$HSE_NBR,' ', data[row[7],col[7]:(col[7]+4)]$STREET, ' ', data[row[7],col[7]:(col[7]+4)]$STTYPE, ', Milwaukee, WI ', data[row[7],col[7]:(col[7]+4)]$ZIP_CODE)
      dest[4] <- paste0(data[row[8],col[8]:(col[8]+4)]$HSE_NBR,' ', data[row[8],col[8]:(col[8]+4)]$STREET, ' ', data[row[8],col[8]:(col[8]+4)]$STTYPE, ', Milwaukee, WI ', data[row[8],col[8]:(col[8]+4)]$ZIP_CODE)
      ori[5] <- paste0(data[row[9],col[9]:(col[9]+4)]$HSE_NBR_home,' ', data[row[9],col[9]:(col[9]+4)]$STREET_home, ' ', data[row[9],col[9]:(col[9]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[9],col[9]:(col[9]+4)]$ZIP_CODE_home)
      dest[5] <- paste0(data[row[10],col[10]:(col[10]+4)]$HSE_NBR,' ', data[row[10],col[10]:(col[10]+4)]$STREET, ' ', data[row[10],col[10]:(col[10]+4)]$STTYPE, ', Milwaukee, WI ', data[row[10],col[10]:(col[10]+4)]$ZIP_CODE)
      ori[6] <- paste0(data[row[11],col[11]:(col[11]+4)]$HSE_NBR_home,' ', data[row[11],col[11]:(col[11]+4)]$STREET_home, ' ', data[row[11],col[11]:(col[11]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[11],col[11]:(col[11]+4)]$ZIP_CODE_home)
      ori[7] <- paste0(data[row[12],col[12]:(col[12]+4)]$HSE_NBR_home,' ', data[row[12],col[12]:(col[12]+4)]$STREET_home, ' ', data[row[12],col[12]:(col[12]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[12],col[12]:(col[12]+4)]$ZIP_CODE_home)
      ori[8] <- paste0(data[row[13],col[13]:(col[13]+4)]$HSE_NBR_home,' ', data[row[13],col[13]:(col[13]+4)]$STREET_home, ' ', data[row[13],col[13]:(col[13]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[13],col[13]:(col[13]+4)]$ZIP_CODE_home)
      dest[6] <- paste0(data[row[14],col[14]:(col[14]+4)]$HSE_NBR,' ', data[row[14],col[14]:(col[14]+4)]$STREET, ' ', data[row[14],col[14]:(col[14]+4)]$STTYPE, ', Milwaukee, WI ', data[row[14],col[14]:(col[14]+4)]$ZIP_CODE)
      dest[7] <- paste0(data[row[15],col[15]:(col[15]+4)]$HSE_NBR,' ', data[row[15],col[15]:(col[15]+4)]$STREET, ' ', data[row[15],col[15]:(col[15]+4)]$STTYPE, ', Milwaukee, WI ', data[row[15],col[15]:(col[15]+4)]$ZIP_CODE)
      dest[8] <- paste0(data[row[16],col[16]:(col[16]+4)]$HSE_NBR,' ', data[row[16],col[16]:(col[16]+4)]$STREET, ' ', data[row[16],col[16]:(col[16]+4)]$STTYPE, ', Milwaukee, WI ', data[row[16],col[16]:(col[16]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = ori[2],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = dest[1],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      doc4 <- mp_directions(
        origin = dest[2],
        destination = ori[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route4 <- mp_get_routes(doc4)
      doc5 <- mp_directions(
        origin = ori[3],
        destination = ori[4],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route5 <- mp_get_routes(doc5)
      doc6 <- mp_directions(
        origin = ori[4],
        destination = dest[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route6 <- mp_get_routes(doc6)
      doc7 <- mp_directions(
        origin = dest[3],
        destination = dest[4],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route7 <- mp_get_routes(doc7)
      doc8 <- mp_directions(
        origin = dest[4],
        destination = ori[5],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route8 <- mp_get_routes(doc8)
      doc9 <- mp_directions(
        origin = ori[5],
        destination = dest[5],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route9 <- mp_get_routes(doc9)
      doc10 <- mp_directions(
        origin = dest[5],
        destination = ori[6],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route10 <- mp_get_routes(doc10)
      doc11 <- mp_directions(
        origin = ori[6],
        destination = ori[7],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route11 <- mp_get_routes(doc11)
      doc12 <- mp_directions(
        origin = ori[7],
        destination = ori[8],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route12 <- mp_get_routes(doc12)
      doc13 <- mp_directions(
        origin = ori[8],
        destination = dest[6],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route13 <- mp_get_routes(doc13)
      doc14 <- mp_directions(
        origin = dest[6],
        destination = dest[7],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route14 <- mp_get_routes(doc14)
      doc15 <- mp_directions(
        origin = dest[7],
        destination = dest[8],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route15 <- mp_get_routes(doc15)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[1,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[1,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[2,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[2,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[3,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[3,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[5,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[5,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[6,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[6,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[7,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[7,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[9,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[9,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[10,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[10,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') %>%
        addPolylines(data=route4, opacity=1, weight=7, color='green') %>%
        addPolylines(data=route5, opacity=1, weight=7, color='brown') %>%
        addPolylines(data=route6, opacity=1, weight=7, color='purple') %>%
        addPolylines(data=route7, opacity=1, weight=7, color='orange') %>%
        addPolylines(data=route8, opacity=1, weight=7, color='black') %>%
        addPolylines(data=route9, opacity=1, weight=7, color='white') %>%
        addPolylines(data=route10, opacity=1, weight=7, color='grey') %>%
        addPolylines(data=route11, opacity=1, weight=7, color='blue') %>%
        addPolylines(data=route12, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route13, opacity=1, weight=7, color='yellow') %>%
        addPolylines(data=route14, opacity=1, weight=7, color='green') %>%
        addPolylines(data=route15, opacity=1, weight=7, color='brown') 
    }
    else if(input$radio==1 & input$integer==0 & input$cluster==3){
      temp <- rhclust(pdis, ptim, 0, data)
      call=substr(temp$merge.route[[6]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[6]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR,' ', data[row[2],col[2]:(col[2]+4)]$STREET, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE)
      ori[2] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR_home,' ', data[row[3],col[3]:(col[3]+4)]$STREET_home, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE_home)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = dest[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = ori[2],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      leaflet() %>%
        addProviderTiles("CartoDB.Positron") %>%
        addCircleMarkers(data = origin_locs[4,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = origin_locs[8,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[4,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = destination_locs[8,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') }
    else if(input$radio==2 & input$integer==0 & input$cluster==3){
      temp <- rhclust(pdis, ptim, 0, data)
      call=substr(temp$merge.route[[5]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[5]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      ori[2] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR_home,' ', data[row[2],col[2]:(col[2]+4)]$STREET_home, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE_home)
      ori[3] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR_home,' ', data[row[3],col[3]:(col[3]+4)]$STREET_home, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      dest[2] <- paste0(data[row[5],col[5]:(col[5]+4)]$HSE_NBR,' ', data[row[5],col[5]:(col[5]+4)]$STREET, ' ', data[row[5],col[5]:(col[5]+4)]$STTYPE, ', Milwaukee, WI ', data[row[5],col[5]:(col[5]+4)]$ZIP_CODE)
      dest[3] <- paste0(data[row[6],col[6]:(col[6]+4)]$HSE_NBR,' ', data[row[6],col[6]:(col[6]+4)]$STREET, ' ', data[row[6],col[6]:(col[6]+4)]$STTYPE, ', Milwaukee, WI ', data[row[6],col[6]:(col[6]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = ori[2],
        destination = ori[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = ori[3],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      doc4 <- mp_directions(
        origin = dest[1],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route4 <- mp_get_routes(doc4)
      doc5 <- mp_directions(
        origin = dest[2],
        destination = dest[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route5 <- mp_get_routes(doc5)
      leaflet() %>%
        addProviderTiles("CartoDB.Positron") %>%
        addCircleMarkers(data = origin_locs[2,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = origin_locs[1,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = origin_locs[10,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[10,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = destination_locs[1,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = destination_locs[2,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') %>%
        addPolylines(data=route4, opacity=1, weight=7, color='green')%>%
        addPolylines(data=route5, opacity=1, weight=7, color='brown')}
    else if(input$radio==3 & input$integer==0 & input$cluster==3){
      temp <- rhclust(pdis, ptim, 0, data)
      call=substr(temp$merge.route[[7]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[7]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      ori[2] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR_home,' ', data[row[2],col[2]:(col[2]+4)]$STREET_home, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR,' ', data[row[3],col[3]:(col[3]+4)]$STREET, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      ori[3] <- paste0(data[row[5],col[5]:(col[5]+4)]$HSE_NBR_home,' ', data[row[5],col[5]:(col[5]+4)]$STREET_home, ' ', data[row[5],col[5]:(col[5]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[5],col[5]:(col[5]+4)]$ZIP_CODE_home)
      ori[4] <- paste0(data[row[6],col[6]:(col[6]+4)]$HSE_NBR_home,' ', data[row[6],col[6]:(col[6]+4)]$STREET_home, ' ', data[row[6],col[6]:(col[6]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[6],col[6]:(col[6]+4)]$ZIP_CODE_home)
      dest[3] <- paste0(data[row[7],col[7]:(col[7]+4)]$HSE_NBR,' ', data[row[7],col[7]:(col[7]+4)]$STREET, ' ', data[row[7],col[7]:(col[7]+4)]$STTYPE, ', Milwaukee, WI ', data[row[7],col[7]:(col[7]+4)]$ZIP_CODE)
      dest[4] <- paste0(data[row[8],col[8]:(col[8]+4)]$HSE_NBR,' ', data[row[8],col[8]:(col[8]+4)]$STREET, ' ', data[row[8],col[8]:(col[8]+4)]$STTYPE, ', Milwaukee, WI ', data[row[8],col[8]:(col[8]+4)]$ZIP_CODE)
      ori[5] <- paste0(data[row[9],col[9]:(col[9]+4)]$HSE_NBR_home,' ', data[row[9],col[9]:(col[9]+4)]$STREET_home, ' ', data[row[9],col[9]:(col[9]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[9],col[9]:(col[9]+4)]$ZIP_CODE_home)
      dest[5] <- paste0(data[row[10],col[10]:(col[10]+4)]$HSE_NBR,' ', data[row[10],col[10]:(col[10]+4)]$STREET, ' ', data[row[10],col[10]:(col[10]+4)]$STTYPE, ', Milwaukee, WI ', data[row[10],col[10]:(col[10]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = ori[2],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = dest[1],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      doc4 <- mp_directions(
        origin = dest[2],
        destination = ori[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route4 <- mp_get_routes(doc4)
      doc5 <- mp_directions(
        origin = ori[3],
        destination = ori[4],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route5 <- mp_get_routes(doc5)
      doc6 <- mp_directions(
        origin = ori[4],
        destination = dest[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route6 <- mp_get_routes(doc6)
      doc7 <- mp_directions(
        origin = dest[3],
        destination = dest[4],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route7 <- mp_get_routes(doc7)
      doc8 <- mp_directions(
        origin = dest[4],
        destination = ori[5],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route8 <- mp_get_routes(doc8)
      doc9 <- mp_directions(
        origin = ori[5],
        destination = dest[5],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route9 <- mp_get_routes(doc9)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[3,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[3,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[5,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[5,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[6,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[6,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[7,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[7,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') %>%
        addPolylines(data=route4, opacity=1, weight=7, color='green') %>%
        addPolylines(data=route5, opacity=1, weight=7, color='brown') %>%
        addPolylines(data=route6, opacity=1, weight=7, color='purple') %>%
        addPolylines(data=route7, opacity=1, weight=7, color='orange') %>%
        addPolylines(data=route8, opacity=1, weight=7, color='black') %>%
        addPolylines(data=route9, opacity=1, weight=7, color='white') 
      ###
    }
    else if(input$radio==1 & input$integer==0 & input$cluster==4){
      temp <- rhclust(pdis, ptim, 0, data)
      call=substr(temp$merge.route[[6]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[6]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR,' ', data[row[2],col[2]:(col[2]+4)]$STREET, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE)
      ori[2] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR_home,' ', data[row[3],col[3]:(col[3]+4)]$STREET_home, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE_home)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = dest[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = ori[2],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      leaflet() %>%
        addProviderTiles("CartoDB.Positron") %>%
        addCircleMarkers(data = origin_locs[4,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = origin_locs[8,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[4,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = destination_locs[8,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') }
    else if(input$radio==2 & input$integer==0 & input$cluster==4){
      temp <- rhclust(pdis, ptim, 0, data)
      call=substr(temp$merge.route[[4]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[4]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      ori[2] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR_home,' ', data[row[2],col[2]:(col[2]+4)]$STREET_home, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR,' ', data[row[3],col[3]:(col[3]+4)]$STREET, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = ori[2],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = dest[1],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      leaflet() %>%
        addProviderTiles("CartoDB.Positron") %>%
        addCircleMarkers(data = origin_locs[5,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = origin_locs[9,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[5,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = destination_locs[9,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') }
    else if(input$radio==3 & input$integer==0 & input$cluster==4){
      temp <- rhclust(pdis, ptim, 0, data)
      call=substr(temp$merge.route[[3]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[3]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      ori[2] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR_home,' ', data[row[2],col[2]:(col[2]+4)]$STREET_home, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR,' ', data[row[3],col[3]:(col[3]+4)]$STREET, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      ori[3] <- paste0(data[row[5],col[5]:(col[5]+4)]$HSE_NBR_home,' ', data[row[5],col[5]:(col[5]+4)]$STREET_home, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[5],col[5]:(col[5]+4)]$ZIP_CODE_home)
      dest[3] <- paste0(data[row[6],col[6]:(col[6]+4)]$HSE_NBR,' ', data[row[6],col[6]:(col[6]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[6],col[6]:(col[6]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = ori[2],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = dest[1],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      doc4 <- mp_directions(
        origin = dest[2],
        destination = ori[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route4 <- mp_get_routes(doc4)
      doc5 <- mp_directions(
        origin = ori[3],
        destination = dest[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route5 <- mp_get_routes(doc5)
      leaflet() %>%
        addProviderTiles("CartoDB.Positron") %>%
        addCircleMarkers(data = origin_locs[6,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = origin_locs[7,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = origin_locs[3,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[7,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = destination_locs[6,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = destination_locs[3,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') %>%
        addPolylines(data=route4, opacity=1, weight=7, color='green')%>%
        addPolylines(data=route5, opacity=1, weight=7, color='brown')}
    else if(input$radio==4 & input$integer==0 & input$cluster==4){
      temp <- rhclust(pdis, ptim, 0, data)
      call=substr(temp$merge.route[[5]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[5]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      ori[2] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR_home,' ', data[row[2],col[2]:(col[2]+4)]$STREET_home, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE_home)
      ori[3] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR_home,' ', data[row[3],col[3]:(col[3]+4)]$STREET_home, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      dest[2] <- paste0(data[row[5],col[5]:(col[5]+4)]$HSE_NBR,' ', data[row[5],col[5]:(col[5]+4)]$STREET, ' ', data[row[5],col[5]:(col[5]+4)]$STTYPE, ', Milwaukee, WI ', data[row[5],col[5]:(col[5]+4)]$ZIP_CODE)
      dest[3] <- paste0(data[row[6],col[6]:(col[6]+4)]$HSE_NBR,' ', data[row[6],col[6]:(col[6]+4)]$STREET, ' ', data[row[6],col[6]:(col[6]+4)]$STTYPE, ', Milwaukee, WI ', data[row[6],col[6]:(col[6]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = ori[2],
        destination = ori[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = ori[3],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      doc4 <- mp_directions(
        origin = dest[1],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route4 <- mp_get_routes(doc4)
      doc5 <- mp_directions(
        origin = dest[2],
        destination = dest[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route5 <- mp_get_routes(doc5)
      leaflet() %>%
        addProviderTiles("CartoDB.Positron") %>%
        addCircleMarkers(data = origin_locs[2,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = origin_locs[1,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = origin_locs[10,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[10,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = destination_locs[1,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = destination_locs[2,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') %>%
        addPolylines(data=route4, opacity=1, weight=7, color='green')%>%
        addPolylines(data=route5, opacity=1, weight=7, color='brown')}
    else if(input$radio==1 & input$integer==0 & input$cluster==5){
      temp <- rhclust(pdis, ptim, 0, data)
      call=substr(temp$merge.route[[6]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[6]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR,' ', data[row[2],col[2]:(col[2]+4)]$STREET, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE)
      ori[2] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR_home,' ', data[row[3],col[3]:(col[3]+4)]$STREET_home, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE_home)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      doc3 <- mp_directions(
        origin = ori[2],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      leaflet() %>%
        addProviderTiles("CartoDB.Positron") %>%
        addCircleMarkers(data = origin_locs[8,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[8,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route3, opacity=1, weight=7) }
    else if(input$radio==2 & input$integer==0 & input$cluster==5){
      temp <- rhclust(pdis, ptim, 0, data)
      call=substr(temp$merge.route[[6]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[6]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR,' ', data[row[2],col[2]:(col[2]+4)]$STREET, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE)
      ori[2] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR_home,' ', data[row[3],col[3]:(col[3]+4)]$STREET_home, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE_home)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      doc3 <- mp_directions(
        origin = ori[1],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      leaflet() %>%
        addProviderTiles("CartoDB.Positron") %>%
        addCircleMarkers(data = origin_locs[4,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[4,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route3, opacity=1, weight=7) }
    else if(input$radio==3 & input$integer==0 & input$cluster==5){
      temp <- rhclust(pdis, ptim, 0, data)
      call=substr(temp$merge.route[[4]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[4]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      ori[2] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR_home,' ', data[row[2],col[2]:(col[2]+4)]$STREET_home, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR,' ', data[row[3],col[3]:(col[3]+4)]$STREET, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = ori[2],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = dest[1],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      leaflet() %>%
        addProviderTiles("CartoDB.Positron") %>%
        addCircleMarkers(data = origin_locs[5,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = origin_locs[9,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[5,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = destination_locs[9,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') }
    else if(input$radio==4 & input$integer==0 & input$cluster==5){
      temp <- rhclust(pdis, ptim, 0, data)
      call=substr(temp$merge.route[[3]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[3]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      ori[2] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR_home,' ', data[row[2],col[2]:(col[2]+4)]$STREET_home, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR,' ', data[row[3],col[3]:(col[3]+4)]$STREET, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      ori[3] <- paste0(data[row[5],col[5]:(col[5]+4)]$HSE_NBR_home,' ', data[row[5],col[5]:(col[5]+4)]$STREET_home, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[5],col[5]:(col[5]+4)]$ZIP_CODE_home)
      dest[3] <- paste0(data[row[6],col[6]:(col[6]+4)]$HSE_NBR,' ', data[row[6],col[6]:(col[6]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[6],col[6]:(col[6]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = ori[2],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = dest[1],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      doc4 <- mp_directions(
        origin = dest[2],
        destination = ori[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route4 <- mp_get_routes(doc4)
      doc5 <- mp_directions(
        origin = ori[3],
        destination = dest[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route5 <- mp_get_routes(doc5)
      leaflet() %>%
        addProviderTiles("CartoDB.Positron") %>%
        addCircleMarkers(data = origin_locs[6,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = origin_locs[7,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = origin_locs[3,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[7,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = destination_locs[6,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = destination_locs[3,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') %>%
        addPolylines(data=route4, opacity=1, weight=7, color='green')%>%
        addPolylines(data=route5, opacity=1, weight=7, color='brown')}
    else if(input$radio==5 & input$integer==0 & input$cluster==5){
      temp <- rhclust(pdis, ptim, 0, data)
      call=substr(temp$merge.route[[5]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[5]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      ori[2] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR_home,' ', data[row[2],col[2]:(col[2]+4)]$STREET_home, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE_home)
      ori[3] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR_home,' ', data[row[3],col[3]:(col[3]+4)]$STREET_home, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      dest[2] <- paste0(data[row[5],col[5]:(col[5]+4)]$HSE_NBR,' ', data[row[5],col[5]:(col[5]+4)]$STREET, ' ', data[row[5],col[5]:(col[5]+4)]$STTYPE, ', Milwaukee, WI ', data[row[5],col[5]:(col[5]+4)]$ZIP_CODE)
      dest[3] <- paste0(data[row[6],col[6]:(col[6]+4)]$HSE_NBR,' ', data[row[6],col[6]:(col[6]+4)]$STREET, ' ', data[row[6],col[6]:(col[6]+4)]$STTYPE, ', Milwaukee, WI ', data[row[6],col[6]:(col[6]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = ori[2],
        destination = ori[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = ori[3],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      doc4 <- mp_directions(
        origin = dest[1],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route4 <- mp_get_routes(doc4)
      doc5 <- mp_directions(
        origin = dest[2],
        destination = dest[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route5 <- mp_get_routes(doc5)
      leaflet() %>%
        addProviderTiles("CartoDB.Positron") %>%
        addCircleMarkers(data = origin_locs[2,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = origin_locs[1,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = origin_locs[10,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[10,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = destination_locs[1,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = destination_locs[2,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') %>%
        addPolylines(data=route4, opacity=1, weight=7, color='green')%>%
        addPolylines(data=route5, opacity=1, weight=7, color='brown')}
    else if(input$radio==1 & input$integer==.5 & input$cluster==2){
      temp <- rhclust(pdis, ptim, 0.5, data)
      call=substr(temp$merge.route[[1]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[1]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR,' ', data[row[2],col[2]:(col[2]+4)]$STREET, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE)
      ori[2] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR_home,' ', data[row[3],col[3]:(col[3]+4)]$STREET_home, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE_home)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = dest[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = ori[2],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[1,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[1,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[7,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[7,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow')
    }
    else if(input$radio==2 & input$integer==.5 & input$cluster==2){
      temp <- rhclust(pdis, ptim, 0.5, data)
      call=substr(temp$merge.route[[8]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[8]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      ori[2] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR_home,' ', data[row[2],col[2]:(col[2]+4)]$STREET_home, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR,' ', data[row[3],col[3]:(col[3]+4)]$STREET, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      ori[3] <- paste0(data[row[5],col[5]:(col[5]+4)]$HSE_NBR_home,' ', data[row[5],col[5]:(col[5]+4)]$STREET_home, ' ', data[row[5],col[5]:(col[5]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[5],col[5]:(col[5]+4)]$ZIP_CODE_home)
      ori[4] <- paste0(data[row[6],col[6]:(col[6]+4)]$HSE_NBR_home,' ', data[row[6],col[6]:(col[6]+4)]$STREET_home, ' ', data[row[6],col[6]:(col[6]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[6],col[6]:(col[6]+4)]$ZIP_CODE_home)
      dest[3] <- paste0(data[row[7],col[7]:(col[7]+4)]$HSE_NBR,' ', data[row[7],col[7]:(col[7]+4)]$STREET, ' ', data[row[7],col[7]:(col[7]+4)]$STTYPE, ', Milwaukee, WI ', data[row[7],col[7]:(col[7]+4)]$ZIP_CODE)
      dest[4] <- paste0(data[row[8],col[8]:(col[8]+4)]$HSE_NBR,' ', data[row[8],col[8]:(col[8]+4)]$STREET, ' ', data[row[8],col[8]:(col[8]+4)]$STTYPE, ', Milwaukee, WI ', data[row[8],col[8]:(col[8]+4)]$ZIP_CODE)
      ori[5] <- paste0(data[row[9],col[9]:(col[9]+4)]$HSE_NBR_home,' ', data[row[9],col[9]:(col[9]+4)]$STREET_home, ' ', data[row[9],col[9]:(col[9]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[9],col[9]:(col[9]+4)]$ZIP_CODE_home)
      dest[5] <- paste0(data[row[10],col[10]:(col[10]+4)]$HSE_NBR,' ', data[row[10],col[10]:(col[10]+4)]$STREET, ' ', data[row[10],col[10]:(col[10]+4)]$STTYPE, ', Milwaukee, WI ', data[row[10],col[10]:(col[10]+4)]$ZIP_CODE)
      ori[6] <- paste0(data[row[11],col[11]:(col[11]+4)]$HSE_NBR_home,' ', data[row[11],col[11]:(col[11]+4)]$STREET_home, ' ', data[row[11],col[11]:(col[11]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[11],col[11]:(col[11]+4)]$ZIP_CODE_home)
      dest[6] <- paste0(data[row[12],col[12]:(col[12]+4)]$HSE_NBR,' ', data[row[12],col[12]:(col[12]+4)]$STREET, ' ', data[row[12],col[12]:(col[12]+4)]$STTYPE, ', Milwaukee, WI ', data[row[12],col[12]:(col[12]+4)]$ZIP_CODE)
      ori[7] <- paste0(data[row[13],col[13]:(col[13]+4)]$HSE_NBR_home,' ', data[row[13],col[13]:(col[13]+4)]$STREET_home, ' ', data[row[13],col[13]:(col[13]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[13],col[13]:(col[13]+4)]$ZIP_CODE_home)
      ori[8] <- paste0(data[row[14],col[14]:(col[14]+4)]$HSE_NBR_home,' ', data[row[14],col[14]:(col[14]+4)]$STREET_home, ' ', data[row[14],col[14]:(col[14]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[14],col[14]:(col[14]+4)]$ZIP_CODE_home)
      dest[7] <- paste0(data[row[15],col[15]:(col[15]+4)]$HSE_NBR,' ', data[row[15],col[15]:(col[15]+4)]$STREET, ' ', data[row[15],col[15]:(col[15]+4)]$STTYPE, ', Milwaukee, WI ', data[row[15],col[15]:(col[15]+4)]$ZIP_CODE)
      dest[8] <- paste0(data[row[16],col[16]:(col[16]+4)]$HSE_NBR,' ', data[row[16],col[16]:(col[16]+4)]$STREET, ' ', data[row[16],col[16]:(col[16]+4)]$STTYPE, ', Milwaukee, WI ', data[row[16],col[16]:(col[16]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = ori[2],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = dest[1],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      doc4 <- mp_directions(
        origin = dest[2],
        destination = ori[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route4 <- mp_get_routes(doc4)
      doc5 <- mp_directions(
        origin = ori[3],
        destination = ori[4],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route5 <- mp_get_routes(doc5)
      doc6 <- mp_directions(
        origin = ori[4],
        destination = dest[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route6 <- mp_get_routes(doc6)
      doc7 <- mp_directions(
        origin = dest[3],
        destination = dest[4],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route7 <- mp_get_routes(doc7)
      doc8 <- mp_directions(
        origin = dest[4],
        destination = ori[5],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route8 <- mp_get_routes(doc8)
      doc9 <- mp_directions(
        origin = ori[5],
        destination = dest[5],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route9 <- mp_get_routes(doc9)
      doc10 <- mp_directions(
        origin = dest[5],
        destination = ori[6],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route10 <- mp_get_routes(doc10)
      doc11 <- mp_directions(
        origin = ori[6],
        destination = dest[6],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route11 <- mp_get_routes(doc11)
      doc12 <- mp_directions(
        origin = dest[6],
        destination = ori[7],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route12 <- mp_get_routes(doc12)
      doc13 <- mp_directions(
        origin = ori[7],
        destination = ori[8],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route13 <- mp_get_routes(doc13)
      doc14 <- mp_directions(
        origin = ori[8],
        destination = dest[7],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route14 <- mp_get_routes(doc14)
      doc15 <- mp_directions(
        origin = dest[7],
        destination = dest[8],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route15 <- mp_get_routes(doc15)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[4,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[4,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[2,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[2,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[3,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[3,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[5,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[5,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[6,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[6,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[8,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[8,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[9,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[9,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[10,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[10,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') %>%
        addPolylines(data=route4, opacity=1, weight=7, color='green') %>%
        addPolylines(data=route5, opacity=1, weight=7, color='brown') %>%
        addPolylines(data=route6, opacity=1, weight=7, color='purple') %>%
        addPolylines(data=route7, opacity=1, weight=7, color='orange') %>%
        addPolylines(data=route8, opacity=1, weight=7, color='black') %>%
        addPolylines(data=route9, opacity=1, weight=7, color='white') %>%
        addPolylines(data=route10, opacity=1, weight=7, color='grey') %>%
        addPolylines(data=route11, opacity=1, weight=7, color='blue') %>%
        addPolylines(data=route12, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route13, opacity=1, weight=7, color='yellow') %>%
        addPolylines(data=route14, opacity=1, weight=7, color='green') %>%
        addPolylines(data=route15, opacity=1, weight=7, color='brown') 
    }
    else if(input$radio==1 & input$integer==.5 & input$cluster==3){
      temp <- rhclust(pdis, ptim, 0.5, data)
      call=substr(temp$merge.route[[1]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[1]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR,' ', data[row[2],col[2]:(col[2]+4)]$STREET, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE)
      ori[2] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR_home,' ', data[row[3],col[3]:(col[3]+4)]$STREET_home, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE_home)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = dest[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = ori[2],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[1,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[1,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[7,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[7,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow')
    }
    else if(input$radio==2 & input$integer==.5 & input$cluster==3){
      temp <- rhclust(pdis, ptim, 0.5, data)
      call=substr(temp$merge.route[[6]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[6]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      ori[2] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR_home,' ', data[row[2],col[2]:(col[2]+4)]$STREET_home, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR,' ', data[row[3],col[3]:(col[3]+4)]$STREET, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      ori[3] <- paste0(data[row[5],col[5]:(col[5]+4)]$HSE_NBR_home,' ', data[row[5],col[5]:(col[5]+4)]$STREET_home, ' ', data[row[5],col[5]:(col[5]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[5],col[5]:(col[5]+4)]$ZIP_CODE_home)
      ori[4] <- paste0(data[row[6],col[6]:(col[6]+4)]$HSE_NBR_home,' ', data[row[6],col[6]:(col[6]+4)]$STREET_home, ' ', data[row[6],col[6]:(col[6]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[6],col[6]:(col[6]+4)]$ZIP_CODE_home)
      dest[3] <- paste0(data[row[7],col[7]:(col[7]+4)]$HSE_NBR,' ', data[row[7],col[7]:(col[7]+4)]$STREET, ' ', data[row[7],col[7]:(col[7]+4)]$STTYPE, ', Milwaukee, WI ', data[row[7],col[7]:(col[7]+4)]$ZIP_CODE)
      dest[4] <- paste0(data[row[8],col[8]:(col[8]+4)]$HSE_NBR,' ', data[row[8],col[8]:(col[8]+4)]$STREET, ' ', data[row[8],col[8]:(col[8]+4)]$STTYPE, ', Milwaukee, WI ', data[row[8],col[8]:(col[8]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = ori[2],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = dest[1],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      doc4 <- mp_directions(
        origin = dest[2],
        destination = ori[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route4 <- mp_get_routes(doc4)
      doc5 <- mp_directions(
        origin = ori[3],
        destination = ori[4],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route5 <- mp_get_routes(doc5)
      doc6 <- mp_directions(
        origin = ori[4],
        destination = dest[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route6 <- mp_get_routes(doc6)
      doc7 <- mp_directions(
        origin = dest[3],
        destination = dest[4],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route7 <- mp_get_routes(doc7)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[2,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[2,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[3,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[3,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[5,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[5,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[8,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[8,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') %>%
        addPolylines(data=route4, opacity=1, weight=7, color='green') %>%
        addPolylines(data=route5, opacity=1, weight=7, color='brown') %>%
        addPolylines(data=route6, opacity=1, weight=7, color='purple') %>%
        addPolylines(data=route7, opacity=1, weight=7, color='orange') 
    }
    else if(input$radio==3 & input$integer==.5 & input$cluster==3){
      temp <- rhclust(pdis, ptim, 0.5, data)
      call=substr(temp$merge.route[[7]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[7]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR,' ', data[row[2],col[2]:(col[2]+4)]$STREET, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE)
      ori[2] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR_home,' ', data[row[3],col[3]:(col[3]+4)]$STREET_home, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE_home)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      ori[3] <- paste0(data[row[5],col[5]:(col[5]+4)]$HSE_NBR_home,' ', data[row[5],col[5]:(col[5]+4)]$STREET_home, ' ', data[row[5],col[5]:(col[5]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[5],col[5]:(col[5]+4)]$ZIP_CODE_home)
      dest[3] <- paste0(data[row[6],col[6]:(col[6]+4)]$HSE_NBR,' ', data[row[6],col[6]:(col[6]+4)]$STREET, ' ', data[row[6],col[6]:(col[6]+4)]$STTYPE, ', Milwaukee, WI ', data[row[6],col[6]:(col[6]+4)]$ZIP_CODE)
      ori[4] <- paste0(data[row[7],col[7]:(col[7]+4)]$HSE_NBR_home,' ', data[row[7],col[7]:(col[7]+4)]$STREET_home, ' ', data[row[7],col[7]:(col[7]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[7],col[7]:(col[7]+4)]$ZIP_CODE_home)
      dest[4] <- paste0(data[row[8],col[8]:(col[8]+4)]$HSE_NBR,' ', data[row[8],col[8]:(col[8]+4)]$STREET, ' ', data[row[8],col[8]:(col[8]+4)]$STTYPE, ', Milwaukee, WI ', data[row[8],col[8]:(col[8]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = dest[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = ori[2],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      doc4 <- mp_directions(
        origin = dest[2],
        destination = ori[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route4 <- mp_get_routes(doc4)
      doc5 <- mp_directions(
        origin = ori[3],
        destination = dest[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route5 <- mp_get_routes(doc5)
      doc6 <- mp_directions(
        origin = dest[3],
        destination = ori[4],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route6 <- mp_get_routes(doc6)
      doc7 <- mp_directions(
        origin = ori[4],
        destination = dest[4],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route7 <- mp_get_routes(doc7)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[4,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[4,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[6,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[6,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[9,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[9,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[10,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[10,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') %>%
        addPolylines(data=route4, opacity=1, weight=7, color='green') %>%
        addPolylines(data=route5, opacity=1, weight=7, color='brown') %>%
        addPolylines(data=route6, opacity=1, weight=7, color='purple') %>%
        addPolylines(data=route7, opacity=1, weight=7, color='orange') 
      ###
    }
    else if(input$radio==1 & input$integer==.5 & input$cluster==4){
      temp <- rhclust(pdis, ptim, 0.5, data)
      call=substr(temp$merge.route[[1]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[1]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR,' ', data[row[2],col[2]:(col[2]+4)]$STREET, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE)
      ori[2] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR_home,' ', data[row[3],col[3]:(col[3]+4)]$STREET_home, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE_home)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = dest[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = ori[2],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[1,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[1,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[7,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[7,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow')
    }
    else if(input$radio==2 & input$integer==.5 & input$cluster==4){
      temp <- rhclust(pdis, ptim, 0.5, data)
      call=substr(temp$merge.route[[6]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[6]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      ori[2] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR_home,' ', data[row[2],col[2]:(col[2]+4)]$STREET_home, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR,' ', data[row[3],col[3]:(col[3]+4)]$STREET, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      ori[3] <- paste0(data[row[5],col[5]:(col[5]+4)]$HSE_NBR_home,' ', data[row[5],col[5]:(col[5]+4)]$STREET_home, ' ', data[row[5],col[5]:(col[5]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[5],col[5]:(col[5]+4)]$ZIP_CODE_home)
      ori[4] <- paste0(data[row[6],col[6]:(col[6]+4)]$HSE_NBR_home,' ', data[row[6],col[6]:(col[6]+4)]$STREET_home, ' ', data[row[6],col[6]:(col[6]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[6],col[6]:(col[6]+4)]$ZIP_CODE_home)
      dest[3] <- paste0(data[row[7],col[7]:(col[7]+4)]$HSE_NBR,' ', data[row[7],col[7]:(col[7]+4)]$STREET, ' ', data[row[7],col[7]:(col[7]+4)]$STTYPE, ', Milwaukee, WI ', data[row[7],col[7]:(col[7]+4)]$ZIP_CODE)
      dest[4] <- paste0(data[row[8],col[8]:(col[8]+4)]$HSE_NBR,' ', data[row[8],col[8]:(col[8]+4)]$STREET, ' ', data[row[8],col[8]:(col[8]+4)]$STTYPE, ', Milwaukee, WI ', data[row[8],col[8]:(col[8]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = ori[2],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = dest[1],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      doc4 <- mp_directions(
        origin = dest[2],
        destination = ori[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route4 <- mp_get_routes(doc4)
      doc5 <- mp_directions(
        origin = ori[3],
        destination = ori[4],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route5 <- mp_get_routes(doc5)
      doc6 <- mp_directions(
        origin = ori[4],
        destination = dest[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route6 <- mp_get_routes(doc6)
      doc7 <- mp_directions(
        origin = dest[3],
        destination = dest[4],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route7 <- mp_get_routes(doc7)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[2,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[2,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[3,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[3,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[5,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[5,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[8,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[8,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') %>%
        addPolylines(data=route4, opacity=1, weight=7, color='green') %>%
        addPolylines(data=route5, opacity=1, weight=7, color='brown') %>%
        addPolylines(data=route6, opacity=1, weight=7, color='purple') %>%
        addPolylines(data=route7, opacity=1, weight=7, color='orange') 
    }
    else if(input$radio==3 & input$integer==.5 & input$cluster==4){
      temp <- rhclust(pdis, ptim, 0.5, data)
      call=substr(temp$merge.route[[2]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[2]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR,' ', data[row[2],col[2]:(col[2]+4)]$STREET, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE)
      ori[2] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR_home,' ', data[row[3],col[3]:(col[3]+4)]$STREET_home, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE_home)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = dest[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = ori[2],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[9,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[9,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[10,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[10,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') 
    }
    else if(input$radio==4 & input$integer==.5 & input$cluster==4){
      temp <- rhclust(pdis, ptim, 0.5, data)
      call=substr(temp$merge.route[[4]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[4]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR,' ', data[row[2],col[2]:(col[2]+4)]$STREET, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE)
      ori[2] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR_home,' ', data[row[3],col[3]:(col[3]+4)]$STREET_home, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE_home)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = dest[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = ori[2],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[4,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[4,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[6,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[6,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') 
    }
    else if(input$radio==1 & input$integer==.5 & input$cluster==5){
      temp <- rhclust(pdis, ptim, 0.5, data)
      call=substr(temp$merge.route[[1]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[1]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR,' ', data[row[2],col[2]:(col[2]+4)]$STREET, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE)
      ori[2] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR_home,' ', data[row[3],col[3]:(col[3]+4)]$STREET_home, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE_home)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = dest[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = ori[2],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[1,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[1,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[7,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[7,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow')
    }
    else if(input$radio==2 & input$integer==.5 & input$cluster==5){
      temp <- rhclust(pdis, ptim, 0.5, data)
      call=substr(temp$merge.route[[5]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[5]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      ori[2] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR_home,' ', data[row[2],col[2]:(col[2]+4)]$STREET_home, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR,' ', data[row[3],col[3]:(col[3]+4)]$STREET, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = ori[2],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = dest[1],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[2,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[2,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[5,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[5,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') 
    }
    else if(input$radio==3 & input$integer==.5 & input$cluster==5){
      temp <- rhclust(pdis, ptim, 0.5, data)
      call=substr(temp$merge.route[[3]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[3]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      ori[2] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR_home,' ', data[row[2],col[2]:(col[2]+4)]$STREET_home, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR,' ', data[row[3],col[3]:(col[3]+4)]$STREET, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = ori[2],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = dest[1],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[3,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[3,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[8,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[8,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') 
    }
    else if(input$radio==4 & input$integer==.5 & input$cluster==5){
      temp <- rhclust(pdis, ptim, 0.5, data)
      call=substr(temp$merge.route[[2]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[2]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR,' ', data[row[2],col[2]:(col[2]+4)]$STREET, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE)
      ori[2] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR_home,' ', data[row[3],col[3]:(col[3]+4)]$STREET_home, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE_home)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = dest[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = ori[2],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[9,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[9,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[10,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[10,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') 
    }
    else if(input$radio==5 & input$integer==.5 & input$cluster==5){
      temp <- rhclust(pdis, ptim, 0.5, data)
      call=substr(temp$merge.route[[4]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[4]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR,' ', data[row[2],col[2]:(col[2]+4)]$STREET, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE)
      ori[2] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR_home,' ', data[row[3],col[3]:(col[3]+4)]$STREET_home, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE_home)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = dest[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = ori[2],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[4,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[4,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[6,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[6,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') 
    }
    else if(input$radio==1 & input$integer==1 & input$cluster==2){
      temp <- rhclust(pdis, ptim, 1, data)
      call=substr(temp$merge.route[[4]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[4]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      ori[2] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR_home,' ', data[row[2],col[2]:(col[2]+4)]$STREET_home, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR,' ', data[row[3],col[3]:(col[3]+4)]$STREET, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = ori[2],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = dest[1],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[5,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[5,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[6,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[6,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') 
    }
    else if(input$radio==2 & input$integer==1 & input$cluster==2){
      temp <- rhclust(pdis, ptim, 1, data)
      call=substr(temp$merge.route[[8]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[8]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR,' ', data[row[2],col[2]:(col[2]+4)]$STREET, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE)
      ori[2] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR_home,' ', data[row[3],col[3]:(col[3]+4)]$STREET_home, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE_home)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      ori[3] <- paste0(data[row[5],col[5]:(col[5]+4)]$HSE_NBR_home,' ', data[row[5],col[5]:(col[5]+4)]$STREET_home, ' ', data[row[5],col[5]:(col[5]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[5],col[5]:(col[5]+4)]$ZIP_CODE_home)
      dest[3] <- paste0(data[row[6],col[6]:(col[6]+4)]$HSE_NBR,' ', data[row[6],col[6]:(col[6]+4)]$STREET, ' ', data[row[6],col[6]:(col[6]+4)]$STTYPE, ', Milwaukee, WI ', data[row[6],col[6]:(col[6]+4)]$ZIP_CODE)
      ori[4] <- paste0(data[row[7],col[7]:(col[7]+4)]$HSE_NBR_home,' ', data[row[7],col[7]:(col[7]+4)]$STREET_home, ' ', data[row[7],col[7]:(col[7]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[7],col[7]:(col[7]+4)]$ZIP_CODE_home)
      dest[4] <- paste0(data[row[8],col[8]:(col[8]+4)]$HSE_NBR,' ', data[row[8],col[8]:(col[8]+4)]$STREET, ' ', data[row[8],col[8]:(col[8]+4)]$STTYPE, ', Milwaukee, WI ', data[row[8],col[8]:(col[8]+4)]$ZIP_CODE)
      ori[5] <- paste0(data[row[9],col[9]:(col[9]+4)]$HSE_NBR_home,' ', data[row[9],col[9]:(col[9]+4)]$STREET_home, ' ', data[row[9],col[9]:(col[9]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[9],col[9]:(col[9]+4)]$ZIP_CODE_home)
      dest[5] <- paste0(data[row[10],col[10]:(col[10]+4)]$HSE_NBR,' ', data[row[10],col[10]:(col[10]+4)]$STREET, ' ', data[row[10],col[10]:(col[10]+4)]$STTYPE, ', Milwaukee, WI ', data[row[10],col[10]:(col[10]+4)]$ZIP_CODE)
      ori[6] <- paste0(data[row[11],col[11]:(col[11]+4)]$HSE_NBR_home,' ', data[row[11],col[11]:(col[11]+4)]$STREET_home, ' ', data[row[11],col[11]:(col[11]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[11],col[11]:(col[11]+4)]$ZIP_CODE_home)
      dest[6] <- paste0(data[row[12],col[12]:(col[12]+4)]$HSE_NBR,' ', data[row[12],col[12]:(col[12]+4)]$STREET, ' ', data[row[12],col[12]:(col[12]+4)]$STTYPE, ', Milwaukee, WI ', data[row[12],col[12]:(col[12]+4)]$ZIP_CODE)
      ori[7] <- paste0(data[row[13],col[13]:(col[13]+4)]$HSE_NBR_home,' ', data[row[13],col[13]:(col[13]+4)]$STREET_home, ' ', data[row[13],col[13]:(col[13]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[13],col[13]:(col[13]+4)]$ZIP_CODE_home)
      ori[8] <- paste0(data[row[14],col[14]:(col[14]+4)]$HSE_NBR_home,' ', data[row[14],col[14]:(col[14]+4)]$STREET_home, ' ', data[row[14],col[14]:(col[14]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[14],col[14]:(col[14]+4)]$ZIP_CODE_home)
      dest[7] <- paste0(data[row[15],col[15]:(col[15]+4)]$HSE_NBR,' ', data[row[15],col[15]:(col[15]+4)]$STREET, ' ', data[row[15],col[15]:(col[15]+4)]$STTYPE, ', Milwaukee, WI ', data[row[15],col[15]:(col[15]+4)]$ZIP_CODE)
      dest[8] <- paste0(data[row[16],col[16]:(col[16]+4)]$HSE_NBR,' ', data[row[16],col[16]:(col[16]+4)]$STREET, ' ', data[row[16],col[16]:(col[16]+4)]$STTYPE, ', Milwaukee, WI ', data[row[16],col[16]:(col[16]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = dest[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = ori[2],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      doc4 <- mp_directions(
        origin = dest[2],
        destination = ori[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route4 <- mp_get_routes(doc4)
      doc5 <- mp_directions(
        origin = ori[3],
        destination = dest[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route5 <- mp_get_routes(doc5)
      doc6 <- mp_directions(
        origin = dest[3],
        destination = ori[4],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route6 <- mp_get_routes(doc6)
      doc7 <- mp_directions(
        origin = ori[4],
        destination = dest[4],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route7 <- mp_get_routes(doc7)
      doc8 <- mp_directions(
        origin = dest[4],
        destination = ori[5],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route8 <- mp_get_routes(doc8)
      doc9 <- mp_directions(
        origin = ori[5],
        destination = dest[5],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route9 <- mp_get_routes(doc9)
      doc10 <- mp_directions(
        origin = dest[5],
        destination = ori[6],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route10 <- mp_get_routes(doc10)
      doc11 <- mp_directions(
        origin = ori[6],
        destination = dest[6],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route11 <- mp_get_routes(doc11)
      doc12 <- mp_directions(
        origin = dest[6],
        destination = ori[7],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route12 <- mp_get_routes(doc12)
      doc13 <- mp_directions(
        origin = ori[7],
        destination = ori[8],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route13 <- mp_get_routes(doc13)
      doc14 <- mp_directions(
        origin = ori[8],
        destination = dest[7],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route14 <- mp_get_routes(doc14)
      doc15 <- mp_directions(
        origin = dest[7],
        destination = dest[8],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route15 <- mp_get_routes(doc15)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[1,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[1,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[2,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[2,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[3,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[3,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[4,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[4,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[7,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[7,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[8,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[8,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[9,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[9,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[10,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[10,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') %>%
        addPolylines(data=route4, opacity=1, weight=7, color='green') %>%
        addPolylines(data=route5, opacity=1, weight=7, color='brown') %>%
        addPolylines(data=route6, opacity=1, weight=7, color='purple') %>%
        addPolylines(data=route7, opacity=1, weight=7, color='orange') %>%
        addPolylines(data=route8, opacity=1, weight=7, color='black') %>%
        addPolylines(data=route9, opacity=1, weight=7, color='white') %>%
        addPolylines(data=route10, opacity=1, weight=7, color='grey') %>%
        addPolylines(data=route11, opacity=1, weight=7, color='blue') %>%
        addPolylines(data=route12, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route13, opacity=1, weight=7, color='yellow') %>%
        addPolylines(data=route14, opacity=1, weight=7, color='green') %>%
        addPolylines(data=route15, opacity=1, weight=7, color='brown') 
    }
    else if(input$radio==1 & input$integer==1 & input$cluster==3){
      temp <- rhclust(pdis, ptim, 1, data)
      call=substr(temp$merge.route[[4]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[4]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      ori[2] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR_home,' ', data[row[2],col[2]:(col[2]+4)]$STREET_home, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR,' ', data[row[3],col[3]:(col[3]+4)]$STREET, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = ori[2],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = dest[1],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[5,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[5,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[6,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[6,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') 
    }
    else if(input$radio==2 & input$integer==1 & input$cluster==3){
      temp <- rhclust(pdis, ptim, 1, data)
      call=substr(temp$merge.route[[6]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[6]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR,' ', data[row[2],col[2]:(col[2]+4)]$STREET, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE)
      ori[2] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR_home,' ', data[row[3],col[3]:(col[3]+4)]$STREET_home, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE_home)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      ori[3] <- paste0(data[row[5],col[5]:(col[5]+4)]$HSE_NBR_home,' ', data[row[5],col[5]:(col[5]+4)]$STREET_home, ' ', data[row[5],col[5]:(col[5]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[5],col[5]:(col[5]+4)]$ZIP_CODE_home)
      dest[3] <- paste0(data[row[6],col[6]:(col[6]+4)]$HSE_NBR,' ', data[row[6],col[6]:(col[6]+4)]$STREET, ' ', data[row[6],col[6]:(col[6]+4)]$STTYPE, ', Milwaukee, WI ', data[row[6],col[6]:(col[6]+4)]$ZIP_CODE)
      ori[4] <- paste0(data[row[7],col[7]:(col[7]+4)]$HSE_NBR_home,' ', data[row[7],col[7]:(col[7]+4)]$STREET_home, ' ', data[row[7],col[7]:(col[7]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[7],col[7]:(col[7]+4)]$ZIP_CODE_home)
      dest[4] <- paste0(data[row[8],col[8]:(col[8]+4)]$HSE_NBR,' ', data[row[8],col[8]:(col[8]+4)]$STREET, ' ', data[row[8],col[8]:(col[8]+4)]$STTYPE, ', Milwaukee, WI ', data[row[8],col[8]:(col[8]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = dest[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = ori[2],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      doc4 <- mp_directions(
        origin = dest[2],
        destination = ori[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route4 <- mp_get_routes(doc4)
      doc5 <- mp_directions(
        origin = ori[3],
        destination = dest[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route5 <- mp_get_routes(doc5)
      doc6 <- mp_directions(
        origin = dest[3],
        destination = ori[4],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route6 <- mp_get_routes(doc6)
      doc7 <- mp_directions(
        origin = ori[4],
        destination = dest[4],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route7 <- mp_get_routes(doc7)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[7,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[7,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[8,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[8,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[9,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[9,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[10,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[10,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') %>%
        addPolylines(data=route4, opacity=1, weight=7, color='green') %>%
        addPolylines(data=route5, opacity=1, weight=7, color='brown') %>%
        addPolylines(data=route6, opacity=1, weight=7, color='purple') %>%
        addPolylines(data=route7, opacity=1, weight=7, color='orange')
    }
    else if(input$radio==3 & input$integer==1 & input$cluster==3){
      temp <- rhclust(pdis, ptim, 1, data)
      call=substr(temp$merge.route[[7]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[7]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR,' ', data[row[2],col[2]:(col[2]+4)]$STREET, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE)
      ori[2] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR_home,' ', data[row[3],col[3]:(col[3]+4)]$STREET_home, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE_home)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      ori[3] <- paste0(data[row[5],col[5]:(col[5]+4)]$HSE_NBR_home,' ', data[row[5],col[5]:(col[5]+4)]$STREET_home, ' ', data[row[5],col[5]:(col[5]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[5],col[5]:(col[5]+4)]$ZIP_CODE_home)
      ori[4] <- paste0(data[row[6],col[6]:(col[6]+4)]$HSE_NBR_home,' ', data[row[6],col[6]:(col[6]+4)]$STREET_home, ' ', data[row[6],col[6]:(col[6]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[6],col[6]:(col[6]+4)]$ZIP_CODE_home)
      dest[3] <- paste0(data[row[7],col[7]:(col[7]+4)]$HSE_NBR,' ', data[row[7],col[7]:(col[7]+4)]$STREET, ' ', data[row[7],col[7]:(col[7]+4)]$STTYPE, ', Milwaukee, WI ', data[row[7],col[7]:(col[7]+4)]$ZIP_CODE)
      dest[4] <- paste0(data[row[8],col[8]:(col[8]+4)]$HSE_NBR,' ', data[row[8],col[8]:(col[8]+4)]$STREET, ' ', data[row[8],col[8]:(col[8]+4)]$STTYPE, ', Milwaukee, WI ', data[row[8],col[8]:(col[8]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = dest[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = ori[2],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      doc4 <- mp_directions(
        origin = dest[2],
        destination = ori[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route4 <- mp_get_routes(doc4)
      doc5 <- mp_directions(
        origin = ori[3],
        destination = ori[4],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route5 <- mp_get_routes(doc5)
      doc6 <- mp_directions(
        origin = ori[4],
        destination = dest[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route6 <- mp_get_routes(doc6)
      doc7 <- mp_directions(
        origin = dest[3],
        destination = dest[4],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route7 <- mp_get_routes(doc7)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[1,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[1,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[2,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[2,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[3,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[3,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[4,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[4,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') %>%
        addPolylines(data=route4, opacity=1, weight=7, color='green') %>%
        addPolylines(data=route5, opacity=1, weight=7, color='brown') %>%
        addPolylines(data=route6, opacity=1, weight=7, color='purple') %>%
        addPolylines(data=route7, opacity=1, weight=7, color='orange')
    }
    else if(input$radio==1 & input$integer==1 & input$cluster==4){
      temp <- rhclust(pdis, ptim, 1, data)
      call=substr(temp$merge.route[[4]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[4]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      ori[2] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR_home,' ', data[row[2],col[2]:(col[2]+4)]$STREET_home, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR,' ', data[row[3],col[3]:(col[3]+4)]$STREET, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = ori[2],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = dest[1],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[5,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[5,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[6,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[6,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') 
    }
    else if(input$radio==2 & input$integer==1 & input$cluster==4){
      temp <- rhclust(pdis, ptim, 1, data)
      call=substr(temp$merge.route[[6]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[6]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR,' ', data[row[2],col[2]:(col[2]+4)]$STREET, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE)
      ori[2] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR_home,' ', data[row[3],col[3]:(col[3]+4)]$STREET_home, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE_home)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      ori[3] <- paste0(data[row[5],col[5]:(col[5]+4)]$HSE_NBR_home,' ', data[row[5],col[5]:(col[5]+4)]$STREET_home, ' ', data[row[5],col[5]:(col[5]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[5],col[5]:(col[5]+4)]$ZIP_CODE_home)
      dest[3] <- paste0(data[row[6],col[6]:(col[6]+4)]$HSE_NBR,' ', data[row[6],col[6]:(col[6]+4)]$STREET, ' ', data[row[6],col[6]:(col[6]+4)]$STTYPE, ', Milwaukee, WI ', data[row[6],col[6]:(col[6]+4)]$ZIP_CODE)
      ori[4] <- paste0(data[row[7],col[7]:(col[7]+4)]$HSE_NBR_home,' ', data[row[7],col[7]:(col[7]+4)]$STREET_home, ' ', data[row[7],col[7]:(col[7]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[7],col[7]:(col[7]+4)]$ZIP_CODE_home)
      dest[4] <- paste0(data[row[8],col[8]:(col[8]+4)]$HSE_NBR,' ', data[row[8],col[8]:(col[8]+4)]$STREET, ' ', data[row[8],col[8]:(col[8]+4)]$STTYPE, ', Milwaukee, WI ', data[row[8],col[8]:(col[8]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = dest[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = ori[2],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      doc4 <- mp_directions(
        origin = dest[2],
        destination = ori[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route4 <- mp_get_routes(doc4)
      doc5 <- mp_directions(
        origin = ori[3],
        destination = dest[3],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route5 <- mp_get_routes(doc5)
      doc6 <- mp_directions(
        origin = dest[3],
        destination = ori[4],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route6 <- mp_get_routes(doc6)
      doc7 <- mp_directions(
        origin = ori[4],
        destination = dest[4],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route7 <- mp_get_routes(doc7)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[7,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[7,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[8,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[8,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[9,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[9,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[10,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[10,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') %>%
        addPolylines(data=route4, opacity=1, weight=7, color='green') %>%
        addPolylines(data=route5, opacity=1, weight=7, color='brown') %>%
        addPolylines(data=route6, opacity=1, weight=7, color='purple') %>%
        addPolylines(data=route7, opacity=1, weight=7, color='orange')
    }
    else if(input$radio==3 & input$integer==1 & input$cluster==4){
      temp <- rhclust(pdis, ptim, 1, data)
      call=substr(temp$merge.route[[5]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[5]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR,' ', data[row[2],col[2]:(col[2]+4)]$STREET, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE)
      ori[2] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR_home,' ', data[row[3],col[3]:(col[3]+4)]$STREET_home, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE_home)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = dest[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = ori[2],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[2,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[2,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[3,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[3,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow')
    }
    else if(input$radio==4 & input$integer==1 & input$cluster==4){
      temp <- rhclust(pdis, ptim, 1, data)
      call=substr(temp$merge.route[[3]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[3]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      ori[2] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR_home,' ', data[row[2],col[2]:(col[2]+4)]$STREET_home, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR,' ', data[row[3],col[3]:(col[3]+4)]$STREET, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = ori[2],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = dest[1],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[1,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[1,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[4,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[4,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow')
    }
    else if(input$radio==1 & input$integer==1 & input$cluster==5){
      temp <- rhclust(pdis, ptim, 1, data)
      call=substr(temp$merge.route[[4]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[4]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      ori[2] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR_home,' ', data[row[2],col[2]:(col[2]+4)]$STREET_home, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR,' ', data[row[3],col[3]:(col[3]+4)]$STREET, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = ori[2],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = dest[1],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[5,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[5,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[6,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[6,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow') 
    }
    else if(input$radio==2 & input$integer==1 & input$cluster==5){
      temp <- rhclust(pdis, ptim, 1, data)
      call=substr(temp$merge.route[[2]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[2]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR,' ', data[row[2],col[2]:(col[2]+4)]$STREET, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE)
      ori[2] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR_home,' ', data[row[3],col[3]:(col[3]+4)]$STREET_home, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE_home)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = dest[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = ori[2],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[7,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[7,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[10,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[10,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow')
    }
    else if(input$radio==3 & input$integer==1 & input$cluster==5){
      temp <- rhclust(pdis, ptim, 1, data)
      call=substr(temp$merge.route[[1]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[1]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR,' ', data[row[2],col[2]:(col[2]+4)]$STREET, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE)
      ori[2] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR_home,' ', data[row[3],col[3]:(col[3]+4)]$STREET_home, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE_home)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = dest[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = ori[2],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[8,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[8,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[9,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[9,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow')
    }
    else if(input$radio==4 & input$integer==1 & input$cluster==5){
      temp <- rhclust(pdis, ptim, 1, data)
      call=substr(temp$merge.route[[5]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[5]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR,' ', data[row[2],col[2]:(col[2]+4)]$STREET, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE)
      ori[2] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR_home,' ', data[row[3],col[3]:(col[3]+4)]$STREET_home, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE_home)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = dest[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = ori[2],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[2,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[2,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[3,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[3,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow')
    }
    else if(input$radio==5 & input$integer==1 & input$cluster==5){
      temp <- rhclust(pdis, ptim, 1, data)
      call=substr(temp$merge.route[[3]],1,1)
      call[call=='p'] <- 5
      call[call=='d'] <- 10
      col=as.numeric(call)
      row=as.numeric(substr(temp$merge.route[[3]],2,4))
      ori <- vector()
      dest <- vector()
      ori[1] <- paste0(data[row[1],col[1]:(col[1]+4)]$HSE_NBR_home,' ', data[row[1],col[1]:(col[1]+4)]$STREET_home, ' ', data[row[1],col[1]:(col[1]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[1],col[1]:(col[1]+4)]$ZIP_CODE_home)
      ori[2] <- paste0(data[row[2],col[2]:(col[2]+4)]$HSE_NBR_home,' ', data[row[2],col[2]:(col[2]+4)]$STREET_home, ' ', data[row[2],col[2]:(col[2]+4)]$STTYPE_home, ', Milwaukee, WI ', data[row[2],col[2]:(col[2]+4)]$ZIP_CODE_home)
      dest[1] <- paste0(data[row[3],col[3]:(col[3]+4)]$HSE_NBR,' ', data[row[3],col[3]:(col[3]+4)]$STREET, ' ', data[row[3],col[3]:(col[3]+4)]$STTYPE, ', Milwaukee, WI ', data[row[3],col[3]:(col[3]+4)]$ZIP_CODE)
      dest[2] <- paste0(data[row[4],col[4]:(col[4]+4)]$HSE_NBR,' ', data[row[4],col[4]:(col[4]+4)]$STREET, ' ', data[row[4],col[4]:(col[4]+4)]$STTYPE, ', Milwaukee, WI ', data[row[4],col[4]:(col[4]+4)]$ZIP_CODE)
      doc <- mp_directions(
        origin = ori[1],
        destination = ori[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route <- mp_get_routes(doc)
      doc2 <- mp_directions(
        origin = ori[2],
        destination = dest[1],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route2 <- mp_get_routes(doc2)
      doc3 <- mp_directions(
        origin = dest[1],
        destination = dest[2],
        alternatives = FALSE,
        key='AIzaSyAwCU0w7-3pLYwtSW_6tA0yRi7B5ENYsGg'
      )
      route3 <- mp_get_routes(doc3)
      
      leaflet() %>%
        addProviderTiles('CartoDB.Positron') %>%
        addCircleMarkers(data = origin_locs[1,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[1,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addCircleMarkers(data = origin_locs[4,], label = ~index,  labelOptions = labelOptions(noHide = T), color="blue") %>%
        addCircleMarkers(data = destination_locs[4,], label = ~index,  labelOptions = labelOptions(noHide = T), color="red") %>%
        addPolylines(data=route, opacity=1, weight=7) %>%
        addPolylines(data=route2, opacity=1, weight=7, color='red') %>%
        addPolylines(data=route3, opacity=1, weight=7, color='yellow')
    }
  })

  
  sliderValues <- reactive({
    data.frame(
      Name = c("Time Factor"),
      Value = as.character(c(input$integer)),
      stringsAsFactors = FALSE)
    
  })
  
  output$values <- renderTable({
    sliderValues()
  })
  
  #OutputPlotfortheDendogram
}
shinyApp(ui, server)
