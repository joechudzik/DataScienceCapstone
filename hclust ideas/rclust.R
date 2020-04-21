library(data.table)

# Calculate the cost of each given path
cost <- function(path, pdist=pdist, ptime=ptime, weight, data) {
  cost <- NULL
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
  if( max(est, na.rm=T) > min(lst, na.rm=T) ){
  #if (sum(est < rev(cummin(rev(lst))), na.rm=T) < length(est)){
    return(100000000)
  }
  else{
    st <- min(lst, max(est), na.rm=T)
    trip_time <- st + dtime
    wait_time <- routeTime$et - trip_time
    for(i in 1:length(wait_time)){ if( is.na(wait_time[i]) ) wait_time[i] <- 0}
    trip_time <- trip_time + wait_time
    new_time <- wait_time + dtime
    for (i in 1:(length(path)-1)) cost[i] <- (1-weight)*pdist[path[i],path[i+1]] + weight*trip_time[i]
    return(sum(cost))
  }
}

# Searching among possible routes by variety of combination of "a" and "b" and pick the optimal one.
rcost <- function(a,b,pdist=pdist,ptime=ptime,weight,data) {
  rcosts <- NULL
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
  for (i in 1:length(routs)) rcosts[i] <- cost(routs[[i]], pdist, ptime, weight, data)
  # May need to add some logic to remove the very high or NA value from rcosts
  indx <- which.min(rcosts)
  return(list(rcost=rcosts[indx], route=routs[[indx]]))
}

# Calculate all of the pairwise costs needed to start the algorithm
pcosts <- function(pdist=pdist,ptime=ptime,weight, data) {
  N <- dim(pdist)[1]/2
  pcosts <- list(rcosts=matrix(NA,nr=N,nc=N), routes=array(list(),dim=c(N,N)))
  for (i in 1:(N-1)) for (j in (i+1):N) {
    pcost.ij <- rcost(i,j,pdist,ptime,weight,data)
    pcosts$rcosts[i,j] <- pcosts$rcosts[j,i] <- pcost.ij$rcost
    pcosts$routes[i,j] <- pcosts$routes[j,i] <- list(pcost.ij$route)
  }
  colnames(pcosts$rcosts) <- rownames(pcosts$rcosts) <- colnames(pcosts$routes) <- rownames(pcosts$routes) <- paste0("-",1:N)
  return(pcosts)
}

# Route hclust function
rhclust <- function(pdist,ptime,weight, data) {
  temp <- list(pcosts(pdist,ptime,weight,data)); N <- dim(pdist)[1]/2;
  rhclust <- list(merge=matrix(0,nr=N-1,nc=2), merge.route=list())
  for (i in 1:(N-2)) {
    indx <- rownames(which(temp[[i]]$rcosts==min(temp[[i]]$rcosts, na.rm = T), arr.ind = T))[1:2]
    print(as.numeric(indx))
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
      cost.j$rcosts[j] <- tmp$rcost; cost.j$routes[[j]] <- tmp$route
    }
    temp[[(i+1)]]$rcosts <- cbind(rbind(temp[[(i+1)]]$rcosts,cost.j$rcosts[-length(cost.j$rcosts)]),cost.j$rcosts)
    temp[[(i+1)]]$routes <- cbind(rbind(temp[[(i+1)]]$routes,cost.j$routes[-length(cost.j$routes)]),cost.j$routes)
    tmp <- colnames(temp[[(i+1)]]$rcosts); tmp[length(tmp)] <- i
    colnames(temp[[(i+1)]]$rcosts) <- rownames(temp[[(i+1)]]$rcosts) <- colnames(temp[[(i+1)]]$routes) <- rownames(temp[[(i+1)]]$routes) <- tmp
  }
  rhclust$merge[(i+1),] <- as.numeric(tmp); 
  rhclust$merge.route[(i+1)] <- temp[[(i+1)]]$routes[1,2]; 
  rhclust$height[(i+1)] <- temp[[(i+1)]]$rcosts[1,2]
  rhclust$labels <- NULL; rhclust$method <- "pool"; rhclust$dist.method <- "pdist"; rhclust$call <- sys.call()
  rhclust$order <- 1:N; class(rhclust) <- "hclust"
  rhclust$order <- unlist(as.dendrogram(rhclust))
  return(rhclust)
}

temp <- rhclust(pdis, ptim, 0.5, data)
cutree(temp,3)
plot(temp)
rect.hclust(temp, k=3, border=2:6)


data <- read.csv('~/Documents/GitHub/DataScienceCapstone/Data/dataForSpecificTimeAlgorithm.csv')
data <- data[1:10,]
pdis <- read.csv('https://raw.githubusercontent.com/joechudzik/DataScienceCapstone/master/Data/pdist.csv')
row.names(pdis) <- pdis[,1]; pdis <- pdis[,-1]
ptim <- read.csv('https://raw.githubusercontent.com/joechudzik/DataScienceCapstone/master/Data/ptime.csv')
row.names(ptim) <- ptim$X; ptim <- ptim[,-1]

test_path <- c('p1','p2','d1','d2')
test_weight = 0.5

