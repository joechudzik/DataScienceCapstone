# Calculate the cost of each given path
cost <- function(path,pdist=pdist) {
  cost <- NULL
  for (i in 1:(length(path)-1)) cost[i] <- pdist[path[i],path[i+1]]
  return(sum(cost))
}

# Searching among possible routes by variety of combination of "a" and "b" and pick the optimal one.
rcost <- function(a,b,pdist=pdist) {
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
  for (i in 1:length(routs)) rcosts[i] <- cost(routs[[i]], pdist)
  indx <- which.min(rcosts)
  return(list(rcost=rcosts[indx], route=routs[[indx]]))
}

# Calculate all of the pairwise costs needed to start the algorithm
pcosts <- function(pdist=pdist) {
  N <- dim(pdist)[1]/2
  pcosts <- list(rcosts=matrix(NA,nr=N,nc=N), routes=array(list(),dim=c(N,N)))
  for (i in 1:(N-1)) for (j in (i+1):N) {
    pcost.ij <- rcost(i,j,pdist)
    pcosts$rcosts[i,j] <- pcosts$rcosts[j,i] <- pcost.ij$rcost
    pcosts$routes[i,j] <- pcosts$routes[j,i] <- list(pcost.ij$route)
  }
  colnames(pcosts$rcosts) <- rownames(pcosts$rcosts) <- colnames(pcosts$routes) <- rownames(pcosts$routes) <- paste0("-",1:N)
  return(pcosts)
}

# Route hclust function
rhclust <- function(pdist) {
  temp <- list(pcosts(pdist)); N <- dim(pdist)[1]/2;
  rhclust <- list(merge=matrix(0,nr=N-1,nc=2), merge.route=list())
  for (i in 1:(N-2)) {
    indx <- rownames(which(temp[[i]]$rcosts==min(temp[[i]]$rcosts, na.rm = T), arr.ind = T))
    rhclust$merge[i,] <- as.numeric(indx); rhclust$merge.route[i] <- temp[[i]]$routes[indx[1],indx[2]]; 
    rhclust$height[i] <- temp[[i]]$rcosts[indx[1],indx[2]]
    ind <- !colnames(temp[[i]]$rcosts)%in%indx
    temp[[(i+1)]] <- list(rcosts=as.matrix(temp[[i]]$rcosts[ind,ind]), routes=as.matrix(temp[[i]]$routes[ind,ind]))
    inds <- as.list(colnames(temp[[i]]$rcosts)[ind])
    if (sum(ind)==1) colnames(temp[[(i+1)]][[1]]) <- rownames(temp[[(i+1)]][[1]]) <- colnames(temp[[(i+1)]][[2]]) <- rownames(temp[[(i+1)]][[2]]) <- inds[[1]]
    cost.j <- list(rcosts=array(NA,length(inds)+1), routes=array(list(),length(inds)+1))
    for (j in 1:length(inds)) {
      if (sign(as.numeric(inds[[j]]))==-1) inds[[j]] <- -as.numeric(inds[[j]]) else inds[[j]] <- rhclust$merge.route[[as.numeric(inds[[j]])]]
      tmp <- rcost(inds[[j]], rhclust$merge.route[[i]], pdist)
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

# Example
pdis <- read.csv('https://raw.githubusercontent.com/joechudzik/DataScienceCapstone/master/Data/pdist.csv')
row.names(pdis) <- pdis[,1]; pdis <- pdis[,-1]

rclust <- rhclust(pdis)
cutree(rclust,3)
rclust$merge.route[5:7]
rclust$height[5:7]

plot(rclust)
rect.hclust(rclust , k = 3, border = 2:6)

library(dendextend)
rclust.dend <- as.dendrogram(rclust)
rclust.col <- color_branches(rclust.dend, k=3)
plot(rclust.col)
