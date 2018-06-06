

# Sub function of cd function below
ss <- function(d) {
  sum(d[lower.tri(d)]^2)
}

# Calculate distance between centroids of two sub-matrices, denoted with labels xlabs and ylabs
cd <- function(m,xlabs,ylabs) {
  x <- as.matrix(m[xlabs,xlabs])
  n.x <- nrow(x)
  y <- as.matrix(m[ylabs,ylabs])
  n.y <- nrow(y)
  xy <- m[c(xlabs,ylabs),c(xlabs,ylabs)]
  return(sqrt( (ss(xy) - (n.x + n.y) * (ss(x)/n.x + ss(y)/n.y)) / (n.x*n.y)))
}

# Calculate distance between matrix centroids for two species.
# Assumes G$s$common has global list of species.
getCdCommon <- function(m,common1,common2) {
  return(cd(m,
            rownames(m)[rownames(m) %in% G$s$SampleID[G$s$common==common1]],
            rownames(m)[rownames(m) %in% G$s$SampleID[G$s$common==common2]]))
}

# Convert distance object to small distance object using centroids of groups
# Vector denotes grouping to use
# Order denotes order for re-ordering matrix before converting to dist
centroidDist <- function(dist, vector, order) {
  distMatrix <- as.matrix(dist)
  newDistMatrix <- matrix(nrow=length(unique(vector)), ncol=length(unique(vector)))
  rownames(newDistMatrix) <- unique(vector)
  colnames(newDistMatrix) <- unique(vector)
  for(i in rownames(newDistMatrix)) {
    for(j in colnames(newDistMatrix)) {
      if(i==j) {
        newDistMatrix[i,j] <- 0
      } else {
        newDistMatrix[i,j] <- getCdCommon(distMatrix,i,j)
      }
    }
  }
  newDistMatrix <-
    newDistMatrix[match(G$speciesTax$common[order],rownames(newDistMatrix)),
                  match(G$speciesTax$common[order],colnames(newDistMatrix))]
  return(as.dist(newDistMatrix))
}

meltDist <- function(dist) {
  tempDf <- as.data.frame(as.matrix(dist))
  tempDf$names <- rownames(tempDf)
  newDf <- melt(tempDf,id.vars = "names")
  newDf$names <- factor(newDf$names,levels=unique(newDf$variable))
  newDf <- newDf[as.numeric(newDf$names) < as.numeric(newDf$variable),]
  return(newDf)
}


split_df <- function(DF, nrows = ceiling(nrow(DF)/3)){
  DF %>%
    rename(xreadCount = readCount) %>%
    mutate(id = rep(1:nrow(.), each = nrows, len = nrow(.))) %>%
    gather(variable, value, -id) %>%
    unite(temp, id, variable) %>%
    group_by(temp) %>%
    mutate(id = 1:n()) %>%
    spread(temp, value) %>%
    select(-id) %>%
    data.frame() %>%
    setNames(rep(names(DF), ceiling(nrow(DF)/nrows))) 
}

string_average <- function(v) {
  return(rank(v) < length(v) / 2)
}


string_pos_cmp <- function(x,y) {
  imax <- min(nchar(x),nchar(y))
  result <- ""
  for(i in seq(1,imax)){
    result <- paste0(result,c(0,1)[(substr(x,i,i) == substr(y,i,i))+1])
  }
  return(result)
}

jaccard_min_max <- function(m) {
  apply(m,1,function(X) {
    apply(m,1,function(Y) {
      return(1-sum(mapply(min,X,Y)) / sum(mapply(max,X,Y)))
    })
  })
}


'%!in%' <- function(x,y)!('%in%'(x,y))

writeTable <- function(x, fileVar="",sepVar=" ",append=FALSE,...) {
  write.table(sepVar,fileVar,append=append,sep="",eol="",row.names = FALSE,col.names = FALSE,quote = FALSE)
  write.table(x,fileVar,sep=sepVar,append=TRUE, quote=FALSE,...)
}


Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


rearrangeCycle <- function(vector) {
  n <- length(vector)
  cycleLength <- floor(sqrt(n))
  newVector <- vector(length=n,mode=class(vector))
  assignedItems <- vector(length=n,mode="numeric")
  nextSlot <- 1
  for(i in seq(1,n)) {
    newVector[nextSlot] <- vector[i]
    assignedItems[nextSlot] <- 1
    nextSlot <- ifelse(nextSlot + cycleLength > n, nextSlot + cycleLength - n, nextSlot + cycleLength)
    while(sum(assignedItems) < n & assignedItems[nextSlot]==1) {
      nextSlot <- ifelse(nextSlot == n, nextSlot + 1, min(which(assignedItems==0)))
    }
  }
  return(newVector)
}

scale_fill_brewer_cycle <- function(limits, ...) {
  cycledLimits <- rearrangeCycle(limits)
  return(scale_fill_brewer(..., limits=cycledLimits, breaks=limits))
}








