

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


# Convert list of character vectors to single character matrix, back-filling NA for sequences of uneven length
charListToMatrix <- function(L) {
  maxLen <- max(sapply(L, function(X) {length(X)}))
  L2 <- lapply(L, function(X) {
    v <- vector(mode="character",length=maxLen)
    v[1:maxLen] <- NA
    v[1:length(X)] <- X
    return(v)
  })
  m <- do.call(rbind, L2)
  return(m)
}

charListToMatrix2 <- function(L) {
  minLen <- min(sapply(L, function(X) {length(X)}))
  L2 <- lapply(L, function(X) {
    return(X[1:minLen])
  })
  m <- do.call(rbind, L2)
  return(m)
}

expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))

vectorToDf <- function(X, colName1, colName2) {
  df <- data.frame(col1=names(X),
                   col2=X)
  colnames(df) <- c(colName1,colName2)
  rownames(df) <- NULL
  return(df)
}

# Tries to load a package silently and installs it if it needs to
get_package <- function(pkg, load = TRUE, silent = FALSE, repos = "http://cran.us.r-project.org") {
  if(!suppressMessages(suppressWarnings(require(pkg, character.only = TRUE, quietly = TRUE)))) {
    try(install.packages(pkg, repos = repos), silent = TRUE)
  }
  if(load) suppressPackageStartupMessages(library(pkg, character.only = TRUE, quietly = TRUE))
  if(load & !silent) message("Loaded ", pkg)
}

# Simultaneously converts relative path to absolute path and checks if file exists
try_filepath <- function(file, desc="unspecified") {
  tryCatch(newfile <- tools::file_path_as_absolute(file),
           error = function(e) stop(paste0("Cannot find ", desc," file: ", file)),
           finally = 1)
  print(paste0("Using ", desc, " file: ", newfile))
  return(newfile)
}

# Simultaneously converts relative path to absolute path and checks if containing directory exists
try_dirpath <- function(file, desc="output") {
  tryCatch(newfile <- paste0(tools::file_path_as_absolute(dirname(file)),"/", basename(file)),
           error = function(e) stop(paste0("Cannot find ", desc, "directory: ", dirname(file))),
           finally = 1)
  print(paste0("Using ", desc, " file: ", newfile))
  return(newfile)
}

try_element <- function(list, name) {
  if(name %in% names(list)) {
    return(list[[name]])
  } else {
    return(NULL)
  }
}

# Render a PDF file using an input RMD script and a target file
render_output <- function(script, pdf, params=NULL) {
  rmarkdown::render(
    input = script,
    output_file = pdf,
    quiet = TRUE,
    clean = TRUE,
    params=params
  )
}

# Reduce a string to a shorter, regex-style string
reduce_string <- function(string) {
  print(paste0("String:",string))
  N <- nchar(string)
  if(N == 0) { return(NULL) }
  if(N <= 3) { return(string) }
  
  get_b_candidates <- function(N) {
    b_candidates <- sapply(seq(1, floor(sqrt(N))), function(X) {
      if(N %% X == 0) { return(c(X, N / X)) } else { return(1) }
    })
    return(sort(unique(unlist(as.vector(b_candidates)))))
  }

  result <- sapply(seq(1, N), function(i) {
    LHS <- substr(string, 1, i)
    N_LHS <- nchar(LHS)
    if(N_LHS <= 3) {
      resultLHS <- LHS
    } else {
      resultLHS <- sapply(get_b_candidates(N_LHS), function(b) {
        snippet <- substr(LHS, 1, N_LHS / b)
        if(paste0(rep(snippet, b), collapse = '') == LHS) {
          return(paste0(snippet, "^", b))
        }
      })
      resultLHS <- unlist(resultLHS)
      resultLHS <- resultLHS[which.min(nchar(resultLHS))]
    }
    
    if(i < N - 3) {
      RHS <- substr(string, i + 1, N)
      resultRHS <- reduce_string(RHS)
      #print(paste0("Returning pair (LHS:", resultLHS, ", RHS:", resultRHS, ")"))
      return(paste0(resultLHS, resultRHS))
    } else if(i < N) {
      resultRHS <- substr(string, i + 1, N)
      #print(paste0("Returning pair (LHS:", resultLHS, ", RHS:", resultRHS, ")"))
      return(paste0(resultLHS, resultRHS))
    } else {
      #print(paste0("Returning ", resultLHS))
      return(resultLHS)
    }
  })
  
  result <- unlist(result)
  result <- result[which.min(nchar(result))]
  
  print(paste0("Returning result ", result, " for string ", string))
  
  return(result)
}

