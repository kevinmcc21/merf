# Various functions to assist in plotting microbiome data in R

# Data reformatting----

alignQiimeData <- function(E) {
  o <- E$o
  s <- E$s
  if(!is.list(o)) stop("Object o is not of type list")
  if(!is.data.frame(s)) stop("Object s is not of type data frame")
  if(!"sample_ids" %in% names(o)) stop("Object o needs element named 'sample_ids'")
  if(!"otu_ids" %in% names(o)) stop("Object o needs element named 'otu_ids'")
  if(!"metadata" %in% names(o)) stop("Object o needs element named 'metadata'")
  if(!"counts" %in% names(o)) stop("Object o needs element named 'counts'")
  if(!"SampleID" %in% colnames(s)) stop("Object s needs column named 'SampleID'")
  if(is.null(colnames(o$counts))) stop("o$counts needs named columns")
  if(is.null(rownames(o$counts))) stop("o$counts needs named rows")
  if(is.null(names(o$metadata))) stop("o$metadata must be a named vector")
  
  s <- s[s$SampleID %in% colnames(o$counts),]
  s <- droplevels(s)
  
  o$counts <- o$counts[o$otu_ids,colnames(o$counts) %in% s$SampleID]
  o$counts <- o$counts[rowSums(o$counts)>0,]
  o$counts <- o$counts[,match(s$SampleID, colnames(o$counts))]
  o$sample_ids <- colnames(o$counts)
  o$metadata <- o$metadata[rownames(o$counts)]
  o$otu_ids <- o$otu_ids[o$otu_ids %in% rownames(o$counts)]
  E$s <- s
  E$o <- o
}

reorderDist <- function(dist, vector) {
  return(as.dist(as.matrix(dist)[match(vector,attr(dist,"Labels")),
                                 match(vector,attr(dist,"Labels"))]))
}

countsToProp <- function(counts) {
  return(apply(counts, 2, function(x)   ifelse (x, x / sum(x),0)))
}

fraction_nonzero <- function (x) sum(x > 0) / length(x)

tidy_lmer <- function(X) {
  sumry <- summary(X)
  coef <- sumry$coefficients
  df <- data.frame(term=rownames(coef),coef,group="fixed",row.names=NULL) %>%
    dplyr::rename(p.value=Pr...t..)
  return(df)
}

# Rank abundance curve plots----

plotRAC <- function(racData, graphLabel="", plotsMin=NULL, plotsMax=NULL, plotLegend=TRUE) {
  numLines <- length(racData)
  
  # get the range for the x and y axis 
  xrange <- c(1,max(sapply(racData,nrow)))
  if(is.null(plotsMin)) {
    plotsMin <- min(sapply(racData, function(x) min(x[,2])))
  }
  if(is.null(plotsMax)) {
    plotsMax <- max(sapply(racData, function(x) max(x[,2])))
  }
  yrange <- c(plotsMin,plotsMax)
  
  
  # set up the plot 
  par(mar=c(4,4,1,8))
  if(!plotLegend) par(mar=c(4,4,1,2))
  plot(xrange, yrange, type="n", xlab="OTU abundance rank",
       ylab="Proportion of sample", yaxt='n',log='y')
  logAxis(2,exponent=FALSE,las=1)
  
  colors <- rainbow(numLines)
  linetype <- c(1:numLines)
  if(!plotLegend) linetype <- rep(1,numLines)
  plotchar <- seq(25,25-numLines)
  
  # add lines
  for (i in 1:numLines) {
    lines(racData[[i]][,1], racData[[i]][,2], type="l", lwd=1.5,
          lty=linetype[i], col=colors[i], pch=plotchar[i])
  } 
  
  # add a title and subtitle 
  title(graphLabel)
  
  # add a legend 
  if(plotLegend) {
    legendX <- par('usr')[2]+0.01*diff(par('usr')[1:2])
    legendY <- 10^par('usr')[4]
    legend(legendX,
           legendY,
           names(racData), cex=0.8, col=colors,
           pch=plotchar, lty=linetype, 
           xjust=0,xpd=NA)
  }
}

plotRACset <- function(racDataDf, graphLabels, graphLabelDict, plotLegend=TRUE) {
  plotsMin <- min(as.numeric(lapply(racDataDf,function(x) min(x[,2]))))
  plotsMax <- max(as.numeric(lapply(racDataDf,function(x) max(x[,2]))))
  for(graphLabel in graphLabels) {
    if(sum(graphLabelDict[,2]==graphLabel)<3) {
      graphLabelDict[,2][graphLabelDict[,2]==graphLabel] <- "Other"
    }
  }
  graphLabels <- unique(graphLabelDict[,2])
  attach(mtcars)
  numGraphs <- length(graphLabels)
  numGraphRows <- max(1,floor(sqrt(numGraphs)))
  numGraphCols <- ceiling(numGraphs/numGraphRows)
  par(mfrow=c(numGraphRows,numGraphCols))
  for(graphLabel in sort(graphLabels)) {
    whichLists <- names(racDataDf)[names(racDataDf) %in% 
                                     graphLabelDict[graphLabelDict[,2]==graphLabel,1]]
    plotRAC(racDataDf[whichLists, drop=FALSE], graphLabel, plotsMin, plotsMax, plotLegend)
  }
}

# PCoA functions----

makeLegendColors <- function(vector) {
  legendColors <- rainbow(length(unique(vector)))
  names(legendColors) <- sort(unique(vector))
  return(legendColors)
}

makeLegendShapes <- function(vector) {
  if(length(unique(vector)) > 8) { stop("Can't make order legend for more than 8 unique values.")}
  vectorOrder <- as.matrix(unique(cbind(as.character(vector),as.numeric(factor(vector)))))
  orderShapeList <- c(15,17,18,19,21,22,23,24,25)
  newVectorOrder <- as.numeric(orderShapeList[as.numeric(vectorOrder[,2])])
  names(newVectorOrder) <- vectorOrder[,1]
  return(newVectorOrder)
}

printPcoaLegends <- function(colors=NULL, color_title="", shapes=NULL, shape_title="") {
  library(dnar)
  if(!is.null(colors)) {
    xpos=ifelse(is.null(shapes),mean(par('usr')[c(1,2)]),mean(par('usr')[c(1,1,2)]))
    legend(
      xpos,
      convertLineToUser(4,1),
      names(colors), pch=15,
      col=colors, xjust=0.5,xpd=NA, y.intersp = 0.75, title=color_title
    )
  }
  if(!is.null(shapes)) {
    xpos=ifelse(is.null(colors),mean(par('usr')[c(1,2)]),mean(par('usr')[c(1,2,2,2,2)]))
    legend(
      xpos,
      convertLineToUser(4,1),
      names(shapes), pch=as.numeric(shapes),
      col="black", xjust=0.5,xpd=NA, y.intersp = 0.75, title=shape_title
    )
  }
}

setPcoaPlotMargins <- function() {
  par(mar=c(12,4,3,1)) # need to leave margin room
}

plotPcoaGroup <- function(E, pcoa, joinVar, fillVar, shapeVar=21, title=NULL, legend=FALSE) {
  if(shapeVar==21) {
    df <- data.frame(Axis.1=pcoa$vectors[,1],Axis.2=pcoa$vectors[,2],
                     ID=rownames(pcoa$vectors)) %>%
      left_join((E$s %>% select(one_of(joinVar), one_of(fillVar))),by=c("ID"=joinVar))
    pcoaPlot <-
      ggplot(df, aes_string(x="Axis.1",y="Axis.2",fill=fillVar)) +
      geom_point(color="black", size=2.5, shape=21)
  } else {
    df <- data.frame(Axis.1=pcoa$vectors[,1],Axis.2=pcoa$vectors[,2],
                     ID=rownames(pcoa$vectors)) %>%
      left_join((E$s %>% select(one_of(joinVar), one_of(fillVar), one_of(shapeVar))),by=c("ID"=joinVar))
    df <- df %>%
      mutate(shapeCol=c("21","22","23","24","25")[as.integer(factor(df %>% select_(shapeVar) %>% pull(1)))])
    pcoaPlot <-
      ggplot(df, aes_string(x="Axis.1",y="Axis.2",fill=fillVar,shape="shapeCol")) +
      geom_point(color="black", size=2.5) +
      scale_shape_manual(values=c(21,22,23,24,25)) +
      guides(fill = guide_legend(override.aes = list(shape = 21)))
  }
  pcoaPlot <- pcoaPlot +
    coord_equal() +
    labs(
      x = sprintf("Axis 1 (%1.1f%%)", pcoa$values[1, 2]*100),
      y = sprintf("Axis 2 (%1.1f%%)", pcoa$values[2, 2]*100)
    ) +
    themePcoa +
    ggtitle(title)
  if(legend) {
    pcoaPlot <- pcoaPlot+theme(legend.position="right")
  }
  return(pcoaPlot)
}

plotUuPcoaGroup <- function(E, joinVar, fillVar, shapeVar=21, title=NULL) {
  plotPcoaGroup(E,E$uu_pcoa,joinVar, fillVar, shapeVar, title)
}

plotWuPcoaGroup <- function(E, joinVar, fillVar, shapeVar=21, title=NULL) {
  plotPcoaGroup(E,E$wu_pcoa, joinVar, fillVar, shapeVar, title)
}

plotJaccardPcoaGroup <- function(E, joinVar, fillVar, shapeVar=21, title=NULL) {
  plotPcoaGroup(E,E$jaccard_pcoa, joinVar, fillVar, shapeVar, title)
}

trioPcoaPlot <- function(E, joinVar="PatientID", fillVar, sampleText="all samples", shapeVar=21) {
  return(plotPcoaGroup(E,E$jaccard_pcoa, joinVar, fillVar, shapeVar=shapeVar, legend=TRUE))
  uuplot <- plotUuPcoaGroup(E, joinVar, fillVar, shapeVar, "Unweighted Unifrac")
  wuplot <- plotWuPcoaGroup(E, joinVar, fillVar, shapeVar, "Weighted Unifrac")
  jplot <- plotJaccardPcoaGroup(E, joinVar, fillVar, shapeVar, "Jaccard")

  cat("\n\n\\pagebreak\n")
  cat(paste0("\n\n# PCOA on ",sampleText,", colored by ",fillVar,"\n"))
  grid.arrange(uuplot,wuplot,jplot, g_legend(plotPcoaGroup(E,E$jaccard_pcoa, joinVar, fillVar, shapeVar=shapeVar, legend=TRUE)), ncol=2)
}

g_legend <- function(a.gplot){ 
  library(ggplot2)
  library(grid)
  library(gridExtra)
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
}


# Diversity functions----

shannonDiv <- function(data){
  H <- numeric(ncol(data))
  mult <- numeric(nrow(data))
  for(i in 1:ncol(data)){
    prop <- data[,i] / sum(data[,i])
    for(j in 1:nrow(data)){
      mult[j] <- prop[j] * log(prop[j])
    }
    H[i] <-- sum(mult, na.rm=T)
  }
  return(cbind(colnames(data),H))
}

evenness <- function(data){
  H <- numeric(ncol(data))
  mult <- numeric(nrow(data))
    for(i in 1:ncol(data)){
    prop <- data[,i] / sum(data[,i])
        for(j in 1:nrow(data)){
      mult[j] <- prop[j] * log(prop[j])
    }
    nonZero <- sum(prop > 0)
    HMax <- -log(1/nonZero)
    H[i] <- -sum(mult, na.rm=T) / HMax
  }
  return(cbind(colnames(data),H))
}

# Prediction trees----

prepareTreeData <- function(df, sampleIdVar, expVar, taxVar, propVar) {
  library(dplyr)
  library(lazyeval)
  treeData <- df %>%
    group_by_(eval(sampleIdVar),eval(expVar), eval(taxVar)) %>%
    summarize_(proportion=interp(~sum(var), var = as.name(propVar))) %>%
    ungroup() %>%
    unite_('idcol',c(sampleIdVar, expVar), sep=" ") %>%
    tidyr::complete_(cols=c(quo(idcol), eval(taxVar)), fill=list(proportion=0)) %>%
    separate_('idcol',c(sampleIdVar, expVar), sep=" ") %>%
    spread(eval(taxVar), proportion) %>%
    mutate_if(names(.)==expVar,as.factor)
  return(treeData)
}

makeDecisionTree <- function(df, sampleIdVar, expVar) {
  library(tree)
  tree <- tree(eval(formula(paste0(expVar, " ~ . - ", sampleIdVar))), 
               data = data.frame(df), model = T)
  return(tree)
}

makePredictionOnTree <- function(df, tree, sampleIdVar, expVar) {
  sp <- df %>%
    select_(sampleIdVar, expVar) %>%
    cbind(predict(tree)) %>%
    gather(Prediction, Probability, -starts_with(sampleIdVar), -starts_with(expVar))
  return(sp)
}

predictionAnalysis <- function(df, sampleIdVar, expVar, taxVar, propVar) {
  data <- prepareTreeData(df, sampleIdVar, expVar, taxVar, propVar)
  tree <- makeDecisionTree(data, sampleIdVar, expVar)
  sp <- makePredictionOnTree(data, tree, sampleIdVar, expVar)
  return(list(data, tree, sp))
}

plotPredictionAnalysis <- function(df) {
  sampleIdVar <- colnames(df)[1]
  expVar <- colnames(df)[2]
  df %>%
    mutate(Prediction = paste("Prediction =", Prediction)) %>%
    mutate(SubjectID = paste(!!names(.[2])," = ", UQ(sym(expVar)))) %>%
    ggplot() +
    geom_point(aes_string(x=sampleIdVar, y="Probability", color=expVar)) +
    facet_grid(as.formula(paste(expVar, "~ Prediction")), space = "free_y", scales="free_y") +
    coord_flip() +
    theme_bw()
}



# Other functions----


getNMostCommonX <- function(df,n,X="Otu",f=filter) {
  return(unique(df %>% f() %>% group_by_(X) %>% summarize(numSamples=n()) %>%
                  top_n(n,numSamples) %>% pull(eval(X))))
}

getNMostAbundantX <- function(df,n,X="Otu",f=filter) {
  return(unique(df %>% f() %>% group_by(SampleID) %>% top_n(n,proportion) %>% pull(eval(X))))
}

getMostAbundantTax <- function(df,n,f=filter) {
  return(unique(df %>% f() %>% group_by(SampleID) %>% top_n(n,proportion) %>% pull(Tax)))
}

alignUnifracData <- function(s, u) {
  return(dist_subset(u,match(s$SampleID,attr(u,"Labels"))))
}

