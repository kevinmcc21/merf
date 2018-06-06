# Functions for modifying and plotting phylogenetic trees in R

sizeOfSubTree <- function(tree,nodeLabel) {
  if(nodeLabel %in% tree$edge[,1]) {
    return(sum(sapply(tree$edge[tree$edge[,1]==nodeLabel,2],
                      function(X) sizeOfSubTree(tree, X))))
  } else {
    return(1)
  }
}

expandTree <- function(tree,replacingLabel,newLabel,newEdgeLength){
  replacingTip <- which(tree$tip.label==replacingLabel)
  replacingEdge <- which(tree$edge[,2]==replacingTip)
  
  numNodes <- max(tree$edge)
  newNode <- numNodes + 2
  newTip <- replacingTip + 1
  numTips <- length(tree$tip.label)
  
  tree$edge[tree$edge > replacingTip] <- tree$edge[tree$edge > replacingTip] + 1
  tree$edge[tree$edge==replacingTip] <- newNode
  tree$edge <- cbind(append(append(tree$edge[,1],newNode,replacingEdge),newNode,replacingEdge),
                     append(append(tree$edge[,2],newTip,replacingEdge),replacingTip,replacingEdge))
  
  tree$edge.length[replacingEdge] <- tree$edge.length[replacingEdge] - newEdgeLength
  tree$edge.length <- append(append(tree$edge.length,newEdgeLength,replacingEdge),newEdgeLength,replacingEdge)
  
  tree$tip.label <- append(tree$tip.label,newLabel,replacingTip)
  tree$Nnode <- tree$Nnode + 1
  tree$node.label <- c(tree$node.label, paste0("'",max(as.numeric(gsub("'","",tree$node.label[-1])))+1,"'"))
  
  return(tree)
}

expandTreeInternal <- function(tree,replacingNode,newLabel,newTipEdgeLength){
  replacingEdge <- which(tree$edge[,2]==replacingNode)
  replacingNodeLeaves <- as(subset(as(tree,"phylo4"),node.subtree=replacingNode),"phylo")$tip.label
  newTip <- max(match(replacingNodeLeaves,tree$tip.label)) + 1
  
  numNodes <- max(tree$edge)
  newNode <- numNodes + 2
  numTips <- length(tree$tip.label)
  replacingNodeHeight <- max(node.depth.edgelength(as(subset(as(tree,"phylo4"),
                                                             node.subtree=replacingNode),"phylo")))
  newInternalEdgeLength <- newTipEdgeLength - replacingNodeHeight
  newReplacingEdgeLength <- tree$edge.length[replacingEdge] - newInternalEdgeLength
  
  
  tree$edge[tree$edge >= newTip] <- tree$edge[tree$edge >= newTip] + 1
  replacingNode <- replacingNode + 1
  tree$edge[tree$edge[,1]==replacingNode,1] <- newNode
  tree$edge <- cbind(append(append(tree$edge[,1],replacingNode,replacingEdge),replacingNode,replacingEdge),
                     append(append(tree$edge[,2],newNode,replacingEdge),newTip,replacingEdge))
  
  tree$edge.length <- append(append(tree$edge.length,newInternalEdgeLength,replacingEdge),
                             newTipEdgeLength,replacingEdge)
  tree$edge.length[replacingEdge] <- newReplacingEdgeLength
  
  tree$tip.label <- append(tree$tip.label,newLabel,newTip - 1)
  tree$Nnode <- tree$Nnode + 1
  tree$node.label <- c(tree$node.label, paste0("'",max(as.numeric(gsub("'","",tree$node.label[-1])))+1,"'"))
  
  return(tree)
}

getNodeDepth <- function(tree, node) {
  if(node %in% tree$edge[,2]) {
    return(1+getNodeDepth(tree,tree$edge[tree$edge[,2]==node,1]))
  } else {
    return(1)
  }
}

getNodeHeight <- function(tree, node) {
  if(node %in% tree$edge[,2]) {
    edgeIndex <- tree$edge[,2]==node
    return(tree$edge.length[edgeIndex]+getNodeHeight(tree,tree$edge[edgeIndex,1]))
  } else {
    return(0)
  }
}

maxEdgeLength <- function(tree, node, prevMax = 0, root=NULL, max = 1) {
  if(is.null(root)) {
    root <- length(tree$tip.label) + 1
  }
  if(node == root) {
    return(prevMax)
  } else {
    edgeIndex <- tree$edge[,2]==node
    thisMax <- max(prevMax, tree$edge.length[edgeIndex])
    if(thisMax > max) {
      return(thisMax)
    } else {
      return(maxEdgeLength(tree, tree$edge[edgeIndex,1], thisMax, root, max))
    }
  }
}

getEdgeLeaves <- function(tree, node, root=NULL) {
  if(is.null(root)) {
    root <- length(tree$tip.label) + 1
  }
  if(node < root) return(node)
  return(unlist(sapply(tree$edge[tree$edge[,1]==node,2],function(X) getEdgeLeaves(tree,X,root=root))))
}

# custom function to replace phyloseq::plot_tree until it is updated on bioconductor
# ape::node_depth_edgelength is no longer compatible for external packages
# Replacing with node.depth.edgelength

dummy_tree_layout <- function (phy, ladderize = FALSE) 
{
  if (inherits(phy, "phyloseq")) {
    phy = phy_tree(phy)
  }
  if (!inherits(phy, "phylo")) {
    stop("tree missing or invalid. Please check `phy` argument and try again.")
  }
  if (is.null(phy$edge.length)) {
    phy$edge.length <- rep(1L, times = nrow(phy$edge))
  }
  if (ladderize != FALSE) {
    if (ladderize == "left") {
      phy <- ladderize(phy, FALSE)
    }
    else if (ladderize == TRUE | ladderize == "right") {
      phy <- ladderize(phy, TRUE)
    }
    else {
      stop("You did not specify a supported option for argument `ladderize`.")
    }
  }
  z = reorder.phylo(phy, order = "postorder")
  Nedge = nrow(phy$edge)[1]
  Nnode = phy$Nnode
  Ntip = length(phy$tip.label)
  ROOT = Ntip + 1
  TIPS = phy$edge[(phy$edge[, 2] <= Ntip), 2]
  NODES = (ROOT):(Ntip + Nnode)
  nodelabels = phy$node.label
  xx = ape::node.depth.edgelength(z)
  yy = ape::node.height(z)
  #xx = dummy_ape_node_depth_edge_length(Ntip, Nnode, z$edge, Nedge, z$edge.length)
  #yy <- numeric(Ntip + Nnode)
  #yy[TIPS] <- 1:Ntip
  #ape_node_height <- function(Ntip, Nnode, edge, Nedge, yy) {
  #  .C(ape:::node_height, PACKAGE = "ape", as.integer(Ntip), 
  #     as.integer(Nnode), as.integer(edge[, 1]), as.integer(edge[, 
  #                                                               2]), as.integer(Nedge), as.double(yy))[[6]]
  #}
  #yy <- ape_node_height(Ntip, Nnode, z$edge, Nedge, yy)
  edgeDT = data.table(phy$edge, edge.length = phy$edge.length, 
                      OTU = NA_character_)
  if (!is.null(phy$tip.label)) {
    edgeDT[, `:=`(OTU, NA_character_)]
    setkey(edgeDT, V2)
    edgeDT[V2 <= Ntip, `:=`(OTU, phy$tip.label)]
  }
  edgeDT[, `:=`(xleft, xx[V1])]
  edgeDT[, `:=`(xright, xx[V2])]
  edgeDT[, `:=`(y, yy[V2])]
  vertDT = edgeDT[, list(x = xleft[1], vmin = min(y), vmax = max(y)), 
                  by = V1, mult = "last"]
  if (!is.null(phy$node.label)) {
    edgeDT[V2 > ROOT, `:=`(x, xright)]
    edgeDT[V2 > ROOT, `:=`(label, phy$node.label[-1])]
    setkey(vertDT, V1)
    vertDT[J(ROOT), `:=`(y, mean(c(vmin, vmax)))]
    vertDT[J(ROOT), `:=`(label, phy$node.label[1])]
  }
  return(list(edgeDT = edgeDT, vertDT = vertDT))
}


dummy_plot_tree <- function (physeq, method = "sampledodge", nodelabf = NULL, color = NULL, 
                        shape = NULL, size = NULL, min.abundance = Inf, label.tips = NULL, 
                        text.size = NULL, sizebase = 5, base.spacing = 0.02, ladderize = FALSE, 
                        plot.margin = 0.2, title = NULL, treetheme = NULL, justify = "jagged") 
{
  fix_reserved_vars = function(aesvar) {
    aesvar <- gsub("^abundance[s]{0,}$", "Abundance", aesvar, 
                   ignore.case = TRUE)
    aesvar <- gsub("^OTU[s]{0,}$", "OTU", aesvar, ignore.case = TRUE)
    aesvar <- gsub("^taxa_name[s]{0,}$", "OTU", aesvar, ignore.case = TRUE)
    aesvar <- gsub("^sample[s]{0,}$", "Sample", aesvar, ignore.case = TRUE)
    return(aesvar)
  }
  if (!is.null(label.tips)) {
    label.tips <- fix_reserved_vars(label.tips)
  }
  if (!is.null(color)) {
    color <- fix_reserved_vars(color)
  }
  if (!is.null(shape)) {
    shape <- fix_reserved_vars(shape)
  }
  if (!is.null(size)) {
    size <- fix_reserved_vars(size)
  }
  if (is.null(phy_tree(physeq, FALSE))) {
    stop("There is no phylogenetic tree in the object you have provided.\n", 
         "Try phy_tree(physeq) to see for yourself.")
  }
  if (!inherits(physeq, "phyloseq")) {
    method <- "treeonly"
  }
  treeSegs <- dummy_tree_layout(phy_tree(physeq), ladderize = ladderize)
  edgeMap = aes(x = xleft, xend = xright, y = y, yend = y)
  vertMap = aes(x = x, xend = x, y = vmin, yend = vmax)
  p = ggplot(data = treeSegs$edgeDT) + geom_segment(edgeMap) + 
    geom_segment(vertMap, data = treeSegs$vertDT)
  if (is.null(text.size)) {
    text.size <- manytextsize(ntaxa(physeq))
  }
  if (!is.null(label.tips) & method != "sampledodge") {
    labelDT = treeSegs$edgeDT[!is.na(OTU), ]
    if (!is.null(tax_table(object = physeq, errorIfNULL = FALSE))) {
      taxDT = data.table(tax_table(physeq), OTU = taxa_names(physeq), 
                         key = "OTU")
      labelDT = merge(x = labelDT, y = taxDT, by = "OTU")
    }
    if (justify == "jagged") {
      labelMap <- aes_string(x = "xright", y = "y", label = label.tips, 
                             color = color)
    }
    else {
      labelMap <- aes_string(x = "max(xright, na.rm=TRUE)", 
                             y = "y", label = label.tips, color = color)
    }
    p <- p + geom_text(labelMap, data = labelDT, size = I(text.size), 
                       hjust = -0.1, na.rm = TRUE)
  }
  if (is.null(nodelabf)) {
    nodelabf = howtolabnodes(physeq)
  }
  p = nodelabf(p, treeSegs$edgeDT[!is.na(label), ])
  p = nodelabf(p, treeSegs$vertDT[!is.na(label), ])
  if (is.null(treetheme)) {
    treetheme <- theme(axis.ticks = element_blank(), axis.title.x = element_blank(), 
                       axis.text.x = element_blank(), axis.title.y = element_blank(), 
                       axis.text.y = element_blank(), panel.background = element_blank(), 
                       panel.grid.minor = element_blank(), panel.grid.major = element_blank())
  }
  if (inherits(treetheme, "theme")) {
    p <- p + treetheme
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  if (method != "sampledodge") {
    return(p)
  }
  dodgeDT = treeSegs$edgeDT[!is.na(OTU), ]
  dodgeDT = merge(x = dodgeDT, y = data.table(psmelt(physeq), 
                                              key = "OTU"), by = "OTU")
  if (justify == "jagged") {
    dodgeDT <- dodgeDT[Abundance > 0, ]
  }
  if (!is.null(color) | !is.null(shape) | !is.null(size)) {
    setkeyv(dodgeDT, cols = c("OTU", color, shape, size))
  }
  else {
    setkey(dodgeDT, OTU, Sample)
  }
  dodgeDT[, `:=`(h.adj.index, 1:length(xright)), by = OTU]
  if (justify == "jagged") {
    dodgeDT[, `:=`(xdodge, (xright + h.adj.index * base.spacing * 
                              max(xright, na.rm = TRUE)))]
  }
  else {
    dodgeDT[, `:=`(xdodge, max(xright, na.rm = TRUE) + h.adj.index * 
                     base.spacing * max(xright, na.rm = TRUE))]
    dodgeDT <- dodgeDT[Abundance > 0, ]
  }
  dodgeMap <- aes_string(x = "xdodge", y = "y", color = color, 
                         fill = color, shape = shape, size = size)
  p <- p + geom_point(dodgeMap, data = dodgeDT, na.rm = TRUE)
  if (!is.null(size)) {
    p <- p + scale_size_continuous(trans = log_trans(sizebase))
  }
  if (any(dodgeDT$Abundance >= min.abundance[1])) {
    pointlabdf = dodgeDT[Abundance >= min.abundance[1], ]
    p <- p + geom_text(mapping = aes(xdodge, y, label = Abundance), 
                       data = pointlabdf, size = text.size, na.rm = TRUE)
  }
  if (!is.null(label.tips)) {
    tiplabDT = dodgeDT
    tiplabDT[, `:=`(xfartiplab, max(xdodge)), by = OTU]
    tiplabDT <- tiplabDT[h.adj.index == 1, .SD, by = OTU]
    if (!is.null(color)) {
      if (color %in% sample_variables(physeq, errorIfNULL = FALSE)) {
        color <- NULL
      }
    }
    labelMap <- NULL
    if (justify == "jagged") {
      labelMap <- aes_string(x = "xfartiplab", y = "y", 
                             label = label.tips, color = color)
    }
    else {
      labelMap <- aes_string(x = "max(xfartiplab, na.rm=TRUE)", 
                             y = "y", label = label.tips, color = color)
    }
    p <- p + geom_text(labelMap, tiplabDT, size = I(text.size), 
                       hjust = -0.1, na.rm = TRUE)
  }
  min.x <- -0.01
  max.x <- dodgeDT[, max(xright, na.rm = TRUE)]
  if ("xdodge" %in% names(dodgeDT)) {
    max.x <- dodgeDT[, max(xright, xdodge, na.rm = TRUE)]
  }
  if (plot.margin > 0) {
    max.x <- max.x * (1 + plot.margin)
  }
  p <- p + scale_x_continuous(limits = c(min.x, max.x))
  return(p)
}

manytextsize <- function(n, mins=0.5, maxs=4, B=6, D=100){
  # empirically selected size-value calculator.
  s <- B * exp(-n/D)
  # enforce a floor.
  s <- ifelse(s > mins, s, mins)
  # enforce a max
  s <- ifelse(s < maxs, s, maxs)
  return(s)
}