
# Making phylogenetic centroids and host range

library(centiserve)

VirPhyloHostRangeMax <- sapply(VirusAssocs, function(a){
  
  if(length(a[a%in%rownames(CytBMatrix)])>0){
    
    b <- a[a%in%rownames(CytBMatrix)]
    
    if(length(b)>1){
      max(CytBMatrix[b, b][upper.tri(CytBMatrix[b, b])])
    } else 0
  } else NA
  
})


VirPhyloHostRangeMean <- sapply(VirusAssocs, function(a){
  
  if(length(a[a%in%rownames(CytBMatrix)])>0){
    
    b <- a[a%in%rownames(CytBMatrix)]
    
    if(length(b)>1){
      mean(CytBMatrix[b, b][upper.tri(CytBMatrix[b, b])])
    } else 0
  } else NA
  
})


VirPhyloHostRangeMedian <- sapply(VirusAssocs, function(a){
  
  if(length(a[a%in%rownames(CytBMatrix)])>0){
    
    b <- a[a%in%rownames(CytBMatrix)]
    
    if(length(b)>1){
      median(CytBMatrix[b, b][upper.tri(CytBMatrix[b, b])])
    } else 0
  } else NA
  
})

Viruses[,c("HostRangeMax",
           "HostRangeMean",
           "HostRangeMedian")] <- cbind(VirPhyloHostRangeMax[Viruses$Sp],
                                        VirPhyloHostRangeMean[Viruses$Sp],
                                        VirPhyloHostRangeMedian[Viruses$Sp])

M2 <- as.matrix(CytBMatrix)

M2[upper.tri(M2)]

PhyloGraph <- graph.adjacency(1-M2, weighted = T, mode = "undirected", diag = F)

SubGraphList <- CentroidList <- list()

for(x in unique(Viruses$Sp)){
  
  SubGraph <- induced_subgraph(PhyloGraph, V(PhyloGraph)$name%in%VirusAssocs[[x]])
  SubGraphList[[x]] <- SubGraph
  
  if(class(try(centroid(SubGraph), silent = T)) != "try-error"){
    SubCentroids <- centroid(SubGraph)
    Centroids <- names(SubCentroids[SubCentroids==min(SubCentroids)])
    CentroidList[[x]] <- Centroids
  } else CentroidList[[x]] <- NA
}

for(x in unique(Viruses$Sp)[which(unique(Viruses$Sp)==x):length(unique(Viruses$Sp))]){
  
  SubGraph <- induced_subgraph(PhyloGraph, V(PhyloGraph)$name%in%VirusAssocs[[x]])
  SubGraphList[[x]] <- SubGraph
  
  if(class(try(centroid(SubGraph), silent = T)) != "try-error"){
    SubCentroids <- centroid(SubGraph)
    Centroids <- names(SubCentroids[SubCentroids==min(SubCentroids)])
    CentroidList[[x]] <- Centroids
  } else CentroidList[[x]] <- NA
  print(x)
}




distances(PhyloGraph, v = V(PhyloGraph)$name%in%VirusAssocs[[x]],
          to = SubCentroid)

