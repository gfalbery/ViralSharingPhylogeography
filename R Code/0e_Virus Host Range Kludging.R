
# Making phylogenetic centroids and host range

library(centiserve); library(tidyverse)

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

# Making centroids ####
# This takes a long time

SubGraphList <- CentroidList <- list()

for(x in unique(Viruses$Sp)){
  
  SubGraph <- induced_subgraph(PhyloGraph, V(PhyloGraph)$name%in%VirusAssocs[[x]])
  SubGraphList[[x]] <- SubGraph
  
  if(class(try(centroid(SubGraph), silent = T)) != "try-error"){
    SubCentroids <- centroid(SubGraph)
    Centroids <- names(SubCentroids[SubCentroids==min(SubCentroids)])
    CentroidList[[x]] <- Centroids
  } else CentroidList[[x]] <- NA
  
  print(x)
}

y = last(names(CentroidList))

for(x in unique(Viruses$Sp)[(which(unique(Viruses$Sp)==y)+1):length(unique(Viruses$Sp))]){
  
  SubGraph <- induced_subgraph(PhyloGraph, V(PhyloGraph)$name%in%VirusAssocs[[x]])
  SubGraphList[[x]] <- SubGraph
  
  if(class(try(centroid(SubGraph), silent = T)) != "try-error"){
    SubCentroids <- centroid(SubGraph)
    Centroids <- names(SubCentroids[SubCentroids==min(SubCentroids)])
    CentroidList[[x]] <- Centroids
  } else CentroidList[[x]] <- NA
  print(x)
}

save(CentroidList, file = "data/CentroidList.Rdata")
load("data/CentroidList.Rdata")

# Creating Viral Host Phylo Matrix ####

VirusHostPD <- matrix(NA, nrow = length(VirusAssocs), ncol = length(VirusAssocs))
dimnames(VirusHostPD) <- list(names(VirusAssocs),names(VirusAssocs))

# 1. Take the centroids of all hosts among viral pairs.
# 2. Average across the sub-CytBmatrix.  

for(i in rownames(VirusHostPD)){
  
  FocalCentroids <- CentroidList[[i]][CentroidList[[i]]%in%rownames(CytBMatrix)]
  
  VirusHostPD[i,] <- sapply(
    names(VirusAssocs), function(x){
      
      if(x %in% intersect(names(CentroidList), 
                          names(VirusAssocs))&
         any(VirusAssocs[[x]] %in% rownames(CytBMatrix))){
        
        CentroidHosts <- CentroidList[[x]][CentroidList[[x]]%in%rownames(CytBMatrix)]
        
        SubMatrix <- CytBMatrix[FocalCentroids, CentroidHosts]
        
        Distance = mean(SubMatrix) 
        
        return(Distance)
        
      } else NA 
    }
  )
  print(i)
}

FVN2 <- reduce(list(SpatialViruses$Sp, 
                    colnames(VirusHostPD)[-which(is.na(VirusHostPD[1,]))],
                    rownames(VirusRangeAdj1)[which(sapply(GridList[unique(SpatialViruses$Sp)], length)>0)]), 
               intersect)

qplot(c(VirusRangeAdj1[FVN2, FVN2]),
      1-c(VirusHostPD[FVN2, FVN2]), geom = c("point", "smooth"))










