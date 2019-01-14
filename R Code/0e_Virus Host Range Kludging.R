
# Making phylogenetic centroids and phylogenetic host range (/distance)

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
           "HostRangeMedian",
           "Hosts")] <- cbind(VirPhyloHostRangeMax[Viruses$Sp],
                              VirPhyloHostRangeMean[Viruses$Sp],
                              VirPhyloHostRangeMedian[Viruses$Sp])

jpeg("Figures/Pairs of records, host range, centrality.jpeg", units = "mm", width = 100, height = 100, res = 300)
Viruses[,c("Records", "HostRangeMean", "vEigenvector")] %>% mutate(Records = log(Records)) %>% ggpairs(lower = list(continuous = "smooth")) +
  ggtitle("Pairs of records, host range, centrality")
dev.off()

ggplot(Viruses, aes(log(Records), HostRangeMean)) + geom_point() + geom_smooth() # More accurate as a GAM process

BarGraph(Viruses, "vFamily", "HostRangeMean") +  # Lot of variation 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Making centroids ####
# This takes a long, long time

M2 <- as.matrix(CytBMatrix)

M2[upper.tri(M2)]

PhyloGraph <- graph.adjacency(1-M2, weighted = T, mode = "undirected", diag = F)

SubGraphList <- CentroidList <- list()

if(file.exists("data/CentroidList.Rdata")) load("data/CentroidList.Rdata") else {
  
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

# Creating Viral Host Phylo Matrix ####

VirusHostPD <- matrix(NA, nrow = length(VirusAssocs), ncol = length(VirusAssocs))
dimnames(VirusHostPD) <- list(names(VirusAssocs),names(VirusAssocs))

# 1. Take the centroids of all hosts among viral pairs, take their interrelationship values.
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

# Doing this just for humans ####

Viruses$Centroid_Human_Distance <- sapply(
  names(VirusAssocs), function(x){
    
    if(x %in% intersect(names(CentroidList), 
                        names(VirusAssocs))&
       any(VirusAssocs[[x]] %in% rownames(CytBMatrix))){
      
      CentroidHosts <- CentroidList[[x]][CentroidList[[x]]%in%rownames(CytBMatrix)]
      
      SubMatrix <- CytBMatrix[CentroidHosts, "Homo_sapiens"]
      
      Distance = mean(SubMatrix) 
      
      return(Distance)
      
    } else NA 
  }
)

# Non-Human Centroid Calculation ####

HumanRemove <- which(rownames(CytBMatrix) == "Homo_sapiens")

M2 <- as.matrix(CytBMatrix[-HumanRemove, -HumanRemove])

PhyloGraph2 <- graph.adjacency(1-M2, weighted = T, mode = "undirected", diag = F)

SubGraphList2 <- CentroidList2 <- list()

if(file.exists("data/CentroidList2.Rdata")) load("data/CentroidList2.Rdata") else {
  
  for(x in unique(Viruses$Sp)){
    
    SubGraph <- induced_subgraph(PhyloGraph2, V(PhyloGraph2)$name%in%VirusAssocs[[x]])
    SubGraphList2[[x]] <- SubGraph
    
    if(class(try(centroid(SubGraph), silent = T)) != "try-error"){
      SubCentroids <- centroid(SubGraph)
      Centroids <- names(SubCentroids[SubCentroids==min(SubCentroids)])
      CentroidList2[[x]] <- Centroids
    } else CentroidList2[[x]] <- NA
    
    print(x)
  }
  
}

y = last(names(CentroidList2))

for(x in unique(Viruses$Sp)[(which(unique(Viruses$Sp)==y)):length(unique(Viruses$Sp))]){
  
  SubGraph <- induced_subgraph(PhyloGraph, V(PhyloGraph)$name%in%VirusAssocs[[x]])
  SubGraphList2[[x]] <- SubGraph
  
  if(class(try(centroid(SubGraph), silent = T)) != "try-error"){
    SubCentroids <- centroid(SubGraph)
    Centroids <- names(SubCentroids[SubCentroids==min(SubCentroids)])
    CentroidList2[[x]] <- Centroids
  } else CentroidList2[[x]] <- NA
  print(x)
}

save(CentroidList2, file = "data/CentroidList2.Rdata")

# Creating Viral Host Phylo Matrix ####

VirusHostPD2 <- matrix(NA, nrow = length(VirusAssocs), ncol = length(VirusAssocs))
dimnames(VirusHostPD2) <- list(names(VirusAssocs),names(VirusAssocs))

# 1. Take the centroids of all hosts among viral pairs, take their interrelationship values.
# 2. Average across the sub-CytBmatrix.  

for(i in rownames(VirusHostPD2)){
  
  FocalCentroids <- CentroidList2[[i]][CentroidList2[[i]]%in%rownames(CytBMatrix)]
  
  VirusHostPD2[i,] <- sapply(
    names(VirusAssocs), function(x){
      
      if(x %in% intersect(names(CentroidList2), 
                          names(VirusAssocs))&
         any(VirusAssocs[[x]] %in% rownames(CytBMatrix))){
        
        CentroidHosts <- CentroidList2[[x]][CentroidList2[[x]]%in%rownames(CytBMatrix)]
        
        SubMatrix <- CytBMatrix[FocalCentroids, CentroidHosts]
        
        Distance = mean(SubMatrix) 
        
        return(Distance)
        
      } else NA 
    }
  )
  print(i)
}

FVN3 <- reduce(list(SpatialViruses$Sp, 
                    colnames(VirusHostPD2)[-which(is.na(VirusHostPD2[1,]))],
                    rownames(VirusRangeAdj1)[which(sapply(GridList[unique(SpatialViruses$Sp)], length)>0)]), 
               intersect)

qplot(c(VirusRangeAdj1[FVN2, FVN2]),
      1-c(VirusHostPD2[FVN2, FVN2]), geom = c("point", "smooth"))

# Doing this just for humans ####

Viruses$Centroid_Human_Distance2 <- sapply(
  names(VirusAssocs), function(x){
    
    if(x %in% intersect(names(CentroidList2), 
                        names(VirusAssocs))&
       any(VirusAssocs[[x]] %in% rownames(CytBMatrix))){
      
      CentroidHosts <- CentroidList2[[x]][CentroidList2[[x]]%in%rownames(CytBMatrix)]
      
      SubMatrix <- CytBMatrix[CentroidHosts, "Homo_sapiens"]
      
      Distance = mean(SubMatrix) 
      
      return(Distance)
      
    } else NA 
  }
)





