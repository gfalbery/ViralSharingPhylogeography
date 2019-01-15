
# Making phylogenetic centroids and phylogenetic host range (/distance)

library(centiserve); library(tidyverse); library(GGally)

human <- which(rownames(CytBMatrix) == "Homo_sapiens")

nhCytBMatrix <- CytBMatrix[-human,-human] # Want this to be non-human host range 

VirPhyloHostRangeMax <- sapply(VirusAssocs, function(a){
  
  if(length(a[a%in%rownames(nhCytBMatrix)])>0){
    
    b <- a[a%in%rownames(nhCytBMatrix)]
    
    if(length(b)>1){
      max(nhCytBMatrix[b, b][upper.tri(nhCytBMatrix[b, b])])
    } else 0
  } else NA
  
})

VirPhyloHostRangeMean <- sapply(VirusAssocs, function(a){
  
  if(length(a[a%in%rownames(nhCytBMatrix)])>0){
    
    b <- a[a%in%rownames(nhCytBMatrix)]
    
    if(length(b)>1){
      mean(nhCytBMatrix[b, b][upper.tri(nhCytBMatrix[b, b])])
    } else 0
  } else NA
  
})

VirPhyloHostRangeMedian <- sapply(VirusAssocs, function(a){
  
  if(length(a[a%in%rownames(nhCytBMatrix)])>0){
    
    b <- a[a%in%rownames(nhCytBMatrix)]
    
    if(length(b)>1){
      median(nhCytBMatrix[b, b][upper.tri(nhCytBMatrix[b, b])])
    } else 0
  } else NA
  
})

Viruses[,c("HostRangeMax",
           "HostRangeMean",
           "HostRangeMedian")] <- cbind(VirPhyloHostRangeMax[as.character(Viruses$Sp)],
                                        VirPhyloHostRangeMean[as.character(Viruses$Sp)],
                                        VirPhyloHostRangeMedian[as.character(Viruses$Sp)])

jpeg("Figures/Pairs of records, host range, centrality.jpeg", units = "mm", width = 100, height = 100, res = 300)
Viruses[,c("Records", "HostRangeMean", "HostRangeMax", "vEigenvector", "vDegree")] %>% 
  mutate(Records = log(Records)) %>% 
  ggpairs(lower = list(continuous = wrap("smooth", method = "loess"))) +
  ggtitle("Pairs plot of records, host range, centrality")
dev.off()

# Making centroids ####
# This takes a long, long time

M2 <- as.matrix(CytBMatrix)

M2[upper.tri(M2)]

PhyloGraph <- graph.adjacency(1-M2, weighted = T, mode = "undirected", diag = F)

SubGraphList <- CentroidList <- list()

for(x in unique(as.character(Viruses$Sp))){
  SubGraph <- induced_subgraph(PhyloGraph, V(PhyloGraph)$name%in%VirusAssocs[[x]])
  SubGraphList[[x]] <- SubGraph
}

if(file.exists("data/CentroidList.Rdata")) load("data/CentroidList.Rdata") else {
    
  for(x in unique(as.character(Viruses$Sp))){
    
    if(class(try(centroid(SubGraph[[x]]), silent = T)) != "try-error"){ # Throws up errors with cliques.
      SubCentroids <- centroid(SubGraph[[x]])
      Centroids <- names(SubCentroids[SubCentroids==min(SubCentroids)])
      CentroidList[[x]] <- Centroids
    } else CentroidList[[x]] <- NA
    
    print(x)
  }
  
  for(y in names(which(sapply(CentroidList, function(a) any(is.na(a)))))){
    
    if(nunique(degree(SubGraphList[[y]]))==1){
      CentroidList[[y]] <- V(SubGraphList[[y]])$name
    }
  }
  
  save(CentroidList, file = "data/CentroidList.Rdata")
  
}


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

# Calculating Distance from Humans ####

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

# Calculating spatial centroid of each virus' host centroid ####

CentroidCentroids <- lapply(
  names(VirusAssocs), function(x){
    
    if(x %in% intersect(names(CentroidList), 
                        names(VirusAssocs))&
       any(VirusAssocs[[x]] %in% HostCentroids$Host)){
      
      CentroidHosts <- HostCentroids[CentroidList[[x]],]
      
      LongMean = mean(CentroidHosts$LongMean, na.rm = T)
      LatMean = mean(CentroidHosts$LatMean, na.rm = T)
      
      return(data.frame(Sp = x, LongMean, LatMean))
      
    } else data.frame(Sp = x, LongMean = NA, LatMean = NA) 
  }
) %>% bind_rows

Viruses <- merge(Viruses, CentroidCentroids, 
                 by = "Sp", all.x = T, 
                 suffixes = c(".Total", ".Centroid"))
