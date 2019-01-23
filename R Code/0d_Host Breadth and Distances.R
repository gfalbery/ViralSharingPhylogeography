
# Making phylogenetic centroids and phylogenetic host range (/distance)

library(centiserve); library(tidyverse); library(GGally); library(igraph)

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

# Making host distances ####

# Making bipartite projections ####

library(igraph)

mFull <- table(AssocsBase[,1:2])
MFull <- as.matrix(mFull)

bipGraphFull <- graph.incidence(MFull, weighted = T)

VirusGraphFull <- bipartite.projection(bipGraphFull)$proj1
HostGraphFull <- bipartite.projection(bipGraphFull)$proj2

hHumanDist <- distances(bipGraphFull, v = V(bipGraphFull)[(dim(Viruses)[1]+1):length(V(bipGraphFull))],
                        to = V(bipGraphFull)[which(names(V(bipGraphFull)) == "Homo_sapiens")])

hHumanDist[,1] <- ifelse(hHumanDist[,1] == Inf, Inf, (hHumanDist[,1])/2)

hHumanDist2 <- data.frame(Sp = rownames(hHumanDist),
                          HumanDist = hHumanDist[,1])

hDomDist <- distances(bipgraph, v = V(bipgraph)[(dim(Viruses)[1]+1):length(V(bipgraph))],
                      to = V(bipgraph)[which(names(V(bipgraph)) %in% Domestics)])

hDomDist <- (hDomDist)/2

hDomDist2 <- data.frame(Sp = rownames(hDomDist),
                        MeanDomDist = apply(hDomDist, 1, function(a) mean(a[!a=="Inf"], na.rm = T)),
                        MinDomDist = apply(hDomDist, 1, function(a) ifelse(all(a==Inf),Inf,min(a[!a==Inf], na.rm = T))),
                        MaxDomDist = apply(hDomDist, 1, function(a) max(a, na.rm = T)))

Hosts <- merge(Hosts, hHumanDist2[,c("Sp","HumanDist")], by = "Sp", all.x = T)
Hosts <- merge(Hosts, hDomDist2[,c("Sp", "MeanDomDist", "MinDomDist", "MaxDomDist")], by = "Sp", all.x = T)

Hosts[Hosts$Sp%in%rownames(CytBMatrix),"TotalPD"] <- rowSums(CytBMatrix[as.character(Hosts[Hosts$Sp%in%rownames(CytBMatrix),"Sp"]),])


# Deriving distances from humans and domestic animals ####

vHumanDist <- distances(bipGraphFull, v = V(bipGraphFull)[1:(dim(Viruses)[1])],
                        to = V(bipGraphFull)[which(names(V(bipGraphFull)) == "Homo_sapiens")])

vHumanDist[,1] <- ifelse(vHumanDist[,1] == Inf, Inf, (vHumanDist[,1]-1)/2)

vHumanDist2 <- data.frame(Sp = rownames(vHumanDist),
                          HumanDist = vHumanDist[,1],
                          Domestic = Viruses$Domestic,
                          Family = Viruses$vFamily)

vDomDist <- distances(bipgraph, v = V(bipgraph)[1:(dim(Viruses)[1])],
                      to = V(bipgraph)[which(names(V(bipgraph)) %in% Domestics)])

vDomDist <- (vDomDist-1)/2

vDomDist2 <- data.frame(Sp = rownames(vDomDist),
                        MeanDomDist = apply(vDomDist, 1, function(a) mean(a[!a=="Inf"], na.rm = T)),
                        MinDomDist = apply(vDomDist, 1, function(a) ifelse(all(a==Inf),Inf,min(a[!a==Inf], na.rm = T))),
                        MaxDomDist = apply(vDomDist, 1, function(a) max(a, na.rm = T)),
                        Human = Viruses$Human,
                        Family = Viruses$vFamily)

# Merging distances and traits ####

Viruses <- merge(Viruses, 
                 vDomDist2[,c("Sp", "MeanDomDist", "MinDomDist", "MaxDomDist")], 
                 by = "Sp", all.x = T)

Viruses <- merge(Viruses, 
                 vHumanDist2[,c("Sp", "HumanDist")], 
                 by = "Sp", all.x = T)

