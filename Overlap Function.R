
# Greg's range extraction function ####

PairsWisely <- function(Rasterstack, Species = "All"){
  
  library(raster); library(tidyverse)
  
  # Converting these to meaningful values ####
  
  FullValuedf <- data.frame(getValues(Rasterstack))
  FullValuedf2 <- reshape2::melt(FullValuedf)
  FullValuedf2$x <- rep(1:Rasterstack[[1]]@ncols, Rasterstack[[1]]@nrows)
  FullValuedf2$y <- rep(Rasterstack[[1]]@nrows:1, each = Rasterstack[[1]]@ncols)
  
  FullRangedf <- FullValuedf %>% 
    dplyr::rename(Host = variable, Presence = value)
  
  if(Species != "All") FullRangedf <- FullRangedf %>% filter(Host %in% Species)
  
  FullRangedf$GridID <- with(FullRangedf, paste(x, y))
  
  Range0 <- levels(FullRangedf$Host)[which(table(FullRangedf$Host)==0)] # Hosts that have no spatial records??
  FullRangedf <- droplevels(FullRangedf) 
  FullRangedf <- FullRangedf[order(FullRangedf$Host),]
  
  # Could use igraph to project it into bipartite host-grid matrix
  # Or could do this bullshit
  
  FullRangeOverlap <- matrix(0, nrow = nlevels(FullRangedf$Host), ncol = nlevels(FullRangedf$Host))
  dimnames(FullRangeOverlap) <- list(levels(FullRangedf$Host),levels(FullRangedf$Host))
  
  for(x in levels(FullRangedf$Host)){
    
    if(x == first(levels(FullRangedf$Host))) t1 <- Sys.time()
    
    Grids <- FullRangedf[FullRangedf$Host==x,"GridID"]
    SubFullRangedf <- FullRangedf[FullRangedf$GridID %in% Grids,]
    
    FullRangeOverlap[x,] <- table(SubFullRangedf$Host)
    
    print(x)
    
    if(x == last(levels(FullRangedf$Host))) t2 <- Sys.time()
    
  }
  
  FullRangeA = matrix(rep(diag(FullRangeOverlap), nrow(FullRangeOverlap)), nrow(FullRangeOverlap))
  FullRangeB = matrix(rep(diag(FullRangeOverlap), each = nrow(FullRangeOverlap)), nrow(FullRangeOverlap))
  
  FullRangeAdj1 <- FullRangeOverlap/(FullRangeA + FullRangeB - FullRangeOverlap) # Weighted evenly
  FullRangeAdj2 <- FullRangeOverlap/(FullRangeA) # Asymmetrical
  
  return(FullRangeAdj1)
  
}