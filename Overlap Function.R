
# Greg's range extraction function ####

PairsWisely <- function(Rasterstack, Species = "All"){
  
  library(raster); library(tidyverse); library(Matrix)
  
  t1 <- Sys.time()
  
  print("Getting the grid values")
  
  Valuedf <- data.frame(getValues(Rasterstack)) %>% as.matrix
  Valuedf[is.na(Valuedf)] <- 0
  
  Valuedf <- Valuedf %>% as("dgCMatrix")
  
  if(Species != "All") Valuedf <- Valuedf[,Species]
  
  RangeOverlap <- matrix(0, nrow = ncol(Valuedf), ncol = ncol(Valuedf))
  dimnames(RangeOverlap) <- list(colnames(Valuedf),colnames(Valuedf))
  
  print("Calculating Overlap")
  
  for(x in colnames(Valuedf)){
    
    print(x)
    
    TrainGrids <- Valuedf[,x]
    
    SubRangedf <- Valuedf[which(TrainGrids==1),]
    
    if(!is.null(dim(SubRangedf))){
      
      RangeOverlap[x,] <- apply(SubRangedf, 2, function(a) length(which(a==1)))
      
    } else   RangeOverlap[x,] <- SubRangedf
    
  }
  
  FullRangeA = matrix(rep(diag(RangeOverlap), nrow(RangeOverlap)), nrow(RangeOverlap))
  FullRangeB = matrix(rep(diag(RangeOverlap), each = nrow(RangeOverlap)), nrow(RangeOverlap))
  
  RangeAdj <- RangeOverlap/(FullRangeA + FullRangeB - RangeOverlap) # Weighted evenly

  TimeTaken = Sys.time() - t1
  
  print(TimeTaken)
  
  return(RangeAdj)
  
}
