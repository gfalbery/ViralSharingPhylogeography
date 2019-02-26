
# Bird Spatial Data Import ####

library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools)

if(!file.exists("data/BirdRanges.Rdata")){
  
  BirdShapes <- st_read("~/AllBirds")
  
  BirdShapes2 <- st_transform(BirdShapes, 
                              "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") # Mollweide projection 
  
  BirdShapes2$Binomial = str_replace(BirdShapes2$SCINAME, " ", "_")
  #BirdShapes2 <- BirdShapes2[BirdShapes2$Binomial%in%Hosts$Sp,]
  BirdShapes2 <- BirdShapes2[order(BirdShapes2$Binomial),]
  BirdRaster <- raster(BirdShapes2, res = 50000) # NB units differ from Mercator!
  
  BirdRanges <- fasterize(BirdShapes2, BirdRaster, by = "Binomial")
  
  save(BirdRanges, file = "data/BirdRanges.Rdata")
  
} else load("data/BirdRanges.Rdata")


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


BirdRangeAdj <- PairsWisely(BirdRanges)

save(BirdRangeAdj, file = "BirdRangeAdj.Rdata")