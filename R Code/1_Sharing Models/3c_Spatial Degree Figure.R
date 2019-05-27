
# Summarising predicted degree by grid 

# Rscript "R Code/1_Sharing Models/3c_Spatial Degree Figure.R"

library(tidyverse); library(raster); library(colorspace)

load("Output Files/Panth1.Rdata")

if(file.exists("Output Files/GridDegree.Rdata")) load("Output Files/GridDegree.Rdata") else{
  
  load("~/LargeFiles/MammalStackFullMercator.Rdata")
  
  RasterListb <- lapply(1:length(AllMammals), function(a){
    
    if(a %% 500==0) print(a)
    
    MammalStackFull[[AllMammals[a]]]
    
  })
  
  ToSkip <- which(sapply(RasterListb, is.null))
  
  names(RasterListb) <- AllMammals
  
  # RasterListb <- RasterListb[-ToSkip]
  # AllMammalsSub <- ToBagSpecies[-ToSkip]
  
  OverlapSums <- rep(0, ncol(RasterListb[[1]])*nrow(RasterListb[[1]]))
  
  for(i in 1:length(RasterListb)){
    
    if(i %% 500 == 0) print(i)
    SubSums <- raster::getValues(RasterListb[[i]]) > 0 %>% as.numeric
    SubSums[is.na(SubSums)] <- 0
    OverlapSums <- OverlapSums + SubSums
    
  }
  
  AllDegree <- rep(0, ncol(RasterListb[[1]])*nrow(RasterListb[[1]]))
  InDegree <- rep(0, ncol(RasterListb[[1]])*nrow(RasterListb[[1]]))
  OutDegree <- rep(0, ncol(RasterListb[[1]])*nrow(RasterListb[[1]]))
  
  for(i in 1:length(RasterListb)){
    
    if(i %% 500 == 0) print(i)
    SubSums <- raster::getValues(RasterListb[[i]]) > 0 %>% as.numeric
    SubSums[is.na(SubSums)] <- 0
    AllDegree <- AllDegree + SubSums*Panth1[Panth1$Sp==AllMammals[i],"AllPredDegree"]
    InDegree <- InDegree + SubSums*Panth1[Panth1$Sp==AllMammals[i],"InDegree"]
    OutDegree <- OutDegree + SubSums*Panth1[Panth1$Sp==AllMammals[i],"OutDegree"]
    
  }
  
  DegreeDF <- data.frame(
    X = rep(1:ncol(MammalStackFull[[1]]), nrow(MammalStackFull[[1]])),
    Y = rep(nrow(MammalStackFull[[1]]):1, each = ncol(MammalStackFull[[1]])),
    
    Richness = OverlapSums,
    AllDegree = AllDegree,
    InDegree = InDegree,
    OutDegree = OutDegree
    
  ) 
  
  UniversalBlank <- raster("Iceberg Input Files/UniversalBlank.tif")
  Land = which(raster::values(UniversalBlank)==0)
  Sea = which(is.na(raster::values(UniversalBlank)))
  
  DegreeDF <- DegreeDF[-Sea,]
  
  GridDegree <- DegreeDF %>% 
    mutate_at(vars(contains("Degree")), function(a) a/DegreeDF$Richness) %>% 
    mutate_at(vars(contains("Degree")), function(a) ifelse(a==0|is.na(a), min(a[a>0], na.rm = T), a)) %>% 
    gather(key = "Metric", value = "Degree", contains("Degree"))
  
  save(GridDegree, file = "Output Files/GridDegree.Rdata")
  
}
