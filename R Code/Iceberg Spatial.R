
# Parsing Cory's maps ####

library(velox);
library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools); library(SpRanger)

IcebergAdjList <- list()

PredReps <- c("Currents", paste0("Futures", 1:4))

x = 1

for(x in 1:length(PredReps)){
  
  print(PredReps[x])
  
  Files <- list.files(paste0("Iceberg Input Files/",PredReps[x]))
  
  VeloxList <- lapply(Files, function(a){
    print(which(Files==a))
    print(a)
    r1 <- velox(paste(paste0("Iceberg Input Files/",PredReps[x]), a, sep = '/'))
  })
  
  RasterLista <- lapply(1:length(VeloxList), function(a){
    
    print(Files[a])
    
    VeloxList[[a]]$as.RasterLayer(band = 1) #%>% rasterToPolygons(dissolve = T)
    
  })
  
  Method = "resample"
  
  if(Method == "resample"){
    
    # Using resample (quicker but maybe doesn't actually work??) ####
    
    ### fix blank
    
    # blank <- matrix(0,360*2,720*2) # proper resolution
    blank <- matrix(0, 360, 720) # quick and dirty resolution
    blank <- raster(blank)
    extent(blank) <- c(-180,180,-90,90)
    projection(blank) <- CRS("+proj=longlat +datum=WGS84")
    
    if(x==1){
      
      RasterListb <- lapply(1:length(RasterLista), function(a){
        
        print(a)
        
        testraster <- RasterLista[[a]]
        testraster <- raster::resample(testraster, blank, method = 'ngb')
        
        return(testraster)
        
      })
      
    } else {
      
      RasterListb <- RasterLista
      
    }
    
    RasterBrick <- raster::brick(RasterListb)
    names(RasterBrick) <- Files %>% str_remove(".tif$") %>% str_replace(" ", "_")
    
    # Colin says don't do this: crs(RasterBrick) <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    
    crs(RasterBrick) <- "+proj=longlat +datum=WGS84"
    
    RasterBrick <- projectRaster(RasterBrick, crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs", 
                                 method = "bilinear")
    
    IcebergAdj <- PairsWisely(RasterBrick)
    
  } else {
    
    # Using rastertopolygons and SF package (which definitely works but takes ages) ####
    
    RasterListb <- lapply(1:length(RasterLista), function(a){
      
      print(Files[a])
      
      RasterLista[[a]] %>% rasterToPolygons(dissolve = T)
      
    })
    
    RasterListc <- st_as_sf(bind(RasterListb))
    
    RasterListc <- st_transform(RasterListc, 
                                "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") # Mollweide projection 
    
    UpperFiles <- Files
    substr(UpperFiles, 1,1) <- substr(UpperFiles, 1,1) %>% toupper()
    
    RasterListc$Binomial <- UpperFiles
    
    IcebergRaster <- raster(RasterListc, res = 25000)
    
    RasterBrick <- fasterize(RasterListc, IcebergRaster, by = "Binomial")
    
    # The next stage ####
    IcebergAdj <- PairsWisely(RasterBrick)
    
  }
  
  save(IcebergAdj, file = paste0(paste0("Iceberg Output Files/",PredReps[x]),"/IcebergRangeAdj.Rdata"))
  save(RasterBrick, file = paste0(paste0("Iceberg Output Files/",PredReps[x]),"/IcebergRasterBrick.Rdata"))
  
  rownames(IcebergAdj) <- colnames(IcebergAdj) <- rownames(IcebergAdj) %>% str_replace('[.]',"_")
  
  IcebergAdjList[[x]] <- IcebergAdj
  
}

save(IcebergAdjList, file = "Iceberg Output Files/IcebergAdjList.Rdata")

# Parsing missing species ####

Removed <- setdiff(rownames(IcebergAdjList[[1]]), rownames(IcebergAdjList[[2]]))

Removed <- intersect(Removed, names(RasterBrick))

sapply(Removed, function(a) RasterBrick[[a]] %>% values %>% na.omit %>% length)
sapply(Removed, function(a) RasterLista[[a]] %>% values %>% na.omit %>% length)

table(sapply(intersect(Removed, names(RasterBrick)), function(a) RasterBrick[[a]] %>% values %>% na.omit %>% length))



