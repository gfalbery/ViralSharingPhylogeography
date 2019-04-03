
# Parsing Cory's maps ####

library(velox);
library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools); library(SpRanger)

Files <- list.files("~/IcebergRasters")

table(table(Files))

ToRemove <- Files[which(sapply(Files, function(a) length(list.files(paste("~/IcebergRasters", a, sep = '/'))))==0)]

file.remove(paste("~/IcebergRasters", ToRemove, sep = '/'))

Files <- setdiff(Files, ToRemove)

#substr(Files, 1,1) <- substr(Files, 1,1) %>% toupper()
#Files <- unique(Files)

VeloxList <- lapply(Files, function(a){
  print(which(Files==a))
  print(a)
  r1 <- velox(list.files(paste("~/IcebergRasters", a, sep = '/'), full.names = TRUE)[1])
})

RasterLista <- lapply(1:length(VeloxList), function(a){
  
  print(Files[a])
  
  VeloxList[[a]]$as.RasterLayer(band = 1) #%>% rasterToPolygons(dissolve = T)
  
})

RasterListb <- lapply(1:length(RasterLista), function(a){
  
  print(Files[a])
  
  RasterLista[[a]] %>% rasterToPolygons(dissolve = T)
  
})

Files <- Files[RasterListb[-which(sapply(RasterListb, length)==0)]]
RasterListb <- RasterListb[-which(sapply(RasterListb, length)==0)]

RasterListc <- st_as_sf(bind(RasterListb))

UpperFiles <- Files
substr(UpperFiles, 1,1) <- substr(UpperFiles, 1,1) %>% toupper()

RasterListc$Binomial <- UpperFiles

IcebergRaster <- raster(RasterListc, res = 0.25)

RasterBrick <- fasterize(RasterListc, IcebergRaster, by = "Binomial")

IcebergAdj <- PairsWisely(RasterBrick)

save(IcebergAdj, file = "IcebergRangeAdj.Rdata")
save(RasterBrick, file = "IcebergRasterBrick.Rdata")

load("data/FullRangeOverlap.Rdata")

IUCN_ENM_Sp <- intersect(rownames(FullRangeAdj), UpperFiles)

SuboverlapIUCN <- FullRangeAdj[IUCN_ENM_Sp,IUCN_ENM_Sp]
SuboverlapENM <- IcebergAdj[IUCN_ENM_Sp,IUCN_ENM_Sp]

IUCN_ENM_DF <- data.frame(
  Sp = rep(IUCN_ENM_Sp, each = length(IUCN_ENM_Sp)),
  Sp2 = rep(IUCN_ENM_Sp, length(IUCN_ENM_Sp)),
  IUCN = c(SuboverlapIUCN),
  ENM = c(SuboverlapENM)
) %>% slice(which(lower.tri(SuboverlapIUCN)))

ggplot(IUCN_ENM_DF, aes(IUCN, ENM)) + geom_point()

ggplot(IUCN_ENM_DF, aes(IUCN, ENM)) + 
  geom_point(alpha = 0.01) + 
  coord_fixed() + 
  geom_smooth() + AlberTheme + 
  ggsave("Iceberg_Correlations2.jpeg", 
         units = "mm", width = 100, height = 100, dpi = 300)



