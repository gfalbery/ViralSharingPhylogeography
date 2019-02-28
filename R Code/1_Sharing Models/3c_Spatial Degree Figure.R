
# Summarising predicted degree by grid 

library(tidyverse); library(raster); library(colorspace)

load("~/Albersnet/data/FullMammalRanges.Rdata")
load("Output Files/Panth1.Rdata")

FullValuedf <- data.frame(getValues(FullMammalRanges)); detach(package:raster)
FullValuedf2 <- reshape2::melt(FullValuedf)
FullValuedf2$x <- rep(1:FullMammalRanges[[1]]@ncols, FullMammalRanges[[1]]@nrows)
FullValuedf2$y <- rep(FullMammalRanges[[1]]@nrows:1, each = FullMammalRanges[[1]]@ncols)

detach(package:raster)

FullRangedf <- FullValuedf2 %>% 
  filter(!is.na(value)) %>% droplevels %>%
  dplyr::rename(Host = variable, Presence = value)

save(FullRangedf, file = "Output Files/FullRangedf.Rdata")
#load("Output Files/FullRangedf.Rdata")

FullRangedf2 <- FullRangedf %>% left_join(Panth1[,c("Sp","AllPredDegree", "InDegree", "OutDegree")], 
                                          by = c("Host" = "Sp")) %>%
  filter(!Host=="Ursus_maritimus")

GridDegree <- FullRangedf2 %>% group_by(x,y) %>% 
  summarise_at(vars(ends_with("Degree")), function(a) mean(a, na.rm = T))

#save(GridDegree, file = "Output Files/GridDegree.Rdata")

GridDegree2 <- gather(GridDegree, key = "Metric", value = "Degree", ends_with("Degree"))

save(GridDegree2, file = "Output Files/GridDegree2.Rdata")

GridDegree %>% mutate_at(vars(ends_with("Degree")), 
                         function(a){ifelse(abs(scale(a))>5,
                                            NA,
                                            a
                         )})

GridDegree3 <- GridDegree %>% 
  mutate(AllPredDegree = ifelse(abs(scale(AllPredDegree)[,1])>5,
                                NA,
                                AllPredDegree
  ),
  InDegree = ifelse(abs(scale(InDegree)[,1])>5,
                    NA,
                    InDegree
  ),
  OutDegree = ifelse(abs(scale(OutDegree)[,1])>5,
                     NA,
                     OutDegree
  )
  )

for(j in c("AllPredDegree","InDegree","OutDegree")){
  
  print(j)
  
  for(i in which(is.na(GridDegree3[,j]))){
    
    print(i/nrow(GridDegree3))
    
    xSurroundings <- GridDegree3[i,"x"]
    ySurroundings <- GridDegree3[i,"y"]
    
    xMeans <- GridDegree3 %>% filter(x == xSurroundings$x) %>% 
      slice(which(y-ySurroundings$y == min(y-ySurroundings$y)))
    
    xMeans <- GridDegree3 %>% 
      slice(-which(is.na(j))) %>%
      filter(x == xSurroundings$x) %>%
      mutate(a = y-ySurroundings$y) 
    
    if(nrow(xMeans)>0){
      xMeans <- xMeans %>%
        mutate(a = a-min(a)) %>%
        filter(a<10)
    } else xMeans <- NULL
    
    yMeans <- GridDegree3 %>% 
      slice(-which(is.na(j))) %>%
      filter(y == ySurroundings$y) %>%
      mutate(a = x-xSurroundings$x)
    
    if(nrow(yMeans)>0){
      yMeans <- yMeans %>%
        mutate(a = a-min(a)) %>%
        filter(a<10)
    } else yMeans <- NULL
    
    if(!(is.null(xMeans)&is.null(yMeans))){
      GridDegree3[i,j] <- mean(c(xMeans[,j],yMeans[,j]))
    }
  }
}

GridDegree4 <- gather(GridDegree3, key = "Metric", value = "Degree", ends_with("Degree"))






