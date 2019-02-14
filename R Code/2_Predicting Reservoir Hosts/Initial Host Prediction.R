
# Predicting the reservoir host(s) of a focal virus for Colin ####

library(tidyverse); library(parallel)

pVirus <- "Crimean-Congo_hemorrhagic_fever_virus" 

"Crimean-Congo"
"Zika"
"Bundibugyo and Tai Forest Ebola"
pVirus <- "Zaire_ebolavirus"

ResHostList <- list()

a = 1

for(a in 1:5){#length(names(VirusAssocs))){
  
  print(names(VirusAssocs)[a])
  
  pHosts <- VirusAssocs[[a]]
  
  pHosts <- intersect(pHosts, AllMammals)
  
  if(length(pHosts)>0){
    
    EstList <- parallel::mclapply(1:length(AllSims), function(x){
      
      FocalNet <- AllSims[[x]] %>% as.matrix
      
      pHosts2 <- intersect(pHosts, rownames(FocalNet))
      pHosts3 <- setdiff(colnames(FocalNet), pHosts2)
      
      Estimates <- FocalNet[pHosts2,pHosts3]
      
      if(is.null(dim(Estimates))) Estimates <- rbind(Estimates, Estimates)/2
      
      Ests <- data.frame(Sp = names(sort(colSums(Estimates))),
                         Count = sort(colSums(Estimates))/nrow(Estimates),
                         Iteration = x)
      
      Ests
      
    }, mc.cores = 20)
    
    Estdf <- EstList %>% 
      bind_rows %>% 
      group_by(Sp) %>% 
      summarise(ResCount = mean(Count)) %>% 
      slice(order(ResCount, decreasing= T))
    
    ResHostList[[names(VirusAssocs)[a]]] <- Estdf
    
  } else ResHostList[[names(VirusAssocs)[a]]] <- NA
  
}

PredHostPolygons <- lapply(ResHostList[2:5], function(a){
  
  FullPolygons %>% filter(Host%in%a$Sp) %>% 
    left_join(a, by = c("Host" = "Sp")) %>%
    filter(ResCount>max(a$ResCount) - sd(a$ResCount))
  
})

PredHostPolygons[[1]] %>% #filter(ResCount>0.9) %>%
  
  ggplot(aes(long, lat, group = paste(Host, group))) + 
  geom_path(data = WorldMap, inherit.aes = F, aes(long, lat, group = group), colour = "black") +
  geom_polygon(aes(fill = Host, alpha = ResCount)) +
  coord_fixed() + ggtitle(names(PredHostPolygons)[1]) #facet_wrap(~Host)

# Doing this to validate rather than predict ####

HostPrediction <- list()

a = 25

for(a in a:length(names(VirusAssocs))){
  
  print(names(VirusAssocs)[a])
  
  pHosts <- VirusAssocs[[a]]
  
  pHosts <- intersect(pHosts, AllMammals)
  
  if(length(pHosts)>0){
    
    EstList <- parallel::mclapply(1:length(AllSims), function(x){
      
      FocalNet <- AllSims[[x]] %>% as.matrix
      
      pHosts2 <- intersect(pHosts, rownames(FocalNet))
      
      ValidEst <- list()
      
      for(b in pHosts2){
        
        pHosts4 <- setdiff(pHosts2, b)
        
        pHosts3 <- setdiff(colnames(FocalNet), pHosts4)
        
        Estimates <- FocalNet[pHosts4, pHosts3]
        
        if(is.null(dim(Estimates))) Estimates <- rbind(Estimates, Estimates)/2
        
        Ests <- data.frame(Sp = names(sort(colSums(Estimates))),
                           Count = sort(colSums(Estimates))/nrow(Estimates),
                           Iteration = x) %>%
          mutate(Focal = ifelse(Sp==b, 1, 0))
        
        rownames(Ests) <- Ests$Sp
        
        ValidEst[[b]] <- Ests
        
      }
      
      ValidEst
      
    }, mc.cores = 20)
    
    HostPrediction[[names(VirusAssocs)[a]]] <- EstList
    
  } else HostPrediction[[names(VirusAssocs)[a]]] <- NA
}

# Trying it out ####

Valid <- HostPrediction %>% lapply(function(a){
  
  if(!is.null(names(a[[1]]))){
    
    b = map(names(a[[1]]), function(b) map(a, b) %>% bind_rows) %>% bind_rows
    
    c = b %>% group_by(Sp, Focal) %>% dplyr::summarise(Count = mean(Count)) %>% slice(order(Count, decreasing = T)) %>%
      mutate(Focal = factor(Focal))
    
  } else c = NA
  
  return(c)
})

save(Valid, file = "ModelValidation.Rdata")

ValidSummary <- data.frame(
  
  NHosts = sapply(VirusAssocs[1:length(HostPrediction)], length)
  
  
)

KeepPredictions <- (1:length(Valid))[-which(sapply(Valid, function(a) any(is.na(a))))]

KeepPredictions %>% 
  lapply(function(a) ggplot(Valid[[a]], aes(Focal, Count, colour = Focal)) + 
           ggforce::geom_sina() + theme(legend.position = "none") +
           ggtitle(names(VirusAssocs)[[a]])) %>% 
  arrange_ggplot2(ncol = 5)




