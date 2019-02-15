
# Doing this to validate rather than predict ####

library(tidyverse); library(parallel)

HostValidation <- list()

a = 1

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
    
    HostValidation[[names(VirusAssocs)[a]]] <- EstList
    
  } else HostValidation[[names(VirusAssocs)[a]]] <- NA
}

# Trying it out ####

Valid <- HostValidation %>% lapply(function(a){
  
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


sapply(Valid[KeepPredictions], function(a) a %>% group_by(Focal) %>% dplyr::summarise(mean(rank(Count))) %>% filter(Focal==1))

FocalRank <- function(x){
  
  y <- x %>% filter(Focal == 1) %>% select(Count)
  z <- x %>% filter(Focal == 0) %>% select(Count)
  
  sapply(y, function(a) rank(c(a,z)))
  
}

