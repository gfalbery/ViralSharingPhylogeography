
# Doing this to validate rather than predict ####

# Rscript "R Code/1_Sharing Models/5d_SubModel Validation.R"

library(tidyverse); library(parallel); library(ggregplot)

source("R Code/00_Master Code.R")

load("SubModelGenData.Rdata")

AllMammals <- intersect(colnames(FullSTMatrix),colnames(FullRangeAdj1))
AllMammals <- AllMammals[order(AllMammals)]

SubResps <- c(#"RNA",
              "Vector","NVector","DNA")

VirusTypeList <- list(which(Viruses$vVectorYNna=="Y"&Viruses$vDNAoRNA=="RNA"),
                      which(Viruses$vVectorYNna=="N"&Viruses$vDNAoRNA=="RNA"),
                      which(Viruses$vDNAoRNA=="DNA"))

for(r in 1:length(SubResps)){
  
  print("Start Validating!")
  
  a = 1
  
  HostValidation <- list()
  
  SubVirusAssocs <- VirusAssocs[VirusTypeList[[r]]]
  
  if(file.exists(paste0(SubResps[r],"ModelValidation.rds"))) Valid <- readRDS(paste0(SubResps[r],"ModelValidation.rds"))
  
  for(a in a:(length(SubVirusAssocs))){
    
    print(names(SubVirusAssocs)[a])
    
    pHosts <- SubVirusAssocs[[a]]
    
    pHosts <- intersect(pHosts, AllMammals)
    
    if(length(pHosts)>0){
      
      EstList <- parallel::mclapply(1:length(SubSims), function(x){
        
        FocalNet <- SubSims[[r]][[x]] %>% as.matrix
        
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
        
      }, mc.cores = 40)
      
      HostValidation[[names(SubVirusAssocs)[a]]] <- EstList
      
    } else HostValidation[[names(SubVirusAssocs)[a]]] <- NA
    
    if(a %% 10 == 0){
      
      Valid <- HostValidation %>% lapply(function(a){
        
        if(!is.null(names(a[[1]]))){
          
          b = map(names(a[[1]]), function(b) map(a, b) %>% bind_rows) %>% bind_rows
          
          c = b %>% group_by(Sp, Focal) %>% dplyr::summarise(Count = mean(Count)) %>% 
            slice(order(Count, decreasing = T)) %>%
            mutate(Focal = factor(Focal))
          
        } else c = NA
        
        return(c)
      })
      
      saveRDS(Valid, file = paste0(SubResps[r],"ModelValidation.rds"))
      
    }
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
  
  saveRDS(Valid, file = paste0(SubResps[r],"ModelValidation.rds"))
  
}
