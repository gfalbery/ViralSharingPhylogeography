
# Comparing outputs for general and specific viruses ####

library(tidyverse)

SubVirusTraits <- VirusTraits %>% filter(vVirusNameCorrected %in% names(VirusAssocs))

VirusTypeList <- list(which(SubVirusTraits$vDNAoRNA=="RNA"),
                      which(SubVirusTraits$vVectorYNna=="Y"&SubVirusTraits$vDNAoRNA=="RNA"),
                      which(SubVirusTraits$vVectorYNna=="N"&SubVirusTraits$vDNAoRNA=="RNA"),
                      which(SubVirusTraits$vDNAoRNA=="DNA"))

Resps <- c("VirusBinary","RNA","DNA","Vector","NVector")

load("Output Files/SubSums.Rdata")

SubGroupValidList <- lapply(1:4, function(a) Validate(VirusAssocs[VirusTypeList[[a]]], SubSums[[a]]))

SubGroupValid <- lapply(SubGroupValidList, function(a){
  
  lapply(a, function(b){
    
    if(!is.null(names(b[[1]]))){
      
      c = b %>% bind_rows %>% as.data.frame()
      
      d = c %>% group_by(Sp, Focal) %>% 
        dplyr::summarise(Count = mean(Count)) %>% 
        slice(order(Count, decreasing = T)) %>%
        mutate(Focal = factor(Focal))
      
    } else d = NA
    
    return(d)
  })
  
})

SubGroupKeeps <- SubGroupValid %>% lapply(function(a){
  
  (1:length(a))[-which(sapply(a, function(b) any(is.na(b))))]
  
})

FocalRank <- function(x){
  
  y <- x[x[,"Focal"]==1,"Count"]
  z <- x[x[,"Focal"]==0,"Count"]
  
  (length(z$Count) + 2) - sapply(y$Count, function(a) rank(c(a,z$Count))[1])
  
}

Valids <- lapply(1:4, function(a){
  
  data.frame(
    
    Virus = names(VirusAssocs)[SubGroupKeeps[[a]]],
    
    NHosts = map(SubGroupValid[[a]][SubGroupKeeps[[a]]], "Focal") %>% 
      sapply(function(a) which(a=="1") %>% length),
    
    No = SubGroupKeeps[[a]],
    
    MeanRank = sapply(SubGroupValid[[a]][SubGroupKeeps[[a]]], function(a) mean(FocalRank(a))),
    
    Resp = Resps[(a+1)]
    
  ) %>% slice(order(MeanRank)) 
  
}) %>% bind_rows()



