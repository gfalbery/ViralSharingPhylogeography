
# Importing Phylopics ####

library(rphylopic); library(tidyverse)

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order)

lapply(levels(Panth1$hOrder), function(a){
  
  Species <- Panth1 %>% filter(hOrder==a) %>% select(Sp)
  
  lapply(Species$Sp, function(b){
    
    pic <- rphylopic::name_search(text = as.character(b), options = "namebankID") %>% return
    
    if(length(pic)>0) return(pic)
    
  })
  
})
