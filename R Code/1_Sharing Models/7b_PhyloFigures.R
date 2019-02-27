
# Importing Phylopics ####

install.packages("rphylopic")

library(rphylopic); library(tidyverse)

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order)

lapply(levels(Panth1$hOrder), function(a){
  
  Species <- Panth1 %>% filter(hOrder==a) %>% select(Species)
  
  lapply(Species, function(b){
    
    rphylopc::name_search(text = b, options = "namebankID")
    
  })
  
})
