
# Importing Phylopics ####

library(rphylopic); library(tidyverse)

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order)

lapply(levels(Panth1$hOrder), function(a){
  
  Species <- Panth1 %>% filter(hOrder==a) %>% select(Sp) 
  
  Species$Sp <- Species$Sp %>% str_replace_all("_" , " ")
  
  lapply(Species$Sp, function(b){
    
    pic <- rphylopic::search_text(text = b, options = "names")
    
    output <- search_images(pic[1], options=c("pngFiles", "credit", "canonicalName"))
    
    pics <- image_data(output[4,], size = "256")
    
    save_png(pics, target = paste0("Figures/",#a,"/",
                                                 b,".png"))
    
  })
  
})



person <- name_search(text = "Homo sapiens", options = "namebankID")
pig <- name_search(text = "Sus scrofa", options = "namebankID")

img <- image_data("27356f15-3cf8-47e8-ab41-71c6260b2724", size = "512")

x <- search_text(text = "Homo sapiens", options = "names")
output <- search_images(x[1], options=c("pngFiles", "credit", "canonicalName"))
pic = image_data(output, size = "256")

BarGraph(Panth1, "hOrder", "AllPredDegree") + add_phylopic(pic[[1]], 1)

