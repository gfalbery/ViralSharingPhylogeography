
# Ecological trait data import ####

library(vegan)

EltonTraits <- read.delim("data/MamFuncDat.txt") %>% na.omit
EltonTraits$Scientific <- EltonTraits$Scientific %>% str_replace(" ","_")

EltonTraits$Carnivore <- ifelse(rowSums(EltonTraits[,c("Diet.Vend","Diet.Vect","Diet.Vfish")])>10,1,0)

DietComp <- EltonTraits %>% select(starts_with("Diet")) %>% select(1:10) %>% as_tibble
Remove <- which(rowSums(DietComp)==0)
DietComp <- DietComp %>% slice(-Remove)

VD <- vegdist(DietComp) %>% as.matrix

colnames(VD) <- rownames(VD) <- EltonTraits %>% slice(-Remove) %>% select(Scientific) %>% unlist

LongDiet <- reshape2::melt(VD) %>%     
  dplyr::rename(Sp = Var1, Sp2 = Var2, DietSim = value)
