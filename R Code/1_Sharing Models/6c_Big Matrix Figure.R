
# Rscript "R Code/1_Sharing Models/6c_Big Matrix Figure.R"

source("R Code/00_Master Code.R")

library(geiger);library(ape);library(picante)

STFull <- read.nexus("data/ele_1307_sm_sa1.tre")[[1]]

library(tidyverse); library(colorspace); library(ggregplot); library(RColorBrewer)

load("Output Files/Panth1.Rdata")
load("Output Files/AllSums.Rdata")

OrderSizeOrder <- Panth1 %>% group_by(hOrder) %>% summarise(A = mean(AllPredDegree)) %>% slice(order(A))

OrderHostOrder <- (Panth1 %>% mutate(hOrder = factor(hOrder, levels = OrderSizeOrder$hOrder)) %>% slice(order(hOrder, Sp)) %>% select(Sp))$Sp %>% unique

m1 = AllSums %>% as.matrix %>% reshape2::melt() %>%
  dplyr::rename(Sp = Var2, Sp2 = Var1, Links = value) %>%
  mutate(Sp = factor(Sp, levels = OrderHostOrder),
         Sp2 = factor(Sp2, levels = OrderHostOrder))

jpeg("SIFigures/BigMatrixPlot.jpeg", units = "mm", width = 350, height = 350, res = 600)

ggplot(m1, aes(Sp, Sp2)) + geom_tile(aes(fill = log10(Links+1))) +
  scale_x_discrete(labels = NULL) + scale_y_discrete(labels = NULL) +
  scale_fill_continuous_sequential(palette = AlberPalettes[2]) +
  coord_fixed() +
  ggtitle("Predicted network links")

dev.off()

PhyloHostOrder <- STFull$tip.label

m1 = AllSums %>% as.matrix %>% reshape2::melt() %>%
  dplyr::rename(Sp = Var2, Sp2 = Var1, Links = value) %>%
  mutate(Sp = factor(Sp, levels = PhyloHostOrder),
         Sp2 = factor(Sp2, levels = PhyloHostOrder))

jpeg("SIFigures/PhyloBigMatrixPlot.jpeg", units = "mm", width = 350, height = 350, res = 600)

ggplot(m1, aes(Sp, Sp2)) + geom_tile(aes(fill = log10(Links+1))) +
  scale_x_discrete(labels = NULL) + scale_y_discrete(labels = NULL) +
  scale_fill_continuous_sequential(palette = AlberPalettes[2]) +
  coord_fixed() +
  ggtitle("Predicted network links")

dev.off()



IntHN = intersect(OrderHostOrder, FHN)

jpeg("SIFigures/SmallMatrixPlot.jpeg", units = "mm", width = 350, height = 350, res = 600)

ggplot(m1[FHN,FHN], aes(Sp, Sp2)) + geom_tile(aes(fill = log10(Links+1))) +
  scale_x_discrete(labels = NULL) + scale_y_discrete(labels = NULL) +
  scale_fill_continuous_sequential(palette = AlberPalettes[2]) +
  coord_fixed() +
  ggtitle("Observed network predicted links")

dev.off()

MLong = HostAdj[IntHN,IntHN] %>% as.matrix %>% reshape2::melt() %>%
  dplyr::rename(Sp = Var2, Sp2 = Var1, Links = value) %>%
  mutate(Sp = factor(Sp, levels = IntHN),
         Sp2 = factor(Sp2, levels = IntHN))

jpeg("SIFigures/ObsMatrixPlot.jpeg", units = "mm", width = 350, height = 350, res = 600)

ggplot(MLong, aes(Sp, Sp2)) + geom_tile(aes(fill = log10(Links+1))) +
  scale_x_discrete(labels = NULL) + scale_y_discrete(labels = NULL) +
  scale_fill_continuous_sequential(palette = AlberPalettes[2]) +
  coord_fixed() +
  ggtitle("Observed network observed links")

dev.off()
