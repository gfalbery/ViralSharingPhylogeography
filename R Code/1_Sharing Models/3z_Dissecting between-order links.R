
# Deconstructing between-order links ####
# Rscript "Dissecting between-order links.R"

library(igraph); library(tidyverse); library(ggregplot); library(parallel); library(SpRanger)

#source("R Code/00_Master Code.R")

load("Output Files/AllSums.Rdata")
load("Output Files/AllSims.Rdata")
load("Output Files/AllSimGs.Rdata")
load("Output Files/Panth1.Rdata")

Combos <- t(combn(levels(Panth1$hOrder),2)) %>% as.data.frame() %>%
  rename(Order1 = V1, Order2 = V2)

hComboList <- list()

for(i in 1:nrow(Combos)){
  
  print(i)
  
  a = Combos[i,] %>% unlist
  
  b1 = Panth1[Panth1$hOrder%in%a[1],"Sp"] %>% intersect(AllMammals)
  b2 = Panth1[Panth1$hOrder%in%a[2],"Sp"] %>% intersect(AllMammals)
  
  if(length(b1)>0&length(b2)>0){
    
    FocalNet <- AllSums[b1,b2] %>% as.matrix
    
    hComboList[[i]] <- data.frame(Degree = c(rowSums(FocalNet), colSums(FocalNet)),
                                  Iteration = i,
                                  Sp = c(b1, b2),
                                  Order = c(rep(a[1], length(b1)), rep(a[2], length(b2))),
                                  Group = rep(1:2, c(length(b1),length(b2))))
  }
}

hComboList <- hComboList %>% bind_rows()

hUniteList <- list()

for(i in levels(Panth1$hOrder)){
  
  print(i)
  
  b1 = Panth1[Panth1$hOrder == i, "Sp"] %>% intersect(AllMammals)
  
  if(length(b1)>0){
    FocalNet <- AllSums[b1,b1] %>% as.matrix
    
    hUniteList[[i]] <- data.frame(Degree = c(rowSums(FocalNet)),
                                  Iteration = i,
                                  Sp = b1,
                                  Order = i)
  }
}

hUniteList <- hUniteList %>% bind_rows()

OutDegrees <- hComboList %>% group_by(Sp) %>% summarise(OutDegree = sum(Degree))
InDegrees <- hUniteList %>% group_by(Sp) %>% summarise(InDegree = sum(Degree))
AllDegrees <- OutDegrees %>% left_join(InDegrees) %>%
  mutate(AllPredDegree = OutDegree + InDegree)

OrderLevelLinks <- Panth1 %>% #dplyr::select(-c("AllPredDegree", "InDegree", "OutDegree")) %>%
  left_join(AllDegrees, by = "Sp") %>%
  group_by(hOrder) %>%
  summarise(HostNumber = n(),
            OutDegree = mean(OutDegree, na.rm = T),
            InDegree = mean(InDegree, na.rm = T),
            AllPredDegree = mean(AllPredDegree, na.rm = T)) %>%
  gather(key = "Metric", value = "Degree", contains("Degree"))

OrderLevelLinks %>% ggplot(aes(log(HostNumber), log(Degree+1))) + geom_point() + 
  coord_fixed() + geom_smooth(method = lm) +
  facet_wrap(~Metric)

hComboList %>% group_by(Iteration, Group) %>%
  summarise(HostNo = n(),
            Degree = sum(Degree)) %>% spread(value = c("HostNo"), key = c("Group")) %>%
  group_by(Iteration) %>%
  summarise(`1` = sum(`1`, na.rm = T), 
            `2` = sum(`2`, na.rm = T),
            Degree = mean(Degree)) %>% bind_cols(Combos) %>%
  ggplot(aes(log(`1`) + log(`2`), log(Degree+1))) + 
  geom_point() + ggrepel::geom_label_repel(aes(label = paste(Order1,Order2))) #+ facet_wrap(~Order2)

hComboList %>% group_by(Iteration, Group) %>%
  summarise(HostNo = n(),
            Degree = sum(Degree)) %>% spread(value = c("HostNo"), key = c("Group")) %>%
  group_by(Iteration) %>%
  summarise(`1` = sum(`1`, na.rm = T), 
            `2` = sum(`2`, na.rm = T),
            Degree = mean(Degree)) %>% bind_cols(Combos) %>%
  ggplot(aes(log(`1`) + log(`2`), log(Degree+1))) + 
  geom_point() + 
  geom_abline(formula = y ~ x)

OrderPairs <- hComboList %>% group_by(Iteration, Group) %>%
  summarise(HostNo = n(),
            Degree = sum(Degree)) %>% spread(value = c("HostNo"), key = c("Group")) %>%
  group_by(Iteration) %>%
  summarise(`1` = sum(`1`, na.rm = T), 
            `2` = sum(`2`, na.rm = T),
            Degree = mean(Degree)) %>% bind_cols(Combos) %>%
  mutate(OrderSizes = log(`1`) + log(`2`))

OrderPairs %>%
  ggplot(aes(log(`1`) + log(`2`), log(Degree+1))) + 
  geom_point() + 
  geom_abline()

OrderPairs %>%
  ggplot(aes(log(`1`) + log(`2`), log(Degree+1), colour = MisspentEuth)) + 
  coord_fixed() +
  geom_point()

OrderPairs %>% filter(Bif1==1&BothEutherians==0)



