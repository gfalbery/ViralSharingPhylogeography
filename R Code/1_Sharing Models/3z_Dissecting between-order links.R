
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
  
  b1 = Panth1$hOrder%in%a[1]
  b2 = Panth1$hOrder%in%a[2]
  
  FocalNet <- (AllSums[b1,b2]/length(AllSims)) %>% as.matrix
  
  hComboList[[i]] <- data.frame(Degree = c(rowSums(FocalNet), colSums(FocalNet)),
                                Iteration = i,
                                Sp = c(Panth1[b1,"Sp"], Panth1[b2,"Sp"]),
                                Order = c(Panth1[b1,"hOrder"], Panth1[b2,"hOrder"]),
                                Group = rep(1:2, c(length(b1[b1==T]),length(b2[b2==T]))))
  
}

hComboList <- hComboList %>% bind_rows()

hUniteList <- list()

for(i in levels(Panth1$hOrder)){
  
  print(i)
  
  b1 = Panth1$hOrder == i
  
  FocalNet <- (AllSums[b1,b1]/length(AllSims)) %>% as.matrix
  
  hUniteList[[i]] <- data.frame(Degree = c(rowSums(FocalNet)),
                                Iteration = i,
                                Sp = c(Panth1[b1,"Sp"]),
                                Order = c(Panth1[b1,"hOrder"]))
  
}

hUniteList <- hUniteList %>% bind_rows()

OutDegrees <- hComboList %>% group_by(Sp) %>% summarise(OutDegree = sum(Degree))
InDegrees <- hUniteList %>% group_by(Sp) %>% summarise(InDegree = sum(Degree))
AllDegrees <- OutDegrees %>% left_join(InDegrees) %>%
  mutate(AllPredDegree = OutDegree + InDegree)

OrderLevelLinks <- Panth1 %>% dplyr::select(-c("AllPredDegree", "InDegree", "OutDegree")) %>%
  left_join(AllDegrees, by = "Sp") %>%
  group_by(hOrder) %>%
  summarise(HostNumber = n(),
            OutDegree = mean(OutDegree),
            InDegree = mean(InDegree),
            AllPredDegree = mean(AllPredDegree)) %>%
  gather(key = "Metric", value = "Degree", contains("Degree"))

OrderLevelLinks %>% ggplot(aes(log(HostNumber), log(Degree+1))) + geom_point() + 
  coord_fixed() + geom_smooth() +
  facet_wrap(~Metric)

hComboList %>% group_by(Iteration, Group) %>%
  summarise(HostNo = n(),
            Degree = sum(Degree)) %>% spread(value = c("HostNo"), key = c("Group")) %>%
  group_by(Iteration) %>%
  summarise(`1` = sum(`1`, na.rm = T), 
            `2` = sum(`2`, na.rm = T),
            Degree = mean(Degree)) %>% bind_cols(Combos) %>%
  ggplot(aes(log(`1`) + log(`2`), log(Degree+1))) + 
  geom_point() + facet_wrap(~Order2)


