
# Deconstructing between-order links ####
# Rscript "Dissecting between-order links.R"

library(igraph); library(tidyverse); library(ggregplot); library(parallel); library(SpRanger)

#source("R Code/00_Master Code.R")

load("Output Files/AllSimGs.Rdata")
load("Output Files/hComboList.Rdata")
load("Output Files/Panth1.Rdata")

Combos <- t(combn(levels(Panth1$hOrder),2)) %>% as.data.frame() %>%
  rename(Order1 = V1, Order2 = V2)


hComboList <- mclapply(1:length(AllSimGs), function(i){
  
  a <- AllSimGs[[i]]
  
  lapply(t(Combos) %>% as.data.frame(), function(x){
    
    distances(a, v = V(a)[Panth1$hOrder == as.character(x[1])], 
              to = V(a)[Panth1$hOrder == as.character(x[2])])
    
  })
  
}, mc.cores = 10)


save(hComboList, file = "Output Files/hComboList.Rdata")

stop()

hComboList <- mclapply(1:length(AllSimGs), function(i){
  lapply(t(Combos) %>% as.data.frame(), function(x) induced_subgraph(AllSimGs[[i]], 
                                                                     as.character(Panth1$hOrder) == x))
}, mc.cores = 10)

Degrees <- lapply(hComboList, function(a) lapply(a, function(b) length(which(b==1))))

Degrees2 <- lapply(1:nrow(Combos), function(a) map(Degrees, a) %>% unlist)

MeanDegrees <- sapply(Degrees2, function(a) mean(colSums(a)))

DegreeScaling <- data.frame(
  
  Order1 = Combos$Order1,
  Order2 = Combos$Order2,
  
  BetweenDegrees = sapply(Degrees2, mean)
  
) %>%
  left_join(Panth1 %>% group_by(hOrder) %>%
              summarise(Number = n()),
            by = c(Order1 = "hOrder")) %>% 
  left_join(Panth1 %>% group_by(hOrder) %>%
              summarise(Number = n()),
            by = c(Order2 = "hOrder"), 
            suffix = c(".1", ".2"))

list(
  
  ggplot(DegreeScaling, aes(log(Number.1), log(BetweenDegrees+1), colour = log(Number.2))) + geom_point(),
  
  ggplot(DegreeScaling, aes(log(Number.2), log(BetweenDegrees+1), colour = log(Number.1))) + geom_point()
  
) %>% arrange_ggplot2

ggplot(DegreeScaling2, aes(log(Number.2) + log(Number.1), log(BetweenDegrees+1))) + geom_point()

DegreeScaling2 = bind_rows(DegreeScaling,
                           DegreeScaling %>% mutate(OrderA = Order2,
                                                    OrderB = Order1) %>%
                             mutate(Order1 = OrderA,
                                    Order2 = OrderB))

ggplot(DegreeScaling2, aes(log(Number.2) + log(Number.1), log(BetweenDegrees+1))) + 
  geom_point() + facet_wrap(~Order2)


