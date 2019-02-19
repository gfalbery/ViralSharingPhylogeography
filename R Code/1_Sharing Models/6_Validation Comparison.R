
# Comparing outputs for general and specific viruses ####

library(tidyverse)

VirusTypeList <- list(which(Viruses$vVectorYNna=="Y"&Viruses$vDNAoRNA=="RNA"),
                      which(Viruses$vVectorYNna=="N"&Viruses$vDNAoRNA=="RNA"),
                      which(Viruses$vDNAoRNA=="DNA"))

load("ModelValidation.Rdata")
VectorValid <- readRDS("VectorModelValidation.rds")
NVectorValid <- readRDS("NVectorModelValidation.rds")
DNAValid <- readRDS("DNAModelValidation.rds")

KeepPredictions2 <- list((1:length(VectorValid))[-which(sapply(VectorValid, function(a) any(is.na(a))))],
                        (1:length(NVectorValid))[-which(sapply(NVectorValid, function(a) any(is.na(a))))],
                        (1:length(DNAValid))[-which(sapply(DNAValid, function(a) any(is.na(a))))])

VectorValidSummary <- data.frame(
  
  Virus = names(VirusAssocs)[VirusTypeList[[1]]][KeepPredictions2[[1]]],
  
  NHosts = map(VectorValid[KeepPredictions2[[1]]], "Focal") %>% 
    sapply(function(a) which(a=="1") %>% length),
  
  No = KeepPredictions2[[1]],
  
  MeanRank = sapply(VectorValid[KeepPredictions2[[1]]], function(a) mean(FocalRank(a)))
  
) %>% slice(order(MeanRank)) %>% merge(ValidSummary, by = "Virus", suffixes = c(".Spec",".Gen"), all.x = T)

NVectorValidSummary <- data.frame(
  
  Virus = names(VirusAssocs)[VirusTypeList[[2]]][KeepPredictions2[[2]]],
  
  NHosts = map(NVectorValid[KeepPredictions2[[2]]], "Focal") %>% sapply(function(a) which(a=="1") %>% length),
  
  No = KeepPredictions2[[2]],
  
  MeanRank = sapply(NVectorValid[KeepPredictions2[[2]]], function(a) mean(FocalRank(a)))
  
) %>% slice(order(MeanRank)) %>% merge(ValidSummary, by = "Virus", suffixes = c(".Spec",".Gen"), all.x = T)


DNAValidSummary <- data.frame(
  
  Virus = names(VirusAssocs)[VirusTypeList[[3]]][KeepPredictions2[[3]]],
  
  NHosts = map(DNAValid[KeepPredictions2[[3]]], "Focal") %>% sapply(function(a) which(a=="1") %>% length),
  
  No = KeepPredictions2[[3]],
  
  MeanRank = sapply(DNAValid[KeepPredictions2[[3]]], function(a) mean(FocalRank(a)))
  
) %>% slice(order(MeanRank))%>% merge(ValidSummary, by = "Virus", suffixes = c(".Spec",".Gen"), all.x = T)







