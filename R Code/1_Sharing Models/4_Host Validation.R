
# Doing this to validate rather than predict ####

# Rscript "R Code/1_Sharing Models/4_Host Validation.R"

library(tidyverse); library(parallel); library(ggregplot); library(ape); library(SpRanger); library(Matrix)

source("R Code/00_Master Code.R")

load("Output Files/AllSims.Rdata")
load("Output Files/AllSums.Rdata")

AllMammals <- intersect(colnames(FullSTMatrix),colnames(FullRangeAdj1))
AllMammals <- AllMammals[order(AllMammals)]

print("Start Validating!")

a = 1

GAMValidation <- Validate(VirusAssocs)

GAMValid <- GAMValidation %>% lapply(function(a){
  
  if(!is.null(names(a[[1]]))){
    
    b = a %>% bind_rows %>% as.data.frame()
    
    c = b %>% group_by(Sp, Focal) %>% 
      dplyr::summarise(Count = mean(Count)) %>% 
      slice(order(Count, decreasing = T)) %>%
      mutate(Focal = factor(Focal))
    
  } else c = NA
  
  return(c)
})

save(GAMValid, file = "Output Files/GAMValidation.Rdata")

{
  load("Output Files/ModelValidation.Rdata")
  load("Output Files/GAMValidation.Rdata")
  
  KeepPredictions <- (1:length(GAMValid))[-which(sapply(GAMValid, function(a) any(is.na(a))))]
  ValidPredictions <- (1:length(Valid))[-which(sapply(Valid, function(a) any(is.na(a))))]
  
  FocalRank <- function(x){
    
    y <- x[x[,"Focal"]==1,"Count"]
    z <- x[x[,"Focal"]==0,"Count"]
    
    (length(z$Count) + 2) - sapply(y$Count, function(a) rank(c(a,z$Count))[1])
    
  }
  
  GAMValidSummary <- data.frame(
    
    Virus = names(GAMValid)[KeepPredictions],
    
    NHosts = map(GAMValid[KeepPredictions], "Focal") %>% sapply(function(a) which(a=="1") %>% length),
    
    No = KeepPredictions,
    
    MeanRank = sapply(GAMValid[KeepPredictions], function(a) mean(FocalRank(a))),
    
    MeanCount = sapply(GAMValid[KeepPredictions], function(a) mean(a$Count)),
    MeanCount1 = sapply(GAMValid[KeepPredictions], function(a) mean(a[a$Focal==1,]$Count)),
    MeanCount0 = sapply(GAMValid[KeepPredictions], function(a) mean(a[a$Focal==0,]$Count))
    
  ) %>% slice(order(MeanRank)) %>%
    mutate(CountDiff = MeanCount1 - MeanCount0)
  
  ValidSummary <- data.frame(
    
    Virus = names(Valid)[ValidPredictions],
    
    NHosts = map(Valid[ValidPredictions], "Focal") %>% sapply(function(a) which(a=="1") %>% length),
    
    No = ValidPredictions,
    
    MeanRank = sapply(Valid[ValidPredictions], function(a) mean(FocalRank(a)))
    
  ) #%>% slice(order(MeanRank))
  
  CompSummary <- left_join(GAMValidSummary, ValidSummary[,c("Virus","MeanRank")], 
                           by = "Virus", suffix = c(".GAM",".GLM"))
  
  CompSummary %>% ggplot(aes(MeanRank.GAM, MeanRank.GLM)) + geom_point() + geom_smooth() + geom_abline() + coord_fixed()
  
  CompSum2 <- gather(CompSummary, key = "key", value = "value", MeanRank.GAM, MeanRank.GLM)
  
  BarGraph(CompSum2, "key", "value")
  
  CompSummary %>% summarise(mean(MeanRank.GAM), 
                            mean(MeanRank.GLM, na.rm = T),
                            median(MeanRank.GAM),
                            median(MeanRank.GLM, na.rm = T))
  
}

