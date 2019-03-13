
# Doing this to validate rather than predict ####

# Rscript "R Code/1_Sharing Models/4_Host Validation.R"

library(tidyverse); library(parallel); library(ggregplot); library(ape); library(SpRanger); library(Matrix)

# source("R Code/00_Master Code.R")

load("Output Files/AllSims.Rdata")
load("Output Files/AllSums.Rdata")

print("Start Validating!")

a = 1

GAMValidation <- Validate(VirusAssocs, AllSums)

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
  
  # load("Output Files/GAMValidation.Rdata")
  
  KeepPredictions <- (1:length(GAMValid))[-which(sapply(GAMValid, function(a) any(is.na(a))))]
  
  FocalRank <- function(x){
    
    y <- x[x[,"Focal"]==1,"Count"]
    z <- x[x[,"Focal"]==0,"Count"]
    
    (length(z$Count) + 2) - sapply(y$Count, function(a) rank(c(a,z$Count))[1])
    
  }
  
  ValidSummary <- data.frame(
    
    Virus = names(GAMValid)[KeepPredictions],
    
    NHosts = map(GAMValid[KeepPredictions], "Focal") %>% sapply(function(a) which(a=="1") %>% length),
    
    No = KeepPredictions,
    
    MeanRank = sapply(GAMValid[KeepPredictions], function(a) mean(FocalRank(a))),
    
    MeanCount = sapply(GAMValid[KeepPredictions], function(a) mean(a$Count)),
    MeanCount1 = sapply(GAMValid[KeepPredictions], function(a) mean(a[a$Focal==1,]$Count)),
    MeanCount0 = sapply(GAMValid[KeepPredictions], function(a) mean(a[a$Focal==0,]$Count))
    
  ) %>% slice(order(MeanRank)) %>%
    mutate(CountDiff = MeanCount1 - MeanCount0)
  
  ValidSummary %>% summarise(mean(MeanRank), 
                             median(MeanRank))
  
}

save(ValidSummary, file = "Output Files/ValidSummary.Rdata")

VirusCovar <- c("IsHoSa","IsHoSa.stringent",
                #"vGenomeMinLength","vGenomeMaxLength","vGenomeAveLength","vWOKcites","vPubMedCites",
                "vCytoReplicTF","vSegmentedTF","vVectorYNna","vSSoDS","vDNAoRNA","vEnvelope",
                "IsZoonotic")

ValidSummary <- VirusTraits %>% dplyr::select(vVirusNameCorrected, VirusCovar) %>%
  dplyr::rename(Virus = vVirusNameCorrected) %>%
  right_join(ValidSummary, by = "Virus")

VirusCovar %>% lapply(function(a){
  
  BarGraph(ValidSummary, "vDNAoRNA", "MeanRank", a, text = "N") + theme(legend.position = "none")
  
}) %>% arrange_ggplot2

