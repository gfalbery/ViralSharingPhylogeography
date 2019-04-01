
# Doing this to validate rather than predict ####

# Rscript "R Code/1_Sharing Models/4_Host Validation.R"

library(tidyverse); library(parallel); library(ggregplot); library(ape); library(SpRanger); library(Matrix)

source("R Code/00_Master Code.R")

load("Output Files/AllSims.Rdata")
load("Output Files/AllSums.Rdata")

print("Start Validating!")

a = 1

GAMValid <- lapply(VirusAssocs, 
                   function(a){
                     NetworkValidate(a, AllSums)
                   })

save(GAMValid, file = "Output Files/GAMValidation.Rdata")

{
  
  load("Output Files/GAMValidation.Rdata")
  
  KeepPredictions <- (1:length(GAMValid))[-which(sapply(GAMValid, function(a) any(is.na(a))))]
  
  FocalRank <- function(x){
    
    y <- x[x[,"Focal"]=="Observed","Count"]
    z <- x[x[,"Focal"]=="Predicted","Count"]
    
    (length(z$Count) + 2) - sapply(y$Count, function(a) rank(c(a,z$Count))[1])
    
  }
  
  ValidSummary <- data.frame(
    
    Virus = names(GAMValid)[KeepPredictions],
    
    NHosts = map(GAMValid[KeepPredictions], "Focal") %>% sapply(function(a) which(a=="Observed") %>% length),
    
    No = KeepPredictions,
    
    MeanRank = sapply(GAMValid[KeepPredictions], function(a) mean(FocalRank(a))),
    
    MeanCount = sapply(GAMValid[KeepPredictions], function(a) mean(a$Count)),
    MeanCount1 = sapply(GAMValid[KeepPredictions], function(a) mean(a[a$Focal=="Observed",]$Count)),
    MeanCount0 = sapply(GAMValid[KeepPredictions], function(a) mean(a[a$Focal=="Predicted",]$Count))
    
  ) %>% slice(order(MeanRank)) %>%
    mutate(CountDiff = MeanCount1 - MeanCount0)
  
  ValidSummary %>% summarise(mean(MeanRank), 
                             median(MeanRank))
  
}

save(ValidSummary, file = "Output Files/ValidSummary.Rdata")

VirusCovar <- c("IsHoSa","IsHoSa.stringent",
                #"vGenomeMinLength","vGenomeMaxLength","vGenomeAveLength","vWOKcites","vPubMedCites",
                "vCytoReplicTF","vSegmentedTF",
                "vVectorYNna","vSSoDS","vDNAoRNA","vEnvelope",
                "IsZoonotic")

VirusTraits <- VirusTraits %>% 
  dplyr::rename(Virus = vVirusNameCorrected) %>%
  mutate(vCytoReplicTF = as.numeric(vCytoReplicTF),
         vSegmentedTF = as.numeric(vSegmentedTF))

ValidSummary <-  ValidSummary %>%
  left_join(VirusTraits, by = "Virus") %>%
  left_join(VirusHostRanges) %>%
  mutate(LogRank = log10(MeanRank),
         LogHosts = log10(NHosts))

Im1 <- INLAModelAdd("LogRank", 1, 
                    c(VirusCovar, "LogHosts", "HostRangeMean"), 
                    Random = "vFamily", "iid", "gaussian", 
                    ValidSummary[!NARows(ValidSummary, c(VirusCovar, "LogHosts", "HostRangeMean", "LogRank")),])



