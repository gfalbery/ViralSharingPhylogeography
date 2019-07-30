
# 5_Host Validation Procedure ####

library(tidyverse); library(parallel); library(ggregplot); library(ape); library(SpRanger); library(Matrix)

# load("Output Files/AllSims.Rdata")
load("Output Files/AllSums.Rdata")

print("Start Validating!")

a = 1

GAMValid <- lapply(VirusAssocs, 
                   function(a){
                     NetworkValidate(a, AllSums, colSums)
                   })

save(GAMValid, file = "Output Files/GAMValid.Rdata")

{
  
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

VirusCovar <- c("IsHoSa", # "IsHoSa.stringent",
                #"vGenomeMinLength","vGenomeMaxLength","vGenomeAveLength","vWOKcites","vPubMedCites",
                "vCytoReplicTF","vSegmentedTF",
                "vVectorYNna","vSSoDS","vDNAoRNA","vEnvelope",
                "IsZoonotic")

VirusTraits <- VirusTraits %>% 
  dplyr::rename(Virus = vVirusNameCorrected) %>%
  mutate(vCytoReplicTF = as.numeric(vCytoReplicTF),
         vSegmentedTF = as.numeric(vSegmentedTF))

ValidDF <-  ValidSummary %>%
  left_join(VirusTraits, by = "Virus") %>%
  left_join(VirusHostRanges) %>%
  mutate(LogRank = log10(MeanRank),
         LogHosts = log10(NHosts)) %>%
  # mutate_if(is.numeric, function(a) scale(a)) %>%
  select(-c(vSubfamily, vGenus, vIsTypeSpecies, vICTVnumber, vGenomeMinLength, vGenomeMaxLength, vGenomeAveLength)) %>%
  na.omit()

Im1 <- INLAModelAdd("LogRank", 1, 
                    c(VirusCovar, "LogHosts", "HostRangeMean"), 
                    Random = "vFamily", "iid", "gaussian", 
                    ValidDF)

list(Im1$AllModels[[1]], Im1$AllModels[[2]]$HostRangeMean, Im1$AllModels[[3]]$vVectorYNna, Im1$AllModels[[4]]$LogHosts) %>% Efxplot

Efxplot(ModelList = Im1["FinalModel"])

LM1 <- lme4::lmer(LogRank ~ HostRangeMean + vVectorYNna + LogHosts + (1|vFamily),
            data = ValidDF)

R2 <- r2glmm::r2beta(LM1, method = 'sgv')
R2
  
# Null predictions using only space and phylogeny ####

GAMValidSpace <- lapply(VirusAssocs, 
                        function(a){
                          NetworkValidate(a, 
                                          FullRangeAdj[rownames(AllSums), 
                                                       rownames(AllSums)])
                        })

GAMValidPhylo <- lapply(VirusAssocs, 
                        function(a){
                          NetworkValidate(a, 
                                          tFullSTMatrix[rownames(AllSums), 
                                                        rownames(AllSums)])
                        })

SpaceValidSummary <- data.frame(
  
  Virus = names(GAMValidSpace)[KeepPredictions],
  
  NHosts = map(GAMValidSpace[KeepPredictions], "Focal") %>% sapply(function(a) which(a=="Observed") %>% length),
  
  No = KeepPredictions,
  
  MeanRank = sapply(GAMValidSpace[KeepPredictions], function(a) mean(FocalRank(a))),
  
  MeanCount = sapply(GAMValidSpace[KeepPredictions], function(a) mean(a$Count)),
  MeanCount1 = sapply(GAMValidSpace[KeepPredictions], function(a) mean(a[a$Focal=="Observed",]$Count)),
  MeanCount0 = sapply(GAMValidSpace[KeepPredictions], function(a) mean(a[a$Focal=="Predicted",]$Count))
  
) %>% slice(order(MeanRank)) %>%
  mutate(CountDiff = MeanCount1 - MeanCount0)

SpaceValidSummary %>% summarise(mean(MeanRank), 
                                median(MeanRank))

PhyloValidSummary <- data.frame(
  
  Virus = names(GAMValidPhylo)[KeepPredictions],
  
  NHosts = map(GAMValidPhylo[KeepPredictions], "Focal") %>% sapply(function(a) which(a=="Observed") %>% length),
  
  No = KeepPredictions,
  
  MeanRank = sapply(GAMValidPhylo[KeepPredictions], function(a) mean(FocalRank(a))),
  
  MeanCount = sapply(GAMValidPhylo[KeepPredictions], function(a) mean(a$Count)),
  MeanCount1 = sapply(GAMValidPhylo[KeepPredictions], function(a) mean(a[a$Focal=="Observed",]$Count)),
  MeanCount0 = sapply(GAMValidPhylo[KeepPredictions], function(a) mean(a[a$Focal=="Predicted",]$Count))
  
) %>% slice(order(MeanRank)) %>%
  mutate(CountDiff = MeanCount1 - MeanCount0)

PhyloValidSummary %>% summarise(mean(MeanRank), 
                                median(MeanRank))

# Trying to use maximum probability ####

SMax <- function(a) apply(a, 2, max)

GAMValidMax <- lapply(VirusAssocs, 
                      function(a){
                        NetworkValidate(a, AllSums, Fun = SMax)
                      })

{
  
  KeepPredictions <- (1:length(GAMValid))[-which(sapply(GAMValid, function(a) any(is.na(a))))]
  
  FocalRank <- function(x){
    
    y <- x[x[,"Focal"]=="Observed","Count"]
    z <- x[x[,"Focal"]=="Predicted","Count"]
    
    (length(z$Count) + 2) - sapply(y$Count, function(a) rank(c(a,z$Count))[1])
    
  }
  
  ValidSummaryMax <- data.frame(
    
    Virus = names(GAMValidMax)[KeepPredictions],
    
    NHosts = map(GAMValidMax[KeepPredictions], "Focal") %>% sapply(function(a) which(a=="Observed") %>% length),
    
    No = KeepPredictions,
    
    MeanRank = sapply(GAMValidMax[KeepPredictions], function(a) mean(FocalRank(a))),
    
    MeanCount = sapply(GAMValidMax[KeepPredictions], function(a) mean(a$Count)),
    MeanCount1 = sapply(GAMValidMax[KeepPredictions], function(a) mean(a[a$Focal=="Observed",]$Count)),
    MeanCount0 = sapply(GAMValidMax[KeepPredictions], function(a) mean(a[a$Focal=="Predicted",]$Count))
    
  ) %>% slice(order(MeanRank)) %>%
    mutate(CountDiff = MeanCount1 - MeanCount0)
  
  ValidSummaryMax %>% summarise(mean(MeanRank), 
                                median(MeanRank))
  
}

# Doing it with EID

EIDVirusAssocs <- lapply(unique(AssocsEID$Cargo), function(a){
  
  unique(AssocsEID[AssocsEID$Cargo==a,"Carrier"])
  
})

names(EIDVirusAssocs) <- unique(AssocsEID$Cargo)

GAMValidEID <- lapply(EIDVirusAssocs, 
                      function(a){
                        NetworkValidate(a, AllSums)
                      })

names(GAMValidEID) <- levels(AssocsEID$Cargo)

{
  
  KeepPredictionsEID <- (1:length(GAMValidEID))[-which(sapply(GAMValidEID, function(a) any(is.na(a))))]
  
  FocalRank <- function(x){
    
    y <- x[x[,"Focal"]=="Observed","Count"]
    z <- x[x[,"Focal"]=="Predicted","Count"]
    
    (length(z$Count) + 2) - sapply(y$Count, function(a) rank(c(a,z$Count))[1])
    
  }
  
  ValidEIDSummary <- data.frame(
    
    Virus = names(GAMValidEID)[KeepPredictionsEID],
    
    NHosts = map(GAMValidEID[KeepPredictionsEID], "Focal") %>% sapply(function(a) which(a=="Observed") %>% length),
    
    No = KeepPredictionsEID,
    
    MeanRank = sapply(GAMValidEID[KeepPredictionsEID], function(a) mean(FocalRank(a))),
    
    MeanCount = sapply(GAMValidEID[KeepPredictionsEID], function(a) mean(a$Count)),
    MeanCount1 = sapply(GAMValidEID[KeepPredictionsEID], function(a) mean(a[a$Focal=="Observed",]$Count)),
    MeanCount0 = sapply(GAMValidEID[KeepPredictionsEID], function(a) mean(a[a$Focal=="Predicted",]$Count))
    
  ) %>% slice(order(MeanRank)) %>%
    mutate(CountDiff = MeanCount1 - MeanCount0)
  
  ValidEIDSummary %>% summarise(mean(MeanRank), 
                                median(MeanRank))
  
}


GAMValidEIDSpace <- lapply(EIDVirusAssocs, 
                           function(a){
                             NetworkValidate(a, FullRangeAdj)
                           })

names(GAMValidEIDSpace) <- levels(AssocsEID$Cargo)

GAMValidEIDPhylo <- lapply(EIDVirusAssocs, 
                           function(a){
                             NetworkValidate(a, tFullSTMatrix)
                           })

names(GAMValidEIDPhylo) <- levels(AssocsEID$Cargo)


{
  
  ValidEIDSummarySpace <- data.frame(
    
    Virus = names(GAMValidEIDSpace)[KeepPredictionsEID],
    
    NHosts = map(GAMValidEIDSpace[KeepPredictionsEID], "Focal") %>% sapply(function(a) which(a=="Observed") %>% length),
    
    No = KeepPredictionsEID,
    
    MeanRank = sapply(GAMValidEIDSpace[KeepPredictionsEID], function(a) mean(FocalRank(a))),
    
    MeanCount = sapply(GAMValidEIDSpace[KeepPredictionsEID], function(a) mean(a$Count)),
    MeanCount1 = sapply(GAMValidEIDSpace[KeepPredictionsEID], function(a) mean(a[a$Focal=="Observed",]$Count)),
    MeanCount0 = sapply(GAMValidEIDSpace[KeepPredictionsEID], function(a) mean(a[a$Focal=="Predicted",]$Count))
    
  ) %>% slice(order(MeanRank)) %>%
    mutate(CountDiff = MeanCount1 - MeanCount0)
  
  ValidEIDSummarySpace %>% summarise(mean(MeanRank), 
                                     median(MeanRank))
  
}

{
  
  ValidEIDSummaryPhylo <- data.frame(
    
    Virus = names(GAMValidEIDPhylo)[KeepPredictionsEID],
    
    NHosts = map(GAMValidEIDPhylo[KeepPredictionsEID], "Focal") %>% sapply(function(a) which(a=="Observed") %>% length),
    
    No = KeepPredictionsEID,
    
    MeanRank = sapply(GAMValidEIDPhylo[KeepPredictionsEID], function(a) mean(FocalRank(a))),
    
    MeanCount = sapply(GAMValidEIDPhylo[KeepPredictionsEID], function(a) mean(a$Count)),
    MeanCount1 = sapply(GAMValidEIDPhylo[KeepPredictionsEID], function(a) mean(a[a$Focal=="Observed",]$Count)),
    MeanCount0 = sapply(GAMValidEIDPhylo[KeepPredictionsEID], function(a) mean(a[a$Focal=="Predicted",]$Count))
    
  ) %>% slice(order(MeanRank)) %>%
    mutate(CountDiff = MeanCount1 - MeanCount0)
  
  ValidEIDSummaryPhylo %>% summarise(mean(MeanRank), 
                                     median(MeanRank))
  
}

# Looking at predicted hosts of viruses in EID2 #### 

BothViruses <- intersect(AssocsEID$Cargo, AssocsBase$Virus)

EIDVirusAssocs

BothAssocs <- merge.list(VirusAssocs, EIDVirusAssocs)

EIDSubAssocs <- EIDVirusAssocs[BothViruses]
HP3SubAssocs <- VirusAssocs[BothViruses]

EIDConfirm <- lapply(BothViruses, function(a){
  
  if(!is.null(VirusAssocs[[a]])){
    
    setdiff(
      EIDVirusAssocs[[a]],
      VirusAssocs[[a]]
    )
  }
})

names(EIDConfirm) <- BothViruses

EIDPredList <- lapply(BothViruses, function(a){
  
  print(a)
  
  pHosts1 <- VirusAssocs[[a]] %>% intersect(AllMammals)
  pHosts2 <- setdiff(AllMammals, pHosts1)
  pHosts3 <- EIDConfirm[[a]] %>% intersect(AllMammals)
  
  if(!is.null(pHosts3)&!is.null(pHosts1)&length(pHosts1)>0){
    
    lapply(pHosts3, function(i){
      
      Subdf <- AllSums[pHosts1,pHosts2] 
      
      if(is.null(dim(Subdf))){
        Subdf <- t(data.frame(Subdf))
      }
      
      F1 <- Subdf %>% colSums()
      F1 <- sort(F1/length(pHosts1), decreasing = T)
      
      which(names(F1)==i) %>% return
      
    }) #%>% median %>% print
    
  }
  
})

# Trying out viral traits' impact on predictability ####

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

Im1 <- INLAModelAdd("LogRank", 1, c(VirusCovar, "LogHosts", "HostRangeMean"), 
                    Random = "vFamily", "iid", "gaussian", 
                    ValidSummary[!NARows(ValidSummary, c(VirusCovar, "LogHosts", "HostRangeMean", "LogRank")),])



