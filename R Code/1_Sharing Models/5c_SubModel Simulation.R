
# Looping STAN Models ####

# Rscript "R Code/1_Sharing Models/5c_SubModel Simulation.R"

library(parallel); library(tidyverse); library(MCMCglmm)

source("R Code/00_Master Code.R")

load("~/Albersnet/RNAModelOutput.Rdata")
load("~/Albersnet/NVectorModelOutput.Rdata")
load("~/Albersnet/VectorModelOutput.Rdata")
load("~/Albersnet/DNAModelOutput.Rdata")

invlogit <- function(a) exp(a)/(1 + exp(a))
logit <- function(a) log(a/(1-a))

tFullSTMatrix <- 1 - (FullSTMatrix - min(FullSTMatrix))/max(FullSTMatrix)

AllMammals <- intersect(colnames(FullSTMatrix),colnames(FullRangeAdj1))
AllMammals <- AllMammals[order(AllMammals)]

AllMammalMatrix <- data.frame(
  Sp = as.character(rep(AllMammals,each = length(AllMammals))),
  Sp2 = as.character(rep(AllMammals,length(AllMammals))),
  Space = c(FullRangeAdj1[AllMammals,AllMammals]),
  Phylo = c(tFullSTMatrix[AllMammals,AllMammals])
)

UpperMammals <- which(upper.tri(FullSTMatrix[AllMammals,AllMammals], diag = T))

AllMammaldf <- AllMammalMatrix[-UpperMammals,]

SubResps <- c(#"RNA",
           "Vector","NVector","DNA")

MX1 = model.matrix( ~ Space + Phylo + Space:Phylo, data = AllMammaldf) %>%
  as.matrix %>% as("dgCMatrix")

N = nrow(AllMammaldf)

SubModels <- list(#RNA = q,
                  Vector = s,
                  NVector = t,
                  DNA = r)

SolList <- map(SubModels, "df")

Estimates <- lapply(SolList, function(a){
  
  Samples <- sample(1:nrow(a), 1000)
  
  apply(a, 2, function(b){
    
    data.frame(Estimate = posterior.mode(as.mcmc(b[Samples])),
               Lower = HPDinterval(as.mcmc(b[Samples]))[1],
               Upper = HPDinterval(as.mcmc(b[Samples]))[2])
    
  }) %>% bind_rows()
  
}) 

SubDFList <- SubPredList <- SubPredDF <- SubSims <- list()

for(z in 1:length(SubResps)){
  
  print(SubResps[z])
  
  # Trying it without random effects ####
  
  ClusterMCMC <- SolList[[z]]
  
  RowsSampled <- sample(1:nrow(ClusterMCMC), 100, replace = F)
  
  XBetas <- c("mu_alpha","beta_space","beta_phylo","beta_inter")
  
  SubPredList[[z]] <- mclapply(1:length(RowsSampled), function(x){
    
    RowSampled <- RowsSampled[x]
    
    XFX <- ClusterMCMC[RowSampled, XBetas] %>% unlist
    XPredictions <- c(XFX %*% t(MX1))
    
    ZPredictions2a <- rnorm(n = N, mean = mean(log(FinalHostMatrix$hDiseaseZACites+1))*ClusterMCMC[RowSampled,"beta_d_cites_s"], sd = ClusterMCMC[RowsSampled[x], "sigma"])
    ZPredictions2b <- rnorm(n = N, mean = mean(log(FinalHostMatrix$hDiseaseZACites+1))*ClusterMCMC[RowSampled,"beta_d_cites_s"], sd = ClusterMCMC[RowsSampled[x], "sigma"])
    
    Predictions <- XPredictions[[1]]@x + ZPredictions2a + ZPredictions2b
    
    PZero <- rbinom(n = N, size = 1, prob = logistic(Predictions))
    
    PZero
    
  }, mc.cores = 10)
  
  # Simulating the network #####
  
  SubSims[[z]] <- parallel::mclapply(1:length(SubPredList[[z]]), function(i){
    
    AssMat <- matrix(NA, 
                     nrow = length(union(AllMammaldf$Sp,AllMammaldf$Sp2)), 
                     ncol = length(union(AllMammaldf$Sp,AllMammaldf$Sp2)))
    
    AssMat[-which(1:length(AssMat)%in%UpperMammals)] <- round(SubPredList[[z]][[i]])
    AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
    diag(AssMat) <- 0
    dimnames(AssMat) <- list(union(AllMammaldf$Sp,AllMammaldf$Sp2),
                             union(AllMammaldf$Sp,AllMammaldf$Sp2))
    
    as(AssMat, "dgCMatrix")
    
  }, mc.cores = 10)
  
}

save(SubDFList, SubPredList, SubPredDF, SubSims, file = "SubModelGenData.Rdata")


