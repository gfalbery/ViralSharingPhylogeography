
# Rscript "R Code/1_Sharing Models/1b_Frequentist GAM Output.R" ####

if(file.exists("Output Files/Finaldf.Rdata")) load("Output Files/Finaldf.Rdata") else source("R Code/00_Master Code.R")

library(mgcv); library(tidyverse); library(ggregplot); library(MASS)

load("Output Files/BAMList.Rdata")
#load("Output Files/BAMList2.Rdata")

Resps <- c("VirusBinary","RNA","DNA","Vector","NVector")

FitList <- PostList <- DrawList <- list()

for(r in 1:length(BAMList)){
  
  Model <- BAMList[[Resps[r]]]
  
  print(Resps[r])
  
  # Model Checking ####
  
  SpCoefNames <- names(Model$coef)[substr(names(Model$coef),1,5)=="SppSp"]
  SpCoef <- Model$coef[SpCoefNames]
  
  # Effects ####
  
  SpaceRange <- seq(from = min(DataList[[Resps[r]]]$Space),
                    to = max(DataList[[Resps[r]]]$Space),
                    length = 100) %>% 
    c(mean(DataList[[Resps[r]]]$Space))
  
  PhyloRange <- seq(from = min(DataList[[Resps[r]]]$Phylo),
                    to = max(DataList[[Resps[r]]]$Phylo),
                    length = 100)  %>% 
    c(mean(DataList[[Resps[r]]]$Phylo))
  
  DietRange <- seq(from = min(DataList[[Resps[r]]]$DietSim),
                   to = max(DataList[[Resps[r]]]$DietSim),
                   length = 10)  %>% 
    c(mean(DataList[[Resps[r]]]$DietSim))
  
  FitList[[Resps[r]]] <- expand.grid(Space = SpaceRange,
                                     Phylo = PhyloRange,
                                     DietSim = DietRange,
                                     MinCites = mean(DataList[[Resps[r]]]$MinCites),
                                     Domestic = 0
  )
  
  FitList[[Resps[r]]]$Spp <- matrix(0 , nrow = nrow(FitList[[Resps[r]]]), ncol = length(SpCoef))# %>% as("dgCMatrix")
  
  FitList[[Resps[r]]] <- FitList[[Resps[r]]] %>% mutate(SpaceQ = cut(Space, quantile(Space, 0:10/10),include.lowest = T, labels = 1:10),
                                                        PhyloQ = cut(Phylo, quantile(Phylo, 0:10/10),include.lowest = T, labels = 1:10))
  
  FitPredictions  <- predict.gam(Model, 
                                 newdata = FitList[[Resps[r]]])
  
  FitList[[Resps[r]]][,"Fit"] <- logistic(FitPredictions)

  print("Getting posterior uncertainty!")
  
  # Posterior Uncertainty Simulation #### https://www.fromthebottomoftheheap.net/2014/06/16/simultaneous-confidence-intervals-for-derivatives/
  
  for(i in c("Space", "Phylo")){
    
    print(i)
    
    PredData = FitList[[Resps[r]]] %>% filter(
      DietSim == dplyr::last(unique(DietSim)))
    
    if(i == "Space") PredData <- PredData %>% filter(Phylo == last(unique(Phylo))) else{
      
      PredData <- PredData %>% filter(Space == last(unique(Space)))
      
    }
    
    lp <- predict(Model, newdata = PredData, 
                  type = "lpmatrix") %>% 
      as.data.frame()
    
    coefs <- coef(Model)
    vc <- vcov(Model)
    
    sim <- mvrnorm(100, mu = coefs, Sigma = vc)
    
    want <- lp %>% colnames
    
    lp <- lp %>% as.matrix #%>% logistic
    
    fits <- lp[, want] %*% t(sim[, want]) %>% as.data.frame() %>%
      mutate(i = PredData[,i])
    
    PostList[[Resps[r]]][[i]] <- gather(fits, key = "Draw", value = "Fit", -i) %>%
      mutate(Fit = logistic(Fit))
    
  }
  
  DrawList[[Resps[r]]] <- list()
  
  for(i in c("Space", "Phylo")){
    
    print(i)
    
    lp = list()
    
    for(j in 1:100){
      
      print(j)
      
      PredData = FitList[[Resps[r]]] %>% filter(
        DietSim == dplyr::last(unique(DietSim)))
      
      if(i == "Space"){
        
        PredData <- PredData %>% filter(Phylo == last(unique(Phylo))) %>%
          mutate(Phylo = sample(DataList[[Resps[r]]]$Phylo, 1))
        
      } else {
        
        PredData <- PredData %>% filter(Space == last(unique(Space))) %>%
          mutate(Space = sample(DataList[[Resps[r]]]$Space, 1))
        
      }
      
      lp[[j]] <- data.frame(Fit = predict(Model, newdata = PredData),
                            Iteration = as.factor(j),
                            i = PredData[,i])
      
    }
    
    DrawList[[Resps[r]]][[i]] <- lp %>% bind_rows() %>% as.data.frame() %>%
      mutate(Fit = logistic(Fit))
    
  }
  
}

save(FitList, PostList, DrawList, file = "Output Files/FitList.Rdata")




