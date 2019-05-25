
# Rscript "R Code/1_Sharing Models/1b_Frequentist GAM Output.R" ####

if(file.exists("Output Files/FitList.Rdata")) load("Output Files/FitList.Rdata") else {
  
  if(file.exists("Output Files/Finaldf.Rdata")) load("Output Files/Finaldf.Rdata") else source("R Code/00_Master Code.R")
  
  library(mgcv); library(tidyverse); library(ggregplot); library(MASS)
  
  load("Output Files/BAMList.Rdata")
  #load("Output Files/BAMList2.Rdata")
  
  Resps <- c("VirusBinary","RNA","DNA","Vector","NVector")
  
  FitList <- PostList <- DrawList <- list()
  
  r = 1
  
  for(r in 1:length(BAMList)){
    
    Model <- BAMList[[Resps[r]]]
    
    print(Resps[r])
    
    # Model Checking ####
    
    SpCoefNames <- names(Model$coef)[substr(names(Model$coef),1,5)=="SppSp"]
    SpCoef <- Model$coef[SpCoefNames]
    
    # Effects ####
    
    SpaceRange <- seq(from = 0,
                      to = 1,
                      length = 101) %>% 
      c(mean(DataList[[Resps[r]]]$Space))
    
    PhyloRange <- seq(from = 0,
                      to = 1,
                      length = 101)  %>% 
      c(mean(DataList[[Resps[r]]]$Phylo))
    
    FitList[[Resps[r]]] <- expand.grid(Space = SpaceRange,
                                       Phylo = PhyloRange,
                                       MinCites = mean(DataList[[Resps[r]]]$MinCites),
                                       Domestic = 0
    ) %>%
      mutate(SpaceQuantile = ifelse(Space == last(unique(Space)), "1.5%",
                                    ifelse(Space == 0, "0%",
                                           ifelse(Space == 0.25, "25%",
                                                  ifelse(Space == 0.5, "50%", NA)))),
             
             PhyloQuantile = ifelse(Phylo == last(unique(Phylo)), "0.1",
                                    ifelse(Phylo == 0, "0",
                                           ifelse(Phylo == 0.25, "0.25",
                                                  ifelse(Phylo == 0.5, "0.5", NA)))),
             Gz = as.numeric(Space==0))
    
    FitList[[Resps[r]]]$Spp <- matrix(0 , nrow = nrow(FitList[[Resps[r]]]), ncol = length(SpCoef))
    
    FitPredictions  <- predict.gam(Model, 
                                   newdata = FitList[[Resps[r]]], 
                                   se.fit = T)
    
    FitList[[Resps[r]]][,c("Fit","Lower", "Upper")] <- logistic(with(FitPredictions, cbind(fit, fit - se.fit, fit + se.fit)))
    
    print("Getting posterior uncertainty!")
    
    # Posterior Uncertainty Simulation #### https://www.fromthebottomoftheheap.net/2014/06/16/simultaneous-confidence-intervals-for-derivatives/
    
    for(i in c("Space", "Phylo")){
      
      print(i)
      
      PredData <- FitList[[Resps[r]]] 
      
      if(i == "Space") PredData <- PredData %>% filter(Phylo == last(unique(Phylo))) else{
        
        PredData <- PredData %>% filter(Space == 0)
        
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
    
    Draw = F
    
    if(Draw){
      
      DrawList[[Resps[r]]] <- list()
      
      for(i in c("Space", "Phylo")){
        
        print(i)
        
        lp = list()
        
        for(j in 1:100){
          
          print(j)
          
          PredData <- FitList[[Resps[r]]]
          
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
    
  }
  
  save(FitList, PostList, DrawList, file = "Output Files/FitList.Rdata")
  
}

plot_grid( FitList[["VirusBinary"]] %>% 
             filter(!is.na(SpaceQuantile)) %>%
             ggplot(aes(Phylo, Fit, colour = SpaceQuantile)) + 
             geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = SpaceQuantile), alpha = 0.2, colour = NA) +
             geom_line(aes(group = as.factor(Space))) +
             labs(y = "Predicted Viral Sharing", x = "Phylogenetic Similarity", 
                  colour = "Overlap", fill = "Overlap") +
             lims(x = c(0,1), y = c(0,1)) +
             coord_fixed() +
             scale_color_discrete_sequential(palette = AlberPalettes[[1]], nmax = 8, order = 5:8)  +
             scale_fill_discrete_sequential(palette = AlberPalettes[[1]], nmax = 8, order = 5:8)  +
             theme(legend.position = c(0.1, 0.8), legend.background = element_rect(colour = "dark grey")) +
             geom_rug(data = DataList[[1]], inherit.aes = F, aes(x = Phylo), alpha = 0.01),
           
           FitList[["VirusBinary"]] %>% 
             filter(!is.na(PhyloQuantile)) %>%
             ggplot(aes(Space, Fit, colour = PhyloQuantile)) + 
             geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = PhyloQuantile), alpha = 0.2, colour = NA) +
             geom_line(aes(group = as.factor(Phylo))) +
             labs(y = "Predicted Viral Sharing", x = "Geographic Overlap", 
                  colour = "Relatedness", fill = "Relatedness") +
             lims(x = c(0,1), y = c(0,1)) +
             coord_fixed() +
             scale_color_discrete_sequential(palette = AlberPalettes[[2]], nmax = 8, order = 5:8)  +
             scale_fill_discrete_sequential(palette = AlberPalettes[[2]], nmax = 8, order = 5:8)  +
             theme(legend.position = c(0.1, 0.8), legend.background = element_rect(colour = "dark grey")) +
             geom_rug(data = DataList[[1]], inherit.aes = F, aes(x = Space), alpha = 0.01),
           
           FitList[["VirusBinary"]] %>% 
             filter(!Phylo == last(unique(Phylo)),
                    !Space == last(unique(Space))) %>%
             ggplot(aes(Space, Phylo)) + 
             geom_tile(aes(fill = Fit)) + 
             labs(x = "Geographic Overlap", 
                  y = "Phylogenetic Similarity",
                  fill = "Estimate") +
             #ggtitle("Tensor Field") +
             lims(x = c(0,1), y = c(0,1)) +
             coord_fixed() +
             theme(legend.position = "bottom") +
             scale_fill_continuous_sequential(palette = "Greens 2", cmax = 20, end = 1),
           
           DataList$VirusBinary %>%
             ggplot(aes(Space, Phylo)) + 
             labs(x = "Geographic Overlap", 
                  y = "Phylogenetic Similarity") +
             #ggtitle("Data Distribution") +
             scale_fill_continuous_sequential(palette = "purp", begin = 0.2) +
             lims(x = c(0,1), y = c(0,1)) +
             coord_fixed() +
             theme(legend.position = "bottom") +
             geom_hex(aes(fill = stat(log(count)))),
           
           nrow = 2, 
           rel_heights = c(1,1.23), 
           labels = "AUTO")