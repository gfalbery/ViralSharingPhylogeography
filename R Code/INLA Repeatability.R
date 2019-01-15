
# INLA ICC function ####

library(tidyverse)

INLARep <- function(Model, Component, Family = "gaussian"){
  
  MySqrt <- function(x) {1 / sqrt(x) }
  
  SigmaList <- list()
  
  for(x in 1:(length(Model$marginals.hyperpar[-which(names(Model$marginals.hyperpar) == paste0("Precision for ", Component))]))){
    subhyperpars <- Model$marginals.hyperpar[-which(names(Model$marginals.hyperpar) == paste0("Precision for ", Component))]
    tau <- subhyperpars[[x]]
    sigma <- inla.emarginal(MySqrt, tau)
    SigmaList[[x]] <- sigma
  }
  
  sigma2 <- inla.emarginal(MySqrt, Model$marginals.hyperpar[[paste0("Precision for ", Component)]])
  
  CI <- 1/sqrt(inla.qmarginal(c(0.025, 0.975),Model$marginals.hyperpar[[paste0("Precision for ", Component)]]))
  
  list(Estimate = sigma2^2/(sum(unlist(SigmaList)^2) + sigma2^2),
       CI = CI^2/(sum(unlist(SigmaList)^2) + CI^2)) %>% return
  
}
