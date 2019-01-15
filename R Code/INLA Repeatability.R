
# INLA ICC function ####

library(tidyverse)

INLARep <- function(Model, Component, Family = "gaussian"){
  
  MySqrt <- function(x) {1 / sqrt(x) }
  tau <- Model$marginals.hyperpar$`Precision for the Gaussian observations`
  sigma <- inla.emarginal(MySqrt, tau)
  sigma
  
  sigma2 <- inla.emarginal(MySqrt, Model$marginals.hyperpar[[paste0("Precision for ", Component)]])

  CI <- 1/sqrt(inla.qmarginal(c(0.025, 0.975),Model$marginals.hyperpar[[paste0("Precision for ", Component)]]))
  
  list(Estimate = sigma2^2/(sigma^2 + sigma2^2),
       CI = CI^2/(sigma^2 + CI^2)) %>% return
       
}
