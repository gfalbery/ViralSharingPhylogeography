
# INLA ICC function ####

library(tidyverse)

INLARep <- function(Model, Family = "gaussian"){
  
  MySqrt <- function(x) {1 / sqrt(x) }
  
  SigmaList <- list()
  
  for(x in 1:(length(Model$marginals.hyperpar))){
    subhyperpars <- Model$marginals.hyperpar
    tau <- subhyperpars[[x]]
    sigma <- inla.emarginal(MySqrt, tau)
    SigmaList[[x]] <- sigma
  }
  
  names(SigmaList) <- names(Model$marginals.hyperpar)
  
  lapply(SigmaList, function(a) a^2/sum(sapply(SigmaList, function(b) b^2)))
  
}
