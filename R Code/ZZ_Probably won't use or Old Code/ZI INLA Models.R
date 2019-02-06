# ZI INLA Models

source("R Code/00_Master Code.R")

library(INLA)

inla.setOption(num.threads = 8)

lapply(1:3, function(i){
  if(i==1) {
    
    saveRDS(inla(data = FinalHostMatrix, # Doesn't fit
                 Virus ~ Space + Phylo2 + Space:Phylo2 + MinCites + DomDom,
                 control.compute = list(dic = TRUE),
                 family = "nbinomial"), file = paste0("ZI INLA Model",i,".Rdata"))}
  if(i==2) {
    
    saveRDS(inla(data = FinalHostMatrix, # Doesn't fit
                 Virus ~ Space + Phylo2 + Space:Phylo2 + MinCites + DomDom,
                 control.compute = list(dic = TRUE),
                 family = "zeroinflatednbinomial1"), file = paste0("ZI INLA Model",i,".Rdata"))}
  
  if(i==3) {
    
    saveRDS(inla(data = FinalHostMatrix, # Doesn't fit
                 Virus ~ Space + Phylo2 + Space:Phylo2 + MinCites + DomDom,
                 control.compute = list(dic = TRUE),
                 family = "zeroinflatednbinomial2"), file = paste0("ZI INLA Model",i,".Rdata"))}
})