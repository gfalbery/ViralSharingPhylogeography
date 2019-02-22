
# Running Frequentist GAMS

library(mgcv)

Resps <- c("VirusBinary","RNA","DNA","Vector","NVector")

BAMList <- DataList <- list()

for(r in 1:length(Resps)){
  
  print(r)
  
  DataList[[Resps[r]]] <- FinalHostMatrix %>% filter(!is.na(Resps[r]))
  
  Formula = as.formula(paste0(Resps[r], "~ t2(Space, scale(Phylo2)) + s(DietSim)"))
  
  BAMList[[Resps[r]]] <- bam(Formula,# + Eaten, # + 
                             #(1 + mmc(dom, dom_2) + mmc(d_cites_s1, d_cites_s2) | mm(Sp, Sp2)),
                             data = DataList[[Resps[r]]], 
                             family = binomial())
  
}






