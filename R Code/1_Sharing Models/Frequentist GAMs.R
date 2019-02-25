
# Running Frequentist GAMS

if(file.exists("Finaldf.Rdata")) load("Finaldf.Rdata") else source("R Code/00_Master Code.R")

FinalHostMatrix$Sp <- factor(FinalHostMatrix$Sp, levels = union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))
FinalHostMatrix$Sp2 <- factor(FinalHostMatrix$Sp2, levels = union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))

library(mgcv)

Resps <- c("VirusBinary","RNA","DNA","Vector","NVector")

BAMList <- DataList <- list()

for(r in 1:length(Resps)){
  
  print(r)
  
  DataList[[Resps[r]]] <- FinalHostMatrix %>% filter(!is.na(Resps[r]))
  
  DataList[[Resps[r]]]$Sp <- factor(DataList[[Resps[r]]]$Sp, levels = union(DataList[[Resps[r]]]$Sp,DataList[[Resps[r]]]$Sp2))
  DataList[[Resps[r]]]$Sp2 <- factor(DataList[[Resps[r]]]$Sp2, levels = union(DataList[[Resps[r]]]$Sp,DataList[[Resps[r]]]$Sp2))
  
  Formula = as.formula(paste0(Resps[r], "~ t2(Space, scale(Phylo2)) + s(DietSim)"))
  
  BAMList[[Resps[r]]] <- bam(Formula,# + Eaten, # + 
                             #(1 + mmc(dom, dom_2) + mmc(d_cites_s1, d_cites_s2) | mm(Sp, Sp2)),
                             data = DataList[[Resps[r]]], 
                             family = binomial())
  
}

# trying random effects ####

r = 1

MZ1 <- model.matrix( ~ Sp - 1, data = DataList[[Resps[r]]]) %>% as.matrix
MZ2 <- model.matrix( ~ Sp2 - 1, data = DataList[[Resps[r]]]) %>% as.matrix

SppMatrix = MZ1 + MZ2

DataList[[Resps[[r]]]]$Spp <- SppMatrix
DataList[[Resps[[r]]]]$Cites <- log(DataList[[Resps[r]]][,c("hDiseaseZACites","hDiseaseZACites.Sp2")] + 1) %>% as.matrix
DataList[[Resps[[r]]]]$Domestic <- cbind(2- FinalHostMatrix$hDom %>% as.factor %>% as.numeric,
                                         2- FinalHostMatrix$hDom.Sp2 %>% as.factor %>% as.numeric) %>% as.matrix

PP <- list(Spp = list(rank = nlevels(DataList[[Resps[r]]]$Sp), diag(nlevels(DataList[[Resps[r]]]$Sp))),
           Cites = list(rank = 2, diag(2)),
           Domestic = list(rank = 2, diag(2)))

Formula2 = as.formula(paste0(Resps[r], 
                             "~ t2(Space, scale(Phylo2)) + s(DietSim) + 
                             Spp + Cites + Domestic"))

m1 <- bam(Formula2,
          data = DataList[[Resps[r]]], 
          family = binomial(),
          paraPen = PP)

