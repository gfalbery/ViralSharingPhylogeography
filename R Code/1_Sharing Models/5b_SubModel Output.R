
# Comparison of effects across different STAN models ####

library(MCMCglmm)

RNABinModel <- readRDS("~/Albersnet/RNABinModel.rds")
DNABinModel <- readRDS("~/Albersnet/DNABinModel.rds")
VectorBinModel <- readRDS("~/Albersnet/VectorBinModel.rds")
NVectorBinModel <- readRDS("~/Albersnet/NVectorBinModel.rds")

list(
  stan_plot(RNABinModel, XBetas[2:4]),
  stan_plot(DNABinModel, XBetas[2:4]),
  stan_plot(VectorBinModel, XBetas[2:4]),
  stan_plot(NVectorBinModel, XBetas[2:4])
) %>% arrange_ggplot2


load("~/Albersnet/Bin Model Output.Rdata")
load("~/Albersnet/NVectorModelOutput.Rdata")
load("~/Albersnet/VectorModelOutput.Rdata")
load("~/Albersnet/DNAModelOutput.Rdata")

SubModels <- list(Full = p,
                  Vector = s,
                  NVector = t,
                  DNA = r)

save(SubModels, file = "SubModels.Rdata")

SolList <- map(SubModels, "df")

Estimates <- lapply(SolList, function(a){
  
  Samples <- sample(1:nrow(a), 1000)
  
  apply(a, 2, function(b){
    
    data.frame(Estimate = posterior.mode(as.mcmc(b[Samples])),
               Lower = HPDinterval(as.mcmc(b[Samples]))[1],
               Upper = HPDinterval(as.mcmc(b[Samples]))[2])
    
  }) %>% bind_rows()
  
}) 

Estdf <- Estimates %>% bind_rows()

Models <- c("Full", "Vector", "NVector", "DNA")

Estdf$Model <- rep(Models, sapply(Estimates, nrow))
Estdf$Var <- SolList %>% lapply(colnames) %>% unlist

ggplot(Estdf, aes(x = Var, y = Estimate, 
                  colour = Model)) + 
  geom_point(position = position_dodge(w = 0.5)) + 
  geom_errorbar(position = position_dodge(w = 0.5), 
                aes(ymin = Lower, 
                    ymax = Upper), 
                size = 0.3, width = 0) + 
  geom_hline(aes(yintercept = 0), lty = 2) + 
  labs(x = NULL) + coord_flip() + 
  scale_x_discrete(limits = c(XBetas,ZBetas))


STANfx <- function(which = 1, var){
  
  model = Estdf %>% filter(Model == Models[which])
  
  paste0(model[model[,"Var"]==var,"Estimate"]%>% round(2),
         " (", 
         model[model[,"Var"]==var,"Lower"]%>% round(2),
         ", ", 
         model[model[,"Var"]==var,"Upper"] %>% round(2),")")
  
}







