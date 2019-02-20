
# Comparison of effects across different STAN models ####

library(MCMCglmm)

BinModel <- readRDS("~/Albersnet/Model Files/RNABinModel.rds")
RNABinModel <- readRDS("~/Albersnet/Model Files/RNABinModel.rds")
DNABinModel <- readRDS("~/Albersnet/Model Files/DNABinModel.rds")
VectorBinModel <- readRDS("~/Albersnet/Model Files/VectorBinModel.rds")
NVectorBinModel <- readRDS("~/Albersnet/Model Files/NVectorBinModel.rds")

list(
  stan_plot(RNABinModel, XBetas[2:4]),
  stan_plot(DNABinModel, XBetas[2:4]),
  stan_plot(VectorBinModel, XBetas[2:4]),
  stan_plot(NVectorBinModel, XBetas[2:4])
) %>% arrange_ggplot2

load("~/Albersnet/Output Files/Bin Model Output.Rdata")
load("~/Albersnet/Output Files/RNAModelOutput.Rdata")
load("~/Albersnet/Output Files/NVectorModelOutput.Rdata")
load("~/Albersnet/Output Files/VectorModelOutput.Rdata")
load("~/Albersnet/Output Files/DNAModelOutput.Rdata")

SubModels <- list(Full = p,
                  RNA = q,
                  Vector = s,
                  NVector = t,
                  DNA = r)

save(SubModels, file = "Output Files/SubModels.Rdata")

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

Models <- c("Full", "RNA", "Vector", "NVector", "DNA")

Estdf$Model <- factor(rep(Models, sapply(Estimates, nrow)), levels = Models)
Estdf$Var <- SolList %>% lapply(colnames) %>% unlist

XBetas <- c("mu_alpha","beta_space","beta_phylo","beta_inter")
ZBetas <- c("beta_d_cites_s","beta_domestic")

ggplot(Estdf, aes(x = Var, y = Estimate, 
                  colour = Model)) + 
  geom_point(position = position_dodge(w = 0.5)) + 
  geom_errorbar(position = position_dodge(w = 0.5), 
                aes(ymin = Lower, 
                    ymax = Upper), 
                size = 0.3, width = 0) + 
  geom_hline(aes(yintercept = 0), lty = 2) + 
  labs(x = NULL) + coord_flip() + 
  scale_x_discrete(limits = c(XBetas, ZBetas))

STANfx <- function(which = 1, var){
  
  model = Estdf %>% filter(Model == Models[which])
  
  paste0(model[model[,"Var"]==var,"Estimate"]%>% round(2),
         " (", 
         model[model[,"Var"]==var,"Lower"]%>% round(2),
         ", ", 
         model[model[,"Var"]==var,"Upper"] %>% round(2),")")
  
}







