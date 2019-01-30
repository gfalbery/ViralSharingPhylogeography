
# Summarising parallel models ####

load("Model Runs/ZI_runs.Rdata")

library(MCMCglmm); library(tidyverse); library(ggregplot)

i = 1

mc <- ZI_runs[1:10 + 10*(i-1)] %>% lapply(function(a) a$Sol[,1:14])
#mc <- ZI_runs[1:10 + 10*(i-1)] %>% lapply(function(a) a$VCV)
mc <- do.call(mcmc.list, mc)
#par(mfrow=c(7,2), mar=c(2,2,1,2), ask = F)
#gelman.plot(mc, auto.layout=F)
#gelman.diag(mc)
par(mfrow=c(7,4), mar=c(2, 1, 1, 1))
plot(mc, ask=F, auto.layout=F)

FullZIModel <- ZI_runs[[1]]
ZINoGModel <- ZI_runs[[11]]
OverlapZIModel <- ZI_runs[[21]]
OverlapNoGZIModel <- ZI_runs[[31]]

ModelList <- list(FullZIModel, ZINoGModel, OverlapZIModel, OverlapNoGZIModel)

for(i in 1:4){
  
  if(i ==3){
    ClusterMCMC <- ZI_runs[1:10 + (i-1)*10][-6] %>% lapply(function(a) as.data.frame(as.matrix(a$Sol[,1:14]))) %>% bind_rows
    ClusterMCMCv <- ZI_runs[1:10 + (i-1)*10][-6] %>% lapply(function(a) as.data.frame(as.matrix(a$VCV))) %>% bind_rows
    
    rows <- sample(1:dim(ClusterMCMC)[1], 1000)
    
    SampledSol <- ClusterMCMC[rows,]
    SampledVCV <- ClusterMCMCv[rows,]
    
    ModelList[[i]]$Sol <- as.mcmc(SampledSol)
    ModelList[[i]]$VCV <- as.mcmc(SampledVCV)
    
  } else {
    ClusterMCMC <- ZI_runs[1:10 + (i-1)*10] %>% lapply(function(a) as.data.frame(as.matrix(a$Sol))) %>% bind_rows
    ClusterMCMCv <- ZI_runs[1:10 + (i-1)*10] %>% lapply(function(a) as.data.frame(as.matrix(a$VCV))) %>% bind_rows
    
    rows <- sample(1:dim(ClusterMCMC)[1], 1000)
    
    SampledSol <- ClusterMCMC[rows,]
    SampledVCV <- ClusterMCMCv[rows,]
    
    ModelList[[i]]$Sol <- as.mcmc(SampledSol)
    ModelList[[i]]$VCV <- as.mcmc(SampledVCV)
  }
}

Efxplot(ModelList)

ZISols <- lapply(1:2, function(a){ summary(ModelList[c(1,2,4)][[a]])$solutions %>% as.data.frame %>% mutate(
  Component = rep(c("Count", "ZI"), dim(summary(ModelList[[1]])$solutions)[1]/2),
  Variable = rep(c("Intercept", "Space", "Phylogeny", "Citations", "DomDom", "DomWild", 
                   "Phylo:Space"),
                 each = 2),
  Name = paste(Component, Variable, sep = ":"),
  Model = c("Full", "No Random")[a]
) %>% rename(Lower = "l-95% CI", Upper = "u-95% CI", Estimate = "post.mean")}) %>% bind_rows


ZISols$Estimate[ZISols$Component=="ZI"] <- -ZISols$Estimate[ZISols$Component=="ZI"]
ZISols[ZISols$Component=="ZI", c("Lower", "Upper")] <- -ZISols[ZISols$Component=="ZI", c("Upper","Lower")]

ggplot(ZISols, aes(x = Variable, y = Estimate, colour = as.factor(Model))) + 
  geom_point(position = position_dodge(w = 0.5)) + 
  geom_errorbar(position = position_dodge(w = 0.5), aes(ymin = Lower, 
                                                        ymax = Upper), size = 0.3, width = 0.2) + 
  geom_hline(aes(yintercept = 0), lty = 2) + labs(x = NULL, colour = "Model") + coord_flip() + AlberTheme +
  facet_wrap(~Component, scales = "free_x")

# Creating replicates of the models as a whole ####

ZI_subrun <- ZI_runs[1:10]

s1 <- lapply(ZI_subrun, function(a) as.data.frame(a$Sol[,1:14])) %>% bind_rows
v1 <- lapply(ZI_subrun, function(a) as.data.frame(a$VCV)) %>% bind_rows

for(x in 1:length(ZI_subrun)){
  ZI_subrun[[x]]$Sol <- as.mcmc(s1[sample(1:dim(s1)[1], 1000),])
  ZI_subrun[[x]]$VCV <- as.mcmc(v1[sample(1:dim(v1)[1], 1000),])
}

ZISols <- lapply(1:10, function(a){ summary(ZI_subrun[[a]])$solutions %>% as.data.frame %>% mutate(
  Component = rep(c("Count", "ZI"), dim(summary(ZI_runs[[1]])$solutions)[1]/2),
  Variable = rep(c("Intercept", "Space", "Phylogeny", "Citations", "DomDom", "DomWild", 
                   "Phylo:Space"),
                 each = 2),
  Name = paste(Component, Variable, sep = ":"),
  Model = a
) %>% rename(Lower = "l-95% CI", Upper = "u-95% CI", Estimate = "post.mean")}) %>% bind_rows

ZISols$Estimate[ZISols$Component=="ZI"] <- -ZISols$Estimate[ZISols$Component=="ZI"]
ZISols[ZISols$Component=="ZI", c("Lower", "Upper")] <- -ZISols[ZISols$Component=="ZI", c("Upper","Lower")]

ggplot(ZISols, aes(x = Variable, y = Estimate, colour = Component)) + 
  geom_point(position = position_dodge(w = 0.5)) + 
  geom_errorbar(position = position_dodge(w = 0.5), aes(ymin = Lower, 
                                                        ymax = Upper), size = 0.3, width = 0.2) + 
  geom_hline(aes(yintercept = 0), lty = 2) + labs(x = NULL) + coord_flip() + AlberTheme +
  facet_wrap(~Model, scales = "free_x")

ggplot(ZISols, aes(x = Variable, y = Estimate, colour = as.factor(Model), shape = as.factor(Component), lty = Component)) + 
  geom_point(position = position_dodge(w = 0.5)) + 
  geom_errorbar(position = position_dodge(w = 0.5), aes(ymin = Lower, 
                                                        ymax = Upper), size = 0.3, width = 0.2) + 
  geom_hline(aes(yintercept = 0), lty = 2) + labs(x = NULL) + coord_flip() + AlberTheme +
  facet_wrap(~Component, scales = "free_x")

# Fourth set ####

ZI_subrun <- ZI_runs[1:10 + 30]

s1 <- lapply(ZI_subrun, function(a) as.data.frame(a$Sol[,1:14])) %>% bind_rows
v1 <- lapply(ZI_subrun, function(a) as.data.frame(a$VCV)) %>% bind_rows

for(x in 1:length(ZI_subrun)){
  ZI_subrun[[x]]$Sol <- as.mcmc(s1[sample(1:dim(s1)[1], 1000),])
  ZI_subrun[[x]]$VCV <- as.mcmc(v1[sample(1:dim(v1)[1], 1000),])
}

ZISols <- lapply(1:10, function(a){ summary(ZI_subrun[[a]])$solutions %>% as.data.frame %>% mutate(
  Component = rep(c("Count", "ZI"), dim(summary(ZI_runs[[1]])$solutions)[1]/2),
  Variable = rep(c("Intercept", "Space", "Phylogeny", "Citations", "DomDom", "DomWild", 
                   "Phylo:Space"),
                 each = 2),
  Name = paste(Component, Variable, sep = ":"),
  Model = a
) %>% rename(Lower = "l-95% CI", Upper = "u-95% CI", Estimate = "post.mean")}) %>% bind_rows

ZISols$Estimate[ZISols$Component=="ZI"] <- -ZISols$Estimate[ZISols$Component=="ZI"]
ZISols[ZISols$Component=="ZI", c("Lower", "Upper")] <- -ZISols[ZISols$Component=="ZI", c("Upper","Lower")]

ggplot(ZISols, aes(x = Variable, y = Estimate, colour = Component)) + 
  geom_point(position = position_dodge(w = 0.5)) + 
  geom_errorbar(position = position_dodge(w = 0.5), aes(ymin = Lower, 
                                                        ymax = Upper), size = 0.3, width = 0.2) + 
  geom_hline(aes(yintercept = 0), lty = 2) + labs(x = NULL) + coord_flip() + AlberTheme +
  facet_wrap(~Model, scales = "free_x")

ggplot(ZISols, aes(x = Variable, y = Estimate, colour = as.factor(Model), shape = as.factor(Component), lty = Component)) + 
  geom_point(position = position_dodge(w = 0.5)) + 
  geom_errorbar(position = position_dodge(w = 0.5), aes(ymin = Lower, 
                                                        ymax = Upper), size = 0.3, width = 0.2) + 
  geom_hline(aes(yintercept = 0), lty = 2) + labs(x = NULL) + coord_flip() + AlberTheme +
  facet_wrap(~Component, scales = "free_x")

