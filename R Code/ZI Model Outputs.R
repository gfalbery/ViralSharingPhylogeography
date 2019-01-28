
# Summarising parallel models ####

load("ZI_runs.Rdata")

mc <- ZI_runs[1:10] %>% lapply(function(a) a$Sol[,1:14])
mc <- do.call(mcmc.list, mc)
par(mfrow=c(7,2), mar=c(2,2,1,2), ask = F)
gelman.plot(mc, auto.layout=F)
gelman.diag(mc)
par(mfrow=c(14,2), mar=c(2, 1, 1, 1))
plot(mc, ask=F, auto.layout=F)

mc <- ZI_runs[1:10+10] %>% lapply(function(a) a$Sol[,1:14])
mc <- do.call(mcmc.list, mc)
par(mfrow=c(7,2), mar=c(2,2,1,2), ask = F)
gelman.plot(mc, auto.layout=F)
gelman.diag(mc)
par(mfrow=c(14,2), mar=c(2, 1, 1, 1))
plot(mc, ask=F, auto.layout=F)

ClusterMCMC <- ZI_10runs %>% lapply(function(a) as.data.frame(as.matrix(a$Sol[,1:14]))) %>% bind_rows
ClusterMCMCv <- ZI_10runs %>% lapply(function(a) as.data.frame(as.matrix(a$VCV))) %:% bind_rows

SampleCluster <- as.mcmc(apply(bind_rows(ClusterMCMC), 2, function(a) as.mcmc(sample(a, 1000))))
SampleClusterv <- as.mcmc(apply(bind_rows(ClusterMCMCv), 2, function(a) as.mcmc(sample(a, 1000))))

FullZIModel <- ZI_runs[[1]]
ZINoGModel <- ZI_runs[[11]]
OverlapZIModel <- ZI_runs[[21]]
OverlapNoGZIModel <- ZI_runs[[31]]

ModelList <- list(FullZIModel, ZINoGModel, OverlapZIModel, OverlapNoGZIModel)

for(i in c(1,2,4)){
  
  ClusterMCMC <- ZI_runs[1:10 + (i-1)*10] %>% lapply(function(a) as.data.frame(as.matrix(a$Sol[,1:14]))) %>% bind_rows
  ClusterMCMCv <- ZI_runs[1:10 + (i-1)*10] %>% lapply(function(a) as.data.frame(as.matrix(a$VCV))) %>% bind_rows
  
  rows <- sample(1:dim(ClusterMCMC)[1], 1000)
  
  SampledSol <- ClusterMCMC[rows,]
  SampledVCV <- ClusterMCMCv[rows,]
  
  ModelList[[i]]$Sol <- as.mcmc(SampledSol)
  ModelList[[i]]$VCV <- as.mcmc(SampledVCV)
  
}

Efxplot(ModelList[c(1,2,4)])

lol1 <- predict(ModelList[[1]])


ZISols <- lapply(1:3, function(a){ summary(ModelList[c(1,2,4)][[a]])$solutions %>% as.data.frame %>% mutate(
  Component = rep(c("Count", "ZI"), dim(summary(ModelList[[1]])$solutions)[1]/2),
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

