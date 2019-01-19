# X_ Exploring Host Similarity Matrix Models

library(MCMCglmm); library(ggregplot)

UpperHosts <- which(upper.tri(HostAdj, diag = T))

HostMatrixdf$PropVirus[HostMatrixdf$PropVirus==0] <- 0.0001
HostMatrixdf$PropVirus[HostMatrixdf$PropVirus==1] <- 0.9999

IM1 <- inla(data = HostMatrixdf[-HostThemselves,], 
            PropVirus ~ Space + Phylo,
            family = "beta")

summary(IM1) # Interesting space is more important here ####

Efxplot(list(IM1))

IM1 <- inla(data = HostMatrixdf[-HostThemselves,], # Doesn't fit
            Virus ~ Space + Phylo,
            control.compute = list(dic = TRUE),
            family = "nbinomial")

IM2 <- inla(data = HostMatrixdf[-HostThemselves,], # Doesn't fit
            Virus ~ Space + Phylo,
            control.compute = list(dic = TRUE),
            family = "zeroinflatednbinomial1")

mf = 5

prior.zi <- list(R = list(V = diag(2), nu = 0, fix = 2))

library(parallel)

setCores <- 6 # use detectCores() by itself if you want all CPUs

# make the cluster
cl <- makeCluster(getOption("cl.cores", setCores))

cl.pkg <- clusterEvalQ(cl, library(MCMCglmm)) # load the MCMCglmm package within the cluster

clusterExport(cl, "prior.zi") # Import each object that's necessary to run the function
clusterExport(cl, "HostMatrixdf")
clusterExport(cl, "UpperHosts")
clusterExport(cl, "mf")

# use parLapply() to execute 10 runs of MCMCglmm(), each with nitt=100000
ClusterMCStart <- Sys.time()
ZI_10runs <- parLapply(cl = cl, 1:10, function(i) {
  
  MCMCglmm(#data = HostMatrixdf[-HostThemselves,], 
    data = HostMatrixdf[-UpperHosts,],
    Virus ~ trait -1 + trait:(Space + Phylo),
    rcov =~ idh(trait):units, 
    prior = prior.zi,
    family = "zipoisson",
    nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
    thin = 10*mf, burnin=3000*mf)
  
})

stopCluster(cl) # Stop running the parallel cluster
ClusterMCEnd <- Sys.time()

ClusterMCMC <- ZI_10runs %>% lapply(function(a) as.data.frame(as.matrix(a$Sol)))
ClusterMCMCv <- ZI_10runs %>% lapply(function(a) as.data.frame(as.matrix(a$VCV)))

SampleCluster <- as.mcmc(apply(bind_rows(ClusterMCMC), 2, function(a) as.mcmc(sample(a, 1000))))
SampleClusterv <- as.mcmc(apply(bind_rows(ClusterMCMCv), 2, function(a) as.mcmc(sample(a, 1000))))

SampleCluster %>% 
  apply(2, function(a) {
    qplot(1:1000, a, geom = "line") + labs(x = "Iteration", y = "Effect", title = colnames(a))
  }) %>% 
  arrange_ggplot2

ZIModel <- ZI_10runs[[1]]

ZIModel$Sol <- SampleCluster
ZIModel$VCV <- SampleClusterv
summary(ZIModel)

Efxplot(list(ZIModel)) + ggtitle("Space and Pylogeny profoundly influence number of viruses shared") +
  scale_x_discrete(limits = c("traitVirus", "traitzi_Virus", 
                              "traitzi_Virus:Space", "traitzi_Virus:Phylo", 
                              "traitVirus:Space", "traitVirus:Phylo"), 
                   labels = c("ZI Intercept", "Poisson Intercept",
                              "Zero:Space", "Zero:Phylo",
                              "Count:Space", "Count:Phylo")
  ) + ggsave("Zero-inflated model output.jpeg", units = "mm", height = 100, width = 100)

mf = 10

prior.zi <- list(R = list(V = diag(2), nu = 0, fix = 2))

library(parallel)

setCores <- 6 # use detectCores() by itself if you want all CPUs

# make the cluster
cl <- makeCluster(getOption("cl.cores", setCores))

cl.pkg <- clusterEvalQ(cl, library(MCMCglmm)) # load the MCMCglmm package within the cluster

clusterExport(cl, "prior.zi") # Import each object that's necessary to run the function
clusterExport(cl, "HostMatrixdf")
clusterExport(cl, "UpperHosts")
clusterExport(cl, "mf")

# use parLapply() to execute 10 runs of MCMCglmm(), each with nitt=100000
ClusterMCStart <- Sys.time()
ZI_10runs2 <- parLapply(cl = cl, 1:10, function(i) {
  
  MCMCglmm(#data = HostMatrixdf[-HostThemselves,], 
    data = HostMatrixdf[-UpperHosts,],
    Virus ~ trait -1 + trait:(Space + Phylo + Space:Phylo + Cites),
    rcov =~ idh(trait):units, 
    prior = prior.zi,
    family = "zipoisson",
    nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
    thin = 10*mf, burnin=3000*mf)
})

stopCluster(cl) # Stop running the parallel cluster
ClusterMCEnd <- Sys.time()

ClusterMCMC2 <- ZI_10runs2 %>% lapply(function(a) as.data.frame(as.matrix(a$Sol)))
ClusterMCMCv2 <- ZI_10runs2 %>% lapply(function(a) as.data.frame(as.matrix(a$VCV)))

library(coda)

mc <- ZI_10runs2 %>% lapply(function(a) a$Sol)
mc <- do.call(mcmc.list, mc)
par(mfrow=c(5,2), mar=c(2,2,1,2))
gelman.plot(mc, auto.layout=F)
gelman.diag(mc)
par(mfrow=c(8,2), mar=c(2, 1, 1, 1))
plot(mc, ask=F, auto.layout=F)

rows <- sample(1:nrow(bind_rows(ClusterMCMC2)), 1000)

SampleCluster2 <- as.mcmc(bind_rows(ClusterMCMC2)[rows,])
SampleClusterv2 <- as.mcmc(bind_rows(ClusterMCMCv2)[rows,])

ZIModel2 <- ZI_10runs2[[1]]

ZIModel2$Sol <- SampleCluster2
ZIModel2$VCV <- SampleClusterv2
summary(ZIModel2)

Efxplot(list(ZIModel)) + ggtitle("Space and Pylogeny profoundly influence number of viruses shared") +
  scale_x_discrete(limits = c("traitVirus", "traitzi_Virus", 
                              "traitzi_Virus:Space", "traitzi_Virus:Phylo", 
                              "traitVirus:Space", "traitVirus:Phylo"), 
                   labels = c("ZI Intercept", "Poisson Intercept",
                              "Zero:Space", "Zero:Phylo",
                              "Count:Space", "Count:Phylo")
  ) + ggsave("Zero-inflated model output.jpeg", units = "mm", height = 100, width = 100)

# Trying with spatial zeroes removed ####

f1 <- as.formula("Virus ~ Space + Phylo")

NoZeroes <- #HostMatrixdf %>% filter(!Sp==Sp2&Space>0) %>%
  MCMCglmm(
    data = HostMatrixdf %>% filter(!Sp==Sp2&Space>0),
    f1,
    #random =~ Sp,
    family = "poisson",
    nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
    thin = 10*mf, burnin=3000*mf)

f2 <- as.formula("Virus ~ Space + Phylo + Space:Phylo + Cites")

NoZeroes2 <- #HostMatrixdf %>% filter(!Sp==Sp2&Space>0) %>%
  MCMCglmm(
    f2,
    data = HostMatrixdf %>% filter(!Sp==Sp2&Space>0),
    family = "poisson",
    nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
    thin = 10*mf, burnin=3000*mf)


f3 <- as.formula("Virus ~ Space + Phylo + Space:Phylo + Cites")

NoZeroes3 <- #HostMatrixdf %>% filter(!Sp==Sp2&Space>0) %>%
  MCMCglmm(
    f3,
    data = HostMatrixdf %>% filter(!Sp==Sp2&Space>0),
    family = "poisson",
    nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
    thin = 10*mf, burnin=3000*mf)

f4 <- as.formula("Virus ~ Space + Phylo + Space:Phylo + Cites")

NoZeroes4 <- #HostMatrixdf %>% filter(!Sp==Sp2&Space>0) %>%
  MCMCglmm(
    f4,
    data = HostMatrixdf %>% filter(!Sp==Sp2&Space>0),
    family = "poisson",
    nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
    thin = 10*mf, burnin=3000*mf)

f5 <- as.formula("Virus ~ Phylo + Cites + Euclid")

OnlyZeroes <- #HostMatrixdf %>% filter(!Sp==Sp2&Space>0) %>%
  MCMCglmm(
    f5,
    data = HostMatrixdf %>% filter(!Sp==Sp2&Space==0),
    family = "poisson",
    nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
    thin = 10*mf, burnin=3000*mf)

# Plotting matrix correlations ####

arrange_ggplot2(list(
  ggplot(HostMatrixdf[-HostThemselves,], aes(Space, Phylo)) + 
    geom_point(colour = AlberColours[1], alpha = 0.3) + geom_smooth(colour = "black") +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2) +
    #coord_fixed(ratio = max(HostMatrixdf$Space)/max(HostMatrixdf$Phylo))
    labs(title = "Space ~ Phylogeny"), #+ ggpubr::stat_cor(method = "spearman"),
  
  ggplot(HostMatrixdf[-HostThemselves,], aes(Space, Virus)) + 
    geom_point(colour = AlberColours[2], alpha = 0.3) + geom_smooth(colour = "black") +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2) +
    labs(title = "Sharing ~ Space"),# + ggpubr::stat_cor(method = "spearman"),
  
  ggplot(HostMatrixdf[-HostThemselves,], aes(Phylo, Virus)) + 
    geom_point(colour = AlberColours[3], alpha = 0.3) + geom_smooth(colour = "black", fill = NA) +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2, col = "black") +
    labs(title = "Sharing ~ Phylogeny")# + ggpubr::stat_cor(method = "spearman")
  
), ncol = 3)


arrange_ggplot2(list(
  
  ggplot(HostMatrixdf[-LowerHosts,], aes(Space, PropVirus)) + 
    geom_point(colour = AlberColours[2], alpha = 0.3) + geom_smooth(colour = "black", fill = NA) +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2, col = "black") +
    labs(title = "Sharing ~ Space") + ggpubr::stat_cor(method = "spearman"),
  
  ggplot(HostMatrixdf[-LowerHosts,], aes(Space, PropVirus2)) + 
    geom_point(colour = AlberColours[3], alpha = 0.3) + geom_smooth(colour = "black", fill = NA) +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2, col = "black") +
    labs(title = "Sharing2 ~ Space") + ggpubr::stat_cor(method = "spearman")
  
), ncol = 2)

# Trying this with no zero-space-sharers

jpeg("Figures/Space and Phylogeny correlate with viral sharing with no spatial zeroes.jpeg",
     units = "mm", height = 150, width = 150, res = 300)

list(
  
  HostMatrixdf %>% filter(!Sp==Sp2) %>%  ggplot(aes(Space, Virus)) + 
    geom_point(colour = AlberColours[2], alpha = 0.3) + geom_smooth(colour = "black") +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2) +
    labs(x = "Space Shared", y = "Viruses Shared", title = "Virus ~ Space"),
  
  HostMatrixdf %>% filter(Space>0, !Sp==Sp2) %>%  ggplot(aes(Space, Virus)) + 
    geom_point(colour = AlberColours[2], alpha = 0.3) + geom_smooth(colour = "black") +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2) +
    lims(x = c(0,1)) +
    labs(x = "Space Shared", y = "Viruses Shared", title = "Virus ~ Space, >0 Space"),
  
  HostMatrixdf %>% filter(!Sp==Sp2) %>%  ggplot(aes(Phylo, Virus)) + 
    geom_point(colour = AlberColours[1], alpha = 0.3) + geom_smooth(colour = "black") +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2) +
    labs(x = "Genetic Similarity", y = "Viruses Shared",title = "Virus ~ Phylogeny"),
  
  HostMatrixdf %>% filter(Space>0, !Sp==Sp2) %>%  ggplot(aes(Phylo, Virus)) + 
    geom_point(colour = AlberColours[1], alpha = 0.3) + geom_smooth(colour = "black") +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2) +
    lims(x = c(0,1)) +
    labs(x = "Genetic Similarity", y = "Viruses Shared", title = "Virus ~ Phylogeny, >0 Space")
  
) %>% arrange_ggplot2(nrow = 2)

dev.off()

arrange_ggplot2(list(
  
  ggplot(HostMatrixdf[-HostThemselves,], aes(Space, PropVirus)) + 
    geom_point(colour = AlberColours[2], alpha = 0.3) + geom_smooth(colour = "black", fill = NA) +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2, col = "black") +
    labs(title = "Sharing ~ Space") + ggpubr::stat_cor(method = "spearman"),
  
  ggplot(HostMatrixdf[-HostThemselves,], aes(Space, PropVirus2)) + 
    geom_point(colour = AlberColours[3], alpha = 0.3) + geom_smooth(colour = "black", fill = NA) +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2, col = "black") +
    labs(title = "Sharing2 ~ Space") + ggpubr::stat_cor(method = "spearman")
  
), ncol = 2)

# For the species that don't overlap, does euclidean distance/min distance describe much? 

CentroidDistance <- Hosts[,c("LongMean", "LatMean")] %>% 
  mutate(LongMean = LongMean - min(LongMean, na.rm = T)) %>% dist %>% as.matrix

dimnames(CentroidDistance) <- list(Hosts$Sp, Hosts$Sp)

CentroidDistance[CentroidDistance>(360*50000)&!is.na(CentroidDistance)] <- 
  (360*50000)*2 - CentroidDistance[CentroidDistance>(360*50000)&!is.na(CentroidDistance)]

HostMatrixdf$Euclid <- c(CentroidDistance[FHN,FHN])
HostMatrixdf$PolyEuclid <- poly(HostMatrixdf$Euclid, 2)[,2]

HostMatrixdf  %>% filter(!Sp==Sp2, Euclid>17000000) %>% 
  ggplot(aes(LongMean, LatMean)) + 
  geom_segment(aes(x = LongMean, xend = LongMean.Sp2, y = LatMean, yend = LatMean.Sp2, group = Sp), alpha = 0.01)

HostMatrixdf  %>% filter(!Sp==Sp2) %>%
  ggplot(aes(Euclid, PropVirus)) + facet_wrap(~Space0, scales = "free") + 
  geom_point(colour = AlberColours[2]) + 
  #geom_smooth(method = lm, formula = y ~ poly(x,2), colour = "black") + 
  geom_smooth(colour = "black") + #ggpubr::stat_cor() +
  ggtitle("Hosts that are far away begin to share the same viruses")

HostMatrixdf %>% filter(!Sp == Sp2) %>% ggplot(aes(Phylo, Virus, colour = SpaceQuantile)) + 
  geom_point(alpha = 0.3) + geom_smooth(method = lm) + facet_wrap(~SpaceQuantile)
