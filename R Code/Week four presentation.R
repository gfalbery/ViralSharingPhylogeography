# Work for week four presentation ####

library(ggplot2); library(igraph); library(ggregplot)

# Checking convergence ####

i=1

mc <- BinModelList[1:10 + 10*(i-1)] %>% lapply(function(a) a$Sol[,1:14])
#mc <- ZI_runs[1:10 + 10*(i-1)] %>% lapply(function(a) a$VCV)
mc <- do.call(mcmc.list, mc)
#par(mfrow=c(7,2), mar=c(2,2,1,2), ask = F)
#gelman.plot(mc, auto.layout=F)
gelman.diag(mc)
par(mfrow=c(7,4), mar=c(2, 1, 1, 1))
plot(mc, ask=F, auto.layout=F)

# Looking at predicted 1's and 0's
list(
  BarGraph(FinalHostMatrix, "PredVirus1Q","VirusBinary", text = "N") + theme(legend.position = "none"), 
  BarGraph(FinalHostMatrix, "PredVirus1bQ","VirusBinary", text = "N") + theme(legend.position = "none"),
  BarGraph(FinalHostMatrix, "PredVirus2Q","VirusBinary", text = "N") + theme(legend.position = "none")
) %>% arrange_ggplot2

df1 <- with(FinalHostMatrix, tapply(VirusBinary, PredVirus1Q, function(a) table(a)/sum(table(a)))) %>%
  bind_rows()

df1b <- with(FinalHostMatrix, tapply(VirusBinary, PredVirus1bQ, function(a) table(a)/sum(table(a)))) %>%
  bind_rows()

df2 <- with(FinalHostMatrix, tapply(VirusBinary, PredVirus2Q, function(a) table(a)/sum(table(a)))) %>%
  bind_rows()

dflong <- rbind(reshape2::melt(df1), reshape2::melt(df1b), reshape2::melt(df2))
dflong$Var2 <- rep(0:1)
dflong$Var3 <- rep(c("1","1b","2"), each = nrow(df1))

ggplot(dflong, aes(variable, value, fill = as.factor(Var2))) + 
  geom_col(position = "stack") + 
  coord_fixed(10) +
  labs(x= "Predicted value", y = "Frequency", fill = "Observed") +
  facet_wrap(~Var3, ncol = 3) +
  ggsave("Quantiles of predicted values.tiff", units = "mm", width = 200, height = 100, dpi = 300)

# Degree predictions ####

list(
  ggplot(Hosts, aes(Degree, PredDegree1)) + geom_point() + geom_smooth() +
    labs(x = "Observed Degree", y = "Predicted Degree") +
    coord_fixed() +
    ggtitle("Full Model"),
  
  ggplot(Hosts, aes(Degree, PredDegree1b)) + geom_point() + geom_smooth() +
    labs(x = "Observed Degree", y = "Predicted Degree") +
    coord_fixed() +
    ggtitle("Full Model (no random effects)"),
  
  ggplot(Hosts, aes(Degree, PredDegree2)) + geom_point() + geom_smooth() +
    labs(x = "Observed Degree", y = "Predicted Degree") +
    coord_fixed() +
    ggtitle("Model without random effects"),
  
  ggplot(Hosts, aes(Degree, PredDegree3)) + geom_point() + geom_smooth() +
    labs(x = "Observed Degree", y = "Predicted Degree") +
    coord_fixed() +
    ggtitle("Full Model, max citations")
) %>% arrange_ggplot2(nrow = 2)

ggplot(Hosts, aes(Degree, PredDegree3b)) + geom_point() + geom_smooth() +
  labs(x = "Observed Degree", y = "Predicted Degree") +
  coord_fixed() +
  ggtitle("Full Model, no random effect, max citations")

# Eigenvector predictions ####

ggplot(Hosts, aes(kader:::cuberoot(Eigenvector), kader:::cuberoot(PredEigen1))) + geom_point() + geom_smooth() +
  labs(x = "Observed Eigenvector", y = "Predicted Eigenvector") +
  coord_fixed() +
  ggtitle("Full Model")

ggplot(Hosts, aes(kader:::cuberoot(Eigenvector), kader:::cuberoot(PredEigen1b))) + geom_point() + geom_smooth() +
  labs(x = "Observed Eigenvector", y = "Predicted Eigenvector") +
  coord_fixed() +
  ggtitle("Full Model (no random effects)")

ggplot(Hosts, aes(kader:::cuberoot(Eigenvector), kader:::cuberoot(PredEigen2))) + geom_point() + geom_smooth() +
  labs(x = "Observed Eigenvector", y = "Predicted Eigenvector") +
  coord_fixed() +
  ggtitle("Model without random effects")

ggplot(Hosts, aes(kader:::cuberoot(Eigenvector), kader:::cuberoot(PredEigen3))) + geom_point() + geom_smooth() +
  labs(x = "Observed Eigenvector", y = "Predicted Eigenvector") +
  coord_fixed() +
  ggtitle("Full Model, max citations")

ggplot(Hosts, aes(kader:::cuberoot(Eigenvector), kader:::cuberoot(PredEigen3b))) + geom_point() + geom_smooth() +
  labs(x = "Observed Eigenvector", y = "Predicted Eigenvector") +
  coord_fixed() +
  ggtitle("Full Model, no random effect, max citations")

ggplot(Hosts, aes(kader:::cuberoot(Eigenvector), kader:::cuberoot(PredEigen3))) + geom_point() + geom_smooth() +
  labs(x = "Observed Eigenvector", y = "Predicted Eigenvector") +
  coord_fixed() +
  ggtitle("Full Model, no random effect, max citations") + facet_wrap(~hOrder, ncol = 8) +
  lims(y = c(0, 1))

ggplot(Hosts, aes(kader:::cuberoot(Eigenvector), kader:::cuberoot(PredEigen3b))) + geom_point() + geom_smooth() +
  labs(x = "Observed Eigenvector", y = "Predicted Eigenvector") +
  coord_fixed() +
  ggtitle("Full Model, no random effect, max citations") + facet_wrap(~hOrder, ncol = 8) +
  lims(y = c(0, 1))

# Getting network-level stats

SimGraphList <- list(SimGraphs1, SimGraphs1b, SimGraphs2, SimGraphs3, SimGraphs3b)

Components <- lapply(SimGraphList, function(a) sapply(a, function(b) components(b)$no))
Degrees <- lapply(SimGraphList, function(a) sapply(a, function(b) mean(degree(b))))

Cluster1 = sapply(SimGraphs1, function(a) transitivity(a)) # all zero, don't bother
Cluster1b = sapply(SimGraphs1b, function(a) transitivity(a))
Cluster2 = sapply(SimGraphs2, function(a) transitivity(a))
Cluster3 = sapply(SimGraphs3, function(a) transitivity(a))
Cluster3b = sapply(SimGraphs3b, function(a) transitivity(a))

Betweenness1 = sapply(SimGraphs1, function(a) mean(betweenness(a)))
Betweenness1b = sapply(SimGraphs1b, function(a) mean(betweenness(a)))
Betweenness2 = sapply(SimGraphs2, function(a) mean(betweenness(a)))
Betweenness3 = sapply(SimGraphs3, function(a) mean(betweenness(a)))
Betweenness3b = sapply(SimGraphs3b, function(a) mean(betweenness(a)))

NetworkCompdf <- data.frame(
  Degree = unlist(Degrees),
  Components = unlist(Components)
  #Between = c(Betweenness1,Betweenness1b,Betweenness2,Betweenness3,Betweenness3b)
) %>% mutate(Model = rep(c("1","1b","2", "3", "3b"), each = 1000))

ggplot(NetworkCompdf, aes(Degree)) + 
  geom_density(aes(colour = Model)) + 
  geom_vline(aes(xintercept = mean(Hosts$Degree)))

ggplot(NetworkCompdf, aes(Components)) + 
  geom_density(aes(colour = Model)) + 
  geom_vline(aes(xintercept = components(Hostgraph)$no))




