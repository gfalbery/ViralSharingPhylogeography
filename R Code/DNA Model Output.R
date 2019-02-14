
# STAN Model Output ####

library(rstan); library(reskew); library(ggregplot); library(parallel)

if(file.exists("~/Albersnet/DNAModelOutput.Rdata")) load("~/Albersnet/DNAModelOutput.Rdata") else{
  
  DNABinModel <- readRDS("~/Albersnet/DNABinModel.rds")
  r <- process_stanfit(DNABinModel, n.pars.to.trim = 3) # Takes a WHILE
  save(r, file = "~/Albersnet/DNAModelOutput.Rdata")
}

f <- FinalHostMatrix %>% filter(!is.na(DNA))

f$Sp <- factor(as.character(f$Sp),
               levels = union(f$Sp, f$Sp2)
)

f$Sp2 <- factor(as.character(f$Sp2),
                levels = union(f$Sp, f$Sp2)
)

UpperDNA <- which(upper.tri(matrix(NA, ncol = nlevels(f$Sp), nrow = nlevels(f$Sp2))))

MCMCSol <- r$df#[rep(501:1000, 8)+rep(0:7*1000, each = 500),]

N = nrow(f)

f$Space_Phylo <- scale(f$Space*f$Phylo2)
f$Space <- scale(f$Space)
f$Space_Phylo <- scale(f$Phylo2)

f <- f %>% mutate(DCites = (log(hDiseaseZACites + 1)),#-mean(species.traits$d_cites))/sd(species.traits$d_cites), 
                  DCites.Sp2 = (log(hDiseaseZACites.Sp2 + 1)))#-mean(species.traits$d_cites))/sd(species.traits$d_cites))

# DNASimulating with specific random effects ####

XMatrix <- cbind(rep(1,N),
                 f[,c("Space","Phylo2","Space_Phylo")]) %>% as.matrix %>% as("dgCMatrix")

MZ1 <- model.matrix( ~ Sp - 1, data = f)
MZ2 <- model.matrix( ~ Sp2 - 1, data = f)

ZMatrixb <- MZ1 + MZ2 %>% as.matrix %>% as("dgCMatrix")
XBetas <- c("mu_alpha","beta_space","beta_phylo","beta_inter")

ZBetas2 <- colnames(r$df)[which(colnames(r$df)=="alpha_species[1]"):
                            which(colnames(r$df)=="alpha_species[164]")]

# Doing the DNASimulating #####

DNAPredList1 <- list()

RowsSampled <- sample(1:nrow(MCMCSol), 1000, replace = F)

DNAPredList1 <- parallel::mclapply(1:length(RowsSampled), function(x){
  
  if(x %% 10 == 0) print(x)
  
  XFX <- MCMCSol[RowsSampled[x], XBetas] %>% unlist
  
  ZFX <- MCMCSol[RowsSampled[x], ZBetas2] %>% unlist
  
  XDNAPredictions <- XFX %*% t(XMatrix)
  
  ZDNAPredictions <- ZFX %*% t(ZMatrixb)
  
  DNAPredictions <- XDNAPredictions + ZDNAPredictions
  
  BinDNAPred <- rbinom(n = N,
                       prob = logistic(DNAPredictions@x),
                       size  = 1)
  
  BinDNAPred
  
}, mc.cores = 10)

DNAPredDF1 <- data.frame(DNAPredList1)
f$DNAPredVirus1 <- apply(DNAPredDF1, 1, mean)
f$DNAPredVirus1Q <- cut(f$DNAPredVirus1,
                        breaks = c(-1:10/10),
                        labels = c(0:10/10))

# DNASimulating without random effects ####

XMatrix <- cbind(rep(1,N),
                 f[,c("Space","Phylo2","Space_Phylo")]) %>% as.matrix %>% as("dgCMatrix")

ZMatrix <- f[,c("DCites", "hDom","DCites.Sp2","hDom.Sp2")] %>% 
  mutate(hDom = ifelse(hDom == "wild", 0, 1), 
         hDom.Sp2 = ifelse(hDom.Sp2 == "wild", 0, 1)) %>%
  as.matrix %>% as("dgCMatrix")

XBetas <- c("mu_alpha","beta_space","beta_phylo","beta_inter")
ZBetas <- c("beta_d_cites_s","beta_domestic")

# Doing the DNASimulating #####

DNAPredList1b <- list()

RowsSampled <- sample(1:nrow(r$df), 1000, replace = F)

DNAPredList1b <- parallel::mclapply(1:length(RowsSampled), function(x){ # to do something non-specific
  
  XFX <- r$df[RowsSampled[x], XBetas] %>% unlist
  
  ZFX <- r$df[RowsSampled[x], ZBetas] %>% unlist
  
  XDNAPredictions <- XFX %*% t(XMatrix)
  
  ZDNAPredictionsa <- ZFX %*% t(ZMatrix[,1:length(ZBetas)])
  ZDNAPredictionsb <- ZFX %*% t(ZMatrix[,(length(ZBetas)+1):(length(ZBetas)*2)])
  
  ZDNAPredictions2a <- rnorm(n = N, mean = ZDNAPredictionsa@x, sd = r$df[RowsSampled[x], "sigma"])
  ZDNAPredictions2b <- rnorm(n = N, mean = ZDNAPredictionsb@x, sd = r$df[RowsSampled[x], "sigma"])
  #ZDNAPredictions2c <- rnorm(n = N, mean = (ZDNAPredictionsa@x + ZDNAPredictionsb@x), sd = r$df[RowsSampled[x], "sigma"])
  
  DNAPredictions <- XDNAPredictions@x + ZDNAPredictions2a + ZDNAPredictions2b
  #DNAPredictions2 <- XDNAPredictions@x + ZDNAPredictions2c
  
  BinDNAPred <- rbinom(n = N,
                       prob = logistic(DNAPredictions),
                       size  = 1)
  
  BinDNAPred
  
}, mc.cores = 10)

DNAPredDF1b <- data.frame(DNAPredList1b)
f$DNAPredVirus1b <- apply(DNAPredDF1b, 1, mean)
f$DNAPredVirus1bQ <- cut(f$DNAPredVirus1b,
                         breaks = c(-1:10/10),
                         labels = c(0:10/10))

# DNASimulating using the random effect ####

DNASimNets1 <- mclapply(1:length(DNAPredList1), function(i){
  
  if(i%%10==0) print(i)
  
  AssMat <- matrix(NA, 
                   nrow = length(union(f$Sp,f$Sp2)), 
                   ncol = length(union(f$Sp,f$Sp2)))
  
  AssMat[-which(1:length(AssMat)%in%UpperDNA)] <- round(DNAPredList1[[i]])
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
  diag(AssMat) <- apply(AssMat,1,function(a) length(a[!is.na(a)&a>0]))
  dimnames(AssMat) <- list(union(f$Sp,f$Sp2),
                           union(f$Sp,f$Sp2))
  
  as(AssMat, "dgCMatrix")
  
}, mc.cores = 10)

DNASimGraphs1 <- mclapply(1:length(DNAPredList1), function(i){
  
  graph.adjacency(DNASimNets1[[i]], mode = "undirected", diag = F)
  
}, mc.cores = 10)

DNADegdf1 <- sapply(DNASimGraphs1, function(a) degree(a))# %>% as.data.frame
DNAPredDegrees1 <- apply(DNADegdf1, 1, mean)
Hosts$DNAPredDegree1 <- DNAPredDegrees1[as.character(Hosts$Sp)]

#Eigendf1 <- sapply(DNASimGraphs1, function(a) eigen_centrality(a)$vector) #%>% as.data.frame
#DNAPredEigen1 <- apply(Eigendf1, 1, mean)
#Hosts$DNAPredEigen1 <- DNAPredEigen1[as.character(Hosts$Sp)]

# DNASimulating without the random effect ####

DNASimNets1b <- mclapply(1:length(DNAPredList1b), function(i){
  
  
  AssMat <- matrix(NA, 
                   nrow = length(union(f$Sp,f$Sp2)), 
                   ncol = length(union(f$Sp,f$Sp2)))
  
  AssMat[-which(1:length(AssMat)%in%UpperDNA)] <- round(DNAPredList1b[[i]])# %>% as.matrix %>% as("dgCMatrix")
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))] #%>% as.matrix %>% as("dgCMatrix")
  diag(AssMat) <- 0
  dimnames(AssMat) <- list(union(f$Sp,f$Sp2),
                           union(f$Sp,f$Sp2))
  
  as(AssMat, "dgCMatrix")
  
}, mc.cores = 10)

DNASimGraphs1b <- mclapply(1:length(DNAPredList1b), function(i){
  
  graph.adjacency(DNASimNets1b[[i]], mode = "undirected")
  
}, mc.cores = 10)

DNADegdf1b <- sapply(DNASimGraphs1b, function(a) degree(a)) %>% as.data.frame
DNAPredDegrees1b <- apply(DNADegdf1b, 1, mean)
Hosts$DNAPredDegree1b <- DNAPredDegrees1b[as.character(Hosts$Sp)]

Eigendf1b <- sapply(DNASimGraphs1b, function(a) eigen_centrality(a)$vector) %>% as.data.frame
DNAPredEigen1b <- apply(Eigendf1b, 1, mean)
Hosts$DNAPredEigen1b <- DNAPredEigen1b[as.character(Hosts$Sp)]

# Getting network-level stats

DNASimGraphList <- list(DNASimGraphs1, DNASimGraphs1b)# , DNASimGraphs2, DNASimGraphs3, DNASimGraphs3b)

Components <- lapply(DNASimGraphList, function(a) sapply(a, function(b) components(b)$no))
ComponentSizes <- lapply(DNASimGraphList, function(a) sapply(a, function(b) components(b)$no))
lapply(DNASimGraphs1b[which(Components[[2]]==2)], function(a) which(components(a)$membership==2))

Degrees <- lapply(DNASimGraphList, function(a) sapply(a, function(b) mean(degree(b))))

Cluster1 = sapply(DNASimGraphs1, function(a) transitivity(a)) # all zero, don't bother
Cluster1b = sapply(DNASimGraphs1b, function(a) transitivity(a))

#Betweenness1 = sapply(DNASimGraphs1, function(a) betweenness(a))
#Betweenness1b = sapply(DNASimGraphs1b, function(a) betweenness(a))

Prev1 = sapply(DNAPredList1, Prev)
Prev1b = sapply(DNAPredList1b, Prev)

#Closeness1 = sapply(DNASimGraphs1, function(a) closeness(a))
#Closeness1b = sapply(DNASimGraphs1b, function(a) closeness(a))

