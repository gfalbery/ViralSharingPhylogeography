
# STAN Model Output ####

library(rstan); library(reskew)

Prev <- function(x){
  
  length(x[x>0])/length(x)
  
}

logistic <- function(a) exp(a)/(1 + exp(a))
logit <- function(a) log(a/(1-a))

traceplot(BinModel,
          pars = c("mu_alpha", 
                   "beta_d_cites_s", 
                   "beta_domestic", 
                   "beta_space",
                   "beta_phylo",
                   "beta_inter",
                   "sigma"))

p <- process_stanfit(BinModel, n.pars.to.trim = 3) # Takes a WHILE
#p <- process_stanfit(BinModel, n.pars.to.trim = 3) # Takes a WHILE

N = nrow(d)

d$Space_Phylo <- d$Space*d$Phylo2

d <- d %>% mutate(DCites = log(hDiseaseZACites + 1), DCites.Sp2 = log(hDiseaseZACites.Sp2 + 1))

# Simulating with specific random effects ####

XMatrix <- cbind(rep(1,N),
                 d[,c("Space","Phylo2","Space_Phylo")]) %>% as.matrix %>% as("dgCMatrix")

MZ1 <- model.matrix( ~ Sp - 1, data = FinalHostMatrix)
MZ2 <- model.matrix( ~ Sp2 - 1, data = FinalHostMatrix)

ZMatrixb <- MZ1 + MZ2 %>% as.matrix %>% as("dgCMatrix")
XBetas <- c("mu_alpha","beta_space","beta_phylo","beta_inter")

ZBetas2 <- colnames(p$df)[which(colnames(p$df)=="alpha_species[1]"):
                   which(colnames(p$df)=="alpha_species[649]")]

# Doing the simulating #####

PredList1 <- list()

RowsSampled <- sample(1:nrow(p$df), 1000, replace = F)

for(x in 1:length(RowsSampled)){ # to do something non-specific
  
  if(x %% 10 == 0) print(x)
  
  XFX <- p$df[RowsSampled[x], XBetas] %>% unlist
  
  ZFX <- p$df[RowsSampled[x], ZBetas2] %>% unlist
  
  XPredictions <- XFX %*% t(XMatrix)
  
  ZPredictions <- ZFX %*% t(ZMatrixb)
  
  Predictions <- XPredictions + ZPredictions
  
  BinPred <- rbinom(n = N,
                    prob = logistic(Predictions@x),
                    size  = 1)
  
  PredList1[[x]] <- BinPred
  
}

PredDF1 <- data.frame(PredList1)

sapply(PredList1, Prev) %>% mean
sapply(PredList1, Prev) %>% qplot

d$PredVirus1 <- apply(PredDF1, 1, mean)

d$PredVirus1Q <- cut(d$PredVirus1,
                     breaks = c(-1:10/10),
                     labels = c(0:10/10))

# Simulating without random effects ####

XMatrix <- cbind(rep(1,N),
                 d[,c("Space","Phylo2","Space_Phylo")]) %>% as.matrix %>% as("dgCMatrix")

ZMatrix <- d[,c("DCites", "hDom","DCites.Sp2","hDom.Sp2")] %>% 
  mutate(hDom = ifelse(hDom == "wild", 0, 1), 
         hDom.Sp2 = ifelse(hDom.Sp2 == "wild", 0, 1)) %>%
  as.matrix %>% as("dgCMatrix")

XBetas <- c("mu_alpha","beta_space","beta_phylo","beta_inter")
ZBetas <- c("beta_d_cites_s","beta_domestic")

# Doing the simulating #####

PredList1b <- PredList1c <- list()

RowsSampled <- sample(1:nrow(p$df), 1000, replace = F)

for(x in 1:length(RowsSampled)){ # to do something non-specific
  
  if(x %% 10 == 0) print(x)
  
  XFX <- p$df[RowsSampled[x], XBetas] %>% unlist
  
  ZFX <- p$df[RowsSampled[x], ZBetas] %>% unlist
  
  XPredictions <- XFX %*% t(XMatrix)
  
  ZPredictionsa <- ZFX %*% t(ZMatrix[,1:length(ZBetas)])
  ZPredictionsb <- ZFX %*% t(ZMatrix[,(length(ZBetas)+1):(length(ZBetas)*2)])
  
  ZPredictions2a <- rnorm(n = N, mean = ZPredictionsa@x, sd = p$df[RowsSampled[x], "sigma"])
  ZPredictions2b <- rnorm(n = N, mean = ZPredictionsb@x, sd = p$df[RowsSampled[x], "sigma"])
  ZPredictions2c <- rnorm(n = N, mean = (ZPredictionsa@x + ZPredictionsb@x), sd = p$df[RowsSampled[x], "sigma"])
  
  Predictions <- XPredictions@x + ZPredictionsa + ZPredictionsb
  Predictions <- XPredictions@x + ZPredictions2a + ZPredictions2b
  Predictions2 <- XPredictions@x + ZPredictions2c
  
  BinPred <- rbinom(n = N,
                    prob = logistic(Predictions),
                    size  = 1)
  
  PredList1b[[x]] <- BinPred
  
  BinPred2 <- rbinom(n = N,
                    prob = logistic(Predictions2),
                    size  = 1)
  
  PredList1c[[x]] <- BinPred2
  
}

PredDF1b <- data.frame(PredList1b)

sapply(PredList1b, Prev) %>% mean
sapply(PredList1b, Prev) %>% qplot

sapply(PredList1c, Prev) %>% mean
sapply(PredList1c, Prev) %>% qplot


d$PredVirus1b <- apply(PredDF1b, 1, mean)

d$PredVirus1bQ <- cut(d$PredVirus1b,
                     breaks = c(-1:10/10),
                     labels = c(0:10/10))

# Simulating using the random effect ####

SimNets1 <- SimGraphs1 <- list()

for(i in 1:length(PredList1)){
  
  if(i%%10==0) print(i)
  
  AssMat <- matrix(NA, 
                   nrow = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)), 
                   ncol = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
  
  AssMat[-which(1:length(AssMat)%in%UpperHosts)] <- round(PredList1[[i]])
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
  diag(AssMat) <- apply(AssMat,1,function(a) length(a[!is.na(a)&a>0]))
  dimnames(AssMat) <- list(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2),
                           union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))
  
  SimNets1[[i]] <- as(AssMat, "dgCMatrix")
  
  SimGraphs1[[i]] <- graph.incidence(AssMat, weighted = TRUE)
  
}

Degdf1 <- sapply(SimGraphs1, function(a) degree(a)) %>% as.data.frame
Eigendf1 <- sapply(SimGraphs1, function(a) eigen_centrality(a)$vector) %>% as.data.frame

PredDegrees1 <- apply(Degdf1, 1, mean)
PredEigen1 <- apply(Eigendf1, 1, mean)

Hosts$PredDegree1 <- PredDegrees1[as.character(Hosts$Sp)]
Hosts$PredEigen1 <- PredEigen1[as.character(Hosts$Sp)]

# Simulating without the random effect ####

SimNets1b <- SimGraphs1b <- list()

for(i in 1:length(PredList1b)){
  
  if(i%%10==0) print(i)
  
  AssMat <- matrix(NA, 
                   nrow = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)), 
                   ncol = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
  
  AssMat[-which(1:length(AssMat)%in%UpperHosts)] <- round(PredList1b[[i]])# %>% as.matrix %>% as("dgCMatrix")
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))] #%>% as.matrix %>% as("dgCMatrix")
  diag(AssMat) <- apply(AssMat,1,function(a) length(a[!is.na(a)&a>0]))
  dimnames(AssMat) <- list(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2),
                           union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))
  
  SimNets1b[[i]] <- as(AssMat, "dgCMatrix")
  
  SimGraphs1b[[i]] <- graph.incidence(AssMat, weighted = TRUE)
  
}

Degdf1b <- sapply(SimGraphs1b, function(a) degree(a)) %>% as.data.frame
Eigendf1b <- sapply(SimGraphs1b, function(a) eigen_centrality(a)$vector) %>% as.data.frame

PredDegrees1b <- apply(Degdf1b, 1, mean)
PredEigen1b <- apply(Eigendf1b, 1, mean)

Hosts$PredDegree1b <- PredDegrees1b[as.character(Hosts$Sp)]
Hosts$PredEigen1b <- PredEigen1b[as.character(Hosts$Sp)]

ggplot(Hosts, aes(Degree, PredDegree1b)) + geom_point() + geom_smooth()

# Packge together species-level varying effect estimates

species.means.df <- data.frame(
  species_mean = species.means,
  species_level = 1:length(species.means)
)

species.traits.mod <- species.traits %>%
  mutate(sp = as.integer(sp))

# Create a predictions data frame

preds <- data.frame(
  sp = d$Sp,
  sp_level = as.integer(d$Sp),
  sp2 = d$Sp2,
  sp2_level = as.integer(d$Sp2),
  # Mean overall intercept fit from the model
  mu_alpha = rep(mean(p$df$mu_alpha), nrow(d))
) %>%
  left_join(., species.means.df, by = c("sp_level" = "species_level")) %>%
  left_join(., species.means.df, by = c("sp2_level" = "species_level")) %>%
  left_join(., species.traits.mod, by = c("sp_level" = "sp")) %>%
  left_join(., species.traits.mod, by = c("sp2_level" = "sp")) %>%
  rename(
    # Mean species-level varying effects fit from the model
    fit_alpha_sp1 = species_mean.x,
    fit_alpha_sp2 = species_mean.y,
    # Species trait data
    domestic_sp1 = domestic.x,
    domestic_sp2 = domestic.y,
    d_cites_standardized_sp1 = d_cites_standardized.x,
    d_cites_standardized_sp2 = d_cites_standardized.y
  ) %>%
  mutate(
    # Manually calculate linear term that feeds into the varying effect
    # for each species. Note: will not match the fitted values, since the
    # fitted values are also subject to residual variation
    linear_term_alpha_sp1 =
      (domestic_sp1 * mean(p$df$beta_domestic) + 
         d_cites_standardized_sp1 * mean(p$df$beta_d_cites_s)),
    linear_term_alpha_sp2 =
      (domestic_sp2 * mean(p$df$beta_domestic) + 
         d_cites_standardized_sp2 * mean(p$df$beta_d_cites_s)),
    # Generate predicted varying effect values for species just like ours
    # (including the residual variation)
    pred_alpha_sp1 = rnorm(linear_term_alpha_sp1, sd = p$df$sigma),
    pred_alpha_sp2 = rnorm(linear_term_alpha_sp2, sd = p$df$sigma)
  )

# Generate predictions

# Predictions using mean fitted varying effect values from the model
preds$pred <- 
  rbinom(nrow(preds), 
         prob = logistic(preds$mu_alpha + preds$fit_alpha_sp1 + preds$fit_alpha_sp2), 
         size = 1
  )

# Predictions using manually calculated linear terms for our species
# (not subject to residual variation)
preds$pred2 <-
  rbinom(nrow(preds),
         prob = logistic(preds$mu_alpha + preds$linear_term_alpha_sp1 + preds$linear_term_alpha_sp2),
         size = 1
  )

# Predictions using manually calculated linear terms for our species
# (and the influence of residual variation)
preds$pred3 <-
  rbinom(nrow(preds),
         prob = logistic(preds$mu_alpha + preds$pred_alpha_sp1 + preds$pred_alpha_sp2),
         size = 1
  )


# How many viral sharing links are observed?

sum(d$VirusBinary)/nrow(d)

sum(preds$pred)/nrow(preds)
sum(preds$pred2)/nrow(preds)
sum(preds$pred3)/nrow(preds)