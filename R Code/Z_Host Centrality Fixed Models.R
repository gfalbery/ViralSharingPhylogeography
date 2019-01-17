
# Modelling Centrality Measures ####

library(INLA); library(tidyverse); library(ggregplot)

# Trying viral records as the response ####

Resp = "Records"

RecordsHostCentCovar = c(
  "hDom",
  "hAllZACites",
  #"hDiseaseZACites",
  "GeogRange",
  "S.Greg1",
  "hArtfclHbttUsrIUCN",
  "hOrder",
  "hMassGrams"
)

TestHosts <- Hosts %>% select(Resp, RecordsHostCentCovar) %>% #, "LongMean", "LatMean") %>%
  
  mutate(hAllZACites = log(hAllZACites + 1),
         #hDiseaseZACites = log(hDiseaseZACites + 1),
         GeogRange = kader:::cuberoot(GeogRange), 
         S.Greg1 = sqrt(S.Greg1),
         Records = c(Records - 1), 
         hMassGrams = log(hMassGrams)) %>%
  
  slice(which(!NARows(Hosts[,c(Resp,RecordsHostCentCovar, "LongMean", "LatMean")])))

IM1 <- INLAModelAdd(Resp, "1", RecordsHostCentCovar, Family = "nbinomial", Data = TestHosts)

RecordsKeptCovar <- RecordsHostCentCovar[-which(RecordsHostCentCovar%in%last(IM1$Removed))]

# Degree is the easiest to incorporate into network generation so:

Resp = "Degree"

DegreeHostCentCovar = c(
  "hAllZACites",
  #"hDiseaseZACites",
  "GeogRange",
  "S.Greg1",
  "hArtfclHbttUsrIUCN",
  "hMassGrams",
  "hDom",
  "hOrder"
)

TestHosts <- Hosts %>% select(Resp, DegreeHostCentCovar, "LongMean", "LatMean") %>%
  
  mutate(hAllZACites = log(hAllZACites + 1),
         #hDiseaseZACites = log(hDiseaseZACites + 1),
         GeogRange = kader:::cuberoot(GeogRange), 
         S.Greg1 = sqrt(S.Greg1),
         hMassGrams = log(hMassGrams)) %>%
  
  slice(which(!NARows(Hosts[,c(Resp,DegreeHostCentCovar, "LongMean", "LatMean")])))

IM2 <- INLAModelAdd(Resp, "1", DegreeHostCentCovar, Family = "nbinomial", Data = TestHosts)

DegreeKeptCovar <- DegreeHostCentCovar[-which(DegreeHostCentCovar%in%last(IM2$Removed))]

# What if records is put as an explanatory variable? ####

Resp = "Degree"

DegreeHostCentCovar2 = c(
  "hAllZACites",
  #"hDiseaseZACites",
  "GeogRange",
  "S.Greg1",
  "hArtfclHbttUsrIUCN",
  "hMassGrams",
  "hDom",
  "hOrder",
  "Records"
)

TestHosts <- Hosts %>% select(Resp, DegreeHostCentCovar2, "LongMean", "LatMean") %>%
  
  mutate(hAllZACites = log(hAllZACites + 1),
         #hDiseaseZACites = log(hDiseaseZACites + 1),
         GeogRange = kader:::cuberoot(GeogRange), 
         S.Greg1 = sqrt(S.Greg1),
         Records = c(log(Records)),
         hMassGrams = log(hMassGrams)) %>%
  
  slice(which(!NARows(Hosts[,c(Resp,DegreeHostCentCovar2, "LongMean", "LatMean")])))

IM3 <- INLAModelAdd(Resp, "1", DegreeHostCentCovar2, Family = "poisson", Data = TestHosts)

DegreeKeptCovar2 <- DegreeHostCentCovar2[-which(DegreeHostCentCovar2%in%last(IM3$Removed))]

# Comparing the output to all these ####

Efxplot(lapply(list(IM1, IM2, IM3), function(a) a$FinalModel))
Efxplot(lapply(list(IM1, IM2, IM3, IM4), function(a) a$FinalModel))

# Now trying Eigenvector j4j

Resp = "Eigenvector"

EigenHostCentCovar = c(
  "hAllZACites",
  #"hDiseaseZACites",
  "GeogRange",
  "S.Greg1",
  "hArtfclHbttUsrIUCN",
  "hMassGrams",
  "hDom",
  "hOrder"
)

TestHosts <- Hosts %>% select(Resp, EigenHostCentCovar, "LongMean", "LatMean") %>%
  
  mutate(hAllZACites = log(hAllZACites + 1),
         #hDiseaseZACites = log(hDiseaseZACites + 1),
         GeogRange = kader:::cuberoot(GeogRange), 
         S.Greg1 = sqrt(S.Greg1),
         hMassGrams = log(hMassGrams)) %>%
  
  slice(which(!NARows(Hosts[,c(Resp,EigenHostCentCovar, "LongMean", "LatMean")])))

IM4 <- INLAModelAdd(Resp, "1", EigenHostCentCovar, Family = "beta", Data = TestHosts)

EigenKeptCovar <- EigenHostCentCovar[-which(EigenHostCentCovar%in%last(IM2$Removed))]

# Adding Phylo or Space Matrices ####

SpaceContacts <- as(solve(RangeAdj1[FHN, FHN]),"dgCMatrix")
GRMatrix <- as(solve(1-CytBMatrix[FHN, FHN]),"dgCMatrix")

# NB phylo matrix has been inverted

FHosts <- Hosts[Hosts$Sp %in% FHN,]

FHosts$IndexSpace = unlist(sapply(FHosts$Sp, function(a) which(FHN==a)))
FHosts$IndexPhylo = unlist(sapply(FHosts$Sp, function(a) which(FHN==a)))

Resp = "Degree"

TestHosts <- FHosts %>% select(Resp, DegreeKeptCovar, "LongMean", "LatMean", 
                              "IndexSpace", "IndexPhylo") %>%
  
  mutate(hAllZACites = log(hAllZACites + 1),
         #hDiseaseZACites = log(hDiseaseZACites + 1),
         GeogRange = kader:::cuberoot(GeogRange), 
         S.Greg1 = sqrt(S.Greg1),
         #hMassGrams = log(hMassGrams),
         Degree = c(Degree)) %>%
  
  slice(which(!NARows(FHosts[,c(Resp,DegreeHostCentCovar, "LongMean", "LatMean")])))

HostMMModelsFixed <- list()

f1 = as.formula(paste("Degree ~ ", paste(DegreeKeptCovar, collapse = " + ")))

HostMMModelsFixed[[1]] <- inla(f1, 
                               data = TestHosts, 
                               family = "gaussian",
                               control.compute = list(dic = T))


f2 =  as.formula(paste("Degree ~ ", paste(DegreeKeptCovar, collapse = " + "), " + ", 'f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
                       constr=TRUE,param = c(0.5, 0.5))'))

HostMMModelsFixed[[2]] <- inla(f2, 
                               data = TestHosts, 
                               family = "gaussian",
                               control.compute = list(dic = T))


f3 =  as.formula(paste("Degree ~ ", paste(DegreeKeptCovar[-which(DegreeKeptCovar == "hOrder")], collapse = " + "), " + ", 'f(IndexPhylo, model="generic0", Cmatrix = GRMatrix,
                       constr=TRUE,param = c(0.5, 0.5))'))

HostMMModelsFixed[[3]] <- inla(f3, 
                               data = TestHosts, 
                               family = "gaussian",
                               control.compute = list(dic = T))


f2.2 =  as.formula(paste("Degree ~ ", paste(DegreeKeptCovar[-which(DegreeKeptCovar == "hOrder")], collapse = " + "), " + ", 'f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
                       constr=TRUE,param = c(0.5, 0.5))'))

HostMMModelsFixed[[4]] <- inla(f2.2, 
                               data = TestHosts, 
                               family = "gaussian",
                               control.compute = list(dic = T))

f4 = hEigenvector ~ f(IndexPhylo, model = "generic0", Cmatrix = GRMatrix,
                      constr = TRUE,param = c(0.5, 0.5)) + 
  f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
    constr=TRUE,param = c(0.5, 0.5))

HostMMModelsFixed[[4]] <- inla(f4, # Takes about an hour ####
                               data = TestHosts, 
                               family = "gaussian",
                               control.compute = list(dic = T))

sapply(HostMMModelsFixed, function(a) a$dic$dic)



