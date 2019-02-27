
# trying to plot GAMs ####

#source("R Code/00_Master Code.R")

library('mgcv')
library('brms')
library('ggplot2')
library('schoenberg')
library('bayesplot')

DNAGAM <- readRDS("Model Files/DNAGAMNoG.rds")

pp_check(DNAGAM, nsamples = 30)
pp_check(DNAGAM, nsamples = 30, type = "ecdf_overlay")

marginal_smooths(DNAGAM)
plot(marginal_smooths(DNAGAM), stype = "raster")

marginal_effects(DNAGAM)

SubResps <- c("RNA","Vector","NVector","DNA")
SubDataList <- StanDataList <- SubGAMModelList <- list()

r = 4

SubDataList[[r]] <- FinalHostMatrix[!is.na(FinalHostMatrix[,SubResps[r]]),]

# Generate species-level trait data
# Get Sp and Sp2 in "d" on the same factor levels

SubDataList[[r]]$Sp <- factor(as.character(SubDataList[[r]]$Sp),
                              levels = union(SubDataList[[r]]$Sp, SubDataList[[r]]$Sp2)
)

SubDataList[[r]]$Sp2 <- factor(as.character(SubDataList[[r]]$Sp2),
                               levels = union(SubDataList[[r]]$Sp, SubDataList[[r]]$Sp2)
)

SubDataList[[r]]$Sharing <- SubDataList[[r]][,SubResps[r]]

SubDataList[[r]]$space_s <- scale(SubDataList[[r]]$Space) %>% c
SubDataList[[r]]$phylo_s <- scale(SubDataList[[r]]$Phylo2) %>% c
SubDataList[[r]]$d_cites_s <- scale(log(SubDataList[[r]]$hDiseaseZACites+1)) %>% c
SubDataList[[r]]$d_cites_s2 <- scale(log(SubDataList[[r]]$hDiseaseZACites.Sp2+1)) %>% c

SubDataList[[r]]$domestic <- ifelse(SubDataList[[r]]$hDom=="domestic",1,0)
SubDataList[[r]]$domestic.Sp2 <- ifelse(SubDataList[[r]]$hDom.Sp2=="domestic",1,0)


SpaceRange <- seq(from = range(SubDataList[[r]]$space_s)[1],
                  to = range(SubDataList[[r]]$space_s)[2],
                  length = 100)

PhyloRange <- seq(from = range(SubDataList[[r]]$phylo_s)[1],
                  to = range(SubDataList[[r]]$phylo_s)[2],
                  length = 100)

GAMPredDF <- expand.grid(space_s = SpaceRange,
                         phylo_s = PhyloRange,
                         d_cites_s = 0,
                         domestic = 0)

GAMPredDF <- GAMPredDF %>% mutate(SpaceQ = cut(space_s, quantile(space_s, 0:10/10),include.lowest = T, labels = 1:10),
                                  PhyloQ = cut(phylo_s, quantile(phylo_s, 0:10/10),include.lowest = T, labels = 1:10))

df = predict(DNAGAM, newdata = GAMPredDF)
df2 = fitted(DNAGAM, newdata = GAMPredDF)

GAMPredDF[,c("Estimate","Lower","Upper")] <- df[,c("Estimate","Q2.5","Q97.5")]
GAMPredDF[,c("Estimate","Lower","Upper")] <- df2[,c("Estimate","Q2.5","Q97.5")]

ggplot(GAMPredDF, aes(space_s, phylo_s, fill = Pred)) + geom_tile()

ggplot(GAMPredDF, aes(space_s, Estimate)) + geom_point() + facet_wrap(~PhyloQ)
ggplot(GAMPredDF, aes(space_s, Estimate)) + geom_line(aes(group = as.factor(phylo_s)), alpha = 0.3)

ggplot(GAMPredDF, aes(space_s, Estimate)) + 
  geom_line(aes(colour = phylo_s, group = as.factor(phylo_s))) +
  labs(y = "Predicted Value") +
  facet_wrap(~PhyloQ)

ggplot(GAMPredDF, aes(phylo_s, Estimate)) + geom_point() + facet_wrap(~SpaceQ)
ggplot(GAMPredDF, aes(phylo_s, Estimate)) + geom_line(aes(group = as.factor(space_s)), alpha = 0.3)

ggplot(GAMPredDF, aes(phylo_s, Estimate)) + 
  geom_line(aes(colour = space_s, group = as.factor(space_s))) +
  facet_wrap(~SpaceQ)



