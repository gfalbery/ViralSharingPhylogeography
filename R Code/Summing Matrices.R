
# Adding together all matrices ####

# Rscript "R Code/Summing Matrices.R"

library(tidyverse); library(Matrix)

load("Output Files/AllPredList.Rdata")

AllPredDF <- AllPredList %>% as.data.frame()

AllPredSums <- apply(AllPredDF,1,sum)

AssMat <- matrix(NA, 
                 nrow = 4276, #length(union(AllMammaldf$Sp,AllMammaldf$Sp2)), 
                 ncol = 4276) #length(union(AllMammaldf$Sp,AllMammaldf$Sp2)))

AssMat[lower.tri(AssMat)] <- AllPredSums
AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
diag(AssMat) <- 0

#dimnames(AssMat) <- list(union(AllMammaldf$Sp,AllMammaldf$Sp2),
#                         union(AllMammaldf$Sp,AllMammaldf$Sp2))

AllSums <- as(AssMat, "dgCMatrix")

save(AllSums, file = "Output Files/AllSums.Rdata")



