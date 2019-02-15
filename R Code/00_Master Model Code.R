
# 00_Running Models and simulating ####

# Master Source Code for Albersnet EcoHealth Alliance Analyses ####

library(dplyr)

# Running data setup scripts ####

CodeRoot2 <- "R Code/1_Sharing Models"

StartTime <- Sys.time()
#source(paste0(CodeRoot2,"/","1a_STAN Model.R"))
source(paste0(CodeRoot2,"/","1b_STAN Simulating.R"))
source(paste0(CodeRoot2,"/","1c_Simulating on All Mammals.R"))

EndTime <- Sys.time()

EndTime - StartTime

