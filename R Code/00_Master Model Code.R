
# 00_Running Models and simulating ####

# Master Output Code for Albersnet EcoHealth Alliance Analyses ####
# Rscript "R Code/00_Master Model Code.R"

library(dplyr)

# Running data setup scripts ####

CodeRoot2 <- "R Code/1_Sharing Models"

StartTime <- Sys.time()

source("R Code/00_Master Code.R")

#if(file.exists("BAMList.Rdata")) load("BAMList.Rdata") else source(paste0(CodeRoot2,"Frequentist GAMs.R"))

#source(paste0(CodeRoot2,"2_Simulating Known Network.R"))
#source(paste0(CodeRoot2,"2b_Known Network Characteristics.R"))

source(paste0(CodeRoot2,"3_Simulating Whole Network.R"))
source(paste0(CodeRoot2,"3b_All Network Characteristics.R"))
source(paste0(CodeRoot2,"3c_Spatial Degree Figure.R"))

#source(paste0(CodeRoot2,"4_Host Validation.R"))

EndTime <- Sys.time()

EndTime - StartTime

