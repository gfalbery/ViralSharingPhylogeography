
# 00_Running Models and simulating ####

# Master Output Code for Albersnet EcoHealth Alliance Analyses ####
# Rscript "R Code/00_Master Model Code.R"

library(dplyr)

# Running data setup scripts ####

source("R Code/00_Master Code.R")

StartTime <- Sys.time()

CodeRoot2 <- "R Code/1_Sharing Models/"

if(file.exists("Output Files/BAMList.Rdata")) load("Output Files/BAMList.Rdata") else 
  source(paste0(CodeRoot2,"1a_Frequentist GAMs.R"))

#print("Known Network Simulation!")
source(paste0(CodeRoot2,"2_Simulating Known Network.R"))
source(paste0(CodeRoot2,"2b_Known Network Characteristics.R"))

print("All Mammal Simulation!")
source(paste0(CodeRoot2,"3_Simulating Whole Network.R"))
source(paste0(CodeRoot2,"3b_All Network Characteristics.R"))
source(paste0(CodeRoot2,"3c_Spatial Degree Figure.R"))

print("Validating!")
source(paste0(CodeRoot2,"4_Host Validation.R"))

EndTime <- Sys.time()

EndTime - StartTime

