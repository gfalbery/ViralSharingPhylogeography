# Master Source Code for Albersnet EcoHealth Alliance Analyses ####

library(dplyr)

setwd("/home/gfalbery/Albersnet")

rm(list = ls())

# Running data setup scripts ####

CodeRoot <- "R Code/0_Data Import"

StartTime <- Sys.time()
print("Import")
source(paste0(CodeRoot,"/","0a_EHA Data Import.R"))
print("Phylo")
source(paste0(CodeRoot,"/","0b_Phylogenetic Data Import.R" ))
print("Space")
source(paste0(CodeRoot,"/","0c_Kludging Spatial Data Import.R"))
print("Ecology")
source(paste0(CodeRoot,"/","0d_Ecological Data Import.R"))
print("Final Dataset")
source(paste0(CodeRoot,"/","0e_Creating Final Host Dataset.R"))
#print("Subsets")
source(paste0(CodeRoot,"/","0g_Creating Viral Subsets.R"))
EndTime <- Sys.time()

EndTime - StartTime

save(FinalHostMatrix, file = "Output Files/Finaldf.Rdata")
