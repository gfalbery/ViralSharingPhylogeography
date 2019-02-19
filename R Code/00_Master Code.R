# Master Source Code for Albersnet EcoHealth Alliance Analyses ####

library(dplyr)

setwd("/home/gfalbery/Albersnet")

rm(list = ls())
#file.remove(LargeFiles)

# Loading files that take a while to obtain ####

#LargeFiles <- paste0("data/",
#                     list.files("data") %>% subset(endsWith(list.files("data"),".Rdata")))#
#
#for(x in LargeFiles) load(x)#
#
#ModelFiles <- paste0("Model Files/",
#                     list.files("Model Files") %>% subset(endsWith(list.files("Model Files"),".Rdata")))

#for(x in ModelFiles) load(x)

# Running data setup scripts ####

CodeRoot <- "R Code/0_Data Import"

StartTime <- Sys.time()
source(paste0(CodeRoot,"/","0a_EHA Data Import.R"))
source(paste0(CodeRoot,"/","0b_Phylogenetic Data Import.R" ))
source(paste0(CodeRoot,"/","0c_Kludging Spatial Data Import.R"))
source(paste0(CodeRoot,"/","0d_Host Breadth and Distances.R"))
source(paste0(CodeRoot,"/","0e_Creating Final Host Dataset.R"))
source(paste0(CodeRoot,"/","0f_Creating Viral Subsets.R"))
EndTime <- Sys.time()

EndTime - StartTime

