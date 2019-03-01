# 0b_Phylogenetic Data

if(!file.exists("data/intermediate/HP3-ST_PDmatrix.csv")){
  
  P <- rprojroot::find_rstudio_root_file
  library(PVR)
  library(ape)
  library(phylogram)
  
  asc <- AssocsBase
  h <- HostTraits
  v <- VirusTraits
  
  #check that host, virus, and asc files all match up with unique viruses and hosts included
  va <- unique(asc$Virus)
  ha <- unique(asc$Host)
  
  setdiff(ha, h$hHostNameFinal) #all spp in asc file are in h file,
  setdiff(va, v$vVirusNameCorrected) #all viruses in asc file are in v file,
  
  #load in phylogenetic trees, cytb and Supertree (ST)
  cytb <- read.tree(P("data/cytb_supertree.tree"))
  ST <- read.tree(P("data/supertree_mammals.tree")) #improved ST, fixed all taxa plus added monotremes 3 april 2014 (NEW ST w/outg)
  
  #drop tips from trees for hosts dropped from most recent host-virus association file
  #cytb$tip.label[sort.list(cytb$tip.label)] #665 spp in original cytb
  cytb_drop <- setdiff(cytb$tip.label, h$hHostNameFinal) #21 spp in cytb tree, not in h file
  cytb2 <- drop.tip(cytb, which(cytb$tip.label %in% cytb_drop)) #new tree, drop all tips in cytb not in h
  
  #ST$tip.label[sort.list(ST$tip.label)] #770 spp in original ST tree
  ST_drop <- setdiff(ST$tip.label, h$hHostNameFinal) #17 spp in ST tree, not in h file
  ST2 <- drop.tip(ST, which(ST$tip.label %in% ST_drop)) #new tree, drop all tips in ST not in h
  
  ## Calculate sp to sp maxtrix of phylo distance (cophenetic)
  vSTphylodist <- as.data.frame(cophenetic(ST2)) #calculate Supertree phylo distance matrix, 753 spp.
  write.csv(vSTphylodist, P("data/intermediate/HP3-ST_PDmatrix.csv")) #write ST PD matrix to file
  vCYTBphylodist <- as.data.frame(cophenetic(cytb2)) #calculate cytb phylo distance matrix, #644 spp.
  write.csv(vCYTBphylodist, P("data/intermediate/HP3-cytb_PDmatrix.csv")) #write cytb PD matrix to file
  
  STMatrix <- as.matrix(vSTphylodist)
  CytBMatrix <- as.matrix(vCYTBphylodist)
  
} else { 
  STMatrix <- read.csv("data/intermediate/HP3-ST_PDmatrix.csv", header = T) 
  CytBMatrix <- read.csv("data/intermediate/HP3-cytb_PDmatrix.csv", header = T) 
  
  rownames(STMatrix) <- STMatrix[,1]
  STMatrix <- STMatrix[,-1]
  STMatrix <- as.matrix(STMatrix)
  
  rownames(CytBMatrix) <- CytBMatrix[,1]
  CytBMatrix <- CytBMatrix[,-1]
  CytBMatrix <- as.matrix(CytBMatrix)
  
}

# Adding in all mammal supertree ####

if(!file.exists("data/intermediate/FullSTMatrix.csv")){
  
  library(geiger);library(ape);library(picante);library(dplyr)
  
  STFull <- read.nexus("data/ele_1307_sm_sa1.tre")[[1]]
  FullSTMatrix <- as.data.frame(cophenetic(STFull)) %>% as.matrix
  
} else{ FullSTMatrix <- as.matrix(read.csv("data/intermediate/FullSTMatrix.csv", header = T)) }
