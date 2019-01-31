
# Trying to add in average host specificity of viruses 
# infecting each mammal species or species pair

HostAssocs <- apply(M, 2, function(a) names(a[a>0]))

FinalHostMatrix[,c("IntersectVirusRel","UnionVirusRel")] <- c(NA,NA)

for(x in 1:nrow(FinalHostMatrix)){
  
  if(x%%100==0) print(x)
  
  Sp1 <- FinalHostMatrix[x,"Sp"] %>% as.character
  Sp2 <- FinalHostMatrix[x,"Sp2"] %>% as.character
  
  SharedViruses1 <-  intersect(AssocsTraits[AssocsTraits$Host==Sp1,"Virus"] %>% as.character,
                               AssocsTraits[AssocsTraits$Host==Sp2,"Virus"] %>% as.character)
  
  OtherHosts1 <- intersect(setdiff(unlist(VirusAssocs[SharedViruses1]), c(Sp1,Sp2)), 
                           colnames(tSTMatrix))
  
  SubPhyloMatrix1 <- tSTMatrix[OtherHosts1, OtherHosts1]
  
  IntersectViralHostRel <- mean(SubPhyloMatrix1)
  
  SharedViruses2 <-  union(AssocsTraits[AssocsTraits$Host==Sp1,"Virus"] %>% as.character,
                           AssocsTraits[AssocsTraits$Host==Sp2,"Virus"] %>% as.character)
  
  OtherHosts2 <- intersect(setdiff(unlist(VirusAssocs[SharedViruses2]), c(Sp1,Sp2)), 
                           colnames(tSTMatrix))
  
  SubPhyloMatrix2 <- tSTMatrix[OtherHosts2, OtherHosts2]
  
  UnionViralHostRel <- mean(SubPhyloMatrix2)
  
  FinalHostMatrix[x,c("IntersectVirusRel","UnionVirusRel")] <- c(IntersectViralHostRel,UnionViralHostRel)
  
}

VirusRels1 <- FinalHostMatrix[,c("IntersectVirusRel","UnionVirusRel")]

# Let's try this using averaging across the associated virus's ranges 

for(x in 1:nrow(FinalHostMatrix)){
  
  if(x%%100==0) print(x)
  
  Sp1 <- FinalHostMatrix[x,"Sp"] %>% as.character
  Sp2 <- FinalHostMatrix[x,"Sp2"] %>% as.character
  
  SharedViruses1 <-  intersect(AssocsTraits[AssocsTraits$Host==Sp1,"Virus"] %>% as.character,
                               AssocsTraits[AssocsTraits$Host==Sp2,"Virus"] %>% as.character)
                          
  IntersectViralHostRel <- mean(VirPhyloHostRangeMean[SharedViruses1])
  
  SharedViruses2 <-  union(AssocsTraits[AssocsTraits$Host==Sp1,"Virus"] %>% as.character,
                           AssocsTraits[AssocsTraits$Host==Sp2,"Virus"] %>% as.character)
  
  UnionViralHostRel <- mean(VirPhyloHostRangeMean[SharedViruses2])
  
  FinalHostMatrix[x,c("IntersectVirusRel2","UnionVirusRel2")] <- c(IntersectViralHostRel,UnionViralHostRel)
  
}

VirusRels2 <- FinalHostMatrix[,c("IntersectVirusRel2","UnionVirusRel2")]

