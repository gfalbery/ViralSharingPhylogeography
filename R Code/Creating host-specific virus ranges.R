
# Trying to add in average host specificity of viruses 
# infecting each mammal species or species pair

HostAssocs <- apply(M, 2, function(a) names(a[a>0]))

FinalHostMatrix[,c("IntersectVirusRel","UnionVirusRel")] <- c(NA,NA)

for(x in nrow(FinalHostMatrix)){
  
  if(x%%100==0) print(x)
  
  Sp1 <- FinalHostMatrix[x,"Sp"] %>% as.character
  Sp2 <- FinalHostMatrix[x,"Sp2"] %>% as.character
  
  SharedViruses1 <-  intersect(AssocsTraits[AssocsTraits$Host==Sp1,"Virus"] %>% as.character,
                              AssocsTraits[AssocsTraits$Host==Sp2,"Virus"] %>% as.character)
  
  OtherHosts1 <- unique(unlist(VirusAssocs[SharedViruses1]))
  OtherHosts1 <- setdiff(OtherHosts1, c(Sp1,Sp2))
  
  SubPhyloMatrix1 <- tSTMatrix[OtherHosts1, OtherHosts1]
  
  IntersectViralHostRel <- mean(SubPhyloMatrix1)
  
  SharedViruses2 <-  union(AssocsTraits[AssocsTraits$Host==Sp1,"Virus"] %>% as.character,
                          AssocsTraits[AssocsTraits$Host==Sp2,"Virus"] %>% as.character)
  
  OtherHosts2 <- unique(unlist(VirusAssocs[SharedViruses2]))
  OtherHosts2 <- setdiff(OtherHosts2, c(Sp1,Sp2))
  
  SubPhyloMatrix2 <- tSTMatrix[OtherHosts2, OtherHosts2]
  
  UnionViralHostRel <- mean(SubPhyloMatrix2)
  
  FinalHostMatrix[x,c("IntersectVirusRel","UnionVirusRel")] <- c(IntersectViralHostRel,UnionViralHostRel)
  
}

# Let's try this using averaging across the associated virus's ranges 

