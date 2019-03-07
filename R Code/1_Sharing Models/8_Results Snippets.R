
# 8_Results snippets

AssocsTraits %>% filter(Host%in%FHN) %>% dim

(AssocsTraits %>% filter(Host%in%FHN) %>% select(Virus))$Virus %>% nunique

length(FHN)

#nrow(FinalHostMatrix)

DevOutput$VirusBinary["Phylo"]
DevOutput$VirusBinary["Space"]

mean(DataList[[1]]$Space, na.rm = T) %>% round(3)*100

DevOutput$VirusBinary["MinCites"]
DevOutput$VirusBinary["Domestic"]

length(AllMammals)
length(intersect(FHN,AllMammals))

table(Panth1$Obs, Panth1$EIDObs)

Panth1 %>% 
  group_by(hOrder) %>% 
  summarise(mean(AllPredDegree)) %>% 
  filter(hOrder == "Rodentia")

Panth1 %>% 
  group_by(hOrder) %>% 
  summarise(mean(AllPredDegree)) %>% 
  filter(hOrder == "Chiroptera")

Panth1 %>% 
  group_by(hOrder) %>% 
  summarise(mean(AllPredDegree)) %>% 
  filter(hOrder == "Carnivora")

Panth1 %>% 
  group_by(hOrder) %>% 
  summarise(mean(AllPredDegree)) %>% 
  filter(hOrder == "Artiodactyla")

Panth1 %>% 
  group_by(hOrder) %>% 
  summarise(mean(OutDegree)) %>% 
  filter(hOrder == "Artiodactyla")

Panth1 %>% 
  group_by(hOrder) %>% 
  summarise(mean(OutDegree)) %>% 
  filter(hOrder == "Carnivora")

Panth1 %>% 
  group_by(hOrder) %>% 
  summarise(mean(OutDegree)) %>% 
  filter(hOrder == "Pholidota")


nrow(ValidSummary)
median(ValidSummary$MeanRank)
length(AllMammals)
ValidSummary %>% filter(MeanRank == 1) %>% nrow
ValidSummary %>% filter(MeanRank <= 10) %>% nrow
ValidSummary %>% 
  filter(MeanRank > length(AllMammals)/2) %>% 
  nrow                        

nlevels(DataList$RNA$Sp)
nlevels(DataList$Vector$Sp)
nlevels(DataList$NVector$Sp)
nlevels(DataList$DNA$Sp)


