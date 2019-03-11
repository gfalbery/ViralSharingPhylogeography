
# Simulating Viral Infection ####

NViruses <- 37000

NHosts <- sapply(VirusAssocs, length)



VirusTraits <- list(
  
  Breadth = rnbinom(NViruses, mu = mean(NHosts), size = 0.45)
  
)

a = 1

InfectedHosts <- list()

SimHosts <- AllSims[[1]]

CaseZero <- lapply(1:NViruses, function(b){
  
  print(b)
  
  PatientZero <- sample(rownames(SimHosts), 1)
  
  Bystanders <- sample(SimHosts[PatientZero,][SimHosts[PatientZero,]==1], VirusTraits$Breadth[1]) %>% names
  
  return(Bystanders)
  
})


InfectedHostTotals <- data.frame(
  VirusTotals = table(unlist(CaseZero)) %>% c,
  Sp = names(table(unlist(CaseZero)))
)

Panth1 %>% 
  left_join(InfectedHostTotals, by = "Sp") %>% 
  BarGraph("hOrder", "VirusTotals", order = T) +
  labs(title = "Simulated viral diversity", x = "Order", y = "Viral richness") +
  ggsave()

