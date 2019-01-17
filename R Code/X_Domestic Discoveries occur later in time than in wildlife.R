
# Exploring temporal sequence of zoonosis emergence

Zoos <- AssocsDiscovery %>% filter(Wildlife.Virus == 1 & Human == 1)

ZooTemps <- reshape2::melt(with(Zoos, tapply(Year, list(Virus, Domestic.Host), min))) %>%
  dplyr::rename(Domestic = Var2)

ggplot(ZooTemps, aes(as.factor(Domestic), value)) + geom_violin() + geom_point(position = "jitter")

ggplot(Hosts, aes(log(DOMYearBP), c(hZoonosisCount))) + geom_point() + geom_smooth(method = lm)

WellsDoms <- c("Bos javanicus, Bos mutus, Bos taurus, Camelus bactrianus, Camelus ferus, Camelus dromedarius, Canis familiaris, Canis lupus, Capra aegagrus, Cavia porcellus, Equus africanus, Equus asinus, Equus caballus, Felis catus, Lama guanicoe, Mus musculus, Oryctolagus cuniculus, Ovis aries, Rattus norvegicus, Rattus rattus, Sus scrofa, Vicugna vicugna")
WellsDoms <- WellsDoms %>% str_split(", ") %>% unlist %>% str_replace(" ","_")

Hosts$WellsDoms <- ifelse(Hosts$Sp%in%WellsDoms, 1, 0)

Hosts %>% filter(Domestic + WellsDoms>0) %>% dplyr::select(Domestic, WellsDoms, Sp)
