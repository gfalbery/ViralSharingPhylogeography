
# Generating Figures ####

# 1.	Panel: Relationships/model outputs from: 
# Space and phylogeny; phylogeny and viral sharing; space and viral sharing for all data



# 2.	Predicted degree centrality with and without random effect for ~
# 700 known viral hosts (possibly supplemental)

# 3.	Mammal order level centrality (bar plots)


# 4.	Gridcell level centrality measure (map), to see if overlaps w species diversity (optional)



# 5.	Panel: viral trait plots (RNA vs. DNA; vector-borne or not; etc)


# 6.	Panel of maps: Illustrative maps of host ranges for mammal species from host predictions (currently unidentified) for specific viruses, e.g. CCHFV; Nipah virus; Ebola; etc. (supplement will list out probability of being a host for all viruses)

PredHostPlot("Crimean-congo_haemorrhagic_fever_virus", focal = c(1,0), facet = T)

list(PredHostPlot("Andes_virus", focal = 0), PredHostPlot("Andes_virus", focal = 1)) %>% 
  arrange_ggplot2




# Supplementary? #####

ggplot(ValidSummary, aes(log10(NHosts), log10(MeanRank))) + geom_smooth() + geom_text(aes(label = Virus))
