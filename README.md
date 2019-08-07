# Albersnet #
## Authors: Greg Albery, Evan Eskew, Noam ross, and Kevin Olival ##

Work for EcoHealth Alliance for a three-month internship, intercalated within Greg's PhD. Internship began fully 07 Jan 2019, and ended March 22nd.

The project examines mammalian viral sharing patterns and their phylogeographic correlates. Using mammal species pairs as the unit of analysis, we constructed Generalised Additive Mixed Models (GAMMs) to untangle the roles of spatial overlap and phylogenetic relatedness when accounting for species-level sampling biases in the network.

Simulating with the resulting estimates produced a neutral network of global viral sharing patterns, which we used to uncover multiple taxonomic and geographic patterns of viral sharing.

## Preprint here: ##

### R scripts are numbered according to the order of use: ###
- scripts beginning with `0_` will import and structure datasets.
- all other scripts will: 
-- run viral sharing GAMMs, 
-- validate the GAMMs, 
-- use the estimates to predict the global mammalian viral sharing network, 
-- summarise traits of the network,
-- use it to predict sharing patterns, and 
-- make figures.

### Data come from three main sources: ###
- The IUCN ( )
- A mammalian supertree ( )
- A mammal-virus association dataset ( )
