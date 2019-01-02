
library(igraph)

V(bipgraph)$color <- 
  rep(c(AlberColours[4], AlberColours[5]), c(dim(Viruses)[1],dim(Hosts)[1]))

V(bipgraph)$frame.color <- NA

E(bipgraph)$color <- adjustcolor("grey", alpha = 0.1)

bp = layout_as_bipartite(bipgraph)
plot(bipgraph, layout = bp)

plot(bipgraph, layout = bp, vertex.label = NA)
plot(bipgraph, layout = layout.circle, vertex.label = NA)

jpeg("bipgraph.jpeg", width = 200, height = 200, units = "mm", res = 1000)
plot(bipgraph, layout = layout.circle, vertex.label = NA)
dev.off()
