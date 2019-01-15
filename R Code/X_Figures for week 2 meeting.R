
# Figures for meeting with KO, NR, EE

library(grid)

jpeg("Domestics do not harbour more zoonoses.jpeg", units = "mm", width = 200, height = 100, res = 300)

arrange_ggplot2(list(
  
  BarGraph(Hosts, "Domestic", "hZoonosisProp", text = "N") + scale_fill_brewer(palette = AlberPalettes[4]) +
    theme(legend.position = "none"),
  
  BarGraph(Hosts, "Domestic", "hZoonosisCount", text = "N") + scale_fill_brewer(palette = AlberPalettes[4]) +
    theme(legend.position = "none")
  
), ncol = 2)

dev.off()
