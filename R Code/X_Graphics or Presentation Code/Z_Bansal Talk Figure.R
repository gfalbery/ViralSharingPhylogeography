
# Figure for Bansal lab meeting ####

Reps <- lapply(CentralityList[["Degree"]][c(1,3,4,7)], INLARep)

df = data.frame(Model = rep(c(1,3,4,7), sapply(Reps, length)),
                Var = c("Resid", "Resid","Space", "Resid", "Phylo", "Resid", "Space", "Phylo"),
                Value = unlist(Reps))

df$Var <- factor(df$Var, levels=c("Resid","Space","Phylo"))
df$Model <- factor(df$Model)

ggplot(df, aes(Model, Value, fill = Var)) + geom_col(position = "stack") + 
  scale_fill_manual(values = c("grey", "#2C7FB8","#DE2D26")) + labs(y = "Variance accounted for", fill = "Component") +
  scale_x_discrete(labels = c(1:4)) +
  ggsave("Figures/Variance Components.jpeg", units = "mm", height = 100, width = 100, dpi = 300)
