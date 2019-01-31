
# X_Matrices ####

Hostadj

HostMatrixdf <- data.frame(Host1 = rep(rownames(Hostadj), each = nrow(Hostadj)),
                           Host2 = rep(rownames(Hostadj), nrow(Hostadj)),
                           Assocs = as.vector(Hostadj)
)

VirusMatrixdf <- data.frame(Virus1 = rep(rownames(Virusadj), each = nrow(Virusadj)),
                           Virus2 = rep(rownames(Virusadj), nrow(Virusadj)),
                           Assocs = as.vector(Virusadj)
)

