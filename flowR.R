library(BiocManager)

library(flowCore)
library(ggcyto)

setwd("D:/data")

path <- file.path(getwd(), "experiments", "4_2_0", "NP190530 - CL4 + JK4-2", "JK4-2_1-0 T L_022.fcs")

x <- read.FCS(path, transformation = "scale")

summary(x)

P_vector <- vector()

for (i in 1:20) {
  P_vector[i] <- paste("$P", i, "E", sep = "")
}

keyword(x, P_vector)

comp <- spillover(x)[[1]]

x_comp <- compensate(x, comp)

autoplot(x_comp,"FSC-H", "SSC-H")

autoplot(x_comp, "SSC-A")

shape_gate <- norm2Filter("FSC-H", "SSC-H", filterId = "shape", scale.factor = 2)

SC_gate <- rectangleGate("FSC-A" = c(0.125, 0.5), "SSC-A" = c(0, 0.25))

sc <- Subset(x_comp, SC_gate)

autoplot(sc, "FSC-A", "Zombie")

zombie_gate <- rectangleGate("[355] 450/50-A" = c(0, 0.025), "FSC-A" = c(0.125, 0.5))

z <- Subset(sc, zombie_gate)

autoplot(z, "CD45", "CD3")

CD3_gate <- rectangleGate("[355] 379/28-A" = c(-0.2, 0.125))

cd3 <- Subset(z, CD3_gate)

autoplot(cd3, "CD45", "CD3")
