rm(list = ls())

library(dplyr)
library(factorMerger)

data("ess")
# Time for FM (this takes quite a bit)
happyFM <- mergeFactors(ess$happy, ess$country,
                        family = "binomial",
                        method = "fast-fixed")

# GIC as AIC
p1 <- plot(happyFM, panel = "GIC",
           title = "",
           panelGrid = FALSE)

# GIC as BIC
p2 <- plot(happyFM, panel = "GIC",
           penalty = log(NROW(ess$country)),
           title = "",
           panelGrid = FALSE)

p3 <- plot(happyFM, panel = "GIC",
           penalty = 500,
           title = "",
           panelGrid = FALSE)

gicComparisons <-
    ggpubr::ggarrange(p1, p2, p3, ncol = 3, labels = c("A", "B", "C"))
ggsave(gicComparisons, file = "./materials/JCGS/examples/ess_gic.pdf",
       width = 17, height = 10)
save(gicComparisons, file = "./materials/JCGS/examples/ess_gic.rda")

p4 <- plot(happyFM, panel = "tree", penalty = 500,
     title = "      European Social Study - happiness proportion",
     panelGrid = FALSE, nodesSpacing = "effects")

essPanels <- ggpubr::ggarrange(gicComparisons, p4, ncol = 1, nrow = 2, labels = c("", "D    "))
library(ggplot2)
ggsave(essPanels, file = "./materials/JCGS/examples/essPanels.pdf",
       width = 15, height = 20)

save(essPanels, file = "./materials/JCGS/examples/essPanels.rda")
