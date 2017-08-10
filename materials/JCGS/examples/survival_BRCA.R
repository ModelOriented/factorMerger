rm(list = ls())
library(dplyr)
library(factorMerger)
library(forcats)
library(survival)

data("BRCA")

BRCA <- BRCA %>%
    filter(!is.na(drugName))

drugName <- fct_lump(BRCA$drugName, prop = 0.05)
brcaSurv <- Surv(time = BRCA$time,
                 event = BRCA$vitalStatus)

drugMerge <- mergeFactors(response = brcaSurv,
                          factor = drugName,
                          family = "survival",
                          method = "adaptive")
print(drugMerge)

survDefault <- plot(drugMerge)
ggsave(survDefault, filename = "./materials/JCGS/examples/survival_default.pdf",
       width = 15, height = 10)

plot(drugMerge,
     nodesSpacing = "effects",
     colorClusters = FALSE,
     showSplit = TRUE,
     title = "BRCA: patient survival vs. drug treatment")

plot(drugMerge,
     nodesSpacing = "effects",
     title = "BRCA: patient survival vs. drug treatment",
     panel = "response",
     palette = "Dark2")



p1 <- plot(drugMerge,
           colorClusters = FALSE,
           showSplit = TRUE,
           nodesSpacing = "effects",
           title = "BRCA: patient survival vs. drug treatment",
           panel = "response",
           palette = "Dark2")

library(ggplot2)
ggsave(p1, file = "./materials/JCGS/examples/survival_no_cluster_colored.pdf",
       width = 11, height = 7)

p2 <- plot(drugMerge,
           nodesSpacing = "effects",
           title = "BRCA: patient survival vs. drug treatment",
           panel = "response",
           palette = "Dark2")

ggsave(p2, file = "./materials/JCGS/examples/survival_clusters_colored.pdf",
       width = 11, height = 7)

ggpubr::ggarrange(p1, p2, ncol = 2)
