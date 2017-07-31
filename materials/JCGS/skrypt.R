library(factorMerger)
library(ggplot2)

pisa2012eur <- pisa2012[pisa2012$CNT %in% c("Austria","Belgium","Germany","Denmark", "Spain","Finland","France","United Kingdom","Italy","Poland","Sweden"), ]
pisa2012eur$CNT <- factor(pisa2012eur$CNT)

fact <- mergeFactors(pisa2012eur$PV1MATH, pisa2012eur$CNT, method = "hclust")

pl1 <- plot(fact,
     panel = "all",
     responsePanel = "tukey",
     palette = "Set2",
     panelGrid = FALSE)


pl0 <- plot(fact,
            panel = "all",
            responsePanel = "tukey",
            palette = "Set2",
            panelGrid = FALSE,
            penalty = 50)
pl0
ggsave(pl0, file = "FM_0.pdf", width = 10, height = 7)

ggsave(pl1, file = "FM_tukey.pdf", width = 10, height = 7)

pl2 <- factorMerger::plotMeansAndConfInt(fact, color = TRUE, clusterSplit = list("GIC", 2), palette = "Set2")
ggsave(pl2, file = "FM_panelCI.pdf", width = 3, height = 5)

pl3 <- factorMerger::plotBoxplot(fact, color = TRUE, clusterSplit = list("GIC", 2), palette = "Set2")
ggsave(pl3, file = "FM_panelBox.pdf", width = 3, height = 5)

pl4 <- factorMerger:::plotTukey(fact, palette = "Set2")
ggsave(pl4, file = "FM_panelTukey.pdf", width = 3, height = 5)

fact2 <- mergeFactors(pisa2012eur$PV1MATH, pisa2012eur$CNT, method = "hclust")


pisa2012eur$highSkills <- ifelse(pisa2012eur$PV1MATH > 545, 1, 0)
pisa2012eur <- na.omit(pisa2012eur)

factB <- mergeFactors(pisa2012eur$highSkills, pisa2012eur$CNT, method = "hclust", family = "binomial")

pl5 <- factorMerger:::plotProportion(factB, color = TRUE, clusterSplit = list("GIC", 2), palette = "Set2")
ggsave(pl5, file = "FM_panelProportion.pdf", width = 3, height = 5)

factMD <- mergeFactors(pisa2012eur[,1:3], pisa2012eur$CNT, method = "hclust")

pl6 <- factorMerger:::plotHeatmap(factMD, color = TRUE, clusterSplit = list("GIC", 2), palette = "Blues")
ggsave(pl6, file = "FM_panelHeatmap.pdf", width = 3, height = 5)

pl7 <- factorMerger::plotProfile(factMD, color = FALSE, clusterSplit = list("GIC", 2), palette = "Blues")
ggsave(pl7, file = "FM_panelProfile.pdf", width = 3, height = 5)

pl8 <- factorMerger::plotFrequency(factMD, color = FALSE, clusterSplit = list("GIC", 2))
ggsave(pl8, file = "FM_panelFrequency.pdf", width = 3, height = 5)


library(survival)
library(survminer)
pisa2012eur <- na.omit(pisa2012eur)
pisa2012eur$one <- 1
factS <- mergeFactors(Surv(pisa2012eur$PV1MATH, pisa2012eur$one),
                      pisa2012eur$CNT, method = "hclust", family = "survival")

pl9 <- factorMerger::plotSurvival(factS, color = FALSE, clusterSplit = list("GIC", 2))
ggsave(pl9, file = "FM_panelSurvival.pdf", width = 3, height = 5)



