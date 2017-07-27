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

ggsave(pl1, file = "FM_tukey.pdf", width = 10, height = 7)

pl2 <- factorMerger::plotMeansAndConfInt(fact, color = TRUE, clusterSplit = list("GIC", 2), palette = "Set2")
ggsave(pl2, file = "FM_panelCI.pdf", width = 3, height = 5)

pl3 <- factorMerger::plotBoxplot(fact, color = TRUE, clusterSplit = list("GIC", 2), palette = "Set2")
ggsave(pl3, file = "FM_panelBox.pdf", width = 3, height = 5)

pl4 <- factorMerger:::plotTukey(fact, palette = "Set2")
ggsave(pl4, file = "FM_panelTukey.pdf", width = 3, height = 5)




