set.seed(123)
sample <- generateSample(100000, 5)
fm <- mergeFactors(sample$response, sample$factor)
p1 <- plot(fm, panel = "tree", title = "")
p2 <- plot(fm, panel = "tree", panelGrid = FALSE, title = "")
p3 <- ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1, labels = c("A", "B"))
ggplot2::ggsave(p3,
       filename = "./materials/JCGS/panelGrid.pdf")
