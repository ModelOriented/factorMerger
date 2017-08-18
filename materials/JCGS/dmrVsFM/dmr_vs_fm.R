library(factorMerger)
set.seed(13)
fact <- rep(1:9, c(2, 10, 2, 2, 10,2, 2, 10, 2))
response <- rnorm(fact, fact/5)
df <- data.frame(fact, response)
df <- df[df$fact %in% c(1,2,4,5,7,8,9), ]

oneDim1 <- mergeFactors(df$response, df$fact, method = "fixed")
t1 <- plot(oneDim1, penalty = 1, panel = "tree", title = '    "fixed" fusing')

oneDim2 <- mergeFactors(df$response, df$fact,  method = "adaptive")
t2 <- plot(oneDim2, penalty = 1, panel = "tree", title = '    "adaptive" fusing')

ggpubr::ggarrange(t1, t2, ncol = 2)

mH1 <- mergingHistory(oneDim1, T)
mH2 <- mergingHistory(oneDim2, T)
mH1$method <- "fixed"
mH2$method <- "adaptive"
mH2$nGroups <- 7:1
mH1$nGroups <- 7:1
bothDims <- as.data.frame(rbind(mH1, mH2))
bothDims %>%
    filter(nGroups > 1) %>%
    ggplot(aes(x = nGroups, y = model, group = method, color = method)) +
    geom_line() + geom_point()

df$fact <- as.factor(oneDim1$factor)
df$fact <- factorMerger:::mergeLevels(df$fact, "(5)", "(8)", groupAB = "(5)(8)")
df$fact <- factorMerger:::mergeLevels(df$fact, "(7)", "(4)", groupAB = "(7)(4)")
df$fact <- factorMerger:::mergeLevels(df$fact, "(7)(4)", "(1)", groupAB = "(7)(4)(1)")

t3 <- plot(mergeFactors(df$response, df$fact, method = "fixed", abbreviate = F), panel = "tree",
           title = '    "fixed" fusing, continued')
t4 <- plot(mergeFactors(df$response, df$fact, method = "adaptive", abbreviate = F),
           panel = "tree",
           title = '    "adaptive" fusing, continued')


differentTrees <- ggpubr::ggarrange(t1, t2, t3, t4, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
ggsave(differentTrees, filename = "./materials/JCGS/dmr_vs_fm/differentTrees.pdf")
