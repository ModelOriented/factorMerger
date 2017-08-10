library(factorMerger)
library(dplyr)

nGroups <- seq(25, 60, 5)
nObs <- seq(10000, 10000, 1000)

grid <- expand.grid(nGroups, nObs)
grid <- grid %>% filter(Var2 == 10000 | Var1 == 25)
res <- data.frame(matrix(nrow = nrow(grid), ncol = 4))

for (i in 1:nrow(grid)) {
    print(i)
    sample <- generateSample(N = grid[i, "Var2"], k = grid[i, "Var1"])
    res[i, 1] <- system.time(mergeFactors(sample$response, sample$factor, method = "adaptive"))[3]
    res[i, 2] <- system.time(mergeFactors(sample$response, sample$factor))[3]
    res[i, 3] <- system.time(mergeFactors(sample$response, sample$factor, method = "fixed"))[3]
    res[i, 4] <- system.time(mergeFactors(sample$response, sample$factor, method = "fast-fixed"))[3]
}

res2 <- res

colnames(res2) <- c("LRT", "LRTSucc", "hClust", "hClustSucc")
res2$nObs <- grid$Var2
res2$nGroups <- grid$Var1
res2 <- res2[complete.cases(res2), ]
res3 <- res2

save(res2, file = "./ignoreFiles/timeComp2.rda")
save(res, file = "./ignoreFiles/timeComp.rda")
library(dplyr)
library(reshape2)
library(ggplot2)

load("./materials/JCGS/timeComparisons/timeComp.rda")
load("./materials/JCGS/timeComparisons/timeComp2.rda")

res <- rbind(res, res2)
res <- as.data.frame(res)

t1 <- res %>% dplyr::filter(nObs == 10000, nGroups > 10) %>%
    select(-nObs) %>% melt(id.vars = "nGroups") %>%
    ggplot(aes(x = nGroups, y = 1 / value * 60, color = variable)) + geom_line() +
    theme_bw() +
    ggtitle("Number of evaluations per 60 sec\n(fixed nObs = 10 000)")

t1 <- t1 +
    scale_color_discrete(guide = guide_legend(title = "method"),
                         labels = c('"adaptive"', '"fast-adaptive"', '"fixed"', '"fast-fixed"')) +
    xlab("Number of groups") + ylab("Number of evaluations")

t2 <- res %>% dplyr::filter(nGroups == 25) %>%
    select(-nGroups) %>% melt(id.vars = "nObs") %>%
    ggplot(aes(x = nObs, y = 1 / value * 60, color = variable)) + geom_line() +
    theme_bw() +
    ggtitle("Number of evaluations per 60 sec\n(fixed nGroups = 25)")

t2 <- t2 +
    scale_color_discrete(guide = guide_legend(title = "method"),
                         labels = c('"adaptive"', '"fast-adaptive"', '"fixed"', '"fast-fixed"')) +
    xlab("Number of observations") + ylab("Number of evaluations")
t2 <- t2 + ggtitle("Number of evaluations per 60 sec\n(fixed number of groups = 25)")
t1 <- t1 + ggtitle("Number of evaluations per 60 sec\n(fixed number of observations = 10 000)")

nEvalPerMinute <- ggpubr::ggarrange(t1, t2, ncol = 2, common.legend = TRUE, labels = c("1", "2"))
ggsave(nEvalPerMinute, filename = "./materials/JCGS/timeComparisons/nEvalsPerMinute.pdf",
       width = 8, height = 5)

colnames(res)[1:4] <- c('"adaptive"', '"fast-adaptive"', '"fixed"', '"fast-fixed"')

p1 <- res %>% dplyr::filter(nGroups < 25) %>% melt(id.vars = c("nGroups", "nObs")) %>%
    ggplot(aes(x = nGroups,
               y = value,
               group = as.factor(nObs),
               col = as.factor(nObs))) + geom_line() +
    theme_bw() + facet_wrap(~variable, scales = "free") +
    ggtitle("Evaluation times for different number of groups (OX axis) and sample sizes (colors)") +
    theme(legend.position = "none") + xlab("Number of groups") + ylab("Evaluation time (sec)")

ggsave(p1, filename = "./materials/JCGS/timeComparisons/nGroupsTrends.pdf")

p2 <- res %>% dplyr::filter(nGroups < 25) %>% melt(id.vars = c("nGroups", "nObs")) %>%
    ggplot(aes(x = nObs,
               y = value,
               group = as.factor(nGroups),
               col = as.factor(nGroups))) + geom_line() +
    theme_bw() + facet_wrap(~variable, scales = "free") +
    ggtitle("Evaluation times for different sample sizes (OX axis) and number of groups (colors)") +
    theme(legend.position = "none") + xlab("Number of observations") + ylab("Evaluation time (sec)")
ggsave(p2, filename = "./materials/JCGS/timeComparisons/nObsTrends.pdf")
