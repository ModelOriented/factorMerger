# ---- How was it done?
library(factorMerger)
library(dplyr)

nGroups <- seq(5, 60, 5)
nObs <- seq(1000, 10000, 1000)

grid <- expand.grid(nGroups, nObs)
res <- data.frame(matrix(ncol = 4))
colnames(res) <- c("nGroups", "nObs", "method", "time")

for (i in 1:nrow(grid)) {
    print(i)
    grid[i, "Var2"] -> N
    grid[i, "Var1"] -> k
    sample <- generateSample(N = N, k = k)

    tmp <- replicate(5, system.time(mergeFactors(sample$response, sample$factor, method = "adaptive"))[3])
    new <- data.frame(nGroups = N, nObs = k, method = "adaptive", time = tmp)
    res <- rbind.data.frame(res, new)

    tmp <- replicate(5, system.time(mergeFactors(sample$response, sample$factor, method = "fast-adaptive"))[3])
    new <- data.frame(nGroups = N, nObs = k, method = "fast-adaptive", time = tmp)
    res <- rbind.data.frame(res, new)

    tmp <- replicate(5, system.time(mergeFactors(sample$response, sample$factor, method = "fixed"))[3])
    new <- data.frame(nGroups = N, nObs = k, method = "fixed", time = tmp)
    res <- rbind.data.frame(res, new)

    tmp <- replicate(5, system.time(mergeFactors(sample$response, sample$factor, method = "fast-fixed"))[3])
    new <- data.frame(nGroups = N, nObs = k, method = "fast-fixed", time = tmp)
    res <- rbind.data.frame(res, new)

}

colnames(res)[1:2] <- colnames(res)[2:1]


# ---- Plots
library(dplyr)
library(reshape2)
library(ggplot2)

load("./materials/JCGS/timeComparisons/timeComp.rda")


t1 <- res %>% dplyr::filter(nObs == 10000, nGroups > 10) %>%
    group_by(nGroups, method) %>% summarize(value = 1 / mean(time) * 60) %>%
    select(-nObs) %>%
    ggplot(aes(x = nGroups, y = value, color = method)) + geom_line() +
    geom_point(size = 0.7) +
    theme_bw() +
    ggtitle("Number of evaluations per 60 sec\n(fixed nObs = 10 000)")

t1 <- t1 +
    scale_color_discrete(guide = guide_legend(title = "method"),
                         labels = c('"adaptive"', '"fast-adaptive"', '"fixed"', '"fast-fixed"')) +
    xlab("Number of groups (k)") + ylab("Number of evaluations")

t2 <-  res %>% dplyr::filter(nGroups == 25) %>%
    group_by(nObs, method) %>% summarize(value = 1 / mean(time) * 60) %>%
    ggplot(aes(x = nObs, y = value, color = method)) + geom_line() +
    theme_bw() + geom_point(size = .7) +
    ggtitle("Number of evaluations per 60 sec\n(fixed nGroups = 25)")

t2 <- t2 +
    scale_color_discrete(guide = guide_legend(title = "method"),
                         labels = c('"adaptive"', '"fast-adaptive"', '"fixed"', '"fast-fixed"')) +
    xlab("Number of observations (n)") + ylab("Number of evaluations")
t2 <- t2 + ggtitle("Number of evaluations per 60 sec\n(fixed number of groups = 25)")
t1 <- t1 + ggtitle("Number of evaluations per 60 sec\n(fixed number of observations = 10 000)")

nEvalPerMinute <- ggpubr::ggarrange(t1, t2, ncol = 2, legend = "bottom", common.legend = TRUE, labels = c("A", "B"))
ggsave(nEvalPerMinute, filename = "./materials/JCGS/timeComparisons/nEvalsPerMinute.pdf",
       width = 8, height = 5)

save(nEvalPerMinute, file = "./materials/JCGS/timeComparisons/nEvalsPerMinute.rda")

p1 <- res %>%
    group_by(nGroups, method, nObs) %>%
    summarise(time = mean(time)) %>%
    ggplot(aes(x = nGroups,
               y = time,
               group = as.factor(nObs),
               col = as.factor(nObs))) + geom_line() + geom_point(size = 0.7) +
    theme_bw() + facet_wrap(~method, scales = "free") +
    ggtitle("Evaluation times for different number of groups (OX axis) and sample sizes (colors)") +
    theme(legend.position = "none") + xlab("Number of groups") + ylab("Evaluation time (sec)")

ggsave(p1, filename = "./materials/JCGS/timeComparisons/nGroupsTrends.pdf")

