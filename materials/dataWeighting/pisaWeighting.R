rm(list = ls())
load("./ignoreFiles/computerStudent2012.rda")
pisa2012 <- computerStudent2012[, c("PV1MATH","PV1READ","PV1SCIE","CNT")]
set.seed(123)
library(dplyr)
idx <- sample.int(NROW(pisa2012), prob = computerStudent2012 %>% pull(W_FSTUWT))
pisa2012 <- pisa2012[idx, ]
rownames(pisa2012) <- NULL
colnames(pisa2012) <- c("math", "reading", "science", "country")
save(pisa2012, file = "./data/pisa2012.rda")
