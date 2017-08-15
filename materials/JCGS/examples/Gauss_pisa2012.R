rm(list = ls())
library(dplyr)
library(factorMerger)

data("pisa2012")

pisa2012 <- pisa2012[complete.cases(pisa2012), ]
save(pisa2012, file = "./data/pisa2012.rda")
oneDimPisa <- mergeFactors(response = pisa2012$math,
                           factor = pisa2012$country,
                           method = "fast-fixed")

save(oneDimPisa, file = "./materials/JCGS/examples/oneDimPisa.rda")

mergingHistory(oneDimPisa, TRUE) %>% head(5)

load("./materials/JCGS/examples/oneDimPisa.rda")
aicPrediction <- cutTree(oneDimPisa,
                         stat = "GIC",
                         value = 2)

aicPrediction %>% table() %>% head()

getOptimalPartitionDf(oneDimPisa,
                      stat = "GIC",
                      value = 2) %>%
    head(5)

