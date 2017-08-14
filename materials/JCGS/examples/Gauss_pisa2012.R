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

# Using factorMerger results in predicting (it is nonsense...) -----------------
load("./materials/JCGS/examples/oneDimPisa.rda")
aicPrediction <- cutTree(oneDimPisa,
                         stat = "GIC",
                         value = 2)

aicPrediction %>% table() %>% head()

getOptimalPartitionDf(oneDimPisa,
                      stat = "GIC",
                      value = 2) %>%
    head(5)

res <- replicate(1000, {
    nTest <- round(0.2 * NROW(pisa2012))
    testIds <- sample(1:NROW(pisa2012), size = nTest)
    test <- pisa2012[testIds, ]
    train <- pisa2012[-testIds, ]

    initialModel <- lm(math ~ country, data = train)
    initialPrediction <- predict(initialModel, test)

    newData <- data.frame(math = pisa2012$math,
                          country = aicPrediction)
    newTest <- newData[testIds, ]
    newTrain <- newData[-testIds, ]
    newModel <- lm(math ~ country, data = newTrain)
    newPrediction <- predict(newModel, newTest)

    return(c(
        Metrics::rmse(test$math, initialPrediction),
        Metrics::rmse(test$math, newPrediction)
    ))
})

res <- res %>% t() %>% as.data.frame()
res <- res[complete.cases(res), ]
mean(res[, 1] < res[, 2])

# ------------------------------------------------------------------------------

