rm(list = ls())
library(dplyr)
library(factorMerger)

data("pisa2012")

oneDimPisa <- mergeFactors(response = pisa2012$math,
                           factor = pisa2012$country,
                           method = "fast-fixed")

save(oneDimPisa, file = "./materials/JCGS/examples/oneDimPisa.rda")

# Using factorMerger results in predicting (it is nonsense...) -----------------
load("./materials/JCGS/examples/oneDimPisa.rda")
bicPrediction <- cutTree(oneDimPisa, value = 10)

nTest <- round(0.2 * NROW(pisa2012))
testIds <- sample(1:NROW(pisa2012), size = nTest)
test <- pisa2012[testIds, ]
train <- pisa2012[-testIds, ]

initialModel <- lm(math ~ country, data = train)
initialPrediction <- predict(initialModel, test)

newData <- data.frame(math = pisa2012$math,
                      country = bicPrediction)
newTest <- newData[testIds, ]
newTrain <- newData[-testIds, ]
newModel <- lm(math ~ country, data = newTrain)
newPrediction <- predict(newModel, newTest)

Metrics::rmse(test$math, initialPrediction)
Metrics::rmse(test$math, newPrediction)

# ------------------------------------------------------------------------------

