rm(list = ls())
library(haven)
library(dplyr)
library(factorMerger)
essData <- haven::read_sav("./ignoreFiles/ESS7e02_1.spss/ESS7e02_1.sav")
essData <- essData %>% as.data.frame()
nLevels <- essData %>% apply(2, function(x) length(unique(x)))
binary <- nLevels == 2
essBinary <- essData[, binary]
happy <- essData[, "happy"]
country <- essData[, "cntry"]
which <- happy %in% 0:10
happy <- (happy %>% as.data.frame())[which, ]
country <- (country %>% as.data.frame())[which, ]
weight <- (essData[, "dweight"] %>% as.data.frame())[which, ]

set.seed(123)
ids <- sample(1:NROW(country), size = 5 * NROW(country), replace = T, prob = weight)
happy <- happy[ids, ]
country <- country[ids, ]

happy <- happy[, 1] > 5
country <- country[, 1] %>% as.factor()

ess <- data.frame(happy = happy,
                  country = country)

save(ess, file = "./data/ess.rda")

library(factorMerger)
data("ess")
# Time for FM
happyFM <- mergeFactors(happy, country, family = "binomial", successive = T)

# GIC as BIC
plot(happyFM, penalty = log(NROW(country)),
     title = "European Social Survey 2014 - HOW HAPPY ARE YOU?",
     panelGrid = FALSE)


# Time for DMR -------------------------------------
happyFM <- mergeFactors(happy, country, family = "binomial",
                        successive = T, method = "hclust")

# GIC as BIC
plot(happyFM, penalty = log(NROW(country)),
     title = "European Social Survey 2014 - HOW HAPPY ARE YOU?",
     panelGrid = FALSE)
