library(factorMerger)

dfWithoutCovariates <- generateMultivariateSample(20,10)
dfWithCovariates <- dfWithoutCovariates
dfWithCovariates$covariates <- data.frame(runif(20), rnorm(20))  


context("Check mergeFactors() function")

test_that("Wrong input",{
  expect_error(mergeFactors())
  expect_error(mergeFactors(dfWithCovariates$response))
  expect_error(mergeFactors(dfWithCovariates$factor))
})

test_that("Output format",{
  expect_is(mergeFactors(dfWithCovariates$response, dfWithCovariates$factor), "factorMerger")
  expect_is(mergeFactors(dfWithCovariates$response, dfWithCovariates$factor, method = "fast-adaptive"), "factorMerger")
})


context("Check groupStats() function")

test_that("Wrong input",{
  expect_error(groupsStats())
  expect_error(groupsStats(dfWithoutCovariates))
})

test_that("Output format",{
  expect_is(groupsStats(mergeFactors(dfWithCovariates$response, dfWithCovariates$factor)), "data.frame")
})