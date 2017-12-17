library(factorMerger)

dfWithoutCovariates <- generateMultivariateSample(20,10)
merged <- mergeFactors(dfWithoutCovariates$response, dfWithoutCovariates$factor) 


context("Check cutTree() function")

test_that("Wrong input",{
  expect_error(cutTree())
})

test_that("Output format",{
  expect_is(cutTree(merged), "factor")
})



context("Check getOptimalPartition() function")

test_that("Wrong input",{
  expect_error(getOptimalPartition())
})

test_that("Output format",{
  expect_is(getOptimalPartition(merged), "character")
})



context("Check getOptimalPartitionDf() function")

test_that("Wrong input",{
  expect_error(getOptimalPartitionDf())
})

test_that("Output format",{
  expect_is(getOptimalPartitionDf(merged), "data.frame")
})


