library(factorMerger)

dfWithoutCovariates <- generateMultivariateSample(20,10)
merged <- mergeFactors(dfWithoutCovariates$response, dfWithoutCovariates$factor) 


context("Check plot() function")

test_that("Wrong input",{
  expect_error(plot())
})

test_that("Output format",{
  expect_is(plot(merged), "gg")
  expect_is(plot(merged), "ggplot")
  expect_is(plot(merged), "ggarrange")
})

context("Check plotHeatmap() function")

test_that("Wrong input",{
  expect_error(plotHeatmap())
})

test_that("Output format",{
  expect_is(plotHeatmap(merged, color=TRUE, clusterSplit=list("pvalue", "GIC")), "gg")
  expect_is(plotHeatmap(merged, color=TRUE, clusterSplit=list("pvalue", "GIC")), "ggplot")
})

# context("Check plotProfile() function")
# 
# test_that("Wrong input",{
#   expect_error(plotProfile())
# })
# 
# test_that("Output format",{
#   expect_is(plotProfile(merged, color=TRUE, clusterSplit=list("pvalue", "GIC")), "gg")
#   expect_is(plotProfile(merged, color=TRUE, clusterSplit=list("pvalue", "GIC")), "ggplot")
# })


# context("Check plotHMeansAndConfInt() function")
# 
# test_that("Wrong input",{
#   expect_error(plotMeansAndConfInt())
# })
# 
# test_that("Output format",{
#   expect_is(plotMeansAndConfInt(merged, color=TRUE, clusterSplit=list("pvalue", "GIC")), "gg")
#   expect_is(plotMeansAndConfInt(merged, color=TRUE, clusterSplit=list("pvalue", "GIC")), "ggplot")
# })

