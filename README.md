# FACTOR MERGER

The aim of this project is to create hierarchical algorithm for post-hoc testing. Detailed description (currently in Polish) and first results may be found here:

https://github.com/geneticsMiNIng/FactorMerger/blob/master/vignettes/factorMerger.Rmd

## Generic functions

* `calculateModel [gaussianFactorMerger, survivalFactorMerge, binomialFactorMerger]`
* `appendProjection [factorMerger, gaussianFactorMerger]`
* `getStatisticName [gaussianFactorMerger, binomialFactorMerger, survivalFactorMerger]`
* `compareModels [lm, coxph, binomglm]`
* `getPvals [lm, mlm]`
* `calculateModelStatistic [coxph]` # logLik
* `logLik [mlm]`
