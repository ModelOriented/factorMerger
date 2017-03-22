# FACTOR MERGER

The aim of this project is to create hierarchical algorithm for post-hoc testing. Detailed description (currently in Polish) and first results may be found here:

https://github.com/geneticsMiNIng/FactorMerger/blob/master/vignettes/factorMerger.Rmd

A short overview in English is available here

https://rawgit.com/geneticsMiNIng/FactorMerger/573948f073a71ababe03dec0a21488118d1f0a48/materials/vignette2.html

## Generic functions

* `calculateModel [gaussianFactorMerger, survivalFactorMerge, binomialFactorMerger]`
* `appendProjection [factorMerger, gaussianFactorMerger]`
* `getStatisticName [gaussianFactorMerger, binomialFactorMerger, survivalFactorMerger]`
* `compareModels [lm, coxph, binomglm]`
* `getPvals [lm, mlm]`
* `calculateModelStatistic [coxph]` # logLik
* `logLik [mlm]`
