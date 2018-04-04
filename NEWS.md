factorMerger 0.3.6
----------------------------------------------------------------
* **Misc**:
    * create a lighter version of Pisa dataset
    * make some plotting functions public


factorMerger 0.3.5
----------------------------------------------------------------
* **Misc**:
    * Weighted regression with the `weights` argument in `mergeFactors` [#63](https://github.com/MI2DataLab/factorMerger/issues/63)
    * Added support for formulas in in `mergeFactors` [#76](https://github.com/MI2DataLab/factorMerger/issues/76)

factorMerger 0.3.4
----------------------------------------------------------------
* **Misc**:
    * Changed intro message (now it's more compact)
    * Added tests [#69](https://github.com/MI2DataLab/factorMerger/issues/69)

* **Fixes**:
    * Vigniette building (#67)
    * Fix in Tukey's plots

factorMerger 0.3.3
----------------------------------------------------------------
* **New features**:
    * In this version arguments: 'method' and 'successive' of the mergeFactors() function are merged in one argument ('method') with its new values:
        * 'fast-fixed' (method = 'hclust', successive = TRUE),
        * 'fixed' (method = 'hclust', successive = FALSE),
        * 'fast-adaptive' (method = 'LRT', successive = TRUE),
        * 'adaptive' (method = 'LRT', successive = FALSE).
* **Fixes**:
    * Updated cheatsheet (#56)
    * Added new dataset pisaEuro (#65)

factorMerger 0.3.2
----------------------------------------------------------------
	
TODO
