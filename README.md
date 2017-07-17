# factorMerger: A Set of Tools to Support Results From Post Hoc Testing

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/factorMerger)](https://cran.r-project.org/package=factorMerger)
[![Build Status](https://api.travis-ci.org/MI2DataLab/factorMerger.png)](https://travis-ci.org/geneticsMiNIng/factorMerger)
[![Pending Pull-Requests](http://githubbadges.herokuapp.com/MI2DataLab/factorMerger/pulls.svg)](https://github.com/geneticsMiNIng/factorMerger/pulls)
[![Github Issues](http://githubbadges.herokuapp.com/MI2DataLab/factorMerger/issues.svg)](https://github.com/MI2DataLab/factorMerger/issues)
[![DOI](https://zenodo.org/badge/70429809.svg)](https://zenodo.org/badge/latestdoi/70429809)

The aim of this project is to create an algorithm of post-hoc testing that would enable to extract hierarchical structure of factors.

A short overview is available here

https://rawgit.com/MI2DataLab/FactorMerger/master/materials/vignette.html

Longer description may be found here

https://github.com/MI2DataLab/FactorMerger/blob/master/materials/useR/user_abstract.pdf


### Installing and loading `factorMerger`

`factorMerger` can be installed from [CRAN](https://cran.r-project.org/package=factorMerger) as follows:

```{r}
install.packages("factorMerger")
```

To install and load the latest version of`factorMerger` from **Github** run:

```{r}
if (!require(devtools)) {
    install.packages("devtools")
}
devtools::install_github("geneticsMiNIng/factorMerger", build_vignettes = FALSE)

library(factorMerger)
```

### Working with `factorMerger`

<img src="https://github.com/MI2DataLab/factorMerger/blob/master/README_workflow.png" alt="fm_workflow" width = '650'/>
