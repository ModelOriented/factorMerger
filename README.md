# factorMerger

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/factorMerger)](https://cran.r-project.org/package=factorMerger)
[![Build Status](https://travis-ci.org/ModelOriented/factorMerger.svg?branch=master)](https://travis-ci.org/ModelOriented/factorMerger)
[![Pending Pull-Requests](http://githubbadges.herokuapp.com/ModelOriented/factorMerger/pulls.svg)](https://github.com/ModelOriented/factorMerger/pulls)
[![Github Issues](http://githubbadges.herokuapp.com/ModelOriented/factorMerger/issues.svg)](https://github.com/ModelOriented/factorMerger/issues)
[![DOI](https://zenodo.org/badge/70429809.svg)](https://zenodo.org/badge/latestdoi/70429809)

## Overview

`factorMerger` is a set of tools to support post-hoc testing that would enable to extract hierarchical structure of factors with respect to a given response.

## Installation

```{r}
# the easiest way to get factorMerger is to install it from CRAN:
install.packages("factorMerger")

# Or the the development version from GitHub:
# install.packages("devtools")
devtools::install_github("ModelOriented/factorMerger", build_vignettes = FALSE)
```

## More

The description of `factorMerger` may be found in the paper *The Merging Path Plot: adaptive fusing of k-groups with likelihood-based model selection* [(here)](https://arxiv.org/abs/1709.04412).


### Workflow

<img src="https://raw.githubusercontent.com/ModelOriented/factorMerger/master/README_workflow.png" alt="fm_workflow" width = '650'/>

### Cheatsheet

[Download cheatsheet](https://raw.githubusercontent.com/ModelOriented/factorMerger/master/materials/factorMerger-cheatsheet.pdf)

![](https://raw.githubusercontent.com/ModelOriented/factorMerger/master/materials/factorMerger-cheatsheet.png)


### Example

#### Survival analysis

```{r}
library(factorMerger)
library(forcats) # distinguish meaningful factors (fct_lump)

data(BRCA)
brcaSurv <- survival::Surv(time = BRCA$time, event = BRCA$vitalStatus)
drugName <- fct_lump(as.factor(BRCA$drugName), prop = 0.05) 

drugNameFM <- mergeFactors(response = brcaSurv[!is.na(drugName)], 
                           factor = drugName[!is.na(drugName)], 
                           family = "survival")

plot(drugNameFM, nodesSpacing = "effects", gicPanelColor = "grey2")

```
