---
title: "factorMerger: a set of tools to support results from post hoc testing"
author: |
   | Agnieszka Sitko^1^ and Przemysław Biecek^1^
   |
   | 1. Faculty of Mathematics, Informatics and Mechanics, University of Warsaw
output: html_document
---

**Keywords**: analysis of variance (*ANOVA*), hierarchical clustering, likelihood ratio test (LRT), post hoc testing

**Webpages**: https://github.com/geneticsMiNIng/FactorMerger

If data is analysed using *ANOVA* it often happens that a more detailed analysis of differences among categorical variable levels might be needed. The traditional approach is to perform *pairwise post hoc tests* - multiple comparisons after *ANOVA*.

However, for a fixed significance level commonly used *single-step post-hoc tests* may result in an unequivocal partition of the data. To deal with this problem, the **factorMerger** package implements stepwise algorithms of merging factor levels. The sequential approach enables to identify factor's hierarchical structure and find non-overlapping groups of its levels with similar distribution of the response. The current implementation of the package handles one-dimensional and multi-dimensional Gaussian models as well as binomial and survival models. 

In this talk we will describe main functionalities of the **factorMerger** package, briefly summarizing the theoretical background, and illustrate its use with the *Merging Path Plots* visualizations created with the **ggplot2** package.  
