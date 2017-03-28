---
output: 
  pdf_document:
    citation_package: natbib
    keep_tex: true
    fig_caption: true
    latex_engine: pdflatex
    template: ./svm-latex-ms.tex
title: "factorMerger: a set of tools to support results from post hoc testing"
author:
- name: Agnieszka Sitko
  affiliation: University of Warsaw
abstract: "`factorMerger` is an *R* package whose purpose is to extend methods of analysing dependencies between groups of a categorical variable after carrying out an analysis of variance (ANOVA). The idea of the package arose from the need to create an algorithm which outputs in a hierachical interpretation of relations between levels of a categorical variable. Thereby, for a given significance level groups may be devided into nonoverlapping clusters. `factorMerger` implements iterative version of post hoc testing based on likelihood ratio test for parametric models: gaussian, binomial and survival. It also provides custom visualizations for each model built on `ggplot2` package. 

\\
\\

*Package webpage*: https://github.com/geneticsMiNIng/FactorMerger"
keywords: "analysis of variance (ANOVA), hierarchical clustering, likelihood ratio test (LRT), post hoc testing"
date: "`r format(Sys.time(), '%B %d, %Y')`"
geometry: margin=1in
fontfamily: mathpazo
fontsize: 11pt
# spacing: double
bibliography: ./factorMerger.bib
biblio-style: apsr
---

# Introduction

If data is analysed using ANOVA a more detailed analysis of the response variable realization differences between categorical variable levels might be needed. The traditional approach will be to perform *post hocs* - multiple comparisons after ANOVA. For each pair of groups we run specific statistical test which output with an information whether response averages in those groups are significantly different. However, if we look from the distant perspective, for a certain significance level, these results may not be consistent. One may consider the case that mean in group A does not differ significantally from the one in group B, similarly with group B and C. In the same time difference between group A and C is detected. Then, data partition is unequivocal, which means impossible to put through. 

The problem of clustering categorical variable into non-overlapping groups has already been present in statistics. J. Tukey proposed an iterative procedure of merging factor levels based on studentized range distribution [@Tukey]. Therefore, this approach was limited to gaussian models. Delete or Merge Regressors algorithm (@Proch), extends partitioning for generalized linear models. *DMR4glm* calculates 



# Algorithm overview

# The *R* package `factorMerger`

# Bibliography

