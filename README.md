# SJD: R package for Structured Joint Decomposition

<!-- badges: start -->
[![Travis build status](https://travis-ci.com/CHuanSite/SJD.svg?branch=main)](https://travis-ci.com/CHuanSite/SJD)
[![R-CMD-check](https://github.com/CHuanSite/SJD/workflows/R-CMD-check/badge.svg)](https://github.com/CHuanSite/SJD/actions)
<!-- badges: end -->

## Overview:

Structured Joint Decomposiiton (SJD) is a R package for visualizing biologically structured gene expression matrix environment based on low rank models. Currently, it provides user 4 different styles of decomposition: (1) separately, (2) concatenately, (3) jointly, (4) sequentially.

To install this package in R, run the following commands:

```R
library(devtools)
install_github("CHuanSite/SJD")
```

This package implements four categories of algorithms to decompose multiple datasets, (1) Separately, (2) Concatenately, (3) Jointly, (4) Two Stage Sequentially. 

For the first three categories, there are three available algorithms: Principal Component Analysis (PCA), Independent Component Analysis (ICA) and Nonnegative Matrix Factorization (NMF). For each method, the algorithm takes three arguments, `dataset`, `group` and `comp_num`, specifying which datasets to be used, what is the structure among the datasets and what's the dimension for each component.

For the last category, there is one algorithm, called two-staged linked component analysis, which is a PCA based statistical model.

## Example Use

Please visit the vignette page [here](https://chuansite.github.io/SJD/)

[Structral Joint Decomposition (SJD)](https://chuansite.github.io/SJD/articles/StructralJointDecomposition.html)

[Two Stage Linked Component Analysis (2s-LCA)](https://chuansite.github.io/SJD/articles/twoStageLCA.html)

## Source

The latest developer version is available in the [Github repository](https://github.com/CHuanSite/SJD).

## Contact

For any improvements and issues that need to be addressed, please contact [Huan Chen](mailto:hchen130@jhu.edu).


