# SJD: R package for Structured Joint Decomposition

## Overview:

Structured Joint Decomposiiton (SJD) is a R package for visualizing biologically structured gene expression matrix environment based on low rank models. Currently, it provides user 4 different styles of decomposition: (1) Separately, (2) Concatenately, (3) Jointly, (4) Sequentially.

To install this package in R, run the following commands:

```R
library(devtools)
install_github("CHuanSite/SJD")
```

This package implements four categories of algorithms to decompose multiple datasets, (1) Separately, (2) Concatenately, (3) Jointly, (4) Two Stage Sequentially. 

For the first three categories, there are three available algorithms: Principal Component Analysis (PCA), Independent Component Analysis (ICA) and Nonnegative Matrix Factorization (NMF). For each method, the algorithm takes three arguments, `dataset`, `group` and `comp_num`, specifying which datasets to be used, what is the structure among the datasets and what's the dimension for each component.

For the last category, there is one algorithm, called two-staged linked component analysis, which is a PCA based statistical model.

## Example usage:

```R
library(SJD)
# Simulation the dataset
dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
               matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
               matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
               matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
               
## Specify the structure among the datasets
group = list(c(1,2,3,4), c(1,2), c(3,4), c(1,3), c(2,4), c(1), c(2), c(3), c(4))
comp_num = c(2,2,2,2,2,2,2,2,2)

## Separate PCA, ICA, NMF
sepPCA_res = sepPCA(dataset, comp_num)
sepICA_res = sepICA(dataset, comp_num)
sepNMF_res = sepNMF(dataset, comp_num)

## Concatenated PCA, ICA, NMF
concatPCA_res = concatPCA(dataset, group, comp_num)
concatICA_res = concatICA(dataset, group, comp_num)
concatNMF_res = concatNMF(dataset, group, comp_num)

## Joint PCA, ICA, NMF
jointPCA_res = jointPCA(dataset, group, comp_num)
jointICA_res = jointICA(dataset, group, comp_num)
jointNMF_res = jointNMF(dataset, group, comp_num)

## twoStageLCA
twoStageLCA_res = twoStageLCA(dataset, group, comp_num)
```

To access the component
```R
concatPCA_res$linked_component_list
```

To access the score
```R
concatPCA_res$score_list
```

## Weighting data sets by different weights

Sometimes it is interesting to weigh data sets differently, to incorporate user's different views on it. The SJD package also has this option.

```R
# Simulation the dataset
dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
               matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
               matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
               matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
               
## Specify the structure among the datasets
group = list(c(1,2,3,4), c(1,2), c(3,4), c(1,3), c(2,4), c(1), c(2), c(3), c(4))
comp_num = c(2,2,2,2,2,2,2,2,2)
weighting = c(2, 1, 4, 3)

## Separate PCA, ICA, NMF
sepPCA_res = sepPCA(dataset, comp_num, weighting)
sepICA_res = sepICA(dataset, comp_num, weighting)
sepNMF_res = sepNMF(dataset, comp_num, weighting)

## Concatenated PCA, ICA, NMF
concatPCA_res = concatPCA(dataset, group, comp_num, weighting)
concatICA_res = concatICA(dataset, group, comp_num, weighting)
concatNMF_res = concatNMF(dataset, group, comp_num, weighting)

## Joint PCA, ICA, NMF
jointPCA_res = jointPCA(dataset, group, comp_num, weighting)
jointICA_res = jointICA(dataset, group, comp_num, weighting)
jointNMF_res = jointNMF(dataset, group, comp_num, weighting)

## twoStageLCA
twoStageLCA_res = twoStageLCA(dataset, group, comp_num, weighting)


```

## Projecting new data sets to extracted components

```R
# Simulation the dataset
dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
               matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
               matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
               matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
               
## Specify the structure among the datasets
group = list(c(1,2,3,4), c(1,2), c(3,4), c(1,3), c(2,4), c(1), c(2), c(3), c(4))
comp_num = c(2,2,2,2,2,2,2,2,2)

## Projected new data sets
proj_dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
proj_group = list(c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE))

## concatenate PCA with projection functionality on
res_concatPCA = concatPCA(dataset, group, comp_num, weighting = NULL, proj_dataset = proj_dataset, proj_group = proj_group)
```

## Single-Cell RNAseq Example

First, install the 'googleDrive' package

```R
install.package('googledrive')
```

Download RNA and explaination data from Google Drive,

```R
library(googledrive)

## Data file 
url_data <- "https://drive.google.com/file/d/1OQovDBPwRX_O2N1GSNY8fzJn-p3-fwQV/view?usp=sharing"
drive_download(url_data, overwrite = TRUE)
unzip('data.zip')

## Explaination file
url_explaination <- 'https://drive.google.com/file/d/1S3HdygRCMvPttmVd9cix4GskWj1VPJaM/view?usp=sharing'
drive_download(url_explaination, overwrite = TRUE)
unzip('data_explaination.zip')
```

Read data into R

```R
library(tidyverse)

## Read in files
inVitro_bulk = read.table('1_inVitro_Bulk_Cortecon.plog2_trimNord.txt', stringsAsFactors = FALSE, header = TRUE) %>% select(-1) %>% as.matrix
inVitro_sc = read.table('2_inVitro_SingleCell_scESCdifBifurc.CelSeq_trimNord.txt', stringsAsFactors = FALSE, header = TRUE) %>% select(-1) %>% as.matrix
inVivo_bulk = read.table('3_inVivo_Bulk_BrainSpan_RNAseq_Gene_DFC_noSVA_plog2_trimNord.txt', stringsAsFactors = FALSE, header = TRUE) %>% select(-1) %>% as.matrix
inVivo_sc = read.table('4_inVivo_SingleCell_CtxDevoSC4kTopoTypoTempo_plog2_trimNord.txt', stringsAsFactors = FALSE, header = TRUE) %>% select(-1) %>% as.matrix

inVitro_bulk_scale = scale(t(scale(t(inVitro_bulk), scale = FALSE)))
inVitro_sc_scale = scale(t(scale(t(inVitro_sc), scale = FALSE)))
inVivo_bulk_scale = scale(t(scale(t(inVivo_bulk), scale = FALSE)))
inVivo_sc_scale= scale(t(scale(t(inVivo_sc), scale = FALSE)))

## legends for the 4 datasets
inVitro_bulk_exp =  read.table("1_inVitro_Bulk_Cortecon.pd.txt",stringsAsFactors = FALSE, header = T)
inVitro_sc_exp = read.table("2_inVitro_SingleCell_scESCdifBifurc.CelSeq.pd.txt", stringsAsFactors = FALSE, header = T)
inVivo_bulk_exp = read.table("3_inVivo_Bulk_BrainSpan.RNAseq.Gene.DFC.pd.txt", stringsAsFactors = FALSE, header = T)
inVivo_sc_exp = read.table("4_inVivo_SingleCell_CtxDevoSC4kTopoTypoTempo.pd.txt", stringsAsFactors = FALSE, header = T)
````

Conduct Two-stage linked component analysis

```R
library(SJD)
## List of datasets and group assignment and number of components
dataset = list(inVitro_bulk_scale, inVitro_sc_scale, inVivo_bulk_scale, inVivo_sc_scale)
group = list(c(1,2,3,4), c(1,2), c(3,4), c(1,3), c(2,4), c(1), c(2), c(3), c(4))
comp_num = c(2,2,2,2,2,2,2,2,2)

## Output result
twoStageLCA_res = twoStageLCA(dataset, group, comp_num)
```
Visualize the result

```R
par(mfrow = c(2,2))

## common component
plot(t(twoStageLCA_res$score_list[[1]][[1]]), col = inVitro_bulk_exp$color, pch = 16, xlab = "PC1", ylab = "PC2", main = "common: inVitro_bulk", cex = 2, cex.axis = 1, cex.lab = 1, cex.main = 1)
plot(t(twoStageLCA_res$score_list[[2]][[1]]), col = inVitro_sc_exp$COLORby.DCX, pch = 16, xlab = "PC1", ylab = "PC2", main = "common: inVitro_sc", cex = 2, cex.axis = 1, cex.lab = 1, cex.main = 1)
plot(t(twoStageLCA_res$score_list[[3]][[1]]), col = inVivo_bulk_exp$color, pch = 16, xlab = "PC1", ylab = "PC2", main = "common: inVivo_bulk", cex = 2, cex.axis = 1, cex.lab = 1, cex.main = 1)
plot(t(twoStageLCA_res$score_list[[4]][[1]]), col = inVivo_sc_exp$COLORby.DCX, pch = 16, xlab = "PC1", ylab = "PC2", main = "common: inVivo_sc", cex = 2, cex.axis = 1, cex.lab = 1, cex.main = 1)


par(mfrow = c(1,2))
## in vitro component
plot(t(twoStageLCA_res$score_list[[1]][[2]]), col = inVitro_bulk_exp$color, pch = 16, xlab = "PC1", ylab = "PC2", main = "common: inVitro_bulk", cex = 2, cex.axis = 1, cex.lab = 1, cex.main = 1)
plot(t(twoStageLCA_res$score_list[[2]][[2]]), col = inVitro_sc_exp$COLORby.DCX, pch = 16, xlab = "PC1", ylab = "PC2", main = "common: inVitro_sc", cex = 2, cex.axis = 1, cex.lab = 1, cex.main = 1)

## in vivo component
plot(t(twoStageLCA_res$score_list[[3]][[3]]), col = inVivo_bulk_exp$color, pch = 16, xlab = "PC1", ylab = "PC2", main = "common: inVivo_bulk", cex = 2, cex.axis = 1, cex.lab = 1, cex.main = 1)
plot(t(twoStageLCA_res$score_list[[4]][[3]]), col = inVivo_sc_exp$COLORby.DCX, pch = 16, xlab = "PC1", ylab = "PC2", main = "common: inVivo_sc", cex = 2, cex.axis = 1, cex.lab = 1, cex.main = 1)

## bulk component
plot(t(twoStageLCA_res$score_list[[1]][[4]]), col = inVitro_bulk_exp$color, pch = 16, xlab = "PC1", ylab = "PC2", main = "common: inVitro_bulk", cex = 2, cex.axis = 1, cex.lab = 1, cex.main = 1)
plot(t(twoStageLCA_res$score_list[[3]][[4]]), col = inVivo_bulk_exp$color, pch = 16, xlab = "PC1", ylab = "PC2", main = "common: inVivo_bulk", cex = 2, cex.axis = 1, cex.lab = 1, cex.main = 1)

## sc component
plot(t(twoStageLCA_res$score_list[[2]][[5]]), col = inVitro_sc_exp$COLORby.DCX, pch = 16, xlab = "PC1", ylab = "PC2", main = "common: inVitro_sc", cex = 2, cex.axis = 1, cex.lab = 1, cex.main = 1)
plot(t(twoStageLCA_res$score_list[[4]][[5]]), col = inVivo_sc_exp$COLORby.DCX, pch = 16, xlab = "PC1", ylab = "PC2", main = "common: inVivo_sc", cex = 2, cex.axis = 1, cex.lab = 1, cex.main = 1)
```

## Source

The latest developer version is available in the [Github repository](https://github.com/CHuanSite/SJD).

## Contact

For any improvements and issues that need to be addressed, please contact [Huan Chen](hchen130@jhu.edu).


