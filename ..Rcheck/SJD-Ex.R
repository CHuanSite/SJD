pkgname <- "SJD"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('SJD')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("BEMA")
### * BEMA

flush(stderr()); flush(stdout())

### Name: BEMA
### Title: BEMA for the standard spiked covariance model
### Aliases: BEMA
### Keywords: Spike numbers

### ** Examples

x = matrix(rnorm(1000, 100), nrow = 1000)
eigen_x = svd(x)
eigen_out = list(eigenvalue = eigen_x$d^2 / 100, p = 1000, n = 100)
BEMA(eigen_out$eigenvalue, p = 1000, n = 100)




cleanEx()
nameEx("NeuroGenesis4")
### * NeuroGenesis4

flush(stderr()); flush(stdout())

### Name: NeuroGenesis4
### Title: NeuroGenesis4
### Aliases: NeuroGenesis4
### Keywords: datasets

### ** Examples


data(NeuroGenesis4)
head(NeuroGenesis4$Meissner.inVitro.bulk.Hs)
head(NeuroGenesis4$LIBD.AZ.inVitro.bulk.Hs)
head(NeuroGenesis4$Geschwind.inVivo.sc.Hs)
head(NeuroGenesis4$Jabaudon.inVivo.sc.Mm)




cleanEx()
nameEx("NeuroGenesis4.afterWrap")
### * NeuroGenesis4.afterWrap

flush(stderr()); flush(stdout())

### Name: NeuroGenesis4.afterWrap
### Title: NeuroGenesis4.afterWrap
### Aliases: NeuroGenesis4.afterWrap
### Keywords: datasets

### ** Examples


data(NeuroGenesis4.afterWrap)
head(NeuroGenesis4.afterWrap$Meissner.inVitro.bulk.Hs)
head(NeuroGenesis4.afterWrap$LIBD.AZ.inVitro.bulk.Hs)
head(NeuroGenesis4.afterWrap$Geschwind.inVivo.sc.Hs)
head(NeuroGenesis4.afterWrap$Jabaudon.inVivo.sc.Mm)




cleanEx()
nameEx("NeuroGenesis4.info")
### * NeuroGenesis4.info

flush(stderr()); flush(stdout())

### Name: NeuroGenesis4.info
### Title: NeuroGenesis4.info
### Aliases: NeuroGenesis4.info
### Keywords: datasets

### ** Examples


data(NeuroGenesis4.info)
head(NeuroGenesis4.info$Meissner.inVitro.bulk.Hs)
head(NeuroGenesis4.info$LIBD.AZ.inVitro.bulk.Hs)
head(NeuroGenesis4.info$Geschwind.inVivo.sc.Hs)
head(NeuroGenesis4.info$Jabaudon.inVivo.sc.Mm)




cleanEx()
nameEx("SJDScorePlotter")
### * SJDScorePlotter

flush(stderr()); flush(stdout())

### Name: SJDScorePlotter
### Title: Plot SJD score
### Aliases: SJDScorePlotter
### Keywords: images

### ** Examples


library(ggplot2)

data(NeuroGenesis4.afterWrap)
data(NeuroGenesis4.info)

SampleMetaNamesTable = data.frame(
   row.names = names(NeuroGenesis4),
   Type = c('Yaxis','Yaxis','2Dscatter','2Dscatter'),
   XaxisColumn = c("X","DAYx","tSNE_1","tsne1:ch1"),
   YaxisColumn = c("PJDscores","PJDscores","tSNE_2","tsne2:ch1"),
   COLaxisColumn = c("color","colorBYlabelsX","PJDscores","PJDscores"),
   PCHColumn = c("","","",""),
   inset = c(TRUE, TRUE, TRUE, TRUE),
   insetLOC = c("topright", "topright", "topright", "topright"),
   insetZoom = c(0.3, 0.3, 0.3, 0.3),
   ordDECREASE=c(FALSE, FALSE, FALSE, FALSE),
   CLRfoldPRB=c(0.5, 0.5, 0.5, 0.5)
)

grp = list(
Shared.All.4 = c(1 : 4),
Shared.bulk.2 = c(1, 2),
Shared.sc.2 = c(3, 4),
Hs.Meisnr.1 = c(1),
Hs.AZ.1 = c(2),
Gesch.1 = c(3),
Telley.1 = c(4)
)

dims = c(2, 2, 2, 2, 2, 2, 2)

twoStageLCA.out = twoStageLCA(dataset = NeuroGenesis4.afterWrap, group = grp, comp_num = dims)

SJDScorePlotter.obj = SJDScorePlotter(
    SJDalg = "twoStageLCA",
    scores = twoStageLCA.out$score_list,
    lbb = "NeuroGenesis4.p2",
    info = NeuroGenesis4.info,
    SampleMetaNamesTable = SampleMetaNamesTable
)




cleanEx()
nameEx("assemble.byComponent")
### * assemble.byComponent

flush(stderr()); flush(stdout())

### Name: assemble.byComponent
### Title: Assemble Files Based on Component
### Aliases: assemble.byComponent
### Keywords: component, images

### ** Examples


library(ggplot2)

data(NeuroGenesis4.afterWrap)
data(NeuroGenesis4.info)

SampleMetaNamesTable = data.frame(
   row.names = names(NeuroGenesis4),
   Type = c('Yaxis','Yaxis','2Dscatter','2Dscatter'),
   XaxisColumn = c("X","DAYx","tSNE_1","tsne1:ch1"),
   YaxisColumn = c("PJDscores","PJDscores","tSNE_2","tsne2:ch1"),
   COLaxisColumn = c("color","colorBYlabelsX","PJDscores","PJDscores"),
   PCHColumn = c("","","","")
)

grp = list(
Shared.All.4 = c(1 : 4),
Shared.bulk.2 = c(1, 2),
Shared.sc.2 = c(3, 4),
Hs.Meisnr.1 = c(1),
Hs.AZ.1 = c(2),
Gesch.1 = c(3),
Telley.1 = c(4)
)
dims = c(2, 2, 2, 2, 2, 2, 2)

lbb = "NeuroGenesis4.p2"

twoStageLCA.out = twoStageLCA(dataset = NeuroGenesis4.afterWrap, group = grp, comp_num = dims)

SJDScorePlotter.obj = SJDScorePlotter(
    SJDalg = "twoStageLCA",
    scores = twoStageLCA.out$score_list,
    lbb = lbb,
    info = NeuroGenesis4.info,
    SampleMetaNamesTable = SampleMetaNamesTable
)

assemble.byComponent.obj = assemble.byComponent(
SJDScorePlotter.obj = SJDScorePlotter.obj,
component = c(1, 2),
SJD_algorithm = "twoStageLCA",
group = 'Shared.All.4')




cleanEx()
nameEx("assemble.byDataset")
### * assemble.byDataset

flush(stderr()); flush(stdout())

### Name: assemble.byDataset
### Title: Assemble files based on Dataset
### Aliases: assemble.byDataset
### Keywords: dataset, images

### ** Examples


library(ggplot2)

data(NeuroGenesis4.afterWrap)
data(NeuroGenesis4.info)

SampleMetaNamesTable = data.frame(
   row.names = names(NeuroGenesis4.afterWrap),
   Type = c('Yaxis','Yaxis','2Dscatter','2Dscatter'),
   XaxisColumn = c("X","DAYx","tSNE_1","tsne1:ch1"),
   YaxisColumn = c("PJDscores","PJDscores","tSNE_2","tsne2:ch1"),
   COLaxisColumn = c("color","colorBYlabelsX","PJDscores","PJDscores"),
   PCHColumn = c("","","","")
)

grp = list(
Shared.All.4 = c(1 : 4),
Shared.bulk.2 = c(1, 2),
Shared.sc.2 = c(3, 4),
Hs.Meisnr.1 = c(1),
Hs.AZ.1 = c(2),
Gesch.1 = c(3),
Telley.1 = c(4)
)
dims = c(2, 2, 2, 2, 2, 2, 2)

lbb = "NeuroGenesis4.p2"

twoStageLCA.out = twoStageLCA(dataset = NeuroGenesis4.afterWrap, group = grp, comp_num = dims)

SJDScorePlotter.obj = SJDScorePlotter(
    SJDalg = "twoStageLCA",
    scores = twoStageLCA.out$score_list,
    lbb = lbb,
    info = NeuroGenesis4.info,
    SampleMetaNamesTable = SampleMetaNamesTable
)

assemble.byDataset.obj = assemble.byDataset(
SJDScorePlotter.obj = SJDScorePlotter.obj,
dataset_name = "Meissner.inVitro.bulk.Hs",
SJD_algorithm = "twoStageLCA",
group = NA)




cleanEx()
nameEx("balanceData")
### * balanceData

flush(stderr()); flush(stdout())

### Name: balanceData
### Title: Balance Data
### Aliases: balanceData
### Keywords: balance, dataset

### ** Examples

dataset = list(matrix(c(1 : 8), nrow = 2), matrix(1 : 6, nrow = 2))
balanceData(dataset)




cleanEx()
nameEx("compNameAssign")
### * compNameAssign

flush(stderr()); flush(stdout())

### Name: compNameAssign
### Title: Components Name Assignment for concatenate, joint and
###   twoStageLCA analysis
### Aliases: compNameAssign
### Keywords: component, name

### ** Examples

linked_component_list = list(matrix(c(1:4), nrow = 2), matrix(c(1:4), nrow = 2))
group_name = c("x", "y")
compNameAssign(linked_component_list, group_name)




cleanEx()
nameEx("compNameAssignSep")
### * compNameAssignSep

flush(stderr()); flush(stdout())

### Name: compNameAssignSep
### Title: Components Name Assignment for Seperate Analysis
### Aliases: compNameAssignSep
### Keywords: component, name

### ** Examples

linked_component_list = list(matrix(c(1:4), nrow = 2), matrix(c(1:4), nrow = 2))
dataset_name = c("x", "y")
compNameAssignSep(linked_component_list, dataset_name)




cleanEx()
nameEx("concatICA")
### * concatICA

flush(stderr()); flush(stdout())

### Name: concatICA
### Title: Concatenated decomposition with Independent Component Analysis
### Aliases: concatICA
### Keywords: ICA pairwise,

### ** Examples

dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
group = list(c(1,2,3,4), c(1,2), c(3,4), c(1,3), c(2,4), c(1), c(2), c(3), c(4))
comp_num = c(2,2,2,2,2,2,2,2,2)
proj_dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
proj_group = list(c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE))
res_concatICA = concatICA(
dataset,
group,
comp_num,
proj_dataset = proj_dataset,
proj_group = proj_group)




cleanEx()
nameEx("concatNMF")
### * concatNMF

flush(stderr()); flush(stdout())

### Name: concatNMF
### Title: Concatenated decomposition with Nonnegative Matrix Factorization
### Aliases: concatNMF
### Keywords: NMF pairwise,

### ** Examples

dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
group = list(c(1,2,3,4), c(1,2), c(3,4), c(1,3), c(2,4), c(1), c(2), c(3), c(4))
comp_num = c(2,2,2,2,2,2,2,2,2)
proj_dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
proj_group = list(c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE))
res_concatNMF = concatNMF(
dataset,
group,
comp_num,
proj_dataset = proj_dataset,
proj_group = proj_group)





cleanEx()
nameEx("concatPCA")
### * concatPCA

flush(stderr()); flush(stdout())

### Name: concatPCA
### Title: Concatenated Decomposition with Principal Component Analysis
### Aliases: concatPCA
### Keywords: PCA pairwise,

### ** Examples

dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
group = list(c(1,2,3,4), c(1,2), c(3,4), c(1,3), c(2,4), c(1), c(2), c(3), c(4))
comp_num = c(2,2,2,2,2,2,2,2,2)
proj_dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
proj_group = list(c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE))
res_concatPCA = concatPCA(
dataset,
group,
comp_num,
weighting = NULL,
proj_dataset = proj_dataset,
proj_group = proj_group)






cleanEx()
nameEx("configuration_setting_generation")
### * configuration_setting_generation

flush(stderr()); flush(stdout())

### Name: configuration_setting_generation
### Title: Configuration for simulated data generation
### Aliases: configuration_setting_generation
### Keywords: configure data for generation of setting the

### ** Examples

configuration_setting_generation()




cleanEx()
nameEx("datasetNameExtractor")
### * datasetNameExtractor

flush(stderr()); flush(stdout())

### Name: datasetNameExtractor
### Title: Dataset Name Extractor
### Aliases: datasetNameExtractor
### Keywords: name

### ** Examples

dataset = list(
x = matrix(c(1 : 4), nrow = 2),
y = matrix(c(1 : 4), nrow = 2))

datasetNameExtractor(dataset)




cleanEx()
nameEx("filterNAValue")
### * filterNAValue

flush(stderr()); flush(stdout())

### Name: filterNAValue
### Title: Filter NA Value
### Aliases: filterNAValue
### Keywords: NA, filter

### ** Examples


x = list(matrix(c(1 : 4), nrow = 2), matrix(c(1 : 4), nrow = 2))
y = list(matrix(c(1 : 4), nrow = 2), matrix(c(1 : 4), nrow = 2))
list_score = list(x, y)
dataset = c(1, 2)
group = list(c(1), c(2))
filterNAValue(list_score, dataset, group)




cleanEx()
nameEx("frameToMatrix")
### * frameToMatrix

flush(stderr()); flush(stdout())

### Name: frameToMatrix
### Title: Transform data.frame to matrix
### Aliases: frameToMatrix
### Keywords: dataframe, matrix

### ** Examples

dataset = list(
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50)
)
frameToMatrix(dataset)




cleanEx()
nameEx("geneNameAssign")
### * geneNameAssign

flush(stderr()); flush(stdout())

### Name: geneNameAssign
### Title: Gene name assignment
### Aliases: geneNameAssign
### Keywords: gene, name

### ** Examples

x = matrix(c(1 : 4), nrow = 2)
y = matrix(c(1 : 4), nrow = 2)
component_list = list(x, y)
gene_name = c("gene.1", "gene.2")
geneNameAssign(component_list, gene_name)




cleanEx()
nameEx("geneNameExtractor")
### * geneNameExtractor

flush(stderr()); flush(stdout())

### Name: geneNameExtractor
### Title: Gene Name Extractor
### Aliases: geneNameExtractor
### Keywords: gene, name

### ** Examples

x = matrix(c(1 : 4), nrow = 2)
rownames(x) = c("row1", "row2")
y = matrix(c(1 : 4), nrow = 2)
rownames(y) = c("row1", "row2")
dataset = list(x, y)
geneNameExtractor(dataset)




cleanEx()
nameEx("geneScreen")
### * geneScreen

flush(stderr()); flush(stdout())

### Name: geneScreen
### Title: Gene Screen
### Aliases: geneScreen
### Keywords: gene screen,

### ** Examples

dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))

screen_prob = c(0.2, 0.2, 0.2, 0.2)

screened_dataset = geneScreen(dataset, screen_prob)




cleanEx()
nameEx("getMatch")
### * getMatch

flush(stderr()); flush(stdout())

### Name: getMatch
### Title: Species matching function
### Aliases: getMatch
### Keywords: genes shared

### ** Examples


data(NeuroGenesis4)
out = getMatch(
rownames(NeuroGenesis4$Meissner.inVitro.bulk.Hs),
inSpecies = 'human',
inType = 'symbol',
newSpecies = 'mouse')




cleanEx()
nameEx("groupNameExtractor")
### * groupNameExtractor

flush(stderr()); flush(stdout())

### Name: groupNameExtractor
### Title: Group Name Extractor
### Aliases: groupNameExtractor
### Keywords: name

### ** Examples

dataset = list(
x = c(1, 2, 3),
y = c(1, 2, 4))




cleanEx()
nameEx("jointICA")
### * jointICA

flush(stderr()); flush(stdout())

### Name: jointICA
### Title: Joint decomposition with Independent Component Analysis
### Aliases: jointICA
### Keywords: ICA Joint,

### ** Examples

dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
group = list(c(1,2,3,4), c(1,2), c(3,4), c(1,3), c(2,4), c(1), c(2), c(3), c(4))
comp_num = c(2,2,2,2,2,2,2,2,2)
proj_dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
proj_group = list(c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE))
res_jointICA = jointICA(
dataset,
group,
comp_num,
proj_dataset = proj_dataset,
proj_group = proj_group)





cleanEx()
nameEx("jointNMF")
### * jointNMF

flush(stderr()); flush(stdout())

### Name: jointNMF
### Title: Joint Decomposition with Nonnegative Matrix Factorization
### Aliases: jointNMF
### Keywords: NMF joint,

### ** Examples

dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
group = list(c(1,2,3,4), c(1,2), c(3,4), c(1,3), c(2,4), c(1), c(2), c(3), c(4))
comp_num = c(2,2,2,2,2,2,2,2,2)
proj_dataset = matrix(runif(5000, 1, 2), nrow = 100, ncol = 50)
proj_group = c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)
res_jointNMF = jointNMF(
dataset,
group,
comp_num,
proj_dataset = proj_dataset,
proj_group = proj_group)





cleanEx()
nameEx("jointPCA")
### * jointPCA

flush(stderr()); flush(stdout())

### Name: jointPCA
### Title: Joint Decomposition with Principal Component Analysis
### Aliases: jointPCA
### Keywords: PCA joint,

### ** Examples

dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
group = list(c(1,2,3,4), c(1,2), c(3,4), c(1,3), c(2,4), c(1), c(2), c(3), c(4))
comp_num = c(2,2,2,2,2,2,2,2,2)
proj_dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
proj_group = list(c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE))
res_jointPCA = jointPCA(
dataset,
group,
comp_num,
proj_dataset = proj_dataset,
proj_group = proj_group)





cleanEx()
nameEx("marchenko_pastur_quantile")
### * marchenko_pastur_quantile

flush(stderr()); flush(stdout())

### Name: marchenko_pastur_quantile
### Title: Quantile for Marchenko Pastur Distribution
### Aliases: marchenko_pastur_quantile
### Keywords: Marchenko-Pastur, quantile

### ** Examples

out = marchenko_pastur_quantile(0.2, 50, 500)




cleanEx()
nameEx("mtx")
### * mtx

flush(stderr()); flush(stdout())

### Name: mtx
### Title: mtx sequence data
### Aliases: mtx
### Keywords: datasets

### ** Examples


data(mtx)
head(mtx$HS.Nico)
head(mtx$Mm.Nico)





cleanEx()
nameEx("normalizeData")
### * normalizeData

flush(stderr()); flush(stdout())

### Name: normalizeData
### Title: Data Normalization
### Aliases: normalizeData
### Keywords: normalize

### ** Examples

dataset = list(
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50)
)
normalizeData(dataset)




cleanEx()
nameEx("projectNMF")
### * projectNMF

flush(stderr()); flush(stdout())

### Name: projectNMF
### Title: Function to estimate sample embeddings for one dataset from a
###   gene loading matrix derived from an NMF decomposition of another
###   dataset.
### Aliases: projectNMF
### Keywords: NMF joint, projection,

### ** Examples

proj_dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
proj_group = c(TRUE, TRUE) # which groupings in the joint decomposition you want to project on.
list_component = jointNMF$linked_component_list # from jointNMF result
res_projNMF = projectNMF(
proj_dataset = proj_dataset,
proj_group = proj_group,
list_component = list_component)

PLEASE MAKE SURE YOUR proj_dataset AND list_component ELEMENTS HAVE MEANINGFUL ROW(GENE) NAMES - they are matched across matrices for the projection. 
#'



cleanEx()
nameEx("projectiLCA")
### * projectiLCA

flush(stderr()); flush(stdout())

### Name: projectiLCA
### Title: Function to estimate sample embeddings for one dataset from a
###   gene loading matrix derived from an iLCA analysis of another dataset.
### Aliases: projectiLCA
### Keywords: iLCA projection, twoStageiLCA,

### ** Examples

proj_dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
proj_group = c(TRUE, TRUE) # which groupings in the twoStageiLCA analysis you want to project on.
list_component = twoStageiLCA_res$linked_component_list # from twoStageiLCA result
ica_score = twoStageiLCA_res$ica_score # from twoStageiLCA result
res_projiLCA = projectiLCA(
proj_dataset = proj_dataset,
proj_group = proj_group,
list_component = list_component,
ica_score = ica_score)

PLEASE MAKE SURE YOUR proj_dataset AND list_component ELEMENTS HAVE MEANINGFUL ROW(GENE) NAMES - they are matched across matrices for the projection. 
#'



cleanEx()
nameEx("pveMultiple")
### * pveMultiple

flush(stderr()); flush(stdout())

### Name: pveMultiple
### Title: Percentage of Variance Explained for Multiple Data sets
### Aliases: pveMultiple
### Keywords: Multiple PVE, analysis

### ** Examples

dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
group = list(c(1,2,3,4), c(1,2), c(3,4), c(1,3), c(2,4), c(1), c(2), c(3), c(4))
comp_num = c(2,2,2,2,2,2,2,2,2)
res_concatPCA = concatPCA(dataset, group, comp_num)
pveMultiple(dataset, group, comp_num, res_concatPCA$score_list, res_concatPCA$linked_component_list)




cleanEx()
nameEx("pveSep")
### * pveSep

flush(stderr()); flush(stdout())

### Name: pveSep
### Title: Percentage of Variance Explained for separate data set
### Aliases: pveSep
### Keywords: PVE, Separate analysis

### ** Examples

dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
comp_num = 2
res_sepPCA = sepPCA(dataset, comp_num)
pveSep(dataset, res_sepPCA$score_list, res_sepPCA$linked_component_list)




cleanEx()
nameEx("rebalanceData")
### * rebalanceData

flush(stderr()); flush(stdout())

### Name: rebalanceData
### Title: Rebalance Data
### Aliases: rebalanceData

### ** Examples

x = list(matrix(c(1 : 4), nrow = 2), matrix(c(1 : 4), nrow = 2))
y = list(matrix(c(1 : 4), nrow = 2), matrix(c(1 : 4), nrow = 2))
list_score = list(x, y)
group = list(c(1), c(2))
dataset = list(matrix(c(1 : 8), nrow = 2), matrix(1 : 6, nrow = 2))
rebalanceData(list_score, group, dataset)




cleanEx()
nameEx("rotate_component")
### * rotate_component

flush(stderr()); flush(stdout())

### Name: rotate_component
### Title: Simulated Component Rotation
### Aliases: rotate_component
### Keywords: component rotation,

### ** Examples

component = svd(matrix(rnorm(100 * 200), nrow = 200))$u[, 1 : 2]
rotate_component(component, pi / 6)




cleanEx()
nameEx("sampleNameAssign")
### * sampleNameAssign

flush(stderr()); flush(stdout())

### Name: sampleNameAssign
### Title: Sample Name Assignment for Concatenate, Joint and TwoStageLCA
###   Analysis
### Aliases: sampleNameAssign
### Keywords: name sample,

### ** Examples

x = list(matrix(c(1 : 4), nrow = 2), matrix(c(1 : 4), nrow = 2))
y = list(matrix(c(1 : 4), nrow = 2), matrix(c(1 : 4), nrow = 2))
score_list = list(x, y)
sample_name = list(c("x.sample.1", "x.sample.2"), c("y.sample.1", "y.sample.2"))
sampleNameAssign(score_list, sample_name)




cleanEx()
nameEx("sampleNameAssignProj")
### * sampleNameAssignProj

flush(stderr()); flush(stdout())

### Name: sampleNameAssignProj
### Title: Sample Name Assignment for projectNMF
### Aliases: sampleNameAssignProj
### Keywords: name sample,

### ** Examples

x = matrix(c(1:4), nrow = 2)
y = matrix(c(1:4), nrow = 2)
score_list = list(x, y)
sample_name = list("x.sample.1", "x.sample.2")
sampleNameAssignProj(score_list, sample_name)




cleanEx()
nameEx("sampleNameAssignSep")
### * sampleNameAssignSep

flush(stderr()); flush(stdout())

### Name: sampleNameAssignSep
### Title: Sample Name Assignment for Seperate Analysis
### Aliases: sampleNameAssignSep
### Keywords: name sample,

### ** Examples

x = matrix(c(1:4), nrow = 2)
y = matrix(c(1:4), nrow = 2)
score_list = list(x, y)
sample_name = list(c("x.sample.1", "x.sample.2"), c("y.sample.1", "y.sample.2"))
sampleNameAssignSep(score_list, sample_name)




cleanEx()
nameEx("sampleNameExtractor")
### * sampleNameExtractor

flush(stderr()); flush(stdout())

### Name: sampleNameExtractor
### Title: Sample name extractor
### Aliases: sampleNameExtractor
### Keywords: extractor name,

### ** Examples

x = matrix(c(1 : 4), nrow = 2)
colnames(x) = c("sp1", "sp2")
y = matrix(c(1 : 4), nrow = 2)
colnames(y) = c("sp3", "sp4")
dataset = list(x, y)
sampleNameExtractor(dataset)




cleanEx()
nameEx("scoreNameAssign")
### * scoreNameAssign

flush(stderr()); flush(stdout())

### Name: scoreNameAssign
### Title: Score Name Assignment for Concatenate, Joint and TwoStageLCA
###   Analysis
### Aliases: scoreNameAssign
### Keywords: name score,

### ** Examples

score_list = list(
list(matrix(c(1 : 4), nrow = 2), matrix(c(1 : 4), nrow = 2)),
list(matrix(c(1 : 4), nrow = 2), matrix(c(1 : 4), nrow = 2)))
dataset_name = c("dat1", "dat2")
group_name = c("comp1", "comp2")
scoreNameAssign(score_list, dataset_name, group_name)




cleanEx()
nameEx("scoreNameAssignProj")
### * scoreNameAssignProj

flush(stderr()); flush(stdout())

### Name: scoreNameAssignProj
### Title: Score Name Assignment for projectNMF
### Aliases: scoreNameAssignProj
### Keywords: name score,

### ** Examples

score_list = list(
list(matrix(c(1 : 4), nrow = 2), matrix(c(1 : 4), nrow = 2)),
list(matrix(c(1 : 4), nrow = 2), matrix(c(1 : 4), nrow = 2)))
group_name = c("comp1", "comp2")
scoreNameAssignProj(score_list, group_name)




cleanEx()
nameEx("scoreNameAssignSep")
### * scoreNameAssignSep

flush(stderr()); flush(stdout())

### Name: scoreNameAssignSep
### Title: Score Name Assign for Seperate Analysis
### Aliases: scoreNameAssignSep
### Keywords: name score,

### ** Examples

score_list = list(matrix(c(1 : 4), nrow = 2), matrix(c(5 : 8), nrow = 2))
dataset_name = c("x", "y")
scoreNameAssignSep(score_list, dataset_name)




cleanEx()
nameEx("score_generation")
### * score_generation

flush(stderr()); flush(stdout())

### Name: score_generation
### Title: Random Score Generation
### Aliases: score_generation
### Keywords: random score

### ** Examples

score_generation(2, 10, c(1,2))




cleanEx()
nameEx("sepICA")
### * sepICA

flush(stderr()); flush(stdout())

### Name: sepICA
### Title: Single Data Set Decomposition with Independent Component
###   Analysis
### Aliases: sepICA
### Keywords: ICA analysis, separate

### ** Examples

dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
comp_num = 2
res_sepICA = sepICA(dataset, comp_num)




cleanEx()
nameEx("sepNMF")
### * sepNMF

flush(stderr()); flush(stdout())

### Name: sepNMF
### Title: Single Data Set Decomposition with Nonnegative Matrix
###   Factorization
### Aliases: sepNMF
### Keywords: NMF analysis, separate

### ** Examples

dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
comp_num = 2
res_sepNMF = sepNMF(dataset, comp_num)




cleanEx()
nameEx("sepPCA")
### * sepPCA

flush(stderr()); flush(stdout())

### Name: sepPCA
### Title: Single Data Set Decomposition with Principal Component Analysis
### Aliases: sepPCA
### Keywords: PCA analysis, separate

### ** Examples

dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
comp_num = 2
res_sepPCA = sepPCA(dataset, comp_num)




cleanEx()
nameEx("simulated_data_generation")
### * simulated_data_generation

flush(stderr()); flush(stdout())

### Name: simulated_data_generation
### Title: Simulated Data Generation
### Aliases: simulated_data_generation
### Keywords: data generation simulated

### ** Examples

configuration_setting = configuration_setting_generation()
simulated_data_generation(configuration_setting)




cleanEx()
nameEx("sjdWrap")
### * sjdWrap

flush(stderr()); flush(stdout())

### Name: sjdWrap
### Title: SJD Wrap
### Aliases: sjdWrap
### Keywords: expression gene shared

### ** Examples


data(NeuroGenesis4)
SJDdataIN = sjdWrap(
data.list = NeuroGenesis4,
species.vector=c("human","human","human","mouse"),
geneType.vector=c("symbol","ensembl","symbol","symbol"),
geneType.out="symbol",
species.out="human")




cleanEx()
nameEx("sjdWrapProjection")
### * sjdWrapProjection

flush(stderr()); flush(stdout())

### Name: sjdWrapProjection
### Title: SJD Wrap for projection
### Aliases: sjdWrapProjection
### Keywords: expression gene shared

### ** Examples


## Load NeuroGenesis4 data into R
data(NeuroGenesis4)

## sjdWrap of the training data sets
SJDdataIN = sjdWrap(
data.list = NeuroGenesis4,
species.vector=c("human","human","human","mouse"),
geneType.vector=c("symbol","ensembl","symbol","symbol"),
geneType.out="symbol",
species.out="human")

## Sample from data, serving as the projection expression matrices
NeuroGenesis4.sample = NeuroGenesis4
NeuroGenesis4.sample[[1]] = NeuroGenesis4.sample[[1]][-5,]
rownames(NeuroGenesis4.sample[[1]])[5] = paste0(rownames(NeuroGenesis4.sample[[1]])[5], ".test")

NeuroGenesis4.sample[[2]] = NeuroGenesis4.sample[[2]][-10,]
rownames(NeuroGenesis4.sample[[2]])[10] = paste0(rownames(NeuroGenesis4.sample[[2]])[10], ".test")

NeuroGenesis4.sample[[3]] = NeuroGenesis4.sample[[3]][-15,]
rownames(NeuroGenesis4.sample[[3]])[15] = paste0(rownames(NeuroGenesis4.sample[[3]])[15], ".test")

NeuroGenesis4.sample[[4]] = NeuroGenesis4.sample[[4]][-20,]
rownames(NeuroGenesis4.sample[[4]])[20] = paste0(rownames(NeuroGenesis4.sample[[4]])[20], ".test")

SJDdataProjection = sjdWrapProjection(
SJDdataIN, NeuroGenesis4.sample, "human", c("human","human","human","mouse"),
"symbol", c("symbol","ensembl","symbol","symbol"))




cleanEx()
nameEx("twoStageLCA")
### * twoStageLCA

flush(stderr()); flush(stdout())

### Name: twoStageLCA
### Title: Two-staged Linked Component Analysis
### Aliases: twoStageLCA
### Keywords: LCA two-staged,

### ** Examples

dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
group = list(c(1,2,3,4), c(1,2), c(3,4), c(1,3), c(2,4), c(1), c(2), c(3), c(4))
comp_num = c(2,2,2,2,2,2,2,2,2)
proj_dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
proj_group = list(c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE))
res_twoStageLCA = twoStageLCA(
dataset,
group,
comp_num,
proj_dataset = proj_dataset,
proj_group = proj_group)



cleanEx()
nameEx("twoStageLCA.rank")
### * twoStageLCA.rank

flush(stderr()); flush(stdout())

### Name: twoStageLCA.rank
### Title: Two-staged LCA and automatic rank selection
### Aliases: twoStageLCA.rank
### Keywords: LCA rank, two-staged,

### ** Examples

dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
group = list(c(1, 2, 3, 4), c(1, 2), c(3, 4), c(1, 3), c(2, 4), c(1), c(2), c(3), c(4))
threshold = c(3, 1.5, 1.5, 1.5, 1.5, 0.5, 0.5, 0.5, 0.5)
res_twoStageLCA.rank = twoStageLCA.rank(
dataset,
group,
threshold = threshold)





cleanEx()
nameEx("twoStageiLCA")
### * twoStageiLCA

flush(stderr()); flush(stdout())

### Name: twoStageiLCA
### Title: Two-staged Independent Linked Component Analysis
### Aliases: twoStageiLCA
### Keywords: LCA independent two-staged,

### ** Examples

dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
group = list(c(1,2,3,4), c(1,2), c(3,4), c(1,3), c(2,4), c(1), c(2), c(3), c(4))
comp_num = c(2,2,2,2,2,2,2,2,2)
proj_dataset = matrix(runif(5000, 1, 2), nrow = 100, ncol = 50)
proj_group = c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)
res_twoStageiLCA = twoStageiLCA(
dataset,
group,
comp_num,
proj_dataset = proj_dataset,
proj_group = proj_group)





cleanEx()
nameEx("twoStageiLCA.rank")
### * twoStageiLCA.rank

flush(stderr()); flush(stdout())

### Name: twoStageiLCA.rank
### Title: Two-staged Independent LCA and automatic rank selection
### Aliases: twoStageiLCA.rank
### Keywords: LCA indepentdent, rank, two-staged,

### ** Examples

dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
group = list(c(1, 2, 3, 4), c(1, 2), c(3, 4), c(1, 3), c(2, 4), c(1), c(2), c(3), c(4))
threshold = c(3, 1.5, 1.5, 1.5, 1.5, 0.5, 0.5, 0.5, 0.5)
res_twoStageiLCA.rank = twoStageiLCA.rank(
dataset,
group,
threshold = threshold)




cleanEx()
nameEx("weightData")
### * weightData

flush(stderr()); flush(stdout())

### Name: weightData
### Title: Weighting Data Set
### Aliases: weightData
### Keywords: weighting

### ** Examples

dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
weighting = c(1, 2)
weighted_dataset = weightData(dataset, weighting)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
