#' Wrapper function for PJD package
#'
#' A wrapper function to align multiple datasets from different species
#'
#' @param data.list Input list of data sets
#' @param species.vector A vector of species types
#' @param geneType.vector A vector of gene types,
#' @param geneType.out Type of genes
#' @param species.out Kind of species output
#'
#' @import biomaRt
#' @importFrom dplyr inner_join
#'
#' @return A list of data sets, after alignment
#'
#' @keywords Wrapper, PJD
#'
#' @examples
#'
#' data(mtx)
#'
#' PJD_example = pjd_wrap(data.list = mtx,
#' species.vector = c("human", "mouse"),
#' geneType.vector = c("symbol", "symbol"),
#' geneType.out = "symbol",
#' species.out = "human")
#'
#' @export


pjd_wrap <- function(data.list, species.vector, geneType.vector, geneType.out="symbol", species.out){
    N = length(data.list) #number of datasets

    if(N != length(species.vector)){
        print("length of data.list != length(species.vector)")
        return(NULL)
    }

    if(N != length(geneType.vector)){
        print("length of data.list != length(geneType.vector)")
        return(NULL)
    }

    genes.tbl.list = list()
    new.data.list = list()

    cat('Using biomaRt to connect gene IDs across', N, 'datasets:\n')

    for (i in 1 : N) {

        cat('Getting biomaRt IDs for dataset',i,'\n')

        species.1 = species.vector[i]

        genes = rownames(data.list[[i]])

        genes.tbl = getMatch(genes, inSpecies=species.1, inType=geneType.vector[i], newSpecies = species.out)
        genes.tbl = na.omit(genes.tbl)

        genes.tbl$ORIGid = genes.tbl[, "genes"]

        if(geneType.out == "symbol") {
            genes.tbl[, "genes"] = genes.tbl[, 5]
            }#geneType.out

        if(geneType.out == "ensembl"){
            genes.tbl[, "genes"] = genes.tbl[, 6]
            }#geneType.out

        indxNoDup = !duplicated(genes.tbl[, "genes"]) & !is.na(genes.tbl[, "genes"]) & genes.tbl[, "genes"] != ""
        genes.tbl.list[[i]] = genes.tbl[indxNoDup, c("genes", "ORIGid")]
        rm(genes.tbl)
    }

    cat('constructed', length(genes.tbl.list), 'tables of cross-species matching genes\n')

    genes.tbl = inner_join(genes.tbl.list[[1]], genes.tbl.list[[2]], by = 'genes', na_matches = "never")#geneType.out

    if(N >= 3){
        for(i in 3 : N){
            genes.tbl = inner_join(genes.tbl, genes.tbl.list[[i]], by = 'genes', na_matches = "never")#geneType.out
        }
    }

    cat('we found', dim(genes.tbl)[1], 'shared genes in', length(genes.tbl.list), 'datasets\n')

    for (i in 1 : N) {
        indxR = match(genes.tbl[,i+1], rownames(data.list[[i]]))#ORIGid
        data = data.list[[i]][indxR, ]
        rownames(data) = genes.tbl[, "genes"]
        new.data.list[[i]] = as.data.frame(data)
        rm(data)
    }

    cat('new data list of', length(new.data.list), 'datasets constructed\n')

    names(new.data.list) = names(data.list)

    return(new.data.list)
}
