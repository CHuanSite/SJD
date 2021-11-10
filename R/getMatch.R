#' Species matching function
#'
#' matching genes from one species (i.e human) to another species (i.e mouse) and output genes
#'
#' @param genes vector of gene names in either ensembl or symbol format
#' @param inSpecies character of species name of input genes i.e 'human'
#' @param inType character of type of input genes i.e 'ensembl'
#' @param newSpecies character of target species for the new vector of genes i.e 'mouse'
#' @param moreAttrIn character of other gene attributes i.e 'species latin'
#' @param moreAttrNew character of other gene attributes i.e 'ensembl.nms'
#'
#' @import biomaRt
#'
#' @return vector of matched genes from the target species
#'
#' @keywords shared genes
#'
#' @examples
#'
#' data(NeuroGenesis4)
#' out = getMatch(
#' rownames(NeuroGenesis4$Meissner.inVitro.bulk.Hs),
#' inSpecies = 'human',
#' inType = 'symbol',
#' newSpecies = 'mouse')
#'
#' @export

getMatch <- function(genes, inSpecies, inType, newSpecies, moreAttrIn = NA, moreAttrNew = NA){
    species = data.frame(
        species.nm = c('human','mouse','roundworm','fruitfly','zebrafish','chicken',
                       'rat','guinea pig','golden hamster','rabbit','pig','sheep',
                       'cow','dog','cat','macaque','bonobo','chimpanzee'),
        species.latin = c('homo sapiens','mus musculus','caenorhabditis elegans','drosphila melanagaster',
                          'danio rerio','gallus gallus','ratus norvegicus','cavia porcellus',
                          'melanochromis auratus','oryctolagus cuniculus','sus scrofa domesticus',
                          'ovis aries','bos taurus','canis lupus familiaris','felis catus',
                          'macaca mulatta','pan paniscus','pan troglodytes'),
        ensembl.nms = c('hsapiens_gene_ensembl','mmusculus_gene_ensembl','celegans_gene_ensembl','dmelanogaster_gene_ensembl',
                        'drerio_gene_ensembl','ggallus_gene_ensembl','rnorvegicus_gene_ensembl',
                        'cporcellus_gene_ensembl','mauratus_gene_ensembl','ocuniculus_gene_ensembl',
                        'sscrofa_gene_ensembl','oaries_gene_ensembl','btaurus_gene_ensembl',
                        'clfamiliaris_gene_ensembl','fcatus_gene_ensembl','mmulatta_gene_ensembl',
                        'ppaniscus_gene_ensembl','ptroglodytes_gene_ensembl'),
        genesymbol.attr = c('hgnc_symbol','mgi_symbol','external_gene_name','external_gene_name','hgnc_symbol',
                            'hgnc_symbol','external_gene_name','hgnc_symbol','mgi_symbol','hgnc_symbol','hgnc_symbol',
                            'hgnc_symbol','hgnc_symbol','hgnc_symbol','hgnc_symbol','hgnc_symbol','hgnc_symbol',
                            'hgnc_symbol'),
        row.names = 1
    )
    # print('d1')
    mart.in = species[inSpecies, ]$ensembl.nms
    mart.new = species[newSpecies, ]$ensembl.nms
    # print('d2')
    df.in = useMart(biomart='ensembl', dataset=mart.in)
    df.new = useMart(biomart='ensembl', dataset=mart.new)
    symbol.in = as.character(species[inSpecies,'genesymbol.attr'])
    symbol.new = as.character(species[newSpecies,'genesymbol.attr'])
    # print('d3')
    if (inType == 'symbol') {
        filter = as.character(species[inSpecies,'genesymbol.attr'])
    } else if (inType == 'ensembl') {
        filter = 'ensembl_gene_id'
    }
    # print('d4')
    cat('You have input ',length(genes),' genes\n')

    atb.in = c(symbol.in, "ensembl_gene_id","description")
    atb.new = c(symbol.new,"ensembl_gene_id","description")
    if (!is.na(moreAttrIn)) {atb.in = c(atb.in, moreAttrIn)}
    if (!is.na(moreAttrNew)) {atb.new = c(atb.new, moreAttrNew)}
    # print(species)
    # print(symbol.in)
    # print(atb.in)
    tbl.match = getLDS(
        mart=df.in, attributes=atb.in, # input genes
        filters=filter, values=genes,
        martL=df.new, attributesL=atb.new
    )
    # print('d5')

    tbl.match[,'Gene.description'] = gsub(' \\[.*','',tbl.match[,'Gene.description'])
    tbl.match[,'Gene.description.1'] = gsub(' \\[.*','',tbl.match[,'Gene.description.1'])

    colnames(tbl.match)[1:length(atb.in)] = paste0(colnames(tbl.match)[1:length(atb.in)],'.',inSpecies)
    colnames(tbl.match)[(length(atb.in)+1):(length(atb.in)+length(atb.new))] = paste0(colnames(tbl.match)[(length(atb.in)+1):(length(atb.in)+length(atb.new))],'.',newSpecies)

    cat('We found ',dim(tbl.match)[1],' matches\n')
    cat(sum(duplicated(tbl.match[,1])),' of those are duplicates and only keeping the 1st of each\n')
    if (inType == 'symbol') {
        tbl.match = tbl.match[match(genes,tbl.match[,1]),]
        rownames(tbl.match) = NULL
        tbl.return = cbind(genes,tbl.match,stringsAsFactors=FALSE)
    } else if (inType == 'ensembl') {
        tbl.match = tbl.match[match(genes,tbl.match[,2]),]
        rownames(tbl.match) = NULL
        tbl.return = cbind(genes,tbl.match,stringsAsFactors=FALSE)
    }
    return(tbl.return)
}
