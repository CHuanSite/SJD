#' Data Normalization
#'
#' Normalize data to have mean zero and std 1
#'
#' @param dataset The input list of data sets matrix
#' @param enable_normalization An argument to decide whether to use normalizaiton or not,  default is TRUE
#'
#' @keywords normalize
#'
#' @examples
#' dataset = list(
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50)
#' )
#' normalizeData(dataset, TRUE)
#'
#' @export

normalizeData <- function(dataset, enable_normalization = TRUE){
    if(enable_normalization){
        dataset = lapply(dataset, FUN = function(x){
            scale(t(scale(t(x), scale = FALSE)))
        })
    }else{
        dataset = lapply(dataset, FUN = function(x){
            t(scale(t(x), scale = FALSE))
        })
    }
    return(dataset)
}
