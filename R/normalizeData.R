#' Data Normalization
#'
#' Normalize data to have mean zero and std 1
#'
#' @param dataset The input list of data sets matrix
#' @param enable_normalization An argument to decide whether to use normalizaiton or not,  default is TRUE
#' @param column_sum_normalization An argument to decide whether to use column sum normalization or not, default it TRUE
#' @param nonnegative_normalization An argument to decide wehther it is nonnegative matrix factorization based or not, default is FALSE
#' @keywords normalize
#'
#' @examples
#' dataset = list(
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50)
#' )
#' normalizeData(dataset)
#'
#' @export

normalizeData <- function(dataset, enable_normalization = TRUE, column_sum_normalization = TRUE, nonnegative_normalization = FALSE){
    if(enable_normalization){
        if(nonnegative_normalization){
            dataset = lapply(dataset, FUN = function(x){
                cs=colSums(x);cs=cs/min(cs);sweep(x,2,cs,"/")
            })
            return(dataset)
        }
        if(column_sum_normalization){
            dataset = lapply(dataset, FUN = function(x){
                t(scale(t(cs=colSums(x);cs=cs/min(cs);sweep(x,2,cs,"/")), scale = FALSE))
            })
        }else{
            dataset = lapply(dataset, FUN = function(x){
                scale(t(scale(t(x), scale = FALSE)))
            })
        }
    }else{
        dataset = lapply(dataset, FUN = function(x){
            t(scale(t(x), scale = FALSE))
        })
    }
    return(dataset)
}
