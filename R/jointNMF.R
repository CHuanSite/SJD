#' Joint Decomposition with Nonnegative Matrix Factorization
#'
#' Joint decomposition of several linked matrices with Nonnegative Matrix Factorization (NMF)
#' It is based on the MSE loss, proposed by Lee, Daniel D., and H. Sebastian Seung. "Learning the parts of objects by non-negative matrix factorization." Nature 401.6755 (1999): 788-791.
#'
#' @param dataset A list of dataset to be analyzed
#' @param group A list of grouping of the datasets, indicating the relationship between datasets
#' @param comp_num A vector indicates the dimension of each compoent
#' @param weighting Weighting of each dataset, initialized to be NULL
#' @param max_ite The maximum number of iterations for the jointNMF algorithms to run, default value is set to 100
#' @param max_err The maximum error of loss between two iterations, or the program will terminate and return, default value is set to be 0.0001
#' @param proj_dataset The datasets to be projected on
#' @param proj_group The grouping of projected data sets
#' @param enable_normalization An argument to decide whether to use normalizaiton or not,  default is TRUE
#' @param column_sum_normalization An argument to decide whether to use column sum normalization or not, default it TRUE
#' @param screen_prob A vector of probabilies for genes to be chosen
#'
#' @return A list contains the component and the score of each dataset on every component after jointNMF algorithm
#'
#' @keywords joint, NMF
#'
#' @examples
#' dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
#' group = list(c(1,2,3,4), c(1,2), c(3,4), c(1,3), c(2,4), c(1), c(2), c(3), c(4))
#' comp_num = c(2,2,2,2,2,2,2,2,2)
#' proj_dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
#' proj_group = list(c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE))
#' res_jointNMF = jointNMF(
#' dataset,
#' group,
#' comp_num,
#' proj_dataset = proj_dataset,
#' proj_group = proj_group)
#'
#' @export

jointNMF <- function(dataset, group, comp_num, weighting = NULL, max_ite = 100, max_err = 0.0001, proj_dataset = NULL, proj_group = NULL, enable_normalization = TRUE, column_sum_normalization = TRUE, screen_prob = NULL){

    ## Obtain names for dataset, gene and samples
    dataset_name = datasetNameExtractor(dataset)
    gene_name = geneNameExtractor(dataset)
    sample_name = sampleNameExtractor(dataset)
    group_name = groupNameExtractor(group)

    ## Preprocess Dataset
    dataset = frameToMatrix(dataset)
    if(!is.null(screen_prob)){
        dataset = geneScreen(dataset, screen_prob)
    }

    dataset = normalizeData(dataset, enable_normalization, column_sum_normalization, nonnegative_normalization = TRUE)
    dataset = balanceData(dataset)
    dataset = weightData(dataset, weighting)

    ## Initialize values for the algorithm
    N = length(dataset)
    K = length(group)
    M = sum(comp_num)
    p = nrow(dataset[[1]])
    N_dataset = unlist(lapply(dataset, ncol))

    ## Initialize the W and H for Nonnegative Matrix Factorization
    max_element = -Inf
    min_element = Inf
    for(i in 1 : N){
        max_element = max(max_element, max(dataset[[i]]))
        min_element = min(min_element, min(dataset[[i]]))
    }

    ## Initialize the values of W and H
    X = c()
    for(i in 1 : N){
        X = cbind(X, dataset[[i]])
    }
    W = matrix(runif(p * M, min_element, max_element), nrow = p)
    H = c()

    for(i in 1 : K){
        H_temp = c()
        for(j in 1 : N){
            if(j %in% group[[i]]){
                H_temp = cbind(H_temp, matrix(runif(N_dataset[j] * comp_num[i], min_element, max_element), nrow = comp_num[i], ncol = N_dataset[j]))
            }else{
                H_temp = cbind(H_temp, matrix(0, nrow = comp_num[i], ncol = N_dataset[j]))
            }
        }
        H = rbind(H, H_temp)
    }

    ## Iteratively estimate the NMF with Euclidean distance
    error_out = c()

    for(ite in 1 : max_ite){
        H = H * (t(W) %*% X) / (t(W) %*% W %*% H)
        W = W * (X %*% t(H)) / (W %*% H %*% t(H))

        H[which(is.na(H))] = 0
        W[which(is.na(W))] = 0

        error_out = c(error_out, sum((X - W %*% H)^2))

        ## Break when the error difference is small
        if(length(error_out) >= 2 && abs(error_out[length(error_out)] - error_out[length(error_out) - 1]) <= max_err){
            break
        }
        # print(ite)
        # print(abs(error_out[length(error_out)] - error_out[length(error_out) - 1]))
    }

    ## Output component and scores
    list_component = list()
    list_score = list()
    for(j in 1 : N){
        list_score[[j]] = list()
    }

    for(i in 1 : K){
        list_component[[i]] = W[, ifelse(i == 1, 1, cumsum(comp_num)[i - 1] + 1) : cumsum(comp_num)[i]]
        for(j in 1 : N){
            list_score[[j]][[i]] = H[ifelse(i == 1, 1, cumsum(comp_num)[i - 1] + 1) : cumsum(comp_num)[i], ifelse(j == 1, 1, cumsum(N_dataset)[j - 1] + 1) : cumsum(N_dataset)[j]]
        }
    }

    ## Assign name for components
    list_component = compNameAssign(list_component, group_name)
    list_component = geneNameAssign(list_component, gene_name)
    list_score = scoreNameAssign(list_score, dataset_name, group_name)
    list_score = sampleNameAssign(list_score, sample_name)
    list_score = filterNAValue(list_score, dataset, group)
    list_score = rebalanceData(list_score, group, dataset)

    ## Project score
    proj_list_score = list()
    if(!is.null(proj_dataset)){
        proj_sample_name = sampleNameExtractor(proj_dataset)
        proj_dataset_name = datasetNameExtractor(proj_dataset)

        proj_dataset = frameToMatrix(proj_dataset)
        proj_dataset = normalizeData(proj_dataset)
        proj_dataset = balanceData(proj_dataset)

        for(i in 1 : length(proj_dataset)){
            proj_list_score[[i]] = list()
            for(j in 1 : length(proj_group[[i]])){
                if(proj_group[[i]][j]){
                    proj_list_score[[i]][[j]] =  t(list_component[[j]]) %*% proj_dataset[[i]]
                }else{
                    proj_list_score[[i]][[j]] = NA
                }
            }
        }


        proj_list_score = scoreNameAssign(proj_list_score, proj_dataset_name, group_name)
        proj_list_score = sampleNameAssign(proj_list_score, proj_sample_name)
    }

    return(list(linked_component_list = list_component, score_list = list_score, proj_score_list = proj_list_score))
}
