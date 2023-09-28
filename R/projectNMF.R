#' Projection function for Joint Decomposition with Nonnegative Matrix Factorization
#' Purpose: finding out the common development patterns between new data and datasets used in the joint decomposition.
#' If you have new dataset(s) and Gene scores produced by JointNMF, this function can produce estimated sample scores for the new datasets.

#' @param proj_dataset The dataset(s) to be projected on. 
#' @param proj_group A listed of boolean combinations indicating which groupings should be used for each projected dataset.The length of proj_group should match the length of proj_dataset, and the length of each concatenated boolean combination should match the length of the parameter group.
#' @param list_component list_component produced from JointNMF decomposition.  
#' @param weighting Weighting of each dataset, initialized to be NULL
#' @param max_ite The maximum number of iterations for the jointNMF algorithms to run, default value is set to 100
#' @param max_err The maximum error of loss between two iterations, or the program will terminate and return, default value is set to be 0.0001
#' @param enable_normalization An argument to decide whether to use normalizaiton or not,  default is TRUE
#' @param column_sum_normalization An argument to decide whether to use column sum normalization or not, default it FALSE
#' 
#' @return A list contains the projected scores of each dataset on every component.
#'
#' @keywords projection, joint, NMF
#'
#' @examples
#' Example 1 (1 matrix in proj_dataset):
#' proj_dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
#' proj_group = list(c(TRUE, TRUE)) # which groupings in the joint decomposition you want to project on.
#' list_component = jointNMF$linked_component_list
#' res_projNMF = projectNMF(
#' proj_dataset = proj_dataset,
#' proj_group = proj_group,
#' list_component = list_component)
#' 
#' Example 2 (2 matrix in proj_dataset):
#' proj_dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50), 
#'                     matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
#' proj_group = list(c(TRUE, TRUE), c(TRUE, TRUE))
#' list_component = jointNMF$linked_component_list
#' res_projNMF = projectNMF(
#' proj_dataset = proj_dataset,
#' proj_group = proj_group,
#' list_component = list_component)



## Projection function
projectNMF <- function(proj_dataset, proj_group, list_component, weighting=NULL, max_ite = 1000, max_err = 0.0001, enable_normalization = TRUE, column_sum_normalization = FALSE){
  p = nrow(proj_dataset[[1]])
  # Initialize the matrix W with zeros
  W <- matrix(0, nrow = p, ncol = sum(comp_num))
  K = length(proj_group[[1]])
  # Reconstruct matrix W using list_component and comp_num
  start_col <- 1
  for (i in 1:K) {
    end_col <- start_col + comp_num[i] - 1
    W[, start_col:end_col] <- list_component[[i]]
    start_col <- end_col + 1
  }
  comp_num = unlist(lapply(list_component, function(x) if (length(x) > 0) ncol(x)))
  
  
  proj_list_score = list()
  if(!is.null(proj_dataset)){
    proj_sample_name = sampleNameExtractor(proj_dataset)
    proj_dataset_name = datasetNameExtractor(proj_dataset)
    group_name = paste0("group", 1:length(proj_group[[1]])) # group name alternative 
    proj_dataset = frameToMatrix(proj_dataset)
    proj_dataset = normalizeData(proj_dataset, enable_normalization, column_sum_normalization, nonnegative_normalization = TRUE)
    proj_dataset = balanceData(proj_dataset)
    proj_dataset = weightData(proj_dataset, weighting)
    
    N = length(proj_dataset)
    M = sum(comp_num) 
    p = nrow(proj_dataset[[1]])
    N_dataset = unlist(lapply(proj_dataset, ncol))
    
    
    ## Initialize the W and H for Nonnegative Matrix Factorization
    max_element = -Inf
    min_element = Inf
    for(i in 1 : N){
      max_element = max(max_element, max(proj_dataset[[i]]))
      min_element = min(min_element, min(proj_dataset[[i]]))
    }
    
    
    ## Initialize the values of W and H
    X = c()
    for(i in 1 : N){
      X = cbind(X, proj_dataset[[i]])
    }
    H = c()
    
    
    for(i in 1 : K){ # K= length(group) = 3 in the example
      H_temp = c()
      for(j in 1 : N){ # K = length(proj_dataset) = 1 in the example
        H_temp = cbind(H_temp, matrix(runif(N_dataset[j] * comp_num[i], min_element, max_element), nrow = comp_num[i], ncol = N_dataset[j]))
      }
      H = rbind(H, H_temp)
    }
    
    
    ## Iteratively estimate the NMF with Euclidean distance
    error_out = c()
    
    for(ite in 1 : max_ite){
      ## H is score (sample matrix S)
      H = H * (t(W) %*% X) / (t(W) %*% W %*% H)
      ## W is list_component (common gene matrix G)
      
      H[which(is.na(H))] = 0
      
      error_out = c(error_out, sum((X - W %*% H)^2))
      
      ## Break when the error difference is small
      if(length(error_out) >= 2 && abs(error_out[length(error_out)] - error_out[length(error_out) - 1]) / abs(error_out[length(error_out) - 1]) <= max_err){
        break
      }
    }
    
    
    for(i in 1 : N){  # i goes to 1 if we only have 1 dataset
      proj_list_score[[i]] = list()
      for(j in 1 : length(proj_group[[i]])){ # j goes from 1 to 3 since we have 3 groups
        if(proj_group[[i]][j]){
          proj_list_score[[i]][[j]] = H[ifelse(j == 1, 1, cumsum(comp_num)[j - 1] + 1) : cumsum(comp_num)[j], ifelse(i == 1, 1, cumsum(N_dataset)[i - 1] + 1) : cumsum(N_dataset)[i]]
        }
      }
    } 
    
    proj_list_score = scoreNameAssign(proj_list_score, proj_dataset_name, group_name)
    proj_list_score = sampleNameAssign(proj_list_score, proj_sample_name)
    for(i in 1 : length(proj_group[[1]])){
      for(j in 1 : N){
        proj_list_score[[j]][[i]] = proj_list_score[[j]][[i]] *  sqrt(ncol(proj_dataset[[j]]))
      }
    }
  }
  
  return(list(proj_score_list = proj_list_score))
}