
## Projection function
projectNMF <- function(proj_dataset, proj_group,p, group, comp_num, list_component, weighting=NULL, max_ite = 1000, max_err = 0.0001, enable_normalization = TRUE, column_sum_normalization = FALSE, screen_prob = NULL){
  
  # Initialize the matrix W with zeros
  W <- matrix(0, nrow = p, ncol = sum(comp_num))
  K = length(group)
  # Reconstruct matrix W using list_component and comp_num
  start_col <- 1
  for (i in 1:K) {
    end_col <- start_col + comp_num[i] - 1
    W[, start_col:end_col] <- list_component[[i]]
    start_col <- end_col + 1
  }

  
  proj_list_score = list()
  if(!is.null(proj_dataset)){
    proj_sample_name = sampleNameExtractor(proj_dataset)
    proj_dataset_name = datasetNameExtractor(proj_dataset)
    group_name = groupNameExtractor(group)
    proj_dataset = frameToMatrix(proj_dataset)
    proj_dataset = normalizeData(proj_dataset, enable_normalization, column_sum_normalization, nonnegative_normalization = TRUE)
    proj_dataset = balanceData(proj_dataset)
    proj_dataset = weightData(proj_dataset, weighting)
    
    N = length(proj_dataset)
    M = sum(comp_num) ## Not sure
    p = nrow(proj_dataset[[1]])
    N_dataset = unlist(lapply(proj_dataset, ncol))
    
    # W = list_component[[1]] ##NOT SURE ABOUT THIS W is still W, unchanged
    
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
    
    
    for(i in 1 : K){ # K= length(group) = 3 in the toy example
      H_temp = c()
      for(j in 1 : N){
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
      # print(ite)
      # print(abs(error_out[length(error_out)] - error_out[length(error_out) - 1]) / abs(error_out[length(error_out) - 1]))
    }
    
    
    for(i in 1 : N){  # i goes to 1 if we only have 1 dataset
      proj_list_score[[i]] = list()
      for(j in 1 : length(proj_group[[i]])){ # j goes from 1 to 3 since we have 3 groups
        if(proj_group[[i]][j]){
          # I switched i and j in this line
          proj_list_score[[i]][[j]] = H[ifelse(j == 1, 1, cumsum(comp_num)[j - 1] + 1) : cumsum(comp_num)[j], ifelse(i == 1, 1, cumsum(N_dataset)[i - 1] + 1) : cumsum(N_dataset)[i]]
        }
      }
    } 
    
    proj_list_score = scoreNameAssign(proj_list_score, proj_dataset_name, group_name)
    proj_list_score = sampleNameAssign(proj_list_score, proj_sample_name)
    proj_list_score = filterNAValue(proj_list_score, proj_dataset, group)
    proj_list_score = rebalanceData(proj_list_score, group, proj_dataset)
  }
  
  return(list(proj_score_list = proj_list_score))
}