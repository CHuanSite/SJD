#' Function to estimate sample embeddings for one dataset from a gene loading matrix derived from an NMF decomposition of another dataset.
#' 
#' projectNMF estimates the embeddings for samples in a new dataset when given a gene loading matrix from an NMF decomposition of another single matrix, 
#' or set of matrices (e.g. the "list_component" from a jointNMF output object)
#'
#' 
#' @param proj_dataset The dataset(s) to be projected on. 
#' @param sample_score a matrix of sample scores 
#' @param list_component a single matrix of gene loadings as a list element, or a list_component produced from a jointNMF() decomposition.
#' @param max_ite The maximum number of iterations for the jointNMF algorithms to run, default value is set to 1000
#' @param max_err The maximum error of loss between two iterations, or the program will terminate and return, default value is set to be 0.0001
#' @param enable_normalization An argument to decide whether to use normalizaiton or not,  default is TRUE
#' @param column_sum_normalization An argument to decide whether to use column sum normalization or not, default it FALSE
#' 
#' @return A list that contains the 1] projected scores of each dataset on every component. and 2] the log of errors as the NMF was iterated.
#'
#' @keywords projection, joint, NMF
#'
#' @examples
#' proj_dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
#' proj_group = c(TRUE, TRUE) # which groupings in the joint decomposition you want to project on.
#' list_component = jointNMF$linked_component_list # from jointNMF result
#' res_projNMF = projectNMF(
#' proj_dataset = proj_dataset,
#' proj_group = proj_group,
#' list_component = list_component)
#' 
#' PLEASE MAKE SURE YOUR proj_dataset AND list_component ELEMENTS HAVE MEANINGFUL ROW(GENE) NAMES - they are matched across matrices for the projection. 
#' #'
#' @export


## Projection function


## Notes for myself: 
## multiple sample scores are produced in projectNMF (because it's multiple matrices), so here how should people input sample scores? 
## Should I make this for only single matrices? Would that be more useful? 
## the values they pass in are optional
## dimension check so that the G = G in data, and S = S in data. 

## DIM OF GENE_SOCRE AND SAMPLE SHOULD num of col < p 
## User should initialize gene score and sample with values they want in columns & rows they want to put in. The rest of the gene score and sample score should be initialize to 0.
## BOTH gene_score and sample score should not be null
newNMF <- function(proj_dataset, gene_score=NULL, sample_score=NULL, normalize = TRUE, max_ite = 1000, max_err = 0.0001, enable_normalization = TRUE, column_sum_normalization = FALSE){
  if(!is.null(proj_dataset)){ 
    if(!is.null(gene_score) && !is.null(sample_score)){
      
      num_patterns_g = ncol(gene_score)
      num_patterns_s = nrow(sample_score)
      
      if(num_patterns_g != num_patterns_s){
        stop("ncol(gene_score) must match nrow(sample_socre)")
      }
      
      proj_dataset = as.matrix(proj_dataset)
      
      p = nrow(proj_dataset)
      col = ncol(proj_dataset)
      
      non_zero_cols_index <- which(colSums(gene_score) > 0)
      non_zero_cols <- gene_score[, non_zero_cols_index]
      
      non_zero_rows_index <- which(rowSums(sample_score) > 0)
      non_zero_rows <- sample_score[non_zero_rows_index, ] 
      
      
      ## Initialize the W and H for Nonnegative Matrix Factorization
      max_element = -Inf
      min_element = Inf
      max_element = max(max_element, max(proj_dataset))
      min_element = min(min_element, min(proj_dataset))
      
      # Initialize the matrix W with random values (gene score)
      W <- matrix(runif(p * num_patterns_g, min(non_zero_cols), max(non_zero_cols)), nrow = p, ncol = num_patterns_g)
      # Initialize the matrix H with random values (sample score)
      H <- matrix(runif(num_patterns_g * col, min(non_zero_rows), max(non_zero_rows)), nrow = num_patterns_g, ncol = col)
      
      #extract non-zero columns and rows
      W[, non_zero_cols_index] = non_zero_cols
      H[non_zero_rows_index, ] = non_zero_rows
      
    
      # create masks for H and W
      mask_W <- matrix(1, nrow = p, ncol = num_patterns_g)
      mask_W[, non_zero_cols_index] = 0
      mask_H <- matrix(1, nrow = num_patterns_g, ncol= col)
      mask_H[non_zero_rows_index, ] = 0
      
      
      proj_sample_name = sampleNameExtractor(proj_dataset)
      
      if (normalize){
        proj_dataset = normalizeData(list(proj_dataset), enable_normalization, column_sum_normalization, nonnegative_normalization = TRUE)[[1]]
      }
      
      X = proj_dataset
      
      ## Iteratively estimate the NMF with Euclidean distance
      error_out = c()
      
      W_temp <- W
      H_temp <- H
      
      for(ite in 1 : max_ite){
        ## H is score (sample matrix S)
        H_temp = H_temp * (t(W_temp) %*% X) / (t(W_temp) %*% W_temp %*% H_temp + mask_H)
        W_temp = W_temp * (X %*% t(H_temp)) / (W_temp %*% H_temp %*% t(H_temp) + mask_W) 
        
        H_temp[non_zero_rows_index, ] = non_zero_rows
        W_temp[, non_zero_cols_index] = non_zero_cols
        
        
        error_out = c(error_out, sum((X - W_temp %*% H_temp)^2))
        
        ## Break when the error difference is small
        if(length(error_out) >= 2 && abs(error_out[length(error_out)] - error_out[length(error_out) - 1]) / abs(error_out[length(error_out) - 1]) <= max_err){
          break
        }
      }
      
      W <- W_temp
      H <- H_temp
      
      H <- t(apply(t(H), 2, min_max_normalization))
      W = W * (X %*% t(H)) / (W %*% H %*% t(H)) 
      
      # proj_list_score = scoreNameAssignProj(proj_list_score, group_name)
      # proj_list_score = sampleNameAssignProj(proj_list_score, proj_sample_name)
      
      # for(i in 1 : length(proj_group)){
      #   if(proj_group[i]){
      #     proj_list_score[[i]] = proj_list_score[[i]] * sqrt(ncol(proj_dataset))
      #   }
      # }
      
      return(list(res_gene_score=W, res_sample_score=H, error_out=error_out))
      
    } else{
      stop("gene_score and sample_score cannot be NULL")
    }
  } else{
    stop("proj_dataset cannot be NULL")
  }
}
