#' Find the best tree among the trees with the number of nodes at most equal to 
#' the sample size
#'
#' @description A function that finds the best tree among all the trees with the number of nodes 
#' at most equal to the sample size
#'
#' @param MSPEall a vector with MSPE values for different trees.
#'
#' @param leng a vector of the same length as \code{MSPEall} indicating the number of 
#' leaf nodes in each tree.
#'
#' @param n the sample size.
#'
#' @return A number indicating the index of the best tree. Note that 1 is returned if 
#' all trees have more than \code{n} nodes.
#'
#' @source Cand{\`e}s, E. and Tao, T. (2007). The Dantzig selector: Statistical estimation when p is much
#' larger than n. Annals of Statistics 35 (6), 2313--2351.
#' 
#' @source Phoa, F. K., Pan, Y. H. and Xu, H. (2009). Analysis of supersaturated 
#' designs via the Dantzig selector. Journal of Statistical Planning and Inference
#' 139 (7), 2362--2372.
#' 
#' 
#'
#' @keywords internal

  find_tree <- function(MSPEall, leng, n){
  # Subset of trees that have more nodes than n
  sub_indexes <- which(leng <  n)
  
  if(length(sub_indexes) == 0){
    # Return the index of the first tree if no trees have less nodes than n
    return(1)
  }
  
  # Subset of MSPE values for trees that have less nodes than n
  sub_MSPEall <- MSPEall[sub_indexes]
  
  # Find the index of the tree with minimum MSPE in the subset
  index_min_subset <- which.min(sub_MSPEall)
  
  # Find the index in the original vector
  index_min <- sub_indexes[index_min_subset]
  
  # Return the index
  return(index_min)
}
