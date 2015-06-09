# distances

#' Calculate tanimoto distance of two bool vectors
#'
#' @param l1 a vector
#' @param l2 a vector
#' @return Tanimoto distance 
tanimoto <- function(l1,l2){
  if( is.vector(l1) & is.vector(l2)){
  return( sum(l1 * l2)/( sum(l1*l1)+sum(l2*l2)-sum(l1 * l2)))
  }else{
    stop("Input are not vectors")
  }
}