# distances

  
tanimoto <- function(l1,l2){
  if( is.vector(l1) & is.vector(l2)){
  return( sum(l1 * l2)/( sum(l1*l1)+sum(l2*l2)-sum(l1 * l2)))
  }else{
    stop("Input are not vectors")
  }
}