

#' Map mols to inchi
#'
#' @param mols A list of molecules
#' @return A list of inchi
mol2Inchi <- function( mols ){
  p_vectorize( mols,  p_mol2Inchi)  
}

#' Map Inchi to  Inchikeys
#'
#' @param mols A list of molecules
#' @return A list of inchikeys
Inchi2InchiKey <- function( mols ){
  p_vectorize( mols,  p_Inchi2InchiKey)  
}
 
#' Compare Inchikeys
#'
#' @param inchikey1 A inchikey
#' @param inchikey2 A inchikey
#' @param check.protonation A check protonation
#' @return Bool
compareInchiKeys<-function( inchikey1, inchikey2 , check.protonation = T ){
  # details on protonation check
  # are in    http://www.inchi-trust.org/wp/wp-content/uploads/2014/06/INCHI-1-API.zip
  # file INCHI_API/inchi_dll/ikey_dll.c
  # function GetStdINCHIKeyFromStdINCHI
  
  if(!check.protonation){
    inchikey1 <- substr(inchikey1,0,24)
    inchikey2 <- substr(inchikey2,0,24)
  }
  return( ( inchikey1 == inchikey2 )[[1]])
}

#' Map mols to Inchikeys
#'
#' @param mols A list of molecules
#' @return A list of inchikeys
mol2InchiKey<-function(mols){
  inchis<-mol2Inchi(mols)
  return(Inchi2InchiKey(inchis))
}
