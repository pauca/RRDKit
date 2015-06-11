# This file is part of RRDKit.
# 
# RRDKit is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Foobar is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Foobar.  If not, see <http://www.gnu.org/licenses/>.


p_vectorize <- function( mols , foo , ... ){
  if(is.list(mols)){
    if( any(sapply(mols, function(m){ is.list(m)}))){
      stop("Invalid input!")
    }
  }
  sapply(mols,  foo)  
}


p_molSupplierApply <- function( molSupplier, fun, ...){
  l1 <- list()
  while(!p_molSupplier_atEnd(molSupplier)){
    mol <- p_molSupplier_next(molSupplier)
    res <- tryCatch( fun(mol,...), error= function(e) NA )
    l1[[length(l1)+1]]<-res
  }
  p_molSupplier_reset(molSupplier)
  return(l1)
}

#' Get unique BRICS fragments for a molecule
#'
#' @param mol a molecule
#' @return Vector of BRICS fragments
getBRICSFragments<-function( mol ){
  f <- fragmentOnBRICSBonds(mol)
  ff <- unique(unlist(strsplit(mol2smiles(f),"\\.")))
  return(ff)
}

#' Get FP defined by SMARTS
#'
#' @param mols A list of molecules
#' @param smarts.mols A list of SMARTS molecules
#' @return A bool matrix with one row per molecule and each column a SMARTS Fingerprint
getSMARTSFP <- function(  mols, smarts.mols ){
  r <- matrix(NA, ncol=length(smarts.mols),nrow=length(mols))
  for( i in 1:length(mols)){
    m<- mols[[i]]
    for( j in 1:length(smarts.mols)){
      sm <- NA
      try({ 
        sm <- SubstructMatch(m, smarts.mols[[j]])  
      },silent=T)
      r[i,j] <- sm
      
    }
  }
  r
}


p_getValidSubrfagemnts <- function(mol,mf){
  indx <- unlist(sapply(mf,function(f){
    SubstructMatch(mol,smarts2mol(f))
  }))
  return(mf[indx])
}

#' Clean a SVG provided by RDKit
#'
#' Cleaning conist on resizing
#'
#' @param svg a svg
#' @param out.w output width
#' @param out.h output height
#' @return a SVG
cleanSVG<-function( svg , out.w, out.h ){
  root<-xmlTreeParse(svg,asText = T)
  w <- xmlGetAttr( root$doc$children[[1]],"width")
  w<-ifelse(is.null(w),"100",w)
  w <- gsub("px","",w)
  h <- xmlGetAttr( root$doc$children[[1]],"height")
  h<-ifelse(is.null(h),"100",h)
  h <- gsub("px","",h)
  
  xmlAttrs(root$doc$children[[1]]) <- c(width = out.w,
                                        height = out.h,
                                        viewBox = paste( " 0 0 ", w, h, sep=" "))
  
  svg2 <- toString(root$doc$children$svg)
  svg2 <- gsub("\n","",svg2)
  svg2 <- gsub("svg:","",svg2)
  svg2 <- gsub("\"","\'",svg2)
  return(svg2)
}

#' Check if object is a "nice" RDKit Molecule
#'
#' @param mol A molecule
is.molecule<-function(mol){
  # check if molecule is sanitized
    result <- F    
    tryCatch({
      kekulize(mol)
      s <- smiles2mol(mol2smiles(mol))      
      result <- T
    }, error = function(e) {})    
    result
}


#' Check if object is a "nice" RDKit smile
#'
#' @param smile A smile
is.nice.smile <-function(smile){
  sapply( smile, function(smi){
  # check if smile
  result <- F    
  tryCatch({
    mol <- smiles2mol(smi)
    kekulize(mol)
    s <- smiles2mol(mol2smiles(mol))      
    result <- T
  }, error = function(e) {})    
  result
  })
}
 
#' unfactor
#'
#' @param f factor
#' @param type type
unfactor <- function (f, type = "n") 
{
  if (!is.factor(f)) {
    return(f)
  }
  else {
    f <- levels(f)[f]
    if (type == "n") {
      return(as.numeric(f))
    }
    if (type == "c") {
      return(as.character(f))
    }
  }
}

#' unfactor shortcut for character factors
#'
#' @param f factor
unfactorc <- function (f) 
{
  unfactor(f, "c")
}

#' unfactor shortcut for numeric factors
#'
#' @param f factor
unfactorn <- function (f) 
{
  unfactor(f, "n")
}