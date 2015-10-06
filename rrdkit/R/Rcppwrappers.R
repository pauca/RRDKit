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
# along with RRDKit.  If not, see <http://www.gnu.org/licenses/>.

library(Rcpp) 

#' Convert smiles to molecules
#'
#' @param smi a smile string or vector of strings
#' @param sanitize toggles H removal and sanitization of the molecule
#' @return List of molecules
#' @examples
#' mol <- smiles2mol("c1ccccc1")
#' mols <- smiles2mol(c("c1ccccc1","CC(=O)OC1=CC=CC=C1C(O)=O"))
smiles2mol <- function( smi  ,  sanitize = TRUE ){
  sapply(smi,p_smile2mol,sanitize )
}

# showmol<-function( ptr , open = T ){
#   svg <- mol2svg(ptr)
#   fileName <- tempfile(pattern = "", tmpdir = tempdir(), fileext = ".svg")
#   file <- file(fileName,"w")
#   cat(svg,file=file)
#   close(file)
#   if(open){ browseURL(paste("file:///",fileName ,sep=""))}
#   return(fileName)
# }

#' Show molecules as 2D in browser
#'
#' @param mols a list of molecules
#' @param open if TRUE show output in browser
#' @param names vector with mols ids
#' @param svg.size  size of pictures
#' @param grid show molecules as grid
#' @return Path to temporary generated html file
#' @examples
#'
#' mols <- smiles2mol(c("CC(=O)NC1=CC=C(O)C=C1","CC(=O)OC1=CC=CC=C1C(O)=O"))
#' names <- c("Paracetamol", "Aspirin")
#' 
#' # showMols(mols, names=names)
showMols<-function( mols , open = T, names = "",  svg.size=200 ,
                    grid = FALSE){
  if(!is.list(mols)){
    mols <- list(mols)
  }

  if(length(names)==1){
    names <- 1:length(mols)
  }  
  
  svgs <- mol2svg( mols )
  svgs <- sapply(svgs, cleanSVG , svg.size,svg.size)
  
  head <- paste( ' <!DOCTYPE html><html>
  <head>
    <style>th {  text-align:left;    } ',
  ifelse( grid ,"
.molbox{
    display: inline-block;
  }",""),'
  </style>
  </head>  
  <body><table>')
  tail <- '</table></body></html>'

  fileName <- tempfile(pattern = "", tmpdir = tempdir(), fileext = ".html")
  file <- file(fileName,"w")
  writeLines(head,file)
  d <- sapply(1:length(svgs),function(i)writeLines(
    paste(
      '<tr class="molbox"><td>',names[i],'</td> <td class="molimg" >',svgs[i],'</td></tr>',sep=""),
      file))
  writeLines(tail,file)
  close(file)
  
 if(open){ browseURL(paste("file:///",fileName ,sep=""))}
  return(fileName)
}

#' Show molecules as 2D in browser
#'
#' @param mols a list of molecules
#' @param open if TRUE show output in browser
#' @param group groups of molecules
#' @param id property id that should be used for naming molecules
#' @param svg.size  size of pictures
#' @return Path to temporary generated html file
#' @examples
#'
#' mols <- smiles2mol(c("CC(=O)NC1=CC=C(O)C=C1","CC(=O)OC1=CC=CC=C1C(O)=O"))
#' names <- c("Paracetamol", "Aspirin")
#' # showMolsGrid(mols)
showMolsGrid<-function( mols , group = 1, id="", open = T ,  svg.size=200){
  if(length(mols)!=length(group)){
    group <- rep(1,length(mols))
  }
  svgs <-  mol2svg ( mols)
  svgs <- sapply(svgs, cleanSVG , svg.size,svg.size)
  if( id == ""){
    ids <- paste( "MOL Nr:", 1:length(mols)," .", sep=" ")
  }else{
    ids <-  sapply( unlist(mols),function(p) RRDKit::molGetProps(p)[[id]])
  }
  head <- ' <!DOCTYPE html><html>
  <head>
  <style>
  .molimg {    
  }

  .molgrid { 
  }

  .molbox{
    display: inline-block;
  }

  .cluster-group{
    border-style: solid;
    border-width: 1px;
  }
  </style>
  </head>
  <body>'
  tail <- '</body></html>'
  fileName <- tempfile(pattern = "", tmpdir = tempdir(), fileext = ".html")
  file <- file(fileName,"w")
  writeLines(head,file)
  
  sapply( sort(unique(group)), function(c){
    svgss<-svgs[group==c]
    idss <- ids[group==c]
    writeLines(paste('<div class="cluster-group"><div class="molgrid"><h1>Group: ',c,'</h1>' ,sep="") ,file)
    d <- sapply(1:length(svgss),function(i)writeLines(paste(
      '<div class="molbox"><h6>',idss[i],'</h6> <div class="molimg" >',svgss[i],'</div></div>',sep=""),
      file))
    writeLines('</div></div><br>',file)
  })
  writeLines(tail,file)
  close(file)
  if(open){ browseURL(paste("file:///",fileName ,sep=""))}
  return(fileName)
}

#' Map molecule properties to a list
#'
#' @param m A molecule
#' @return A list of properties
molGetProps <- function( m ){
  if(is.list(m)){ stop("Invalid input (needs a single molecule)")}
  l <- list()
  for( p in p_molGetPropList(m)){
     l[[p]]<- tryCatch( p_molGetProp(m,p), error=function(e){return (NA)})
  }
  return(l)
}

#' Get Properties of molecules
#'
#' @param mols a list of molecules
#' @return data.frame with properties
molsGetProps <- function( mols ){
  if(!is.list(mols)){ 
    mols <- list(mols)
  }
  
  lprops<- lapply(mols,function( m ){
    if(is.list(m)){ stop("Invalid input (needs a single molecule)")}
    l <- list()
    for( p in p_molGetPropList(m)){
      l[[p]]<- tryCatch( p_molGetProp(m,p), error=function(e){return (NA)})
    }
    return(l)
  })
  n <- sort(unique(unlist(lapply(lprops,function(l)names(l)))))
  m <- lapply(lprops,function(x) { u <- unlist(x)[n];names(u)<-n;return(u)})
  as.data.frame(do.call(rbind,m))
}

# molSetProp <- function( m  ,key , value){
#   p_molSetProp(m, as.character(value), key)
# }

#' Set Properties to Molecules
#'
#' @param ms A molecule or list of molecules
#' @param key A key to be added or modified
#' @param v A value or vector 
#' @param df A dataframe with columns to be included in the sdf
molsSetProps <- function( ms, key=NULL, v=NULL ,df=NULL  ){
  if(!is.list(ms)){ 
    ms <- list(ms)
  }
  
  if(!is.null(df)){
    for(col in colnames(df)){
      molsSetProps( ms, key=col, v=df[,col])
    }
    
  }else{
    for(i in 1:length(ms)){
      #molSetProp(ms[[i]],key,v[i])
      p_molSetProp(ms[[i]], as.character(v[i]), key)
    }
  }
}

#' Read a SDF file and create molecules
#'
#' @param file file name
#' @param sanitize  - if true sanitize the molecule before returning it
#' @param removeHs  - if true remove Hs from the molecule before returning it (triggers sanitization)
#' @param strictParsing  - if not set, the parser is more lax about correctness of the contents.
#' @return A list of molecules
read.sdf<-function(file,sanitize=T, removeHs=T,  strictParsing =T){
  ms <- p_molSupplier(file,sanitize,  removeHs, strictParsing )
  obj<- unlist(p_molSupplierApply(ms,function(m)return(m)))
  obj  
} 

#' Write a SDF file
#'
#' @param file File name
#' @param mols A list of molecules
#' @param setForceV3000 force V3000 
#' @examples
#'
#' mols <- smiles2mol(c("CC(=O)NC1=CC=C(O)C=C1","CC(=O)OC1=CC=CC=C1C(O)=O"))
#' write.sdf("mols.sdf",mols)
#' 
#' compute2D(mols)
#' write.sdf("mols.sdf",mols)
write.sdf <- function(file,mols,setForceV3000=F){
  p_writeSdf(file,mols,setForceV3000)
}

#' Computes 2D coordinates for molecules
#'
#' @param mols A list of molecules
compute2D <- function(mols){
  if(!is.list(mols)){ 
    mols <- list(mols)
  }
  
  for( i in 1:length(mols)){
    p_molCompute2DCoords( mols[[i]])
  }
}

#' Reads a smi file into a data.frame
#'
#' @param file File name
#' @param colnames Vector with colnames
#' @return A data frame
read.smi <- function( file, colnames="" ){
  df <- read.table(file,header=F, comment.char="",stringsAsFactors=F,
             quote = "\"",
             sep="\t" )
  if(length(colnames)>1)
  colnames(df)<-colnames
  df
}



#' Map molecules to SVG
#'
#' @param mols A list of molecules
#' @return A list of svg
mol2svg <- function( mols ){
  p_vectorize( mols,  p_mol2svg)  
}

#' Calculate molecular weight
#'
#' @param mols A list of molecules
#' @return A list of molecular weights
mol2mw <- function( mols ){
  p_vectorize( mols,  p_mol2mw)  
}

#' Get Smiles of molecules
#'
#' @param mols A list of molecules
#' @return A list of molecular weights
mol2smiles <- function( mols ){
  p_vectorize( mols,  p_mol2smiles)  
}

#' Calculate TPSA
#'
#' @param mols A list of molecules
#' @return TPSA
mol2TPSA <- function( mols ){
  p_vectorize( mols,  p_mol2TPSA)  
} 

#' Calculate LogP
#'
#' @param mols A list of molecules
#' @return LogP
mol2LogP <- function( mols ){
  p_vectorize( mols,  p_mol2LogP)  
} 

#' Calculate CalcNumRotatableBonds
#'
#' @param mols A list of molecules
#' @return  NumRotatableBonds
mol2NumRotableBonds <- function( mols ){
  p_vectorize( mols,  p_mol2NumRotatableBonds)  
} 

#' Calculate CalcNum Hidrogen bond acceptors
#'
#' @param mols A list of molecules
#' @return  NumRotatableBonds
mol2NumHBA <- function( mols ){
  p_vectorize( mols,  p_mol2NumHBA)  
} 


#' Calculate CalcNum Hidrogen bond acceptors
#'
#' @param mols A list of molecules
#' @return  NumRotatableBonds
mol2NumHBD <- function( mols ){
  p_vectorize( mols,  p_mol2NumHBD)  
} 




#' Calculate murcko scaffold
#'
#' @param mols A list of molecules
#' @return murcko scaffolds
mol2murcko <- function( mols ){
  p_vectorize( mols,  p_mol2murcko)  
} 

#' computeGasteigerCharges
#'
#' @param mols A list of molecules
#' @return computeGasteigerCharges
computeGasteigerCharges <- function( mols ){
  p_vectorize( mols,  p_computeGasteigerCharges)  
} 

#' Calculate kekulize
#'
#' @param mols A list of molecules
#' @return kekulize
kekulize <- function( mols ){
  p_vectorize( mols,  p_kekulize)  
} 

#' map a smarts to a molecule
#'
#' @param smarts a smarts string
#' @return a molecule
smarts2mol <- function(smarts){
  p_vectorize( smarts,  p_smarts2mol)  
}
