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

  
molSupplierApply <- function( molSupplier, fun, ...){
  l1 <- list()
  while(!molSupplier_atEnd(molSupplier)){
    mol <- molSupplier_next(molSupplier)
    res <- tryCatch( fun(mol,...), error= function(e) NA )
    l1[[length(l1)+1]]<-res
  }
  molSupplier_reset(molSupplier)
  return(l1)
}

getBRICSFragments<-function( mol ){
  f <- fragmentOnBRICSBonds(mol)
  ff <- unique(unlist(strsplit(mol2smiles(f),"\\.")))
  return(ff)
}

p_getValidSubrfagemnts <- function(mol,mf){
  indx <- unlist(sapply(mf,function(f){
    SubstructMatch(mol,smarts2mol(f))
  }))
  return(mf[indx])
}

 

