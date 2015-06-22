/*# This file is part of RRDKitInchi.
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
# along with RRDKitInchi.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <Rcpp.h>
using namespace Rcpp;
#include "inchi.h"
 
using namespace RDKit;

// [[Rcpp::export]]
std::string p_mol2Inchi( SEXP xp ){
    RWMol *mol=(RWMol*)(R_ExternalPtrAddr(xp));
    ExtraInchiReturnValues tmp;
    std::string inchi = MolToInchi(*mol,tmp);
    return (inchi); 
}

// [[Rcpp::export]]
std::string p_Inchi2InchiKey(  SEXP x ){
  std::string inchi = as<std::string>(x);
  return( InchiToInchiKey(  inchi));
}


//' Map Inchi to molecule
//'
//' @param x a  molecule
//' @param sanitize bool
//' @param removeHs bool
//' @return a molecule
// [[Rcpp::export]]
SEXP  Inchi2Mol(  SEXP x  , bool sanitize = false, bool removeHs = false){
  ExtraInchiReturnValues rv ;
  std::string inchi = as<std::string>(x);
  RWMol *m =  InchiToMol(  inchi, rv , sanitize, removeHs);
  if(m){
    SEXP ptr =  R_MakeExternalPtr(m, R_NilValue ,R_NilValue);
   return (ptr);
  }else{
    throw std::invalid_argument("Inchi Parse Error");
  }
}
 
