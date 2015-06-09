/*    This file is part of RRDKit.

    RRDKit is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Foobar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <RcppCommon.h>
#include <Rcpp.h>
#include <Rcpp/vector/Vector.h>
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/MolDrawing/MolDrawing.h>
#include <GraphMol/MolDrawing/DrawingToSVG.h> 
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/PartialCharges/GasteigerCharges.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/Dict.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/Fingerprints/MACCS.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <DataStructs/ExplicitBitVect.h> 

#include <vector>
#include <algorithm>
#include <iostream> 
 
 
using namespace Rcpp;
using namespace RDKit;
#include <R.h>
#include <Rdefines.h>



 
// [[Rcpp::export]] 
SEXP   p_smile2mol(std::string smi  ,  bool sanitize = true ){
  int sanitize_ = 0;
  if( sanitize ){
    sanitize_ = 1;
  }

  RWMol *m = SmilesToMol( smi ,0 , sanitize_  );
  if(m){
    SEXP ptr =  R_MakeExternalPtr(m, R_NilValue ,R_NilValue);
    return (ptr);
  }else{
    throw std::invalid_argument("Smiles Parse Error" );
  }
}

//' map a smarts to a molecule
//'
//' @param smarts a smarts string
//' @return a molecule
// [[Rcpp::export]] 
SEXP   smarts2mol(std::string smarts ){
  RWMol *m = SmartsToMol( smarts );
  if(m){
    SEXP ptr =  R_MakeExternalPtr(m, R_NilValue ,R_NilValue);
    return (ptr);
  }else{
    //::Rf_error(""); 
    throw std::invalid_argument("Smarts Parse Error");
  }
}

// [[Rcpp::export]] 
bool p_molIsNil(SEXP xp){
  RWMol *mol=(RWMol*)(R_ExternalPtrAddr(xp));
  return (mol == NULL);
}

RWMol*  p_getMol(SEXP xp){
  // get mol in a saver way...
  if( TYPEOF(xp) != EXTPTRSXP){
    ::Rf_error("Not a pointer");    
  }
  
  RWMol *mol=(RWMol*)(R_ExternalPtrAddr(xp));
    
  bool rwmol = false;
  if(dynamic_cast<RWMol*>(mol)){
    rwmol = true;
  }
  
  if(rwmol){
    //any is true
  }else{
    ::Rf_error("not molecule"); 
  }

  if( !mol ){ 
    ::Rf_error("null molecule"); 
  }

  return mol;
}



SDMolSupplier*  p_getMolSDMolSupplier(SEXP xp){
  if( TYPEOF(xp) != EXTPTRSXP){
    ::Rf_error("Not a pointer");    
  }

  SDMolSupplier *mol=(SDMolSupplier*)(R_ExternalPtrAddr(xp));

  if(!dynamic_cast<SDMolSupplier*>(mol)){
    ::Rf_error("not SDMolSupplier"); 
  }

  if( !mol ){ 
    ::Rf_error("null SDMolSupplier"); 
  }
  return mol;
}

// [[Rcpp::export]]
std::string p_mol2smiles( SEXP xp ){
    RWMol *mol = p_getMol(xp);  
    return MolToSmiles(*mol);
}

// [[Rcpp::export]]
std::string p_mol2svg( SEXP xp ){
    ROMol *mol = (ROMol*)(R_ExternalPtrAddr(xp));
    RDDepict::compute2DCoords(*mol);
    std::vector<int> drawing=RDKit::Drawing::MolToDrawing(*mol);
    std::string svg=RDKit::Drawing::DrawingToSVG(drawing);
    return svg;
}



//' map a molecule to a maccs fingerprints
//'
//' @param xp a molecule
//' @return a vector
// [[Rcpp::export]]
IntegerVector mol2maccs(  SEXP xp){
    int i;
    RWMol *mol = p_getMol(xp);  
    
    if( mol ){
     ExplicitBitVect *bv = MACCSFingerprints::getFingerprintAsBitVect ( *mol  ); 
      
     if(bv){
      IntegerVector v = IntegerVector( bv->getNumBits ());
      for( i = 0;  i< bv->getNumBits ();i++ ){
        v(i) = bv->getBit(i)?1:0;
      }      
      return v;
      }
    }
    return IntegerVector(0); 
}

//// [[Rcpp::export]]
//IntegerVector mol2topologicalFingerprints(  SEXP xp){
//    int i;
//    RWMol *mol = p_getMol(xp);  
//    
//    if( mol ){
//     ExplicitBitVect *bv = RDKit::RDKFingerprintMol ( *mol  ); 
//      
//     if(bv){
//      IntegerVector v = IntegerVector( bv->getNumBits ());
//      for( i = 0;  i< bv->getNumBits ();i++ ){
//        v(i) = bv->getBit(i)?1:0;
//      }      
//      return v;
//      }
//    }
//    return IntegerVector(0); 
//}


//' map a molecule to a morgan fingerprints (atom environement)
//'
//' @param xp a molecule
//' @param radius radius of scaffold
//' @param nBits final length
//' @param useFeatures useFeatures
//' @return a vector
// [[Rcpp::export]]
IntegerVector mol2morgan(  SEXP xp , unsigned int radius=2,unsigned int nBits=2048, 
    bool useFeatures=false){
  // radius = 2 => ecfp4 2^2
  // radius = 2 & useFeatures=True => fcfp4
    int i;
    RWMol *mol = p_getMol(xp);  
    
    if( mol ){
     ExplicitBitVect *bv = MorganFingerprints::getFingerprintAsBitVect ( *mol  , radius,nBits,
        0,0,false,useFeatures,false,0); 

     if(bv){
      IntegerVector v = IntegerVector( bv->getNumBits ());
      for( i = 0;  i< bv->getNumBits ();i++ ){
        v(i) = bv->getBit(i)?1:0;
      }   
      delete bv;
      return v;
      }
    }
    return IntegerVector(0); 
}


//// [[Rcpp::export]]
//void p_molSupplierWrite(  std::string file, SEXP xp ){
//    
//    SDMolSupplier * molsupp = p_getMolSDMolSupplier(xp);  
//    SDWriter::SDWriter  sdw  = SDWriter::SDWriter(file);
//    
//    while(! molsupp->atEnd()){
//      ROMol  * mol = molsupp->next();
//      sdw.write( *mol , RDKit::defaultConfId);
//    }
//    sdw.close();     
//} 

// [[Rcpp::export]]
void p_writeSdf(  std::string file,  SEXP  pv, bool setForceV3000 = false){
    Rcpp::List v(pv);
    int i =0;
    SDWriter::SDWriter  sdw  = SDWriter::SDWriter(file); 
    sdw.setForceV3000(setForceV3000);
    for (unsigned i = 0; i < v.size (); ++ i){
        RWMol *mol=(RWMol*)(R_ExternalPtrAddr(v[i]));
         sdw.write(  *mol , RDKit::defaultConfId);
         //sdw.write(  *mol , 0);
    }
    sdw.close();
} 
 

// [[Rcpp::export]] 
SEXP  p_molSupplier(std::string file,bool sanitize=true, 
              bool removeHs=true, bool strictParsing =true){
  SDMolSupplier * sp =  new SDMolSupplier( file, sanitize, removeHs, strictParsing );
  if(sp){
    //XPtr<RWMol> ptr(m, true);
    SEXP ptr =  R_MakeExternalPtr(sp,R_NilValue,R_NilValue);
    return (ptr);
  }else{
    ::Rf_error("Error"); 
  }
}

// [[Rcpp::export]] 
bool  p_molSupplier_atEnd(SEXP xp ){
  SDMolSupplier * molsupp = p_getMolSDMolSupplier(xp);        
  bool b = molsupp->atEnd(); 
  return b;
}

// [[Rcpp::export]] 
SEXP  p_molSupplier_next(SEXP xp ){
  SDMolSupplier * molsupp = p_getMolSDMolSupplier(xp);  
  bool b = molsupp->atEnd(); 
  if( !b){
    ROMol  * mol = molsupp->next();
    SEXP ptr =  R_MakeExternalPtr(mol,R_NilValue,R_NilValue);
    return (ptr);
  }
   ::Rf_error("Out of bonds"); 
}

// [[Rcpp::export]] 
void  p_molSupplier_reset(SEXP xp ){
  SDMolSupplier * molsupp = p_getMolSDMolSupplier(xp);  
  molsupp->reset();   
}

 
// [[Rcpp::export]]
std::vector<std::string> p_molGetPropList(  SEXP xp ){    
    ROMol *mol=(ROMol*)(R_ExternalPtrAddr(xp));
    std::vector<std::string> prop = mol->getPropList(false,false);
    return prop;
} 
 
// [[Rcpp::export]]
std::string p_molGetProp(  SEXP xp , SEXP property ){
    std::string property_ = as<std::string>(property);
    ROMol *mol=(ROMol*)(R_ExternalPtrAddr(xp));
    std::string  p_s = mol->getProp<std::string> ( property_ ) ;
    return p_s;
} 

// [[Rcpp::export]]
void p_molSetProp(  SEXP xp , SEXP property, SEXP key ){
    std::string property_ = as<std::string>(property);
    std::string key_ = as<std::string>(key);
    ROMol *mol=(ROMol*)(R_ExternalPtrAddr(xp));
    mol->setProp<std::string> ( key_,property_ ) ;
} 

//' split molecule to brics fragment
//'
//' @param xp a molecule
//' @return a BRICS fragment
// [[Rcpp::export]]
SEXP fragmentOnBRICSBonds(  SEXP xp ){    
    ROMol *mol=(ROMol*)(R_ExternalPtrAddr(xp));
    //std::vector<std::string> prop = mol->getPropList(false,false);
    ROMol* molf = MolFragmenter::fragmentOnBRICSBonds(*mol);
    SEXP ptr =  R_MakeExternalPtr(molf, R_NilValue ,R_NilValue);
    return (ptr);
} 

//' SubstructMatch
//'
//' @param xp_mol a molecule
//' @param xp_query a template molecule (from a smart)
//' @return bool
// [[Rcpp::export]]
bool SubstructMatch(  SEXP xp_mol , SEXP xp_query  ){    
    MatchVectType match;
    ROMol *mol  =(ROMol*)(R_ExternalPtrAddr(xp_mol));
    ROMol *query=(ROMol*)(R_ExternalPtrAddr(xp_query));
    return  SubstructMatch(*mol,*query,match);    
} 


// [[Rcpp::export]]
double p_mol2mw(  SEXP xp ){    
    RWMol * mol =  p_getMol(xp);        
    return  Descriptors::calcExactMW( *mol  );
}

// [[Rcpp::export]]
double p_mol2TPSA(  SEXP xp ){    
    RWMol * mol =  p_getMol(xp);        
    return  Descriptors::calcTPSA( *mol  );
}


// [[Rcpp::export]]
double p_mol2LogP(  SEXP xp ){    
    RWMol * mol =  p_getMol(xp);    
    double logp=0 ;
    double mw =0  ;
    Descriptors::calcCrippenDescriptors( *mol , logp , mw);
    return logp;
}

 
// [[Rcpp::export]]
std::string p_mol2murcko( SEXP xp ){
    RWMol *mol  =  p_getMol(xp);  
    ROMol *molm = MurckoDecompose ( *mol  ); 
    std::string smi = MolToSmiles(*molm);
    free(molm);
    return (smi);
}

//' Compute GasteigerCharges
//'
//' @param xp a molecule
//' @return the GasteigerCharges
// [[Rcpp::export]]
std::vector< double > p_computeGasteigerCharges( SEXP xp ){
  //computeGasteigerCharges(smiles2mol("c1cccc2c1CCCC2"))
    ROMol *mol  =  p_getMol(xp); 
    std::vector< double >  charges( mol->getNumAtoms());
    RDKit::computeGasteigerCharges (  *mol,  charges, 12, false);
    return  charges;    
}

//' Kekulize molecule
//'
//' @param xp a molecule
//' @param  markAtomsBonds  markAtomsBonds
//' @param  maxBackTracks  maxBackTracks
// [[Rcpp::export]]
void p_kekulize( SEXP xp,bool    markAtomsBonds = true,
		unsigned int  	maxBackTracks = 100  ){
    RWMol *mol  =  p_getMol(xp);  
    MolOps::Kekulize(*mol,markAtomsBonds,maxBackTracks);    
}

// [[Rcpp::export]]
void p_molCompute2DCoords( SEXP xp ){
    ROMol *mol = (ROMol*)(R_ExternalPtrAddr(xp));
    RDDepict::compute2DCoords(*mol);
}
