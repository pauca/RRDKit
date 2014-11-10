RRDKit
======

A pragmatic interface to RDKit (C++ API) in R.

RRDKit package provides a pragmatic interface to some of the RDKit functions in R. It is intended to work smoothly with R. RRDKit aims to be a tool to perform
basic operations from RDKit. If you are looking for a more richer tool check RDKit web site.


## Prerequisites

* R >= 3.1.0 

* An RDKit installation and $RDBASE configured. ( You can follow the
  instructions in [http://www.rdkit.org/docs/Install.html](http://www.rdkit.org/docs/Install.html)). Note: no Python bindings are needed. Check "Building the RDKit" section.
  
## Installation

* Download RRDKit_X.X.tar.gz
* Run R CMD INSTALL RRDKit_X.X.tar.gz.
  
  
## Functions

### Read and Write
read.sdf( file )  
write.sdf( file , mols )  

smiles2mol( smile )  
smarts2mol( smart )  

mol2smiles( mol )  

### Changing Properties
molGetProps( mol )  
molsGetProps( mols )  

molSetProp( mol  ,key , value)  
molsSetProp ( mols  ,key , values )  

### Molecule viewers 

Next functions open a browser with a 2D representation of the molecules.

showmol(mol)  
showmols(mol)  
showmols.grid(mols)  
mol2svg(mol)   

### Descriptors
mol2maccs(mol)  
mol2morgan(mol)  
mol2mw(mol)  
mol2TPSA(mol)  
mol2LogP(mol)  
mol2murcko(mol)  
computeGasteigerCharges(mol)  

### Others
SubstructMatch(  mol , query )  

