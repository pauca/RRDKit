RRDKit
======

A pragmatic interface to RDKit (C++ API) in R.

RRDKit package provides a pragmatic interface to some of the RDKit functions in R. It is intended to work smoothly with R. RRDKit aims to be a tool to perform
basic operations from RDKit. If you are looking for a more richer tool check RDKit web site.


## Prerequisites

* R >= 3.1.0

* R Packages: Rcpp, testthat.

* An RDKit installation and $RDBASE configured. Follow the
  instuctions in [http://www.rdkit.org/docs/Install.html](http://www.rdkit.org/docs/Install.html)). Check "Building the RDKit" section.
  
* $LD_LIBRARY_PATH should contain path folder with libboost_python.so.1.XX.XX and $RDBASE/lib .  
  
## Installation

* Download latest RRDKit_X.X.tar.gz
* Run R CMD INSTALL RRDKit_X.X.tar.gz.
  
## Usage

```
library(RRDKit)  
mols1 <- read.sdf(system.file("data/aspirine.sdf", package="RRDKit"))  
mols2 <- read.sdf(system.file("data/clozapine.sdf", package="RRDKit"))  
mols <- c(mols1,mols2)
sapply( mols, mol2mw )  
showmols(mols)  
```

## Functions

### Read and Write

```
read.sdf( file )  
write.sdf( file , mols )  

smiles2mol( smile )  
smarts2mol( smart )  

mol2smiles( mol )  
```
### Changing Properties
```
molGetProps( mol )  
molsGetProps( mols )  

molSetProp( mol  ,key , value)  
molsSetProp ( mols  ,key , values )  
```
### Molecule viewers 

Next functions open a browser with a 2D representation of the molecules.
```
showmol(mol)  
showmols(mols)  
showmols.grid(mols)  
mol2svg(mol)   
```
### Descriptors
```
mol2maccs(mol)  
mol2morgan(mol)  
mol2mw(mol)  
mol2TPSA(mol)  
mol2LogP(mol)  
mol2murcko(mol)  
computeGasteigerCharges(mol)  
```
### Others
```
SubstructMatch(  mol , query )  
```
