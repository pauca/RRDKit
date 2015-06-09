RRDKit
======

A pragmatic interface to RDKit (C++ API) in R.

RRDKit package provides a pragmatic interface to some of the RDKit functions in R. It is intended to work smoothly with R. RRDKit aims to be a tool to perform
basic operations from RDKit. If you are looking for a more richer tool check RDKit web site.


## Prerequisites

* R >= 3.1.0

* R Packages: Rcpp, testthat.

* A RDKit installation. Preferably use RDKit version RDKit_2014_09_1. Follow the
  instuctions in [http://www.rdkit.org/docs/Install.html](http://www.rdkit.org/docs/Install.html)). Check "Building the RDKit" section. 
  
* Note that Python wrappers can be disabled and INCHI support enabled (check RDKitInchi):
```
cmake    -D RDK_BUILD_PYTHON_WRAPPERS= -D RDK_BUILD_INCHI_SUPPORT=ON ..
                                          
```

* RDBASE (the root directory of the RDKit distribution  e.g. ~/RDKit  ) configured. 
  
* LD_LIBRARY_PATH must include $RDBASE/lib.
  
## Installation

* Download latest RRDKit and Install:
```
library(devtools)
install_github("pauca/RRDKit/RRDKit")
```

  
## Usage

```
library(RRDKit)  
mols1 <- read.sdf(system.file("extdata/aspirine.sdf", package="RRDKit"))  
mols2 <- read.sdf(system.file("extdata/clozapine.sdf", package="RRDKit"))  
mols <- c(mols1,mols2)
mol2mw(mols)
showMols(mols)  
```

## Functions

### Read and Write

```
read.sdf( file )  
write.sdf( file , mols )  

#Read a smi file as data frame:
read.smi(file)


smiles2mol( smile )  
smarts2mol( smart )  

mol2smiles( mol )  
```
### Changing Properties
```
molsGetProps( mols )  

molsSetProp ( mols  ,key , values )  
```
### Molecule viewers 

Next functions open a browser with a 2D representation of the molecules.
```

showMols(mols)  
showMolsGrid(mols)  
mol2svg(mol)  
molCompute2DCoords
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
