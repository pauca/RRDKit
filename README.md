rrdkit
======

A pragmatic interface to RDKit (C++ API) from R.

rrdkit package provides a pragmatic interface to some of the RDKit functions in R. It is intended to work smoothly with R. rrdkit aims to be a tool to perform
basic operations from RDKit. If you are looking for a more richer tool check RDKit web site.


## Prerequisites

* R >= 3.2.0

* R Packages: Rcpp, testthat, XML.

* A RDKit installation. Preferably use latest RDKit version ( Latest Tested rdkit-Release_2016_03_1 ). Follow the
  instuctions in [http://www.rdkit.org/docs/Install.html](http://www.rdkit.org/docs/Install.html)). Check "Building the RDKit" section. 
  
* Note that Python wrappers can be disabled (optional) and INCHI support must be enabled:
```
cmake -D RDK_BUILD_PYTHON_WRAPPERS= -D RDK_BUILD_INCHI_SUPPORT=ON ..
                                          
```

* RDBASE (the root directory of the RDKit distribution  e.g. ~/RDKit  ) configured. 
  
* LD_LIBRARY_PATH must include $RDBASE/build/lib.
  
## Installation

* with devtools
```
library(devtools)
install_github("pauca/rrdkit/rrdkit")
```

  
## Examples of usage

```
library(rrdkit)  
mols1 <- read.sdf(system.file("extdata/aspirine.sdf", package="rrdkit"))  
mols2 <- read.sdf(system.file("extdata/clozapine.sdf", package="rrdkit"))  
mols <- c(mols1,mols2)
mol2mw(mols)
showMols(mols)

inchi <- mol2Inchi(mols)  
Inchi2InchiKey(inchi)
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
mol2svg(mols)  
compute2D(mols)

mols <- smiles2mol(c("CC(=O)NC1=CC=C(O)C=C1","CC(=O)OC1=CC=CC=C1C(O)=O"))
df <- data.frame(id=c(1,2),names=c("Paracetamol", "Aspirin"),mw=mol2mw(mols))
showMolsDF(mols,df)

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
