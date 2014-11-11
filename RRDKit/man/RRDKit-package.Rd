\name{RRDKit}
\alias{RRDKit}
\alias{RRDKit}
\docType{package}
\title{
A pragmatic interface to RDKit in R
}
\description{
This package provides a pragmatic interface to some of the RDKit functions in R. It is intended to work smoothly with R workflows.

RRDKit uses the RDKit installation in $RDBASE.

RRDKit provides some basic functions to visualize the molecules.
}
\details{
\tabular{ll}{
Package: \tab RRDKit\cr
Type: \tab Package\cr
Version: \tab 0.6\cr
Date: \tab 2014-11-05\cr
License: \tab GPL-3\cr
}

# Available functions
# Read and Write
read.sdf( file )\cr
write.sdf( file , mols )\cr

smiles2mol( smile )\cr
smarts2mol( smart )\cr

mol2smiles( mol )\cr

# Changing Properties
molGetProps( mol )\cr
molsGetProps( mols )\cr

molSetProp( mol  ,key , value)\cr
molsSetProp ( mols  ,key , values )\cr

# molecule viewers
showmol(mol)\cr
showmols(mol)\cr
showmols.grid(mols)\cr

# descriptors
mol2svg(mol) \cr
mol2maccs(mol)\cr
mol2morgan(mol)\cr
mol2mw(mol)\cr
mol2TPSA(mol)\cr
mol2LogP(mol)\cr
mol2murcko(mol)\cr
computeGasteigerCharges(mol)\cr

# others
SubstructMatch(  mol , query )\cr

}
\author{
Pau Carrio

Maintainer: Pau Carrio <paucarrio@gmail.com>
}
\references{
www.rdkit.org
}

\keyword{ rrdkit, molecule }
 
\examples{
%% ~~ simple examples of the most important functions ~~
}
