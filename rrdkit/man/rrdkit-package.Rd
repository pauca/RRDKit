\name{rrdkit}
\alias{rrdkit}
\alias{rrdkit}
\docType{package}
\title{
A pragmatic interface to RDKit in R
}
\description{
This package provides a pragmatic interface to some of the RDKit functions in R. It is intended to work smoothly with R workflows.

rrdkit uses the RDKit installation in $RDBASE.

rrdkit provides some basic functions to visualize the molecules.
}
\details{
\tabular{ll}{
Package: \tab rrdkit\cr
Type: \tab Package\cr
Version: \tab 0.11\cr
Date: \tab 2015-06-09\cr
License: \tab GPL-3\cr
}


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
library(rrdkit)  
mols1 <- read.sdf(system.file("extdata/aspirine.sdf", package="rrdkit"))  
mols2 <- read.sdf(system.file("extdata/clozapine.sdf", package="rrdkit"))  
mols <- c(mols1,mols2)
mol2mw(mols)
\dontrun{
showMols(mols)  
}
}
