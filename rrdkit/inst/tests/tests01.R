library(Rcpp)
library(testthat)
library(rrdkit)
context("smiles conversions and operations")

test_that("Smiles to mol", {
  expect_that( mol2smiles(smiles2mol("C"))[[1]] , equals("C"))
  expect_that( mol2smiles(smiles2mol("c1ccccc1"))[[1]] , equals("c1ccccc1"))
  expect_that( mol2smiles(smiles2mol("Cl.c1ccccc1"))[[1]] , equals("Cl.c1ccccc1"))
  expect_that( smiles2mol("C1") , throws_error())
  expect_that( smiles2mol("XYZ") , throws_error())
})


test_that("Smiles to mol and operations", {
  expect_that( mol2mw(smiles2mol("C"))[[1]] , equals(16.0313))
  expect_that( mol2mw(smiles2mol("c1ccc1"))[[1]] , equals(52.031300128))
  
})

context("SVG functions")

test_that("SVG functions", {
  expect_that( mol2svg(smiles2mol("C"))[[1]] , equals("<?xml version='1.0' encoding='iso-8859-1'?>\n<svg:svg version='1.1' baseProfile='full'\n              xmlns:svg='http://www.w3.org/2000/svg'\n                        xmlns:xlink='http://www.w3.org/1999/xlink'\n                  xml:space='preserve'\nwidth='40px' height='50px' >\n<svg:rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='40' height='50' x='0' y='0'> </svg:rect>\n<svg:g transform='translate(2,2.5) scale(.85,.85)'><svg:g transform='translate(20,30)'><svg:rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='150' height='50' x='-75' y='-25'> </svg:rect>\n<svg:text style='font-size:50px;font-style:normal;font-weight:normal;line-height:125%;letter-spacing:0px;word-spacing:0px;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:middle;fill:#000000' y='18.75'><svg:tspan>CH4</svg:tspan></svg:text></svg:g>\n</svg:g></svg:svg>"))  
})

context("Inchi functions")
test_that("Inchi to mol and operations", {
  mol1 <- read.sdf(system.file("extdata/aspirine.sdf", package="rrdkit"))
  mol2 <- read.sdf(system.file("extdata/clozapine.sdf", package="rrdkit"))
  expect_that(  mol2Inchi(mol1)  , equals("InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"))
  expect_that(  mol2Inchi(c(mol1,mol2))[[1]]  , equals("InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"))
  expect_that(  mol2InchiKey(c(mol1,mol2))[[1]], equals("BSYNRYMUTXBXSQ-UHFFFAOYSA-N"))
  
  inchikeys <- mol2InchiKey(c(mol1,mol2))
  expect_that( compareInchiKeys( inchikeys[1] ,inchikeys[2]), equals(FALSE))
  expect_that( compareInchiKeys( inchikeys[1] ,inchikeys[2],check.protonation = TRUE), equals(FALSE))
  expect_that( compareInchiKeys( inchikeys[1] ,inchikeys[1],check.protonation = TRUE), equals(TRUE))
})