//g++ -o svgcheck.o svgcheck.cpp -ggdb -I$RDBASE/Code  -L$RDBASE/build/lib -lAlignment         -lFileParsers          -lOptimizer  -lCatalogs          -lFingerprints         -lPartialCharges  -lChemicalFeatures  -lForceFieldHelpers    -lRDBoost  -lChemReactions     -lForceField           -lRDGeneral  -lChemTransforms    -lFragCatalog          -lRDGeometryLib  -lDataStructs       -lGraphMol             -lShapeHelpers  -lDepictor          -lhc                   -lSimDivPickers  -lDescriptors       -lMolAlign             -lSLNParse  -lDistGeometry      -lMolCatalog           -lSmilesParse  -lDistGeomHelpers   -lMolChemicalFeatures  -lSubgraphs  -lEigenSolvers      -lMolTransforms        -lSubstructMatch 


#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/MolDrawing/MolDrawing.h>
#include <GraphMol/MolDrawing/DrawingToSVG.h>

using namespace RDKit;
std::string MolToSVG(const ROMol &mol){
std::vector<int> drawing=RDKit::Drawing::MolToDrawing(mol);
std::string svg=RDKit::Drawing::DrawingToSVG(drawing);
return svg;
}

int main(void){
RWMol *mol=SmilesToMol("CCC");
std::string svg=MolToSVG(*mol);
std::cout <<svg<<std::endl;
delete mol;
return(0);
}