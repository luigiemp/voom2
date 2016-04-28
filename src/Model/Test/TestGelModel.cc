// Test for Model
  
#include "GelModel.h"
#include "GelMesh.h"
#include "EigenResult.h"
#include "Spring.h"
#include "AngleSpring.h"

using namespace voom;

int main(int argc, char** argv) {

  
  
  //Test Mechanics model
  {
    cout << " ------------------------------- " << endl;
    cout << " TEST OF GEL MODEL " << endl << endl;
            
    GelMesh TestGelMesh("../../Mesh/Test/Filament2.node", "../../Mesh/Test/Filament2.ele");
    
    // Initialize Model
    uint NodeDoF = TestGelMesh.getDimension();
    uint NumFil = TestGelMesh.getNumberOfFilaments();
    uint NumNod = TestGelMesh.getNumberOfNodes();
    uint PbDoF = NumNod*NodeDoF;

    cout << "Number of elements is " << NumFil << endl;
    cout << "Number of degrees of freedom is " << PbDoF << endl;

    // Initlize filament materials
    vector<FilamentMaterial * > springs;
    springs.reserve(NumFil);
    
    vector<FilamentMaterial * > angleSprings;
    angleSprings.reserve(NumFil);
    
    Real k = 1.0;
    Real kappa = 1.0;
    for (int nFil = 0; nFil < NumFil; nFil++) {        
      Spring* Spr = new Spring(nFil,k); 
      springs.push_back(Spr);
      AngleSpring* Ang = new AngleSpring(nFil,kappa);
      angleSprings.push_back(Ang);
    }
    
    // Construct Gel
    GelModel myModel(&TestGelMesh, springs, angleSprings, NodeDoF,0,1);
    
    EigenResult R(PbDoF, 1);
   
    R.setRequest(3);
      
    // Isotropic deformation
    myModel.initializeField(1.1);    
       
    // Compute Energy & residual
    myModel.compute(R);


    cout << "Total energy is " << R.getEnergy() << endl;
    cout << endl << " END OF TEST OF GEL MODEL " << endl;
    cout << " ------------------------------ " << endl << endl;
  }


}
