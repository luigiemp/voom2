// Test for Model

#include "GelModel.h"
#include "GelMesh.h"
#include "EigenResult.h"
#include "Spring.h"

using namespace voom;

int main(int argc, char** argv) {

  
  
  //Test Mechanics model
  {
    cout << " ------------------------------- " << endl;
    cout << " TEST OF GEL MODEL " << endl << endl;
    
    
    GelMesh TestGelMesh("../../Mesh/Test/Filament2.node", "../../Mesh/Test/Filament2.ele");
    
    // Initialize Model
    uint NodeDoF = 3;

    uint NumFil = TestGelMesh.getNumberOfFilaments();
    
    uint PbDoF = NumFil*NodeDoF;

    std::cout << "Number of elements is " << NumFil << std::endl;
    
    vector<FilamentMaterial * > materials;
    materials.reserve(NumFil);
    
    Real k;
    Vector3d d0;
        
    for (int nFil = 0; nFil < NumFil; nFil++) {        
      Spring* Mat = new Spring(nFil,k,d0);
      materials.push_back(Mat);
    }
    
    GelModel myModel(&TestGelMesh, materials, NodeDoF,0,1);
    
    EigenResult R(PbDoF, 1);
   
    R.setRequest(7);
     
    myModel.initializeField(1.1);    
    
    myModel.compute(R);

    
    cout << endl << " END OF TEST OF GEL MODEL " << endl;
    cout << " ------------------------------ " << endl << endl;
  }


}
