// Test for PreferredNormal Body

#include "PreferredNormalBody.h"
#include "FEMesh.h"
#include "EigenResult.h"

using namespace voom;

int main(int argc, char** argv) {

  //Test PreferredNormal body
  {
    srand(time(NULL));
    cout << " ----------------------------------- " << endl;
    cout << " TEST OF PREFERREDNORMAL BODY " << endl << endl;
  
    State myState;
    bool CheckOverlap = false;
    int dofPerNode = 3;
    // FEMesh myFEmesh("../../Mesh/Test/Cube.node", "../../Mesh/Test/SurfCube.ele", &myState, dofPerNode, CheckOverlap);
    FEMesh myFEmesh("../../Mesh/Test/CubeQuad.node", "../../Mesh/Test/SurfCubeQuad.ele", &myState, dofPerNode, CheckOverlap);   
    uint NumQP = 1;

    
    Real Krot = double(rand())/RAND_MAX;
    PreferredNormalBody myBody(&myFEmesh, &myState, Krot);
    int PbDoF = myState.getDOFcount();
    cout << "PbDoF = " << PbDoF << endl;
    myState.printPhi();
    myBody.writeOutputVTK("PreferredNormalBody", 0);



    EigenResult myResult(PbDoF, 0);
    myResult.initializeResults(PbDoF*4);



    Real perturbationFactor = 0.1;
    int myRequest = FORCE | STIFFNESS; // Check both Forces and Stiffness
    Real myH = 1e-6;
    Real myTol = 1e-7;
    myResult.setRequest(ENERGY);
    myResult.resetResults(ENERGY);
    myBody.compute(&myResult);
    cout << endl << "Body energy is  = " << myResult.getEnergy() << endl;
    // Change field and recompute energy
    for (int i=0; i<PbDoF; i++) {
      myState.linearizedUpdate(i, double(0.1*rand())/RAND_MAX);
    }
    myResult.resetResults(ENERGY);
    myBody.compute(&myResult);
    cout << endl << "Body energy is  = " << myResult.getEnergy() << endl;

    myBody.checkConsistency(&myResult, perturbationFactor, myRequest, myH, myTol);
   
    

    cout << endl << " END OF TEST OF PREFERREDNORMAL BODY " << endl;
    cout <<         " ---------------------------------- " << endl << endl;
  }


}
