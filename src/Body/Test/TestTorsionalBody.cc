// Test for Torsional Body

#include "TorsionalBody.h"
#include "MechanicsBody.h"
#include "FEMesh.h"
#include "EigenResult.h"
#include "CompNeoHookean.h"

using namespace voom;

int main(int argc, char** argv) {

  //Test Torsional Body
  {
    srand(time(NULL));
    cout << " ----------------------------------- " << endl;
    cout << " TEST OF TORSIONAL BODY " << endl << endl;

    State myState;
    bool CheckOverlap = false;
    int dofPerNode = 3;
    FEMesh myFEmesh("../../../../voom2/src/Mesh/Test/Cube.node", "../../../../voom2/src/Mesh/Test/Cube.ele", &myState, dofPerNode, CheckOverlap);   
    uint NumQP = 1;

    int NumMat = myFEmesh.getNumberOfElements()*NumQP;
    vector<MechanicsMaterial * > materials;
    materials.reserve(NumMat);
    for (int k = 0; k < NumMat; k++) {
      CompNeoHookean* Mat = new CompNeoHookean(k, 1.0+double(rand())/RAND_MAX, 3.0+double(rand())/RAND_MAX);
      materials.push_back(Mat);
    }
    
    MechanicsBody MechBody(&myFEmesh, &myState, materials);



    vector<int > BodyNodes;
    for (int i = 18; i < 27; i++) {
      BodyNodes.push_back(i);
    }
    Real TorK = double(rand())/RAND_MAX;
    Vector3d Center; Center << double(rand())/RAND_MAX, double(rand())/RAND_MAX, double(rand())/RAND_MAX;
    TorsionalBody TorBody(&myFEmesh, &myState, BodyNodes, TorK, Center);
 

 


    // Test consistency
    int PbDoF = myState.getDOFcount();
    EigenResult myResult(PbDoF, 0);
    myResult.initializeResults(PbDoF*4);
 

    Real perturbationFactor = 0.1;
    int myRequest = FORCE | STIFFNESS; // Check both Forces and Stiffness
    Real myH = 1e-6;
    Real myTol = 1e-7;
    myResult.setRequest(ENERGY);
    myResult.resetResults(ENERGY);
    TorBody.compute(&myResult);
    cout << endl << "TorsionalBody energy is  = " << myResult.getEnergy() << endl;
    // Change field and recompute energy
    for (int i=0; i<PbDoF; i++) {
      myState.linearizedUpdate(i, double(0.1*rand())/RAND_MAX);
    }
    myResult.resetResults(ENERGY);
    TorBody.compute(&myResult);
    cout << endl << "TorsionalBody energy is  = " << myResult.getEnergy() << endl;

    TorBody.checkConsistency(&myResult, perturbationFactor, myRequest, myH, myTol);
   

    
    cout << endl << " END OF TEST OF 1st TORSIONAL BODY " << endl;
    cout <<         " ---------------------------------- " << endl << endl;
  }



}
