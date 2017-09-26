// Test for Pressure Body

#include "PotentialBody.h"
#include "MechanicsBody.h"
#include "FEMesh.h"
#include "EigenResult.h"
#include "CompNeoHookean.h"
#include "Harmonic.h"

using namespace voom;

int main(int argc, char** argv) {

  //Test Potential model
  {
    srand(time(NULL));
    cout << " ----------------------------------- " << endl;
    cout << " TEST OF POTENTIAL BODY " << endl << endl;
  
    State myState, AuxState;
    bool CheckOverlap = false;
    int  dofPerNode   = 3;
    FEMesh CubeMesh("../../Mesh/Test/Cube.node", "../../Mesh/Test/Cube.ele", &myState, dofPerNode, CheckOverlap);
    FEMesh MembraneMesh("../../Mesh/Test/CubeTZ.node", "../../Mesh/Test/SurfCube.ele", &AuxState, dofPerNode, CheckOverlap);
    // FEMesh myFEmesh("../../Mesh/Test/CubeQuad.node", "../../Mesh/Test/SurfCubeQuad.ele", &myState, dofPerNode, CheckOverlap);   
    uint NumQP = 1;
    

    // Itialize Auxiliary Mech Body
    int NumMat = CubeMesh.getNumberOfElements()*NumQP;
    vector<MechanicsMaterial * > materials;
    materials.reserve(NumMat);
    for (int k = 0; k < NumMat; k++) { 
      CompNeoHookean* Mat = new CompNeoHookean(k, 1.0+double(rand())/RAND_MAX, 3.0+double(rand())/RAND_MAX);
      materials.push_back(Mat);
    }
    
    MechanicsBody MecBody(&CubeMesh, &myState, materials);
    int PbDoF = myState.getDOFcount();
    MecBody.writeOutputVTK("MecBody", 0);


    // Itialize Potential Body
    Real k  = double(rand())/RAND_MAX;
    Real r0 = double(rand())/RAND_MAX;
    Harmonic PairMat( k, r0 );
    Real SearchR = 3.0;
    vector<int > BodyNodes;
    for (int i = 19; i < 27; i++)
      BodyNodes.push_back(i);
    PotentialBody PotBody(&MembraneMesh, &myState, &PairMat, BodyNodes, SearchR);
    PotBody.writeOutputVTK("PotBody", 0);


    // Test consistency
    EigenResult myResult(PbDoF, 0);
    myResult.initializeResults(PbDoF*4);


    Real perturbationFactor = 0.1;
    int myRequest = FORCE | STIFFNESS; // Check both Forces and Stiffness
    Real myH = 1e-6;
    Real myTol = 1e-7;
    myResult.setRequest(ENERGY);
    myResult.resetResults(ENERGY);
    PotBody.compute(&myResult);
    cout << endl << "PotentialBody energy is  = " << myResult.getEnergy() << endl;
    // Change field and recompute energy
    for (int i=0; i<PbDoF; i++) {
      myState.linearizedUpdate(i, double(0.1*rand())/RAND_MAX);
    }
    myResult.resetResults(ENERGY);
    PotBody.compute(&myResult);
    cout << endl << "PotentialBody energy is  = " << myResult.getEnergy() << endl;

    PotBody.checkConsistency(&myResult, perturbationFactor, myRequest, myH, myTol);
   

    
    cout << endl << " END OF TEST OF 1st POTENTIAL BODY " << endl;
    cout <<         " ---------------------------------- " << endl << endl;
  }



}
