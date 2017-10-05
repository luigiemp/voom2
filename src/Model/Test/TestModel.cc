// Test for Model

#include "Model.h"
#include "EigenResult.h"
#include "MechanicsBody.h"
#include "PressureBody.h"
#include "FEMesh.h"
#include "PassMyoA.h"
#include "CompNeoHookean.h"



using namespace voom;

int main(int argc, char** argv) {

  //Test Mechanics model
  {
    srand(time(NULL));
    cout << " ----------------------------------- " << endl;
    cout << " TEST MODEL CLASS USING SINGLE MECHANICS BODY " << endl << endl;

    State myState;
    bool CheckOverlap = false;
    int dofPerNode = 3;
    FEMesh myFEmesh("../../../../voom2/src/Mesh/Test/Cube.node", "../../../../voom2/src/Mesh/Test/Cube.ele", &myState, dofPerNode, CheckOverlap);   
    int NumQP = 1;
    

    int NumMat = myFEmesh.getNumberOfElements()*NumQP;
    vector<MechanicsMaterial * > materials;
    materials.reserve(NumMat);
    for (int k = 0; k < NumMat; k++) {
      vector<Vector3d > Fibers;
      Vector3d N; N << 1.0, 0.0, 0.0;
      Fibers.push_back(N);
      // PassMyoA* Mat = new PassMyoA(k, 1.0+double(rand())/RAND_MAX, 3.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX,  1.0+double(rand())/RAND_MAX, 2.0+double(rand())/RAND_MAX, Fibers);

      CompNeoHookean* Mat = new CompNeoHookean(k, 1.0+double(rand())/RAND_MAX, 3.0+double(rand())/RAND_MAX);
      materials.push_back(Mat);
    }
    
    MechanicsBody myBody(&myFEmesh, &myState, materials);
    int PbDoF = myState.getDOFcount();
    vector<Body *> Bodies;
    Bodies.push_back(&myBody);

    Model myModel(Bodies, &myState);

    EigenResult myResult(PbDoF, NumMat*(materials[0]->getMaterialParameters()).size() );
    cout << endl << "Material parameters per element = " << (materials[0]->getMaterialParameters()).size() << endl;
    myResult.initializeResults(PbDoF*4);



    Real perturbationFactor = 0.1;
    int myRequest = FORCE | STIFFNESS; // Check both Forces and Stiffness
    Real myH = 1e-6;
    Real myTol = 1e-7;
    myResult.setRequest(ENERGY);
    myResult.resetResults(ENERGY);
    myModel.compute(&myResult);
    cout << endl << "Model energy is  = " << myResult.getEnergy() << endl;

    myModel.checkConsistency(&myResult, perturbationFactor, myRequest, myH, myTol);


    
    cout << endl << " END OF TEST OF 1st MODEL WITH 1 BODY" << endl;
    cout <<         " ---------------------------------- " << endl << endl;
  }

  {
    srand(time(NULL));
    cout << " ----------------------------------- " << endl;
    cout << " TEST MODEL CLASS USING MULTIPLE BODIES " << endl << endl;

    // Initialize State and Mesh
    State myState;
    bool CheckOverlap = true;
    int dofPerNode = 3;
    FEMesh CubeMeshOne("../../../../voom2/src/Mesh/Test/Cube.node", "../../../../voom2/src/Mesh/Test/Cube.ele", &myState, dofPerNode, CheckOverlap);
    FEMesh CubeMeshTwo("../../../../voom2/src/Mesh/Test/CubeTZ.node", "../../../../voom2/src/Mesh/Test/Cube.ele", &myState, dofPerNode, CheckOverlap);
    FEMesh SurfFEmesh("../../../../voom2/src/Mesh/Test/Cube.node", "../../../../voom2/src/Mesh/Test/SurfCube.ele", &myState, dofPerNode, CheckOverlap);   

  

    // Initialize Materials   
    // Cube 1
    vector<MechanicsMaterial * > materialsOne;
    vector<GeomElement* > ElementsOne = CubeMeshOne.getElements();
    int indMat = 0;
    for (int e = 0; e < ElementsOne.size(); e++) {
      for (int q = 0; q < ElementsOne[e]->getNumberOfQuadPoints(); q++) { // One material object per QP
	// PassMyoA* Mat = new PassMyoA(k, 39.02, 7.31, 100.0, 0.01, 2.85, 2.82);
	// Vector3d N; N << 1.0/sqrt(3.0), 1.0/sqrt(3.0), 1.0/sqrt(3.0);
	// Vector3d N; N << 0.0, 0.0, 1.0;
	// Mat->setN(N);
	// materials.push_back(Mat);
	materialsOne.push_back(new CompNeoHookean(indMat, 10.0, 3.0) );
	indMat++;
      }
    }
    
    // Initialize Mech Body One
    MechanicsBody CubeBodyOne(&CubeMeshOne, &myState, materialsOne);
    
    // Cube 2
    vector<MechanicsMaterial * > materialsTwo;
    vector<GeomElement* > ElementsTwo = CubeMeshTwo.getElements();
    for (int e = 0; e < ElementsTwo.size(); e++) {
      for (int q = 0; q < ElementsTwo[e]->getNumberOfQuadPoints(); q++) { // One material object per QP
	// PassMyoA* Mat = new PassMyoA(k, 39.02, 7.31, 100.0, 0.01, 2.85, 2.82);
	// Vector3d N; N << 1.0/sqrt(3.0), 1.0/sqrt(3.0), 1.0/sqrt(3.0);
	// Vector3d N; N << 0.0, 0.0, 1.0;
	// Mat->setN(N);
	// materials.push_back(Mat);
	materialsTwo.push_back(new CompNeoHookean(indMat, 10.0, 3.0) );
	indMat++;
      }
    }
    
    // Initialize Mech Body Two
    MechanicsBody CubeBodyTwo(&CubeMeshTwo, &myState, materialsTwo);


    // Initialize Pressure Body
    Real Pressure = 1.0;
    PressureBody PressureBody(&SurfFEmesh, &myState, Pressure);



    // Initialize Model
    vector<Body *> Bodies;
    Bodies.push_back(&CubeBodyOne);
    Bodies.push_back(&CubeBodyTwo);
    Bodies.push_back(&PressureBody);
    Model myModel(Bodies, &myState);


    int PbDoF = myState.getDOFcount();
    EigenResult myResult(PbDoF, (indMat-1)*(materialsOne[0]->getMaterialParameters()).size() );
    myResult.initializeResults(PbDoF*4);



    Real perturbationFactor = 0.1;
    int myRequest = FORCE | STIFFNESS; // Check both Forces and Stiffness
    Real myH = 1e-6;
    Real myTol = 1e-7;
    myResult.setRequest(ENERGY);
    myResult.resetResults(ENERGY);
    myModel.compute(&myResult);
    cout << endl << "Model energy is  = " << myResult.getEnergy() << endl;

    myModel.checkConsistency(&myResult, perturbationFactor, myRequest, myH, myTol);


  
    cout << endl << " END OF TEST OF MODEL WITH MULTIPLE BODIES" << endl;
    cout <<         " ---------------------------------- " << endl << endl;
  }
}
