#include "State.h"
#include "FEMesh.h"
#include "Model.h"
#include "PotentialBody.h"
#include "MechanicsBody.h"
#include "Harmonic.h"
#include "ImposedKinematics.h"
#include "EigenNRsolver.h"
#include "LBFGSB.h"

using namespace voom;

int main(int argc, char** argv)
{

  // Timing
  time_t start, end;
  time(&start);
  
  // Initialize State and Mesh
  State myState, AuxState;
  bool CheckOverlap = false;
  int dofPerNode = 3;
  FEMesh CubeMesh("../../Mesh/Test/Cube.node", "../../Mesh/Test/Cube.ele", &myState, dofPerNode, CheckOverlap);
  // Membrane meshes
  FEMesh MembraneMeshA("../../Mesh/Test/Cube.node", "../../Mesh/Test/SurfCubeZ.ele", &AuxState, dofPerNode, CheckOverlap);
  CheckOverlap = true;
  FEMesh MembraneMeshB("../../Mesh/Test/Cube.node", "../../Mesh/Test/SurfCubeX.ele", &AuxState, dofPerNode, CheckOverlap);
  FEMesh MembraneMeshC("../../Mesh/Test/Cube.node", "../../Mesh/Test/SurfCubeY.ele", &AuxState, dofPerNode, CheckOverlap);
 
 
  

  // Initialize bodies
  // Initialize Materials   
  // Cube 
  vector<MechanicsMaterial * > materials;
  vector<GeomElement* > Elements = CubeMesh.getElements();
  int indMat = 0;
  vector<Real > Alphas(3, 0.0), Stretches(3, 0.0); 
  Alphas[0]    = 1.0;    Alphas[1]    = 1.0;    Alphas[2]    = 1.0;
  Stretches[0] = 1.2;    Stretches[1] = 1.2;    Stretches[2] = 0.5;
  vector<Vector3d > Directions(3, Vector3d::Zero());
  Vector3d N1, N2, N3;
  N1 << 1.0, 0.0, 0.0; N2 << 0.0, 1.0, 0.0; N3 << 0.0, 0.0, 1.0;
  Directions[0] = N1; Directions[1] = N2; Directions[2] = N3;
  Real Beta = 1.0;
  for (int e = 0; e < Elements.size(); e++) {
    for (int q = 0; q < Elements[e]->getNumberOfQuadPoints(); q++) { // One material object per QP
      materials.push_back(new ImposedKinematics(indMat, Alphas, Stretches, Directions, Beta));
      indMat++;
    }
  }
    
  // Initialize Mech Body One
  MechanicsBody Cube(&CubeMesh, &myState, materials);

       int PbDoF = myState.getDOFcount();
       EigenResult TestResult(PbDoF, 0);
       TestResult.initializeResults(PbDoF*4);
 
       TestResult.setRequest(ENERGY);
       TestResult.resetResults(ENERGY);
       Cube.compute(&TestResult);
       cout << endl << "Cube energy is  = " << TestResult.getEnergy() << endl;
  

  Real kLJ = 1.0;
  Real r0  = 0.0;
  Real SearchR = 1.1;
  Harmonic PairMat( kLJ, r0 );
  vector<int > BodyNodesA(9, 0), BodyNodesB(9, 0), BodyNodesC(9, 0);
  BodyNodesA[0] = 0; BodyNodesB[0] = 0;  BodyNodesC[0] = 0;
  BodyNodesA[1] = 1; BodyNodesB[1] = 3;  BodyNodesC[1] = 1;
  BodyNodesA[2] = 2; BodyNodesB[2] = 6;  BodyNodesC[2] = 2;
  BodyNodesA[3] = 3; BodyNodesB[3] = 9;  BodyNodesC[3] = 9;
  BodyNodesA[4] = 4; BodyNodesB[4] = 12; BodyNodesC[4] = 10;
  BodyNodesA[5] = 5; BodyNodesB[5] = 15; BodyNodesC[5] = 11;
  BodyNodesA[6] = 6; BodyNodesB[6] = 18; BodyNodesC[6] = 18;
  BodyNodesA[7] = 7; BodyNodesB[7] = 21; BodyNodesC[7] = 19;
  BodyNodesA[8] = 8; BodyNodesB[8] = 24; BodyNodesC[8] = 20;

  PotentialBody MemA(&MembraneMeshA, &myState, &PairMat, BodyNodesA, SearchR);
  TestResult.resetResults(ENERGY);
  MemA.compute(&TestResult);
  cout << endl << "MemA energy is  = " << TestResult.getEnergy() << endl;
  PotentialBody MemB(&MembraneMeshB, &myState, &PairMat, BodyNodesB, SearchR);
  TestResult.resetResults(ENERGY);
  MemB.compute(&TestResult);
  cout << endl << "MemB energy is  = " << TestResult.getEnergy() << endl;
  PotentialBody MemC(&MembraneMeshC, &myState, &PairMat, BodyNodesC, SearchR);
  TestResult.resetResults(ENERGY);
  MemC.compute(&TestResult);
  cout << endl << "MemC energy is  = " << TestResult.getEnergy() << endl;




  
  // Initialize Model
  vector<Body *> Bodies;
  Bodies.push_back(&Cube);
  Bodies.push_back(&MemA);
  Bodies.push_back(&MemB);
  Bodies.push_back(&MemC);
  Model myModel(Bodies, &myState);
  




  // Initialize Solver
  vector<int > DoFid;
  vector<Real > DoFvalues;

  EigenNRsolver mySolver(&myState, &myModel, 
			 DoFid, DoFvalues,
			 CHOL, 1.0e-12, 10);

  Cube.writeOutputVTK("Cube", 0);
  MemA.writeOutputVTK("MemA", 0);
  MemB.writeOutputVTK("MemB", 0);
  MemC.writeOutputVTK("MemC", 0);
  mySolver.solve(DISP);
  Cube.writeOutputVTK("Cube", 1);
  MemA.writeOutputVTK("MemA", 1);
  MemB.writeOutputVTK("MemB", 1);
  MemC.writeOutputVTK("MemC", 1);
  
  EigenResult myR(myState.getDOFcount(), 0);
  Real ResX = 0.0, ResY = 0.0, ResZ = 0.0;
  myR.setRequest(FORCE);
  myModel.compute(&myR);
  for (int n = 0; n < myState.getXsize(); n++) {
    ResX +=  myR.getResidual(n*3);
    ResY +=  myR.getResidual(n*3+1);
    ResZ +=  myR.getResidual(n*3+2);
  }
  cout << endl << "ResX = " << ResX << " - ResY = " << ResY << " - ResZ = " << ResZ << endl;

  

  // Timing
  time (&end);
  cout << endl << "BVP solved in " << difftime(end,start) << " s" << endl;



  return 0;
}
