#include "PassMyoA.h"
#include "CompNeoHookean.h"
#include "FEMesh.h"
#include "EigenEllipticResult.h"
#include "MechanicsModel.h"
#include "EigenNRsolver.h"
#include "LBFGSB.h"

using namespace voom;

int main(int argc, char** argv)
{
  // Initialize Mesh
  // FEMesh BodyMesh("./OneTet.node", "./OneTet.ele");
  // FEMesh SurfMesh("./OneTet.node", "./SurfOneTet.ele");
  // FEMesh BodyMesh("./QuadTet.node", "./QuadTet.ele");
  // FEMesh SurfMesh("./QuadTet.node", "./QuadTet.surf");
  // EMesh BodyMesh("../../Mesh/Test/Cube.node", "../../Mesh/Test/Cube.ele");
  // FEMesh SurfMesh("../../Mesh/Test/Cube.node", "../../Mesh/Test/SurfCube.ele");
  // FEMesh BodyMesh("../../Mesh/Test/CubeQuad.node", "../../Mesh/Test/CubeQuad.ele");
  // FEMesh SurfMesh("../../Mesh/Test/CubeQuad.node", "../../Mesh/Test/SurfCubeQuad.ele");
  FEMesh BodyMesh("../../Mesh/Test/CoarseLV.node", "../../Mesh/Test/CoarseLV.ele");
  FEMesh SurfMesh("../../Mesh/Test/CoarseLV.node", "../../Mesh/Test/CoarseLV.surf");
 
  
    
  // Initialize Material
  uint NumMat =  BodyMesh.getNumberOfElements();
  vector<MechanicsMaterial * > materials;
  vector<Vector3d > Fibers;
  
  materials.reserve(NumMat);
  for (int k = 0; k < NumMat; k++) {
    PassMyoA* Mat = new PassMyoA(k, 30.48, 7.25, 100.0, 1.0, 2.60, 2.82);
    Vector3d N; N << 0.0, 0.0, 1.0;
    Fibers.push_back(N);
    materials.push_back(Mat);
  }
  
  BodyMesh.setFibers(Fibers);

  // Initialize Model
  int NodeDoF = 3;
  int PressureFlag = 1;
  int NodalForcesFlag = 0;
  Real Pressure = 1.0;
  MechanicsModel myModel(&BodyMesh, materials, NodeDoF, PressureFlag, &SurfMesh, NodalForcesFlag);
  myModel.updatePressure(Pressure);
  
 
  // Initialize Result
  uint PbDoF = (BodyMesh.getNumberOfNodes())*myModel.getDoFperNode();
  EigenEllipticResult myResults(PbDoF, NumMat);

  // Check on applied pressure
  int myRequest = 2;
  myResults.setRequest(myRequest);
  myModel.compute(myResults);
  VectorXd Fstart = *(myResults._residual);
  Real pX = 0.0, pY = 0.0, pZ = 0.0;
  for (int i = 0; i < Fstart.size(); i += 3) {
    pX += Fstart(i);
    pY += Fstart(i+1);
    pZ += Fstart(i+2);
  }
  cout << "Pressure = " << pX << " " << pY << " " << pZ << endl;



  return 0;
}
