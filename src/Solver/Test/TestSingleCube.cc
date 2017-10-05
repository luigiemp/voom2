#include "State.h"
#include "FEMesh.h"
#include "Model.h"
#include "PressureBody.h"
#include "PassMyoA.h"
#include "CompNeoHookean.h"
#include "EigenNRsolver.h"
#include "LBFGSB.h"

using namespace voom;

int main(int argc, char** argv)
{

  // Timing
  time_t start, end;
  time(&start);
  
  // Initialize State and Mesh
  State myState;
  bool CheckOverlap = false;
  int dofPerNode = 3;
  // FEMesh CubeMesh("../../../../voom2/src/Mesh/Test/CubeQuad.node", "../../../../voom2/src/Mesh/Test/CubeQuad.ele", &myState, dofPerNode, CheckOverlap);
  FEMesh CubeMesh("../../../../voom2/src/Mesh/Test/CubeQuad.node", "../../../../voom2/src/Mesh/Test/CubeQuad.ele", &myState, dofPerNode, CheckOverlap);
  CheckOverlap = true;
  FEMesh SurfFEmesh("../../../../voom2/src/Mesh/Test/CubeQuad.node", "../../../../voom2/src/Mesh/Test/SurfCubeQuad.ele", &myState, dofPerNode, CheckOverlap);   


  // Initialize Materials   
  vector<MechanicsMaterial * > materials;
  vector<GeomElement* > Elements = CubeMesh.getElements();
  int indMat = 0;
  for (int e = 0; e < Elements.size(); e++) {
    for (int q = 0; q < Elements[e]->getNumberOfQuadPoints(); q++) { // One material object per QP
      // PassMyoA* Mat = new PassMyoA(k, 39.02, 7.31, 100.0, 0.01, 2.85, 2.82);
      // Vector3d N; N << 1.0/sqrt(3.0), 1.0/sqrt(3.0), 1.0/sqrt(3.0);
      // Vector3d N; N << 0.0, 0.0, 1.0;
      // Mat->setN(N);
      // materials.push_back(Mat);
      materials.push_back(new CompNeoHookean(indMat, 10.0, 3.0) );
      indMat++;
    }
  }

  // Initialize Body
  MechanicsBody CubeBody(&CubeMesh, &myState, materials);

  // Initialize Pressure Body
  Real Pressure = 2.0;
  PressureBody PBody(&SurfFEmesh, &myState, Pressure);
  


  // Initialize Model
  vector<Body *> Bodies;
  Bodies.push_back(&CubeBody);
  Bodies.push_back(&PBody);
  Model myModel(Bodies, &myState);

  

  // Initialize Solver
  vector<int > DoFid;
  vector<Real > DoFvalues;
  for (int n = 0; n < myState.getXsize(); n++) {
    Vector3d Xtemp = myState.getX(n);
    if ( Xtemp(2) < 1.0e-8 ) {
      if ( Xtemp(0) < 1.0e-8 && Xtemp(1) < 1.0e-8 ) {
	for (int i = 0; i < 3; i++) {
	  DoFid.push_back(n*3 + i);   DoFvalues.push_back(Xtemp(i));
	}
      }
      else if ( Xtemp(0) < 1.0e-8 ) {
	DoFid.push_back(n*3    );   DoFvalues.push_back(Xtemp(0));
	DoFid.push_back(n*3 + 2);   DoFvalues.push_back(Xtemp(2));
      }
      else if ( Xtemp(1) < 1.0e-8 ) {
	DoFid.push_back(n*3 + 1);   DoFvalues.push_back(Xtemp(1));
	DoFid.push_back(n*3 + 2);   DoFvalues.push_back(Xtemp(2));
      }
      else {
	DoFid.push_back(n*3 + 2);   DoFvalues.push_back(Xtemp(2));
      }
    }
    // else if ( Xtemp(2) > 2-1.0e-8 ) {
    //   DoFid.push_back(n*3 + 2);
    //   DoFvalues.push_back(Xtemp(2)-0.2);
    // }
  }

         // EigenResult myR(myState.getDOFcount(), 0);
	 // for (int i = 0; i < DoFid.size(); i++) {
	 //   myState.setPhi(DoFid[i], DoFvalues[i]); // Set known displacements
	 // }
	 // myR.setRequest(ENERGY);
	 // PBody.compute(&myR);
	 // cout << "SurfEneryg is = " << myR.getEnergy() << endl;
	 // CubeBody.compute(&myR);
	 // cout << "CubeEneryg is = " << myR.getEnergy() << endl;

            // Real ResX = 0.0, ResY = 0.0, ResZ = 0.0;
	    // myR.setRequest(FORCE);
	    // myModel.compute(&myR);
	    // for (int n = 0; n < myState.getXsize(); n++) {
	    //   cout << myState.getX(n, 0) << " " << myState.getX(n, 1) << " " << myState.getX(n, 2) << " - ";
	    //   for (int i = 0; i<3; i++) {
	    // 	cout << myR.getResidual(n*3+i) << " ";
	    //   }
	    //   cout << endl;
	    //   ResX +=  myR.getResidual(n*3);
	    //   ResY +=  myR.getResidual(n*3+1);
	    //   ResZ +=  myR.getResidual(n*3+2);
	    // }
	    // cout << ResX << " " << ResY << " " << ResZ << endl;



  EigenNRsolver mySolver(&myState, &myModel, 
			 DoFid, DoFvalues,
			 CHOL, 1.0e-8, 100);

  // CubeBody.writeOutputVTK("SingleCube_NBC", 0);
  mySolver.solve(DISP);
  CubeBody.writeOutputVTK("SingleCube_NBC", 1);



            EigenResult myR(myState.getDOFcount(), 0);
            Real ResX = 0.0, ResY = 0.0, ResZ = 0.0;
	    myR.setRequest(FORCE);
	    myModel.compute(&myR);
	    for (int n = 0; n < myState.getXsize(); n++) {
	      // cout << myState.getX(n, 0) << " " << myState.getX(n, 1) << " " << myState.getX(n, 2) << " - ";
	      // for (int i = 0; i<3; i++) {
	      // 	cout << myR.getResidual(n*3+i) << " ";
	      // }
	      // cout << endl;
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
