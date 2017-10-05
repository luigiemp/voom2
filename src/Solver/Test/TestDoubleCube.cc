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
  // FEMesh CubeMeshOne("../../../../voom2/src/Mesh/Test/Cube.node", "../../../../voom2/src/Mesh/Test/Cube.ele", &myState, dofPerNode, CheckOverlap);
  // CheckOverlap = true;
  // FEMesh CubeMeshTwo("../../../../voom2/src/Mesh/Test/CubeTZ.node", "../../../../voom2/src/Mesh/Test/Test/Cube.ele", &myState, dofPerNode, CheckOverlap);
  // FEMesh SurfFEmesh("../../../../voom2/src/Mesh/Test/Cube.node", "../../../../voom2/src/Mesh/Test/SurfCube.ele", &myState, dofPerNode, CheckOverlap);
  FEMesh CubeMeshOne("../../../../voom2/src/Mesh/Test/CubeTet938.node", "../../../../voom2/src/Mesh/Test/CubeTet938.ele", &myState, dofPerNode, CheckOverlap);
  CheckOverlap = true;
  FEMesh CubeMeshTwo("../../../../voom2/src/Mesh/Test/CubeTet938TZ.node", "../../../../voom2/src/Mesh/Test/CubeTet938.ele", &myState, dofPerNode, CheckOverlap);
  FEMesh SurfFEmesh("../../../../voom2/src/Mesh/Test/CubeTet938.node", "../../../../voom2/src/Mesh/Test/SurfCubeTet938.ele", &myState, dofPerNode, CheckOverlap);   
 
 

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
  Real Pressure = 2.0;
  PressureBody PBody(&SurfFEmesh, &myState, Pressure);
 


  // Initialize Model
  vector<Body *> Bodies;
  Bodies.push_back(&CubeBodyOne);
  Bodies.push_back(&CubeBodyTwo);
  Bodies.push_back(&PBody);
  Model myModel(Bodies, &myState);
  


  // Initialize Solver
  vector<int > DoFid;
  vector<Real > DoFvalues;
  // for (int n = 0; n < myState.getXsize(); n++) {
  //   Vector3d Xtemp = myState.getX(n);
  //   if ( Xtemp(2) < 1.0e-8 ) {
  //     if ( Xtemp(0) > 1-1.0e-8 && Xtemp(0) < 1+1.0e-8 && 
  // 	   Xtemp(1) > 1-1.0e-8 && Xtemp(1) < 1+1.0e-8 ) {
  // 	for (int i = 0; i < 3; i++) {
  // 	  DoFid.push_back(n*3 + i);   DoFvalues.push_back(Xtemp(i));
  // 	}
  //     }
  //     else if ( Xtemp(0) > 1-1.0e-8 && Xtemp(0) < 1+1.0e-8 ) {
  // 	DoFid.push_back(n*3    );   DoFvalues.push_back(Xtemp(0));
  // 	DoFid.push_back(n*3 + 2);   DoFvalues.push_back(Xtemp(2));
  //     }
  //     else if ( Xtemp(1) > 1-1.0e-8 && Xtemp(1) < 1+1.0e-8 ) {
  // 	DoFid.push_back(n*3 + 1);   DoFvalues.push_back(Xtemp(1));
  // 	DoFid.push_back(n*3 + 2);   DoFvalues.push_back(Xtemp(2));
  //     }
  //     else {
  // 	DoFid.push_back(n*3 + 2);   DoFvalues.push_back(Xtemp(2));
  //     }
  //   }
  //   // else if ( Xtemp(2) > (2.0 - 1.0e-5) ) {
  //   // 	DoFid.push_back(n*3 + 2);   DoFvalues.push_back(Xtemp(2)-0.2);
  //   // }
  // }

  // for (int n = 0; n < myState.getXsize(); n++) {
  //   Vector3d Xtemp = myState.getX(n);
  //   if ( Xtemp(2) < 1.0e-8 ) {
  //     if ( Xtemp(0) < 1.0e-8 && Xtemp(1) < 1.0e-8 ) {
  // 	for (int i = 0; i < 3; i++) {
  // 	  DoFid.push_back(n*3 + i);   DoFvalues.push_back(Xtemp(i));
  // 	}
  //     }
  //     else if ( Xtemp(0) < 1.0e-8 ) {
  // 	DoFid.push_back(n*3    );   DoFvalues.push_back(Xtemp(0));
  // 	DoFid.push_back(n*3 + 2);   DoFvalues.push_back(Xtemp(2));
  //     }
  //     else if ( Xtemp(1) < 1.0e-8 ) {
  // 	DoFid.push_back(n*3 + 1);   DoFvalues.push_back(Xtemp(1));
  // 	DoFid.push_back(n*3 + 2);   DoFvalues.push_back(Xtemp(2));
  //     }
  //     else {
  // 	DoFid.push_back(n*3 + 2);   DoFvalues.push_back(Xtemp(2));
  //     }
  //   }
  //   // else if ( Xtemp(2) > (2.0 - 1.0e-5) ) {
  //   //   DoFid.push_back(n*3 + 2);   DoFvalues.push_back(Xtemp(2)-0.2);
  //   // }
  // }

  Real Xb = -0.5, Yb = -0.5, Zb = 0.0;
  for (int n = 0; n < myState.getXsize(); n++) {
    Vector3d Xtemp = myState.getX(n);
    if ( Xtemp(2) < Zb + 1.0e-8 ) {
      if ( Xtemp(0) < Xb + 1.0e-8 && Xtemp(1) < Yb + 1.0e-8 ) {
  	for (int i = 0; i < 3; i++) {
  	  DoFid.push_back(n*3 + i);   DoFvalues.push_back(Xtemp(i));
  	}
      }
      else if ( Xtemp(0) < Xb + 1.0e-8 ) {
  	DoFid.push_back(n*3    );   DoFvalues.push_back(Xtemp(0));
  	DoFid.push_back(n*3 + 2);   DoFvalues.push_back(Xtemp(2));
      }
      else if ( Xtemp(1) < Yb + 1.0e-8 ) {
  	DoFid.push_back(n*3 + 1);   DoFvalues.push_back(Xtemp(1));
  	DoFid.push_back(n*3 + 2);   DoFvalues.push_back(Xtemp(2));
      }
      else {
  	DoFid.push_back(n*3 + 2);   DoFvalues.push_back(Xtemp(2));
      }
    }
    // else if ( Xtemp(2) > (1.0 - 1.0e-5) && Xtemp(2) < (1.0 + 1.0e-5) ) {
    //   DoFid.push_back(n*3 + 2);   DoFvalues.push_back(Xtemp(2)-0.05);
    // }
  }

         // myState.printPhi();
	 // for (int i = 0; i<DoFid.size(); i++) {
	 //   cout << DoFid[i] << " " << DoFvalues[i] << endl;
	 // }
	 
	 // cout << myState.getDOFcount() << endl;
	 // for (int i = 0; i < myState.getXsize(); i++) {
	 //   cout << myState.getX(i,0)   << " " << myState.getX(i,1)   << " " << myState.getX(i,2) << " ";
	 //   cout << myState.getPhi(i,0) << " " << myState.getPhi(i,1) << " " << myState.getPhi(i,2) << " ";
	 //   cout << myState.getNodeDof(i) << " " << myState.getGdof(i) << endl; 
	 // }



  EigenNRsolver mySolver(&myState, &myModel, 
			 DoFid, DoFvalues,
			 CHOL, 1.0e-12, 100);

  CubeBodyOne.writeOutputVTK("DoubleCubeOne_NBC", 0);
  CubeBodyTwo.writeOutputVTK("DoubleCubeTwo_NBC", 0);
  PBody.writeOutputVTK("DoubleCubeSurf_NBC", 0);
  mySolver.solve(DISP);
  CubeBodyOne.writeOutputVTK("DoubleCubeOne_NBC", 1);
  CubeBodyTwo.writeOutputVTK("DoubleCubeTwo_NBC", 1);
  PBody.writeOutputVTK("DoubleCubeSurf_NBC", 1);



  EigenResult myR(myState.getDOFcount(), 0);
  Real ResX = 0.0, ResY = 0.0, ResZ = 0.0;
  myR.setRequest(FORCE);
  myModel.compute(&myR);
  for (int n = 0; n < myState.getXsize(); n++) {
    // cout << myState.getX(n, 0) << " " << myState.getX(n, 1) << " " << myState.getX(n, 2) <<;
    // for (int i = 0; i<3; i++) {
    //   cout << myR.getResidual(n*3+i) << " ";
    // }
    // cout << endl;
    ResX +=  myR.getResidual(n*3);
    ResY +=  myR.getResidual(n*3+1);
    ResZ +=  myR.getResidual(n*3+2);
  }
  cout << endl << "ResX = " << ResX << " - ResY = " << ResY << " - ResZ = " << ResZ << endl;
 
  // myState.printPhi();


  // Timing
  time (&end);
  cout << endl << "BVP solved in " << difftime(end,start) << " s" << endl;



  return 0;
}
