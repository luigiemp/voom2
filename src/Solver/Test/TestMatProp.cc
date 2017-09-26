#include "State.h"
#include "FEMesh.h"
#include "Model.h"
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
  FEMesh CubeOneMesh("../../Mesh/Test/Cube.node", "../../Mesh/Test/Cube.ele", &myState, dofPerNode, CheckOverlap);
  CheckOverlap = true;
  FEMesh CubeTwoMesh("../../Mesh/Test/CubeTranslatedZ.node", "../../Mesh/Test/Cube.ele", &myState, dofPerNode, CheckOverlap);

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
    
    // Initialize Body
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
    
    // Initialize Body
    MechanicsBody CubeBodyTwo(&CubeMeshTwo, &myState, materialsTwo);



  // Initialize Model
  vector<Body *> Bodies;
  Bodies.push_back(&CubeBodyOne);
  Bodies.push_back(&CubeBodyTwo);
  Model myModel(Bodies, &myState);

  
  // Initialize Solver
  vector<int > DoFid;
  vector<Real > DoFvalues;
  for (int n = 0; n < myState.getXsize(); n++) {
    Vector3d Xtemp = myState.getX(n);
    if ( Xtemp(2) < 1.0e-8 ) {
      for (int i = 0; i < 3; i++) {
	DoFid.push_back(n*3 + i);
	DoFvalues.push_back(Xtemp(i));
      }
    } else if ( Xtemp(2) > 4.0-1.0e-8 ) {
      DoFid.push_back(n*3 + 2);
      DoFvalues.push_back(Xtemp(2)-0.2);
    }
  }

  EigenNRsolver mySolver(&myState, &myModel, 
			 DoFid, DoFvalues,
			 CHOL, 1.0e-8, 100);

  CubeBodyOne.writeOutputVTK("BVPtestOne", 0);
  CubeBodyTwo.writeOutputVTK("BVPtestTwo", 0);
  mySolver.solve(DISP);
  CubeBodyOne.writeOutputVTK("BVPtestOne", 1);
  CubeBodyTwo.writeOutputVTK("BVPtestTwo", 1);


 
    /*
  // Initialize Result
  uint PbDoF = (Cube.getNumberOfNodes())*myModel.getDoFperNode();
  EigenEllipticResult myResults(PbDoF, NumMat*2);

  // Run Consistency check
  Real perturbationFactor = 0.1;
  uint myRequest = 7; // Check both Forces and Stiffness
  Real myH = 1e-6;
  Real myTol = 1e-7;
  myModel.checkConsistency(myResults, perturbationFactor, myRequest, myH, myTol);
  myModel.checkDmat(myResults, perturbationFactor, myH, myTol);
  
  {
    // Solve for displacement first
    // Print initial configuration
    myModel.writeOutputVTK("CubeQuad_", 0);

    // Check on applied pressure
    myRequest = 2;
    myResults.setRequest(myRequest);
    myModel.compute(myResults);
    VectorXd Fstart = *(myResults._residual);
    Real pX = 0.0, pY = 0.0, pZ = 0.0;
    for (int i = 0; i < Fstart.size(); i += 3) {
      pX += Fstart(i);
      pY += Fstart(i+1);
      pZ += Fstart(i+2);
    }
    cout << endl << "Pressure = " << pX << " " << pY << " " << pZ << endl << endl;
    
    
    
    // EBC
    
    // Quad tet cube applied pressure //
    vector<int > DoFid(13,0);
    vector<Real > DoFvalues(13, 0.0);
    DoFid[0] = 0;  DoFid[1] = 1;  DoFid[2] = 2;  
    DoFid[3] = 5;  DoFid[4] = 8;  DoFid[5] = 11;  
    DoFid[6] = 26; DoFid[7] = 29; DoFid[8] = 35;          
    DoFid[9] = 38; DoFid[10] = 47;
    DoFid[11] = 3; DoFid[12] = 10;
    
    // Lin tet cube applied pressure //
    // vector<int > DoFid(13,0);
    // vector<Real > DoFvalues(13, 0.0);  
    // for (uint i = 0; i < 9; i++) DoFid[i] = i*3 + 2; // Constrained base in z direction
    // DoFid[9] = 0; DoFid[10] = 1;
    // DoFid[11] = 7; DoFid[12] = 18;
 
    for(int i = 0; i < DoFid.size(); i++) {
      DoFvalues[i] = Cube.getX(floor(double(DoFid[i])/3.0) ,DoFid[i]%3);
      // cout << DoFid[i] << " " <<  DoFvalues[i] << endl;
    }

    // Lin tet cube - applied displacements //
    // for (uint i = 18; i < 27; i++) { 
    //   DoFid[i-5] = i*3 + 2; // Lower top in z direction by 0.1
    //   DoFvalues[i-5] = 1.95;
    // }

    // One Tet test //
    // vector<int > DoFid(8,0);
    // DoFid[0] = 0; DoFid[1] = 1; DoFid[2] = 2;
    // DoFid[3] = 4; DoFid[4] = 5; DoFid[5] = 6;
    // DoFid[6] = 8;
    // // DoFid[7] = 9; DoFid[8] = 10;
    // DoFid[7] = 11;
    // vector<Real > DoFvalues(8, 0.0);
    // DoFvalues[7] = 2.20839;
    
    // Solver
    Real NRtol = 1.0e-12;
    uint NRmaxIter = 100;
    EigenNRsolver mySolver(&myModel, DoFid, DoFvalues, CHOL, NRtol, NRmaxIter);
    mySolver.solve(DISP); 

    myModel.writeOutputVTK("CubeQuad_", 1);
    myModel.writeField(DispFile);


    // int m = 5; 
    // double factr = 1.0e+1;
    // double pgtol = 1.0e-5;
    // int iprint = 0;
    // int maxIterations = -1;
    
    // vector<int > nbd(PbDoF, 0);
    // vector<double > lowerBound(PbDoF, 0.0);
    // vector<double > upperBound(PbDoF, 0.0);
    // for(int i = 0; i < DoFid.size(); i++) {
    //   nbd[DoFid[i]] = 2;
    //   lowerBound[DoFid[i]] = DoFvalues[i];
    //   upperBound[DoFid[i]] = DoFvalues[i];
    // }
    // cout << "PbDoF = " << PbDoF << endl;
    // for (int j = 0; j < PbDoF; j++) {
    //   cout << nbd[j] << " " << lowerBound[j] << " " << upperBound[j] << endl;
    // }
    
    // LBFGSB mySolver(& myModel, & myResults, m, factr, pgtol, iprint, maxIterations);
    // mySolver.setBounds(nbd, lowerBound, upperBound);
    // mySolver.solve();
    
    // // Print final configuration
    // myModel.writeOutputVTK("CubeQuad_", 1);
    
    // // Print _field
    // cout << endl;
    // myModel.printField();
    // cout << endl;
    
    // Compute Residual
    myRequest = 2;
    myResults.setRequest(myRequest);
    myModel.compute(myResults);
    VectorXd F = *(myResults._residual);
    for (int i = 0; i < DoFid.size(); i ++) {
      ForcesID.push_back( DoFid[i] );
      Forces.push_back( -F(DoFid[i]) );
    }
  } // End of displacement solution

  // Next recompute material properties using EMFO
  {
    ifstream BCinp(DispFile.c_str());
    int DoF;

    BCinp >> DoF;
    vector<int > DoFid(DoF, 0);
    vector<Real > DoFvalues(DoF, 0);
    for(int i = 0; i < DoF; i++) {
      DoFid[i] = i;
      BCinp >> DoFvalues[i];
      // cout <<  DoFid[i] << " " <<  DoFvalues[i] << endl;
    }
   
    // Change Material properties
    // Find unique material parameters
    set<MechanicsMaterial *> UNIQUEmaterials;
    for (uint i = 0; i < materials.size(); i++) 
      UNIQUEmaterials.insert(materials[i]);

    srand (time(NULL));
    for (set<MechanicsMaterial *>::iterator MatIt = UNIQUEmaterials.begin(); MatIt != UNIQUEmaterials.end(); MatIt++) {
      int MatID = (*MatIt)->getMatID();
      vector<Real > MatProp = (*MatIt)->getMaterialParameters();
      int NumPropPerMat = MatProp.size();
      for (int m = 0; m < NumPropPerMat; m++) {
    	MatProp[m] *= double(rand())/double(RAND_MAX);
	cout << MatProp[m] << " ";
      }
      cout << endl;
      (*MatIt)->setMaterialParameters(MatProp);
    }

    // Solve for correct material properties
    Real NRtol = 1.0e-12;
    uint NRmaxIter = 10;
    EigenNRsolver mySolver(&myModel, DoFid, DoFvalues, CHOL, NRtol, NRmaxIter);
    mySolver.solve(MAT); 
 
    // Print found material properties
    for (set<MechanicsMaterial *>::iterator MatIt = UNIQUEmaterials.begin(); MatIt != UNIQUEmaterials.end(); MatIt++) {
      int MatID = (*MatIt)->getMatID();
      vector<Real > MatProp = (*MatIt)->getMaterialParameters();
      int NumPropPerMat = MatProp.size();
      for (int m = 0; m < NumPropPerMat; m++) {
	cout << MatProp[m] << " ";
      }
      cout << endl;
    }

  } // End of material properties solution
    */

  // Timing
  time (&end);
  cout << endl << "BVP solved in " << difftime(end,start) << " s" << endl;



  return 0;
}
