#include "PassMyoA.h"
#include "CompNeoHookean.h"
#include "Guccione.h"
#include "FEMesh.h"
#include "EigenEllipticResult.h"
#include "MechanicsModel.h"
#include "EigenNRsolver.h"

#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include <Eigen/Eigenvalues> 

using namespace voom;

void writeNeumannBC(vector<int > BCid, vector<Real > Forces, string OutputFile, int step);

int main(int argc, char* argv[])
{
  // Timing
  time_t start, end;
  time(&start);

  // Read in input parameters from command line
  Real a1 = 0.0, a2 = 0.0;
  if (argc != 3) {
    cout << "Input missing, a1 and a2 required!" << endl;
    return -1;
  } else {
    a1 = atof(argv[1]);
    a2 = atof(argv[2]);
    cout << "a1 = " << a1 << ";   a2 = " << a2 << endl;
  }

  // Initialize Mesh
  FEMesh LVmesh("/u/home/l/luigiemp/project-cardio/voom2/src/Applications/PassiveLV/Small_A/Small_A.node", 
		"/u/home/l/luigiemp/project-cardio/voom2/src/Applications/PassiveLV/Small_A/Small_A.ele");
  FEMesh LVsurf("/u/home/l/luigiemp/project-cardio/voom2/src/Applications/PassiveLV/Small_A/Small_A.node", 
		"/u/home/l/luigiemp/project-cardio/voom2/src/Applications/PassiveLV/Small_A/Small_A.surfEle");
  string BCfile = "/u/home/l/luigiemp/project-cardio/voom2/src/Applications/PassiveLV/Small_A/Small_A.bc";
 
  cout << endl;
  cout << "Number Of Nodes   : " << LVmesh.getNumberOfNodes()    << endl;
  cout << "Number Of Element : " << LVmesh.getNumberOfElements() << endl;
  cout << "Mesh Dimension    : " << LVmesh.getDimension()        << endl << endl;
  
  // Initialize Materials and Fibers
  uint NumEl =  LVmesh.getNumberOfElements();
  vector<MechanicsMaterial * > materials;
  materials.reserve(NumEl);

  string FiberFile = "/u/home/l/luigiemp/project-cardio/voom2/src/Applications/PassiveLV/Small_A/Small_A.fiber";
  ifstream FiberInp(FiberFile.c_str());
  int NumMat = 1;
  int NumMatProp = 2;
  
  // CompNeoHookean* Mat = new CompNeoHookean(0, 1.0, 1.0);
  // Guccione* Mat = new Guccione(0, 3.0, 11.1, 1.76, 10.0);
  // PassMyoA* Mat = new PassMyoA(0, 35.19, 7.06, 100.0, 0.025, 2.87, 2.82);
  // PassMyoA* Mat = new PassMyoA(0, 35.19, 7.06, 0.0, 0.025, 2.87, 2.82);
  


  for (int k = 0; k < NumEl; k++) {
    vector<Vector3d > Fibers(3, Vector3d::Zero());
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
	FiberInp >> Fibers[i](j);
      }
    }

    // PassMyoA* Mat = new PassMyoA(0, 35.19, 7.06, 0.0, 0.025, 2.87, 2.82, Fibers);
    // PassMyoA* Mat = new PassMyoA(0, 3.10202, 13.5451, 100.0, 0.025, a1, a2, Fibers);
    Guccione* Mat = new Guccione(0, 3.0, 11.1, 1.76, 10.0, Fibers);
    materials.push_back(Mat);
  }



  // Initialize Model
  int NodeDoF = 3;
  int PressureFlag = 1;
  int NodalForcesFlag = 0;
  MechanicsModel myModel(&LVmesh, materials, NodeDoF, PressureFlag, &LVsurf,
			 NodalForcesFlag);
 
  // Initialize Result
  uint PbDoF = (LVmesh.getNumberOfNodes())*myModel.getDoFperNode();
  EigenEllipticResult myResults(PbDoF, NumMat*NumMatProp);
 
  // Print initial configuration
  myModel.writeOutputVTK("/u/home/l/luigiemp/project-cardio/voom2/src/Applications/PassiveLV/TestInVivoDir/LV_F_", 0);



  // Solve for displacement first
  // Check on applied pressure
  myModel.updatePressure(1.0);
  int myRequest = 2;
  myResults.setRequest(myRequest);
  myModel.compute(myResults);
  VectorXd F = *(myResults._residual);
  Real pX = 0.0, pY = 0.0, pZ = 0.0;
  for (int i = 0; i < F.size(); i += 3) {
    pX += F(i);
    pY += F(i+1);
    pZ += F(i+2);
  }
  cout << endl << "Pressure = " << pX << " " << pY << " " << pZ << endl << endl;

  // EBC
  int NumBC = 0, node = 0, ind = 0;;
  vector<int > BCid;
  vector<Real > BCvalues;
  ifstream BCinp(BCfile.c_str());

  BCinp >> NumBC;
  BCid.reserve(NumBC*3);
  BCvalues.reserve(NumBC*3);
  for(int i = 0; i < NumBC; i++) {
    BCinp >> node;
    for (int j = 0; j < 3; j++) {
      BCid.push_back(node*3 + j);
      BCvalues.push_back(LVmesh.getX(node, j));
      // cout << BCid[ind] << " " <<  BCvalues[ind] << endl;
      ind++;
    }
  }
  
 // Loop through pressure
  int NumPsteps = 10;
  Real DeltaP = -0.1;
    
  {
    // Solve for displacement first
    // Solver
    Real NRtol = 1.0e-10;
    uint NRmaxIter = 100;
    EigenNRsolver mySolver(&myModel, BCid, BCvalues, CHOL, NRtol, NRmaxIter);  

    for (int s = 1; s <= NumPsteps; s++) {
      myModel.updatePressure(Real(s)*DeltaP);
      mySolver.solve(DISP); 

      myModel.writeOutputVTK("/u/home/l/luigiemp/project-cardio/voom2/src/Applications/PassiveLV/TestInVivoDir/LV_Fcomp_", s);
      myModel.writeField("/u/home/l/luigiemp/project-cardio/voom2/src/Applications/PassiveLV/TestInVivoDir/LV_FdispComp_", s);
      
      // Compute Residual
      myResults.setRequest(FORCE);
      myModel.compute(myResults);
      VectorXd F = *(myResults._residual);
      vector<Real > Forces(BCid.size(), 0.0);
      for (int i = 0; i < BCid.size(); i ++) {
  	Forces[i] = -F(BCid[i]);
      }
      writeNeumannBC(BCid, Forces, "/u/home/l/luigiemp/project-cardio/voom2/src/Applications/PassiveLV/TestInVivoDir/LV_FforcesComp_", s);
    }
  }
  // End of displacement solution
  cout << endl << "End of displacement solution" << endl << endl << endl;

 
  /* 
  // Next compute material properties using EMS
  {
    cout << endl << "Material properties identification." << endl;
    // Need to apply all external forces: pressure + any Neumann BC
    myModel.setNodalForcesFlag(1);

    // Initialize solver parameters
    Real NRtol = 1.0e-12;
    uint NRmaxIter = 10;
    uint DeltaSteps = 1;
    vector<vector<double > > MatPropList(NumPsteps, vector<double>(NumMatProp, 0.0));

    // Considering each loading step separately
    // for ( int s = 1; s <= NumPsteps; s+=DeltaSteps ) { \\ here
      // NR loop
      Real error = 1.0;
      uint iter = 0;
      while (iter < NRmaxIter && error > NRtol)
      {
	myResults.setRequest(8);
	myModel.setResetFlag(1);

	for ( int s = 1; s <= NumPsteps; s+=DeltaSteps ) { // here
	cout << "Step considered = " << s << endl;

	int DoF;

	// Read in displacements
	stringstream DispFileNameStream;
	DispFileNameStream << "/u/home/l/luigiemp/project-cardio/voom2/src/Applications/PassiveLV/TestInVivoDir/LV_Fdisp_" << s << ".dat";
	ifstream DBCinp;
	DBCinp.open( (DispFileNameStream.str()).c_str() );

	DBCinp >> DoF;
	vector<int > DoFid(DoF, 0);
	vector<Real > DoFvalues(DoF, 0);
	for(int i = 0; i < DoF; i++) {
	  DoFid[i] = i;
	  DBCinp >> DoFvalues[i];
	  myModel.setField(i, DoFvalues[i]);
	}
   
	// Read in forces
	stringstream ForceFileNameStream;
	ForceFileNameStream << "/u/home/l/luigiemp/project-cardio/voom2/src/Applications/PassiveLV/TestInVivoDir/LV_Fforces_" << s << ".dat";
	ifstream NBCinp;
	NBCinp.open( (ForceFileNameStream.str()).c_str() );

	NBCinp >> DoF;
	vector<int > ForceDoF(DoF, 0);
	vector<Real > ForceValues(DoF, 0);
	for(int i = 0; i < DoF; i++) {
	  NBCinp >> ForceDoF[i];
	  NBCinp >> ForceValues[i];
	}
	myModel.updateNodalForces( &ForceDoF, &(ForceValues) );

	// Update pressure
	myModel.updatePressure(Real(s)*DeltaP);

	// Compute Hg and Gradg
	myModel.compute(myResults);
	myModel.setResetFlag(0);
	} // here

	VectorXd DeltaAlpha;
	SimplicialCholesky<SparseMatrix<Real > > chol(*(myResults._Hg));  // performs a Cholesky factorization of _stiffness
	DeltaAlpha = chol.solve(*(myResults._Gradg));                     // use the factorization to solve for the given right hand side

	// Update material properties in material
	for (int k = 0; k < NumEl; k++) {
	  int MatID = materials[k]->getMatID();
	  vector<Real > MatProp = materials[k]->getMaterialParameters();
	  int NumPropPerMat = MatProp.size();
	  for (int m = 0; m < NumPropPerMat; m++) {
	    MatProp[m] -= DeltaAlpha(MatID*NumPropPerMat + m);
	  }
	  materials[k]->setMaterialParameters(MatProp);
	}

	// Compute condition number
	MatrixXd dHg;
	dHg = MatrixXd(*(myResults._Hg));
	cout << "Hessian = " << dHg << endl;
	VectorXcd evals = dHg.eigenvalues();
	cout << "Max eigenvalues = " << evals.real().maxCoeff() << "; Min eigenvalues = " << evals.real().minCoeff() 
	     << "; Condition Number = " << evals.real().maxCoeff()/evals.real().minCoeff() << endl;
	cout << "Max imag eigenvalues = " << evals.imag().maxCoeff() << "; Min imag eigenvalues = " << evals.imag().minCoeff() << endl;
	
	// Update iter and error
	iter++;
	error = DeltaAlpha.norm();
	myResults.setRequest(ENERGY);
	myModel.compute(myResults);
	cout << "Energy = " << myResults.getEnergy() << "   - NR iter = " << iter << "   -  NR error = " << error << endl;
      
      } // while loop  
   
      cout << endl << "Identified material properties" << endl;

      // Print found material properties  
      int matid = 0;
      int MatID = materials[matid]->getMatID();
      vector<Real > MatProp = materials[matid]->getMaterialParameters();
      for (int m = 0; m < NumMatProp; m++) {
	cout << MatProp[m] << " ";           // here
	// MatPropList[s-1][m] = MatProp[m]; // here
      }
      cout << endl;

      // } // end for loop // here
   
     //   here
    // // Print all identified mat prop 
    // ofstream alphaFile ("/u/home/l/luigiemp/project-cardio/voom2/src/Applications/PassiveLV/TestInVivoDir/alphas.dat");
    // for ( int s = 1; s <= NumPsteps; s+=DeltaSteps ) {
    //   for (int m = 0; m < NumMatProp; m++) {
    // 	cout      << MatPropList[s-1][m] << " ";
    // 	alphaFile << MatPropList[s-1][m] << " ";
    //   }
    //   cout << endl;
    //   alphaFile << endl;
    // }
    // alphaFile.close();
      
  } // End of material properties solution
*/

  // Timing
  time (&end);
  cout << endl << "BVP solved in " << difftime(end,start) << " s" << endl;


    
  return 0;
}



void writeNeumannBC(vector<int > BCid, vector<Real > Forces, string OutputFile, int step) {
  // Create outputFile name
  stringstream FileNameStream;
  FileNameStream << OutputFile << step << ".dat";
  ofstream out;
  out.open( (FileNameStream.str()).c_str() );
  
  out << BCid.size() << endl;
  for (uint i = 0; i < BCid.size(); i++) {
    out << BCid[i] << " " << setprecision(15) << Forces[i] << endl;
  }
  out.close();
}

