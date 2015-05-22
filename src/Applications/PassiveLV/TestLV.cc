#include "PassMyoA.h"
#include "CompNeoHookean.h"
#include "FEMesh.h"
#include "EigenEllipticResult.h"
#include "MechanicsModel.h"
#include "EigenNRsolver.h"

using namespace voom;

void writeNeumannBC(vector<int > BCid, vector<Real > Forces, string OutputFile, int step);

int main(int argc, char** argv)
{
  // Timing
  time_t start, end;
  time(&start);

  // Initialize Mesh
  FEMesh LVmesh("../../Mesh/Test/CoarseLV.node", "../../Mesh/Test/CoarseLV.ele");
  FEMesh LVsurf("../../Mesh/Test/CoarseLV.node", "../../Mesh/Test/CoarseLV.surf");
  string BCfile = "CoarseLV.BaseSurfBC";
 
  cout << endl;
  cout << "Number Of Nodes   : " << LVmesh.getNumberOfNodes() << endl;
  cout << "Number Of Element : " << LVmesh.getNumberOfElements() << endl;
  cout << "Mesh Dimension    : " << LVmesh.getDimension() << endl << endl;
  
  // Initialize Materials and Fibers
  uint NumEl =  LVmesh.getNumberOfElements();
  vector<MechanicsMaterial * > materials;
  materials.reserve(NumEl);

  string FiberFile = "CoarseLV.fiber";
  ifstream FiberInp(FiberFile.c_str());
  int NumMat = 1;
  
  PassMyoA* Mat = new PassMyoA(0, 35.19, 7.06, 100.0, 0.025, 2.87, 2.82);

  vector<Vector3d > Fibers;
  Fibers.reserve(NumEl);
  for (int k = 0; k < NumEl; k++) {
    materials.push_back(Mat);
    Vector3d N = Vector3d::Zero();
    FiberInp >> N(0);
    FiberInp >> N(1);
    FiberInp >> N(2);
    Fibers.push_back(N);

    // materials.push_back(new CompNeoHookean(k, 10.0, 10.0) );
  }
  LVmesh.setFibers(Fibers);



  // Initialize Model
  int NodeDoF = 3;
  int PressureFlag = 1;
  int NodalForcesFlag = 0;
  MechanicsModel myModel(&LVmesh, materials, NodeDoF, PressureFlag, &LVsurf,
			 NodalForcesFlag);
 
  // Initialize Result
  uint PbDoF = (LVmesh.getNumberOfNodes())*myModel.getDoFperNode();
  EigenEllipticResult myResults(PbDoF, NumMat*2);
 
  // Print initial configuration
  myModel.writeOutputVTK("LV_F_", 0);



  // Solve for displacement first
  // Check on applied pressure
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
  int NumPsteps = 30;
  Real DeltaP = 0.05;
  /*  
  {
    // Solve for displacement first
    // Solver
    Real NRtol = 1.0e-10;
    uint NRmaxIter = 100;
    EigenNRsolver mySolver(&myModel, BCid, BCvalues, CHOL, NRtol, NRmaxIter);  

    for (int s = 1; s <= NumPsteps; s++) {
      myModel.updatePressure(Real(s)*DeltaP);
      mySolver.solve(DISP); 

      myModel.writeOutputVTK("LV_F_", s);
      myModel.writeField("LV_Fdisp_", s);
      
      // Compute Residual
      myResults.setRequest(FORCE);
      myModel.compute(myResults);
      VectorXd F = *(myResults._residual);
      vector<Real > Forces(BCid.size(), 0.0);
      for (int i = 0; i < BCid.size(); i ++) {
	Forces[i] = -F(BCid[i]);
      }
      writeNeumannBC(BCid, Forces, "LV_Fforces_", s);
    }
  }
  // End of displacement solution
  cout << endl << "End of displacement solution" << endl << endl << endl;
*/
 

  // Next recompute material properties using EMFO
  {
    // Change Material properties
    // Find unique material parameters
    set<MechanicsMaterial *> UNIQUEmaterials;
    for (uint i = 0; i < materials.size(); i++) 
      UNIQUEmaterials.insert(materials[i]);

    vector<Real > ErrorInMatProp;
    cout << endl << "Original material properties" << endl;
    for (int k = 0; k < NumMat; k++) {
      vector<Real > MatProp = materials[k]->getMaterialParameters();
      int NumPropPerMat = MatProp.size();
      for (int m = 0; m < NumPropPerMat; m++) {
	cout << MatProp[m] << " ";
	ErrorInMatProp.push_back(MatProp[m]);
      }
      cout << endl;
    }

    cout << endl << endl << "Initial guess for material properties" << endl;
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
    cout << endl << "Material properties identification." << endl;

    // Need to apply all external forces: pressure + any Neumann BC
    myModel.setNodalForcesFlag(1);

    // Initialize solver parameters
    Real NRtol = 1.0e-12;
    uint NRmaxIter = 10;
    Real error = 1.0;
    uint iter = 0;
    uint DeltaSteps = 3;

    // NR loop
    while (iter < NRmaxIter && error > NRtol)
    {
      myResults.setRequest(8);
      myModel.setResetFlag(1);

      for ( int s = 1; s <= NumPsteps; s+=DeltaSteps ) {
	
	cout << "Step considered = " << s << endl;

	int DoF;

	// Read in displacements
	stringstream DispFileNameStream;
	DispFileNameStream << "LV_Fdisp_" << s << ".dat";
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
	ForceFileNameStream << "LV_Fforces_" << s << ".dat";
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
      }

      VectorXd DeltaAlpha;
      SimplicialCholesky<SparseMatrix<Real > > chol(*(myResults._Hg));  // performs a Cholesky factorization of _stiffness
      DeltaAlpha = chol.solve(*(myResults._Gradg));                             // use the factorization to solve for the given right hand side

      // Update material properties in material
      for (set<MechanicsMaterial *>::iterator MatIt = UNIQUEmaterials.begin(); MatIt != UNIQUEmaterials.end(); MatIt++) {
	int MatID = (*MatIt)->getMatID();
	vector<Real > MatProp = (*MatIt)->getMaterialParameters();
	int NumPropPerMat = MatProp.size();
	for (int m = 0; m < NumPropPerMat; m++) {
	  MatProp[m] -= DeltaAlpha(MatID*NumPropPerMat + m);
	}
	(*MatIt)->setMaterialParameters(MatProp);
      }

      // Update iter and error
      iter++;
      error = DeltaAlpha.norm();
      myResults.setRequest(ENERGY);
      myModel.compute(myResults);
      cout << "Energy = " << myResults.getEnergy() << "   - NR iter = " << iter << "   -  NR error = " << error << endl;
      
    } // while loop  
   
    cout << endl << "Identified material properties" << endl;

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
    
    int ind = 0;
    for (int k = 0; k < NumMat; k++) {
      vector<Real > MatProp = materials[k]->getMaterialParameters();
      int NumPropPerMat = MatProp.size();
      for (int m = 0; m < NumPropPerMat; m++) {
	ErrorInMatProp[ind] -= MatProp[m];
	ind++;
      }
    }

    Real MatPropError = 0.0;
    for (int i = 0; i < ErrorInMatProp.size(); i++) {
      MatPropError += ErrorInMatProp[i]*ErrorInMatProp[i];
    }

    cout << "L2 norm of Mat Prop Error = " << pow(MatPropError, 0.5) << endl;
  } // End of material properties solution



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

