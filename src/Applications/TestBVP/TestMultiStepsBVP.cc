#include "PassMyoA.h"
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
  FEMesh Cube("../../Mesh/Test/CubeQuad.node", "../../Mesh/Test/CubeQuad.ele");
  FEMesh surfMesh("../../Mesh/Test/CubeQuad.node", "../../Mesh/Test/SurfCubeQuad.ele");
 
  cout << endl;
  cout << "Number Of Nodes   : " << Cube.getNumberOfNodes() << endl;
  cout << "Number Of Element : " << Cube.getNumberOfElements() << endl;
  cout << "Mesh Dimension    : " << Cube.getDimension() << endl << endl;
  
    
  // Initialize Material
  uint NumMat =  Cube.getNumberOfElements();
  vector<MechanicsMaterial * > materials;
  materials.reserve(NumMat);
  vector<Vector3d > Fibers; 
  Fibers.reserve(NumMat);

  srand (time(NULL));
  cout << "Original material properties" << endl;
  for (int k = 0; k < NumMat; k++) {
    Real alpha1 = 35.19 *( 0.9 + 0.2*double(rand())/double(RAND_MAX) ),
         alpha2 = 7.06 *( 0.9 + 0.2*double(rand())/double(RAND_MAX) );
    cout << alpha1 << " " << alpha2 << endl;
    PassMyoA* Mat = new PassMyoA(k, alpha1, alpha2, 100.0, 0.025, 2.87, 2.82);
    materials.push_back(Mat);

    Vector3d N; N << 1.0/sqrt(3.0), 1.0/sqrt(3.0), 1.0/sqrt(3.0);
    Fibers.push_back(N);
  }
  cout << endl;

  // Load fibers into elements
  Cube.setFibers(Fibers);

  //// EBC
  vector<int > BCid(13,0);
  vector<Real > BCvalues(13, 0.0);
  BCid[0] = 0;  BCid[1] = 1;  BCid[2] = 2;  
  BCid[3] = 5;  BCid[4] = 8;  BCid[5] = 11;  
  BCid[6] = 26; BCid[7] = 29; BCid[8] = 35;          
  BCid[9] = 38; BCid[10] = 47;
  BCid[11] = 3; BCid[12] = 10;
 
  for(int i = 0; i < BCid.size(); i++) {
    BCvalues[i] = Cube.getX(floor(double(BCid[i])/3.0) ,BCid[i]%3);
  }

  // Initialize Model
  int NodeDoF = 3;
  int PressureFlag = 1;
  int NodalForcesFlag = 0;
  vector<vector< Real> > ForcesStorage;
  MechanicsModel myModel(&Cube, materials, NodeDoF, PressureFlag, &surfMesh,
			 NodalForcesFlag);
 
  // Initialize Result
  uint PbDoF = (Cube.getNumberOfNodes())*myModel.getDoFperNode();
  EigenEllipticResult myResults(PbDoF, NumMat*2);

  // Loop through pressure
  int NumPsteps = 10;
  Real DeltaP = 0.2;
  
  {
    // Solve for displacement first
    // Print initial configuration
    myModel.writeOutputVTK("CubeQuad_", 0);
      
    // Solver
    Real NRtol = 1.0e-12;
    uint NRmaxIter = 100;
    EigenNRsolver mySolver(&myModel, BCid, BCvalues, CHOL, NRtol, NRmaxIter);
    
    for (int s = 1; s <= NumPsteps; s++) {
      myModel.updatePressure(Real(s)*DeltaP);
      mySolver.solve(DISP); 

      myModel.writeOutputVTK("CubeQuad_", s);
      myModel.writeField("CubeQuad_", s);

      // Compute Residual
      myResults.setRequest(FORCE);
      myModel.compute(myResults);
      VectorXd F = *(myResults._residual);
      vector<Real > Forces;
      for (int i = 0; i < BCid.size(); i ++) {
	Forces.push_back( -F(BCid[i]) );
      }
      ForcesStorage.push_back(Forces);
      // writeNeumannBC(BCid, Forces, "LV_Fforces_", s);
    }
   
  }
  // End of displacement solution
  cout << endl << "End of displacement solution" << endl << endl << endl;

  // Next recompute material properties using EMFO
  {

    // Solver
    Real NRtol = 1.0e-12;
    uint NRmaxIter = 10;

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
    Real error = 1.0;
    uint iter = 0;
       
    // NR loop
    while (iter < NRmaxIter && error > NRtol)
    {
      myResults.setRequest(8);
      myModel.setResetFlag(1);

      for ( int s = 1; s <= NumPsteps; s+=3 ) {
	cout << "Step considered = " << s << endl;
	myModel.updateNodalForces(&BCid, &(ForcesStorage[s-1]));
	
	stringstream FileNameStream;
	FileNameStream << "CubeQuad_" << s << ".dat";
	ifstream BCinp;
	BCinp.open( (FileNameStream.str()).c_str() );
	
	int DoF;
	
	BCinp >> DoF;
	vector<int > DoFid(DoF, 0);
	vector<Real > DoFvalues(DoF, 0);
	for(int i = 0; i < DoF; i++) {
	  DoFid[i] = i;
	  BCinp >> DoFvalues[i];
	  // cout <<  DoFid[i] << " " <<  DoFvalues[i] << endl;
	}
	
	// Initialize field vector
	for (int i = 0; i < DoFid.size(); i++) {
	  myModel.setField(DoFid[i], DoFvalues[i]); // Set known displacements
	}
	
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
      myResults.setRequest(1);
      myModel.compute(myResults);
      cout << "Energy = " << myResults.getEnergy() << "   - NR iter = " << iter << "   -  NR error = " << error << endl;
      
    } // while loop  
      
    cout << endl << "Identified material properties" << endl;
    // Print found material properties
    // for (set<MechanicsMaterial *>::iterator MatIt = UNIQUEmaterials.begin(); MatIt != UNIQUEmaterials.end(); MatIt++) {
    //   int MatID = (*MatIt)->getMatID();
    //   vector<Real > MatProp = (*MatIt)->getMaterialParameters();
    //   int NumPropPerMat = MatProp.size();
    //   for (int m = 0; m < NumPropPerMat; m++) {
    // 	cout << MatProp[m] << " ";
    //   }
    //   cout << endl;
    // }
    int ind = 0;
    for (int k = 0; k < NumMat; k++) {
      vector<Real > MatProp = materials[k]->getMaterialParameters();
      int NumPropPerMat = MatProp.size();
      for (int m = 0; m < NumPropPerMat; m++) {
	cout << MatProp[m] << " ";
	ErrorInMatProp[ind] -= MatProp[m];
	ind++;
      }
      cout << endl;
    }
    // myModel.writeOutputVTK("CubeQuad_", 3);

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
