#include "PassMyoA.h"
#include "CompNeoHookean.h"
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

int main(int argc, char** argv)
{
  // Timing
  time_t start, end;
  time(&start);

  // Initialize Mesh
  FEMesh LVmesh("../CoarseLV.node", "../CoarseLV.ele");
  FEMesh LVsurf("../CoarseLV.node", "../CoarseLV.surf");
  string BCfile = "../CoarseLV.BaseSurfBC";
 
  cout << endl;
  cout << "Number Of Nodes   : " << LVmesh.getNumberOfNodes() << endl;
  cout << "Number Of Element : " << LVmesh.getNumberOfElements() << endl;
  cout << "Mesh Dimension    : " << LVmesh.getDimension() << endl << endl;
  
  // Initialize Materials and Fibers
  uint NumEl =  LVmesh.getNumberOfElements();
  vector<MechanicsMaterial * > materials;
  materials.reserve(NumEl);

  string FiberFile = "../CoarseLV.fiber";
  ifstream FiberInp(FiberFile.c_str());
  int NumMat = 1;
  
  // PassMyoA* Mat = new PassMyoA(0, 35.19, 7.06, 100.0, 0.025, 2.87, 2.82);
  PassMyoA* Mat = new PassMyoA(0, 35.19, 7.06, 0.0, 0.025, 2.87, 2.82);

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







    // Generate noise according to a specified SNR
    Real SNR = 10000.0;
    boost::mt19937 eng;
    boost::normal_distribution<double> dist(0.0, 1.0);
    boost::variate_generator< boost::mt19937, boost::normal_distribution<double> > gen(eng, dist);

    // Store undeformed nodal position
    int NumDoF = LVmesh.getNumberOfNodes() * 3;
    VectorXd X0(NumDoF);
    vector<VectorXd > noise( 10, VectorXd(NumDoF) );
    for (int n = 0; n < LVmesh.getNumberOfNodes(); n++) {
      X0(n*3) = LVmesh.getX(n, 0);
      X0(n*3 + 1) = LVmesh.getX(n, 1);
      X0(n*3 + 2) = LVmesh.getX(n, 2);
    }
    
    // Loop over time steps
    ind = 0;
    Real IM_Max = 1.2;
    for ( int s = 1; s <= NumPsteps; s+=DeltaSteps ) {
     
      // Read in displacements
      stringstream DispFileNameStream;
      DispFileNameStream << "LV_Fdisp_" << s << ".dat";
      ifstream DBCinp;
      DBCinp.open( (DispFileNameStream.str()).c_str() );
      Real temp = 0.0;
      DBCinp >> temp;

      // Add noise
      for (int n = 0; n < LVmesh.getNumberOfNodes()*3; n++) {
	// Compute displacement values
	DBCinp >> temp;
	// cout << temp << "  ";
	// Decompose displacement in real and imaginary component
	Real u = (temp - X0(n));
	// cout << u << " ";
	Real phi = M_PI * fabs(u/IM_Max);
	// cout << phi << " ";
	Real a = cos(phi), b = sin(phi);
	// cout << a << " " << b << endl;
	// Compute randomly distributed noise consisten with given SNR
	Real r1 = gen(), r2 = gen();
	if (fabs(r1) > 1.96) {r1 = 1.96;}; if (fabs(r2) > 1.96) {r2 = 1.96;}; 
	Real c = (1.0/sqrt(SNR*SNR-1.0))*r1;
	Real d = (1.0/sqrt(SNR*SNR-1.0))*r2;
	// Sum noise to signal and compute new phi
	a += c; b += d;
	// cout << a << " " << b << endl;
	Real phiNoise = atan(b/a);
	// cout << phiNoise <<  endl;
	if ( a < 0 ) { phiNoise += 0.5*M_PI; };
	
	// Transform back noise to displacement
	temp = phiNoise*IM_Max*double( (u/IM_Max > 0) - (u/IM_Max < 0) )/M_PI;
	// cout << X0(n) + temp << endl;
	noise[ind](n) = X0(n) + temp;
      }
      ind++;
      
    }
    






    // NR loop
    while (iter < NRmaxIter && error > NRtol)
    {
      myResults.setRequest(8);
      myModel.setResetFlag(1);

      ind = 0;
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
	Real temp = 0.0;
	for(int i = 0; i < DoF; i++) {
	  DoFid[i] = i;
	  DBCinp >> DoFvalues[i];
	  DoFvalues[i] +=  abs(DoFvalues[i] - noise[ind](i)); // 0.001*(double(rand())/double(RAND_MAX));
	  // if ( abs(DoFvalues[i] - noise[ind](i)) > 1.0e-3 ) {
 	  // DoFvalues[i] = noise[ind](i);
	  // cout << DoFvalues[i] << " " << noise[ind](i) << endl; }
	  myModel.setField(i, DoFvalues[i]);
	}
	ind++;
   
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

      // Compute condition number
      MatrixXd dHg;
      dHg = MatrixXd(*(myResults._Hg));
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

    // Print final configuration
    // myModel.writeOutputVTK("LV_F_", 31);
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

