#include "PassMyoA.h"
#include "CompNeoHookean.h"
#include "PlasticMaterial.h"
#include "HillForceVelPotential.h"
#include "BlankPotential.h"
#include "APForceVelPotential.h"
#include "Potential.h"
#include "ViscousPotential.h"
#include "BlankViscousPotential.h"
#include "NewtonianViscousPotential.h"
#include "FEMesh.h"
#include "EigenEllipticResult.h"
#include "MechanicsModel.h"
#include "EigenNRsolver.h"
#include "LBFGSB.h"

using namespace voom;

int main(int argc, char** argv)
{

  // Timing
  time_t start, end;
  time(&start);

  // Conduction Velocity
  double cv = 0.05;

  // OutputString
  string outputString = "StripLVFlipped_attempt2";

  // Fibers Plotting String
  string fiberPlotString = "StripFiberLV_attempt2.vtk";


  // Initialize Mesh
  // Assumptions to use this main as is: strip has a face at z=0; tetrahedral mesh
  // FEMesh Cube("Cube6.node", "Cube6.ele");
  // FEMesh surfMesh("Cube6.node", "Cube6Surf.ele");
  // Real xmax = 1.0;
  // FEMesh Cube("Strip36.node", "Strip36.ele");
  // FEMesh surfMesh("Strip36.node", "Strip36Surf.ele");
  FEMesh Cube("Strip144.node", "Strip144.ele");
  FEMesh surfMesh("Strip144.node", "Strip144Surf.ele");
  Real xmax = 6.0;
  
  // Theta Values
  double theta_min = 60.0;
  double theta_max = -60.0;

  cout << endl;
  cout << "Number Of Nodes   : " << Cube.getNumberOfNodes() << endl;
  cout << "Number Of Element : " << Cube.getNumberOfElements() << endl;
  cout << "Mesh Dimension    : " << Cube.getDimension() << endl << endl;
  

    
  // Initialize Material
  uint NumMat =  Cube.getNumberOfElements();
  
  vector<MechanicsMaterial * > PLmaterials;
  PLmaterials.reserve(NumMat);
  
  vector<GeomElement*> cubeElements = Cube.getElements();
  cubeElements[0]->getNumberOfQuadPoints();

    ofstream out;

  // Fibers VTK File Headers:
  out.open(fiberPlotString.c_str());
  out << "# vtk DataFile Version 3.1" << endl;
  out << "Fiber vector representation" << endl;
  out << "ASCII" << endl;
  out << "DATASET UNSTRUCTURED_GRID" << endl;
  out << "POINTS " << NumMat << " FLOAT" << endl;

  vector <Vector3d> fiberVectors;
  vector <Vector3d> sheetVectors;
  vector <Vector3d> sheetNormalVectors;
  

  CompNeoHookean PassiveMat(0, 0.4, 0.04);
  CompNeoHookean ActiveMat(0, 40.0, 4.0);
  APForceVelPotential TestPotential(1.0, 500.0);
  // HillForceVelPotential TestPotential(4.4*pow(10,-3), .01*0.59, 25);
  // BlankPotential TestPotential;

  // BlankViscousPotential ViscPotential;
  NewtonianViscousPotential ViscPotential(0., 0.);

  vector <Vector3d> dirvec(3, Vector3d::Zero(3,1));
  dirvec[0] << 1., 0., 0.;
  dirvec[1] << 0., 1., 0.;
  dirvec[2] << 0., 0., 1.;

  Vector3d HardParam(0.,0.,0.);

  // Read in Activation File:
  ifstream myfile;
  myfile.open ("ActivationFunction_LargeTimestep.dat");
  vector <double> Time(795, 0.0);
  vector <double> ActivationFactor(795, 0.0);
  for (int i = 0; i < 795; i++)
  {
    myfile >> Time[i];
    myfile >> ActivationFactor[i];
  }
  myfile.close();


  for (int el_iter = 0; el_iter < cubeElements.size(); el_iter++)
  {
    int el_numQuadPoints = cubeElements[el_iter]->getNumberOfQuadPoints();
    int el_numNodes = cubeElements[el_iter]->getNodesPerElement();
    vector<int> el_nodeIds = cubeElements[el_iter]->getNodesID();
    Real quadPointX = 0.0; Real quadPointY = 0.0; Real quadPointZ = 0.0;

    // TODO: Works only for 1 Quad Point. Need to fix this.
    for (int el_node_iter = 0; el_node_iter < el_numNodes; el_node_iter++) {
      quadPointX += cubeElements[el_iter]->getN(0, el_node_iter) * Cube.getX(el_nodeIds[el_node_iter], 0);
      quadPointY += cubeElements[el_iter]->getN(0, el_node_iter) * Cube.getX(el_nodeIds[el_node_iter], 1);
      quadPointZ += cubeElements[el_iter]->getN(0, el_node_iter) * Cube.getX(el_nodeIds[el_node_iter], 2);
    }

    out << quadPointX << " " << quadPointY << " " << quadPointZ << endl;
    
    // Compute Theta based on x-position
    double el_theta = (theta_max - theta_min)/xmax * quadPointX + theta_min;
  
    // Set three directional vectors based on Theta.
    vector <Vector3d> el_vectors(3, Vector3d::Zero(3,1));
    el_vectors[0] << 0., cos(el_theta*3.14159 / 180.), sin(el_theta*3.14159 / 180.);
    el_vectors[1] << 0., -sin(el_theta*3.14159 / 180.), cos(el_theta*3.14159 / 180.);
    el_vectors[2] << 1., 0., 0.;

    fiberVectors.push_back(el_vectors[0]);
    sheetVectors.push_back(el_vectors[1]);
    sheetNormalVectors.push_back(el_vectors[2]);

    PlasticMaterial* PlMat = new PlasticMaterial(el_iter, &ActiveMat, &PassiveMat, &TestPotential, &ViscPotential);
    PlMat->setDirectionVectors(el_vectors);
    PlMat->setHardeningParameters(HardParam);
    PlMat->setActiveDeformationGradient(Matrix3d::Identity(3,3));
    PlMat->setTotalDeformationGradient(Matrix3d::Identity(3,3));

    PlMat->setTimestep(0.01);
    PlMat->setActivationMultiplier(0.0);

    PLmaterials.push_back(PlMat);
    // PLmaterials.push_back(&PassiveMat);
  }

  out << "POINT_DATA " << cubeElements.size() << endl;
  out << "VECTORS FiberDirection FLOAT" << endl;
  for (int el_iter = 0; el_iter < cubeElements.size(); el_iter++)
    out << fiberVectors[el_iter][0] << " " << fiberVectors[el_iter][1] << " " << fiberVectors[el_iter][2] << endl;

  out.close();

  // Initialize Model
  int NodeDoF = 3;
  int PressureFlag = 0;
  Real Pressure = 0.0;
  int NodalForcesFlag = 0;
  vector<int > ForcesID;
  vector<Real > Forces;
  MechanicsModel myModel(&Cube, PLmaterials, NodeDoF, PressureFlag, &surfMesh,
			 NodalForcesFlag);
  myModel.updatePressure(Pressure);
  myModel.updateNodalForces(&ForcesID, &Forces);
 
  // Initialize Result
  uint myRequest;
  uint PbDoF = (Cube.getNumberOfNodes())*myModel.getDoFperNode();
  EigenEllipticResult myResults(PbDoF, NumMat*2);

  // Run Consistency check

  Real perturbationFactor = 0.1;
  myRequest = 7; // Check both Forces and Stiffness
  Real myH = 1e-6;
  Real myTol = 1e-7;
  // myModel.checkConsistency(myResults, perturbationFactor, myRequest, myH, myTol);
  // myModel.checkDmat(myResults, perturbationFactor, myH, myTol);





    // Print initial configuration
    myModel.writeOutputVTK("CubeSmall_", 0);

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
    vector<int > DoFid;
    vector<Real > DoFvalues;
    vector<int > DoFxmax;
    int ind = 0;
    for(int i = 0; i < Cube.getNumberOfNodes(); i++) {
      // All nodes at x = xmax/2
      if ( Cube.getX(i, 0) < xmax/2.0 + 1.0e-12 && Cube.getX(i,0) > xmax/2.0 - 1.0e-12) {
	DoFid.push_back(i*3);
	ind++;
	// All nodes at y = 0
	if ( Cube.getX(i, 1) < 1.0e-12 ) {
	  DoFid.push_back(i*3 + 1);
	  ind++;
	}
	// All nodes at z = 0
	if ( Cube.getX(i, 2) < 1.0e-12 ) {
	  DoFid.push_back(i*3 + 2);
	  ind++;
	}
      }
      /*
      // All nodes at x = xmax
      if ( Cube.getX(i, 0) > xmax - 1.0e-12 ) {
	DoFid.push_back(i*3);
	DoFxmax.push_back(ind);
	ind++;
      }
     */
    }
 
    for(int i = 0; i < DoFid.size(); i++) {
      DoFvalues.push_back( Cube.getX(floor(double(DoFid[i])/3.0) ,DoFid[i]%3) );
      // cout << DoFid[i] << " " <<  DoFvalues[i] << endl;
    }

    // Lin tet cube - applied displacements //
    // for (uint i = 18; i < 27; i++) { 
    //   DoFid[i-5] = i*3 + 2; // Lower top in z direction by 0.1
    //   DoFvalues[i-5] = 1.95;
    // }
    
    // Solver
    Real NRtol = 1.0e-12;
    uint NRmaxIter = 100;
    EigenNRsolver mySolver(&myModel, DoFid, DoFvalues, CHOL, NRtol, NRmaxIter);

    ind = 0;
    int NumPassiveInc = 0;
    Real DeltaX = 0.2*xmax/NumPassiveInc;
    // cout << DoFxmax.size() << endl;
    for (int s = 0; s < NumPassiveInc; s++) {
      // Impose BC
      for (int i = 0; i < DoFxmax.size(); i++) {
	DoFvalues[DoFxmax[i]] += DeltaX;
      }
      ind++;
      mySolver.solve(DISP);
      for (int k = 0; k < NumMat; k++)
      {
	(PLmaterials[k])->updateStateVariables();
      }
      myModel.writeOutputVTK(outputString, ind);
    }

    for (int s = 0; s < 900; s++) {
      cout << "Step " << s << endl;
      cout << "Activation Factor: " << ActivationFactor[s] << endl;
      
      if (s == 0)  // Free Contraction = 0, Fixed = 200
      {
	for (int i = 0; i < DoFxmax.size(); i++)
	{
	  DoFid.erase(DoFid.begin() + DoFxmax[i] - i);
	  DoFvalues.erase(DoFvalues.begin() + DoFxmax[i] - i);
	  ViscPotential.setViscosity(1.0E-5);
	  // ViscPotential.setViscosity(0.0);
	  // High  = 1E-4, Low = 1E-5
	}
      }

      for (int k = 0; k < NumMat; k++)
      {
    	(PLmaterials[k])->setTimestep(0.5/1000);
	if(s >= 795 || ActivationFactor[s] < 0.0001)
    		(PLmaterials[k])->setActivationMultiplier(0.0);
	else
	    	(PLmaterials[k])->setActivationMultiplier(ActivationFactor[s]);
      }
      
      ind++;
      mySolver.solve(DISP);
     
      for (int k = 0; k < NumMat; k++)
      {
    	(PLmaterials[k])->updateStateVariables();
      }
      myModel.writeOutputVTK(outputString, ind);
      // myModel.writeField("CubeSmall_", 1);
    }

  // Timing
  time (&end);
  cout << endl << "BVP solved in " << difftime(end,start) << " s" << endl;



  return 0;
}
