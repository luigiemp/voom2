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
#include "MechanicsModel.h"
#include "EigenNRsolver.h"

using namespace voom;

int main(int argc, char** argv)
{
  cout << string(50, '\n'); // Clear Screen

  // Timing
  time_t start, end;
  time(&start);

  // Activation Sequence? Set True to make activation a function of x
  bool activationSequence = false;

  // Fibers across the wall thickness in Strip? Set True if non-unidirectional fibers
  bool fibersAcrossWall = false;
  double theta_min = 60.0;
  double theta_max = -60.0;

  // Conduction Velocity (cm/ms)
  double cv = 0.06;
  
  // Simulation Time (in ms)
  double simTime = 10;  
  
  // Time Step (in ms)
  double deltaT = 0.01;

  // writeOutputFlag set to true for output.
  bool writeOutputFlag = false;

  // OutputString
  string outputString = "Output/Cube/NoViscosity/Cube_NoViscosity_";

  // string outputString = "SingleElementTest_";

  // Fiber Visualization String
  string fiberPlotString = outputString + "Fiber.vtk";
  // string fiberPlotString = "SingElementTest_Fiber.vtk";

  // Initialize Mesh
  // Assumptions to use this main as is: strip has a face at z=0; tetrahedral mesh
  FEMesh Cube("Mesh/Cube1.node", "Mesh/Cube1.ele");
  FEMesh surfMesh("Mesh/Cube1.node", "Mesh/Cube1Surf.ele");
  // FEMesh Cube("Mesh/Cube6.node", "Mesh/Cube6.ele");
  // FEMesh surfMesh("Mesh/Cube6.node", "Mesh/Cube6Surf.ele");
  Real xmax = 1.0;
  // FEMesh Cube("Mesh/Strip36.node", "Mesh/Strip36.ele");
  // FEMesh surfMesh("Mesh/Strip36.node", "Mesh/Strip36Surf.ele");
  // Real xmax = 1.0;
  // FEMesh Cube("Mesh/Strip144.node", "Mesh/Strip144.ele");
  // FEMesh surfMesh("Mesh/Strip144.node", "Mesh/Strip144Surf.ele");
  // Real xmax = 6.0; 


  // ***************************************** //
  // **************Begin Computing************ //
  // ***************************************** //
  cout << endl;
  cout << "\033[1;32mNumber Of Nodes \t : \033[0m" << Cube.getNumberOfNodes() << endl;
  cout << "\033[1;32mNumber Of Element \t : \033[0m" << Cube.getNumberOfElements() << endl;
  cout << "\033[1;32mMesh Dimension \t\t : \033[0m" << Cube.getDimension() << endl << endl;
  
  
   
  // Initialize Material
  uint NumMat =  Cube.getNumberOfElements();
  vector<GeomElement*> meshElements = Cube.getElements();
  
  vector<MechanicsMaterial * > PLmaterials;
  PLmaterials.reserve(NumMat);

  CompNeoHookean PassiveMat(0, 4.0, 0.4);
  CompNeoHookean ActiveMat(0, 4.0, 0.4);
  // APForceVelPotential TestPotential(1.0, 500.0);
  // HillForceVelPotential TestPotential(4.4*pow(10,-3), .01*0.59, 25);
  // BlankPotential TestPotential;

  BlankViscousPotential ViscPotential;
  // NewtonianViscousPotential ViscPotential(0.5, 0.5);

  // Visualize Fiber directions
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

  Vector3d HardParam(0.,0.,0.);

  // Read in Activation File:
  ifstream myfile;
  myfile.open ("ActivationFunction_1ms.dat");
  vector <double> Time(398, 0.0);
  vector <double> ActivationFactor(398, 0.0);
  for (int i = 0; i < 398; i++)
  {
    myfile >> Time[i];
    myfile >> ActivationFactor[i];
  }
  myfile.close();

  deltaT = Time[1] - Time[0];
  
  // Cycle Length (in ms)
  double cycleLength = Time[397];

  // Calculate the activation time for each quadrature point in each element:
  vector <double> activationTimesQP(NumMat, 0.0);

  for (int el_iter = 0; el_iter < meshElements.size(); el_iter++) {
    int el_numQuadPoints = meshElements[el_iter]->getNumberOfQuadPoints();
    int el_numNodes = meshElements[el_iter]->getNodesPerElement();
    vector<int> el_nodeIds = meshElements[el_iter]->getNodesID();
    Real quadPointX = 0.0; Real quadPointY = 0.0; Real quadPointZ = 0.0;
      
    // TODO: Works only for 1 Quad Point. Need to fix this.
    for (int el_node_iter = 0; el_node_iter < el_numNodes; el_node_iter++) {
      quadPointX += meshElements[el_iter]->getN(0, el_node_iter) * Cube.getX(el_nodeIds[el_node_iter], 0);
      quadPointY += meshElements[el_iter]->getN(0, el_node_iter) * Cube.getX(el_nodeIds[el_node_iter], 1);
      quadPointZ += meshElements[el_iter]->getN(0, el_node_iter) * Cube.getX(el_nodeIds[el_node_iter], 2);
    }
    
    // Activation Time Calculations:
    if(activationSequence)
      activationTimesQP[el_iter] = quadPointX/cv;
    else
      activationTimesQP[el_iter] = 0;  // Everything gets activated right away
      
    out << quadPointX << " " << quadPointY << " " << quadPointZ << endl;

    vector <Vector3d> el_vectors(3, Vector3d::Zero(3,1));
      
    if (fibersAcrossWall)
    {
      // Compute Theta based on x-position
      double el_theta = (theta_max - theta_min)/xmax * quadPointX + theta_min;
      
      // Set three directional vectors based on Theta.
      el_vectors[0] << 0., cos(el_theta*3.14159 / 180.), sin(el_theta*3.14159 / 180.);
      el_vectors[1] << 0., -sin(el_theta*3.14159 / 180.), cos(el_theta*3.14159 / 180.);
      el_vectors[2] << 1., 0., 0.;
    }
    else
    {
      el_vectors[0] << 1., 0., 0.;
      el_vectors[1] << 0., 1., 0.;
      el_vectors[2] << 0., 0., 1.;
    }
    fiberVectors.push_back(el_vectors[0]);
    sheetVectors.push_back(el_vectors[1]);
    sheetNormalVectors.push_back(el_vectors[2]);

    // PLmaterials.push_back(&PassiveMat);
    APForceVelPotential* TestPotential = new APForceVelPotential(1.0, 500.0);
    PlasticMaterial* PlMat = new PlasticMaterial(el_iter, &ActiveMat, &PassiveMat, TestPotential, &ViscPotential);
    PlMat->setDirectionVectors(el_vectors);
    PlMat->setHardeningParameters(HardParam);
    PlMat->setActiveDeformationGradient(Matrix3d::Identity(3,3));
    PlMat->setTotalDeformationGradient(Matrix3d::Identity(3,3));

    PlMat->setTimestep(deltaT);
    PlMat->setActivationMultiplier(0.0);

    PLmaterials.push_back(PlMat);
  }

  // Finish vtk file for plotting fiber data:
  out << "POINT_DATA " << meshElements.size() << endl;
  out << "VECTORS FiberDirection FLOAT" << endl;
  for (int el_iter = 0; el_iter < meshElements.size(); el_iter++)
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
  EigenResult myResults(PbDoF, NumMat*2);


  // Print initial configuration
  if (writeOutputFlag)
  	myModel.writeOutputVTK(outputString, 0);

  // Check on applied pressure
  myRequest = 2;
  myResults.setRequest(myRequest);
  myModel.compute(&myResults);
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
    // All nodes at x = 0
    if ( Cube.getX(i, 0) < 1.0e-12 ) {
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

    // All nodes at x = xmax
    if ( Cube.getX(i, 0) > xmax - 1.0e-12 ) {
      DoFid.push_back(i*3);
      DoFxmax.push_back(ind);
      ind++;
    }
    
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
  int NumPassiveInc = 00;
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
    if (writeOutputFlag)
      myModel.writeOutputVTK(outputString, ind);
  }

  for (int s = 0; s < simTime/deltaT; s++) {
    cout << endl << "============================" << endl;
    cout << "Step " << s << endl;
    // cout << "Activation Factor: " << ActivationFactor[s] << endl;
      
    if (s == 0)  // Free Contraction = 0, Fixed = 200
      {
	for (int i = 0; i < DoFxmax.size(); i++)
	  {
	    DoFid.erase(DoFid.begin() + DoFxmax[i] - i);
	    DoFvalues.erase(DoFvalues.begin() + DoFxmax[i] - i);
	    // ViscPotential.setViscosity(0.0E-5);
	    // High  = 1E-4, Low = 1E-5
	  }
      }

    for (int k = 0; k < NumMat; k++)
      {
    	(PLmaterials[k])->setTimestep(deltaT/1000);
	
	if (s * deltaT < activationTimesQP[k] || s * deltaT > activationTimesQP[k] + cycleLength)
	  (PLmaterials[k])->setActivationMultiplier(0.0);
	else
	  {
	    // Figure out time in cycle
	    double tempNormalizedTime = s * deltaT - activationTimesQP[k];
	    (PLmaterials[k])->setActivationMultiplier(ActivationFactor[tempNormalizedTime/deltaT]);
	  }
      }
      
    ind++;
    mySolver.solve(DISP);
     
    for (int k = 0; k < NumMat; k++)
      {
    	(PLmaterials[k])->updateStateVariables();
      }
    if (writeOutputFlag)
    	myModel.writeOutputVTK(outputString, ind);
    // myModel.writeField("CubeSmall_", 1);
  }

  // Timing
  time (&end);
  cout << endl << "BVP solved in " << difftime(end,start) << " s" << endl;



  return 0;
}
