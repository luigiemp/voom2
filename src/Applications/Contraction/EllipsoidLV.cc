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
//#include "EigenResult.h"
#include "MechanicsModel.h"
#include "EigenNRsolver.h"
// #include "LBFGSB.h"

using namespace voom;

int main(int argc, char** argv)
{

  cout << string(50, '\n'); // Clear Screen


  // Timing
  time_t start, end;
  time(&start);

  // Activation Sequence? Set True to make activation a function of z
  bool activationSequence = true;

  // Fibers across the wall thickness?
  bool fibersAcrossWall = true;

  // Conduction Velocity (cm/ms)
  double cv = 0.06;
  
  // Simulation Time (in ms)
  double simTime = 2000;  
  
  // Time Step (in ms)
  double deltaT = 0.01;

  // OutputString
  string outputString = "Ellipsoid_CV_3_";

  // Fiber Visualization String
  string fiberPlotString = "Ellipsoid_CV_Fiber.vtk";



  // Initialize Mesh
  // Assumptions to use this main as is: strip has a face at z=0; tetrahedral mesh
  FEMesh Cube("Mesh/Small_A.node", "Mesh/Small_A.ele");
  FEMesh surfMesh("Mesh/Small_A.node", "Mesh/Small_A.surfEle");
  string FiberFile = "Mesh/Small_A.fiber";
  string BCfile = "Mesh/Small_A.bc";
  ifstream FiberInp(FiberFile.c_str());

  double z_min = 0.0;
  for (int node_iter = 0; node_iter < Cube.getNumberOfNodes(); node_iter++)
  {
    VectorXd tempNodalPos = Cube.getX(node_iter);
    if(tempNodalPos[2] < z_min)
      z_min = tempNodalPos[2];
  }

  cout << endl;
  cout << "Number Of Nodes   : " << Cube.getNumberOfNodes() << endl;
  cout << "Number Of Element : " << Cube.getNumberOfElements() << endl;
  cout << "Mesh Dimension    : " << Cube.getDimension() << endl << endl;
  
  // Initialize Material
  uint NumMat =  Cube.getNumberOfElements();
  vector<GeomElement*> meshElements = Cube.getElements();
  
  vector<MechanicsMaterial * > PLmaterials;
  PLmaterials.reserve(NumMat);

  CompNeoHookean PassiveMat(0, 0.4, 0.04);
  CompNeoHookean ActiveMat(0, 0.4, 0.04);
  APForceVelPotential TestPotential(1.0, 500.0);
  // HillForceVelPotential TestPotential(4.4*pow(10,-3), .01*0.59, 25);
  // BlankPotential TestPotential;

  BlankViscousPotential ViscPotential;
  Vector3d HardParam(0.,0.,0.);

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

  // Read in Activation File:
  ifstream myfile;
  myfile.open ("ActivationFunction_1ms.dat");
  vector <double> Time(795, 0.0);
  vector <double> ActivationFactor(795, 0.0);
  for (int i = 0; i < 795; i++)
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
      activationTimesQP[el_iter] = (quadPointZ - z_min)/cv;
    else
      activationTimesQP[el_iter] = 0;  // Everything gets activated right away
      
    out << quadPointX << " " << quadPointY << " " << quadPointZ << endl;
    // cout << activationTimesQP[el_iter] << endl;

    vector <Vector3d> el_vectors(3, Vector3d::Zero(3,1));
      
    if (fibersAcrossWall)
    {
      // Load Fiber Data:
      FiberInp >> el_vectors[0][0]; FiberInp >> el_vectors[0][1]; FiberInp >> el_vectors[0][2];
      FiberInp >> el_vectors[1][0]; FiberInp >> el_vectors[1][1]; FiberInp >> el_vectors[1][2];
      FiberInp >> el_vectors[2][0]; FiberInp >> el_vectors[2][1]; FiberInp >> el_vectors[2][2];
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

  // Run Consistency check

  Real perturbationFactor = 0.1;
  myRequest = 7; // Check both Forces and Stiffness
  Real myH = 1e-6;
  Real myTol = 1e-7;

  // Before checking consistency, the perturbed deformation state must be 
  // set to the current deformation state.

  // myModel.checkConsistency(myResults, perturbationFactor, myRequest, myH, myTol);
  // myModel.checkDmat(myResults, perturbationFactor, myH, myTol);


    // Print initial configuration
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
    cout << "********" << " Setting up EBCs " << "********" << endl;
    int NumBC = 0, node = 0, ind = 0;;
    vector<int > BCnodes;
    vector<int > BCid;
    vector<Real > BCvalues;
    ifstream BCinp(BCfile.c_str());
   
    if (BCinp.is_open())
	cout << "BC File Opened Successfully!" << endl; 
    else
	cout << "ERROR: BC File Failed to Open!" << endl;
   
    BCinp >> NumBC;
    cout << "Number of Nodes with EBC: " << NumBC << endl;
    BCid.reserve(NumBC*3);
    BCvalues.reserve(NumBC*3);
    for(int i = 0; i < NumBC; i++) {
      BCinp >> node;
      BCnodes.push_back(node);
      for (int j = 0; j < 3; j++) {
        BCid.push_back(node*3 + j);
        BCvalues.push_back(Cube.getX(node, j));
        cout << BCid[ind] << " " <<  BCvalues[ind] << endl;
        ind++;
      }
    }
    

    // Lin tet cube - applied displacements //
    // for (uint i = 18; i < 27; i++) { 
    //   DoFid[i-5] = i*3 + 2; // Lower top in z direction by 0.1
    //   DoFvalues[i-5] = 1.95;
    // }
    
    // Solver
    Real NRtol = 1.0e-12;
    uint NRmaxIter = 100;
    EigenNRsolver mySolver(&myModel, BCid, BCvalues, CHOL, NRtol, NRmaxIter);

    ind = 0;
    for (int s = 0; s < simTime/deltaT; s++) {
      cout << "Step " << s << endl;
      // cout << "Activation Factor: " << ActivationFactor[s] << endl;
      
      if (s == 0)  // Free Contraction = 0, Fixed = 200
      {
	  ViscPotential.setViscosity(0.0E-5);
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
      myModel.writeOutputVTK(outputString, ind);
      // myModel.writeField("CubeSmall_", 1);
      
      /*
      // Check EBCs after Solve:
      // Get the Field from model
      vector <double > tempField(Cube.getNumberOfNodes() * 3, 0.0);
      myModel.getField(tempField);
      for (int node_iter = 0; node_iter < NumBC; node_iter++)
      {
        for (int dim_iter = 0; dim_iter < 3; dim_iter++)
          cout << BCnodes[node_iter] * 3 + dim_iter << "\t" << tempField[BCnodes[node_iter] * 3 + dim_iter] - Cube.getX(BCnodes[node_iter], dim_iter) << endl;
      }
      */
    }

  // Timing
  time (&end);
  cout << endl << "BVP solved in " << difftime(end,start) << " s" << endl;

  return 0;
}
