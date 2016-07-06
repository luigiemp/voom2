#include "PassMyoA.h"
#include "CompNeoHookean.h"
#include "PlasticMaterial.h"
#include "HillForceVelPotential.h"
#include "BlankPotential.h"
#include "APForceVelPotential.h"
#include "Humphrey_Compressible.h"
#include "LinYinActive_Compressible.h"
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

  // Activation Sequence? Set True to make activation occur with activation Time
  bool activationSequence = true;

  bool fibersAcrossWall = true;

  // Simulation Time (in ms)
  double simTime = 2000;  
  
  // Time Step (in ms)
  double deltaT = 0.01;

  // OutputString
  string outputString = "Output/RabbitLV_Test/RabbitLV";

  // Fiber Visualization String
  string fiberPlotString = outputString + "Fiber.vtk";

  // Spring BC:
  int SpringBCflag = 1;

  // Spring Stiffnes
  Real SpringK = 0.0e-6;

  // Initialize Mesh
  // Assumptions to use this main as is: strip has a face at z=0; tetrahedral mesh
  FEMesh Cube("Mesh/RabbitLV/RabbitLV.node", "Mesh/RabbitLV/RabbitLV.ele");
  FEMesh surfMesh("Mesh/RabbitLV/RabbitLV.node", "Mesh/RabbitLV/RabbitLV.springSurfEle");
  string FiberFile = "Mesh/RabbitLV/RabbitLV.fiber";
  string BCfile = "Mesh/RabbitLV/RabbitLV.apexNodeBC";
  string ActivationTimeFile = "Mesh/RabbitLV/RabbitLV.activationTime";   // This is the Element Activation Time File
  string ActivationFile = "ActivationFunction_1ms.dat";  // This is the Calcium Transient
  // string BCfile = "Mesh/Small_A.bc";
  ifstream FiberInp(FiberFile.c_str());
  ifstream ActTime(ActivationTimeFile.c_str());

  cout << endl;
  cout << "Number Of Nodes   : " << Cube.getNumberOfNodes() << endl;
  cout << "Number Of Element : " << Cube.getNumberOfElements() << endl;
  cout << "Mesh Dimension    : " << Cube.getDimension() << endl << endl;
  
  // Initialize Material
  uint NumMat =  Cube.getNumberOfElements();
  vector<GeomElement*> meshElements = Cube.getElements();
  
  vector<MechanicsMaterial * > PLmaterials;
  PLmaterials.reserve(NumMat);

  Humphrey_Compressible PassiveMat(0, 15.98, 55.85, 0.0, -33.27, 30.21, 3.590, 64.62);
  LinYinActive_Compressible ActiveMat(0, -38.70, 40.83, 25.12, 9.51, 171.18);

  BlankViscousPotential ViscPotential;
  // NewtonianViscousPotential ViscPotential(0.005, 0.5);
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
  myfile.open (ActivationFile.c_str());
  double tempTime = 0.0;
  double tempActivationFactor = 0.0;

  vector <double> Time;
  vector <double> ActivationFactor;

  while(myfile >> tempTime >> tempActivationFactor)
  {
    Time.push_back(tempTime);
    ActivationFactor.push_back(tempActivationFactor);
  }

  double cycleLength = Time[Time.size() - 1];
  myfile.close();
  deltaT = Time[1] - Time[0];

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
      ActTime >> activationTimesQP[el_iter];
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
    APForceVelPotential* TestPotential = new APForceVelPotential(1.0, 50000.0);
    PlasticMaterial* PlMat = new PlasticMaterial(el_iter, &ActiveMat, &PassiveMat, TestPotential, &ViscPotential);
    PlMat->setDirectionVectors(el_vectors);
    PlMat->setHardeningParameters(HardParam);
    PlMat->setActiveDeformationGradient(Matrix3d::Identity(3,3));
    PlMat->setTotalDeformationGradient(Matrix3d::Identity(3,3));

    PlMat->setTimestep(deltaT/1000.);
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
			 NodalForcesFlag, SpringBCflag);
  myModel.updatePressure(Pressure);
  myModel.updateNodalForces(&ForcesID, &Forces);

  myModel.initSpringBC("Mesh/RabbitLV/RabbitLV.springSurf", &surfMesh, SpringK);
 
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

  // EBC
  cout << "********" << " Setting up EBCs " << "********" << endl;
  int NumBC = 0, BCidtemp  = 0, ind = 0;;
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
    BCinp >> BCidtemp;
    // BCnodes.push_back(node);
    // for (int j = 0; j < 3; j++) {
    BCid.push_back(BCidtemp);
    BCvalues.push_back(Cube.getX(floor(BCidtemp/3), BCidtemp%3));
    cout << BCid[ind] << " " <<  BCvalues[ind] << endl;
    ind++;
    // }
  }

  // Solver
  Real NRtol = 1.0e-8;
  uint NRmaxIter = 100;
  EigenNRsolver mySolver(&myModel, BCid, BCvalues, CHOL, NRtol, NRmaxIter);

  ind = 0;
  for (int s = 0; s < simTime/deltaT; s++) {
    cout << "Step " << s << endl;

    for (int k = 0; k < NumMat; k++)
      {
    	(PLmaterials[k])->setTimestep(deltaT/1000);
	
	if (s * deltaT < activationTimesQP[k] || s * deltaT > activationTimesQP[k] + cycleLength)
	  (PLmaterials[k])->setActivationMultiplier(0.0);
	else
	{
	  // cout << "El: " << k << "\t" << s*deltaT << "<" <<activationTimesQP[k] << endl;
	  // Figure out time in cycle
	  double tempNormalizedTime = s * deltaT - activationTimesQP[k];
	  (PLmaterials[k])->setActivationMultiplier(ActivationFactor[tempNormalizedTime/deltaT]);
	  // (PLmaterials[k])->setActivationMultiplier(0.0);
	}
      }
      
      ind++;
      mySolver.solve(DISP);
      myModel.setPrevField();
     
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
