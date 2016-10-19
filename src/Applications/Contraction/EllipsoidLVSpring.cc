#include "PassMyoA.h"
#include "CompNeoHookean.h"
#include "PlasticMaterial.h"
#include "HillForceVelPotential.h"
#include "BlankPotential.h"
#include "APForceVelPotential.h"
#include "APForceVelPotential_Test.h"
#include "Humphrey_Compressible.h"
#include "LinYinActive_Compressible.h"
#include "Potential.h"
#include "ViscousPotential.h"
#include "BlankViscousPotential.h"
#include "NewtonianViscousPotential.h"
#include "FEMesh.h"
#include "Jacobian.h"
#include "MechanicsModel.h"
#include "EigenNRsolver.h"

using namespace voom;

double calculateEjectionFraction(MechanicsModel* cavityModel, const MechanicsModel* myocardiumModel, const vector<int> surfaceNodes, const vector<double> currentMyocardiumField);

int main(int argc, char** argv)
{
  cout << string(50, '\n'); // Clear Screen

  // Timing
  time_t start, end;
  time(&start);

  // Activation Sequence? Set True to make activation a function of z
  bool activationSequence = true;

  // Constant Activation?
  bool constantActivation = false;
  double activationValue = 0.0;

  // Fibers across the wall thickness?
  bool fibersAcrossWall = true;

  // Calculate Ejection Fraction?
  bool calculateEjectionFractionFlag = false;

  // Use Conduction Velocity (if false - uses activation Time file)
  bool useConductionVelocity = false;

  // Conduction Velocity (cm/ms)
  double cv = 0.06;

  // Pressure File
  bool pressureFlag = true;
  string pressureFile = "Pressure_1ms.dat";
  
  // Simulation Time (in ms)
  double simTime = 2000;  
  
  // Time Step (in ms)
  double deltaT = 0.01;

  // OutputString
  string outputString = "/u/project/cardio/adityapo/ScratchResults/EllipsoidSpring2/Ellipsoid";

  // Fiber Visualization String
  string fiberPlotString = outputString + "_Fiber.vtk";

  // Printing Volume to File
  string volumeFile = outputString + "_Volume.txt";

  // Spring BC:
  int SpringBCflag = 0;

  // Spring Stiffnes
  Real SpringK = 0.0; // 1.0e4;

  // Initialize Mesh
  // Assumptions to use this main as is: strip has a face at z=0; tetrahedral mesh
  FEMesh Cube("Mesh/EllipsoidMesh/Small_A.node", "Mesh/EllipsoidMesh/Small_A.ele");
  FEMesh surfMesh("Mesh/EllipsoidMesh/Small_A.node", "Mesh/EllipsoidMesh/Small_A.outerSurfEle");
  FEMesh innerSurfMesh("Mesh/EllipsoidMesh/Small_A.node", "Mesh/EllipsoidMesh/Small_A.innerSurfEle");
  string FiberFile = "Mesh/EllipsoidMesh/Small_A.fiber";
  string BCfile = "Mesh/EllipsoidMesh/Small_A.bc";
  string ActivationTimeFile = "Mesh/EllipsoidMeshFiner/Small_B.activationTime";   // This is the Element Activation Time File
  string ActivationFile = "ActFunc_600ms_1msInterval.dat";  // This is the Calcium Transient
  ifstream FiberInp(FiberFile.c_str());
  ifstream ActTime(ActivationTimeFile.c_str());

  // Code for calculating Ejection Fraction
  MechanicsModel* cavityModel = NULL;
  vector<MechanicsMaterial * > cavityMaterials;
  FEMesh* cavityMesh;
  vector <int> surfaceNodes;

  if (calculateEjectionFractionFlag)
  {
    // Setup mesh for calculating ejection fraction
    cavityMesh = new FEMesh("Mesh/EllipsoidMesh/Small_A_Cavity.node", 
      		            "Mesh/EllipsoidMesh/Small_A_Cavity.ele");
    string outerSurfFile = "Mesh/EllipsoidMesh/Small_A.innerSurfNode";

    ifstream cavityNodeFileStream (outerSurfFile.c_str());

    if (cavityNodeFileStream.is_open())
	cout << "** Opened Cavity Surface File Successfully." << endl;
    else
    {
	cout << "** Cannot Open Cavity Surface File." << endl;
	exit(1);
    }

    double numSurfaceNodes;
    cavityNodeFileStream >> numSurfaceNodes;

    for (int i = 0; i < numSurfaceNodes; i++)
    {
      double tempSurfaceCavityNode;
      cavityNodeFileStream >> tempSurfaceCavityNode;
      surfaceNodes.push_back(tempSurfaceCavityNode);
    }
    cavityNodeFileStream.close();

    // Initialize Cavity Model
    uint cavityNodeDoF = 3;
    uint cavityNumQP = 1;
  
    uint NumMatCavity = cavityMesh->getNumberOfElements()*cavityNumQP;
    cavityMaterials.reserve(NumMatCavity);
  
    for (int k = 0; k < NumMatCavity; k++) {
      Jacobian* Mat = new Jacobian(k);
      cavityMaterials.push_back(Mat);
    }
  
    cavityModel = new MechanicsModel(cavityMesh, cavityMaterials, cavityNodeDoF);
  }

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


  APForceVelPotential TestPotential(4.0, 0.0, 3.0);	// 50.0 for 2nd parameter, force
  // APForceVelPotential_Test TestPotential(4.0, 500000.0);

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

  // Read in Pressure File:
  double tempTime = 0.0;
  double tempPressure = 0.0;

  vector <vector <double> > pressureData;

  if (pressureFlag)
  {
    ifstream pressureFileStream;
    pressureFileStream.open(pressureFile.c_str());
    
    if (pressureFileStream.is_open())
    {
      cout << "** Opened Pressure File Successfuly." << endl;
      while (pressureFileStream >> tempTime >> tempPressure) {
	vector <double> tempTimePressure;
	tempTimePressure.push_back(tempTime);
	tempTimePressure.push_back(tempPressure);
	pressureData.push_back(tempTimePressure);
      }
    }
    else {
      cout << "** Failed to Open Pressure File." << endl;
      exit(1);
    }
    
  }

  // Read in Activation File:
  ifstream myfile;
  myfile.open (ActivationFile.c_str());
  tempTime = 0.0;
  double tempActivationFactor = 0.0;
  
  vector <double> Time;
  vector <double> ActivationFactor;

  while(myfile >> tempTime >> tempActivationFactor)
  {
    Time.push_back(tempTime);
    if (!constantActivation)
      ActivationFactor.push_back(tempActivationFactor);
    else
      ActivationFactor.push_back(activationValue);
  }
  
  double cycleLength = Time[Time.size() - 1];
  myfile.close();
  deltaT = Time[1] - Time[0];
  

  // Calculate the activation time for each quadrature point in each element:
  vector <double> activationTimesQP(NumMat, 0.0);

  if (!useConductionVelocity)
  {
    if (ActTime.is_open())
        cout << "** Opened Time Activation File Successfully." << endl;
    else
    {
        cout << "** Cannot Open Time Activation File." << endl;
        exit(1);
    }
  }



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
    if(activationSequence) {
      if (useConductionVelocity)
        activationTimesQP[el_iter] = (quadPointZ - z_min)/cv;
      else
        ActTime >> activationTimesQP[el_iter];
    }
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

    Humphrey_Compressible* PassiveMat = new Humphrey_Compressible(0, 15.98, 55.85, 0.0, -33.27, 30.21, 3.590, 64.62, el_vectors);
    LinYinActive_Compressible* ActiveMat = new LinYinActive_Compressible(0, -38.70, 40.83, 25.12, 9.51, 171.18, el_vectors);

    PlasticMaterial* PlMat = new PlasticMaterial(el_iter, ActiveMat, PassiveMat, &TestPotential, &ViscPotential);
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
  ActTime.close();

  // Initialize Model
  int NodeDoF = 3;
  int PressureFlag = 1;
  Real Pressure = 0.0;
  int NodalForcesFlag = 0;
  vector<int > ForcesID;
  vector<Real > Forces;
  MechanicsModel myModel(&Cube, PLmaterials, NodeDoF, PressureFlag, &innerSurfMesh,
			 NodalForcesFlag, SpringBCflag);
  myModel.updatePressure(Pressure);
  myModel.updateNodalForces(&ForcesID, &Forces);

  myModel.initSpringBC("Mesh/EllipsoidMesh/Small_A.outerSurf", &surfMesh, SpringK);
 
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

  // Solver
  Real NRtol = 1.0e-5;
  uint NRmaxIter = 100;
  EigenNRsolver mySolver(&myModel, BCid, BCvalues, CHOL, NRtol, NRmaxIter);

  
  // SOLVE:

  ind = 0;

  ofstream outVolume;
  outVolume.open(volumeFile.c_str());

  for (int s = 0; s < simTime/deltaT; s++)
  {
    cout << "Step " << s << endl;
    
    // Update pressure	
    myModel.updatePressure(-1.0/20 * pressureData[s][1]);
    // myModel.updatePressure(-1.0 * Real(s)/15);

    for (int k = 0; k < NumMat; k++)
    {
      (PLmaterials[k])->setTimestep(deltaT/1000);
	
      if (s * deltaT < activationTimesQP[k] || s * deltaT > activationTimesQP[k] + cycleLength)
        (PLmaterials[k])->setActivationMultiplier(0.025);
      else
      {
        // Figure out time in cycle
        double tempNormalizedTime = s * deltaT - activationTimesQP[k];
        double tempActivationMultiplier = ActivationFactor[tempNormalizedTime/deltaT];
        if (tempActivationMultiplier < 0.025)
          tempActivationMultiplier = 0.025;
	(PLmaterials[k])->setActivationMultiplier(tempActivationMultiplier);
      }
    }
      
    ind++;
    mySolver.solve(DISP);
    myModel.setPrevField();

    // Print out myocardium volume and ejection fraction (if requested)
    outVolume << s * deltaT << "\t" << myModel.computeCurrentVolume();

    if (calculateEjectionFractionFlag)
    {
      vector <double> currentMyocardiumField(Cube.getNumberOfNodes() * Cube.getDimension(), 0.0);
      myModel.getField(currentMyocardiumField);
      outVolume << "\t" << calculateEjectionFraction(cavityModel, &myModel, surfaceNodes, currentMyocardiumField);
      cavityModel->writeOutputVTK(outputString + "Cavity_", ind);
    }
    outVolume << endl;

    cout << "Volume data written." << endl;

    // Update State Variables:     
    for (int k = 0; k < NumMat; k++)
    {
      (PLmaterials[k])->updateStateVariables();
    }    

    cout << "State Variables Updated." << endl;
    
    // Write Output
    myModel.writeOutputVTK(outputString, ind);
    cout << "Output Written for step." << endl;
  }
  outVolume.close();

  // Timing
  time (&end);
  cout << endl << "BVP solved in " << difftime(end,start) << " s" << endl;

  return 0;
}

double calculateEjectionFraction(MechanicsModel* cavityModel, const MechanicsModel* myocardiumModel, const vector<int> surfaceNodes, const vector<double> currentMyocardiumField)
{
  if (cavityModel == NULL)
  {
    cout << "** Trying to calculate ejection fraction without turning the flag on." << endl;
    exit(1);
  }

  double referenceVolume = cavityModel->computeRefVolume();
 
  // Set the field of the cavity surface nodes 
  for (int i = 0; i < surfaceNodes.size(); i++)
  {
    for (int dim_iter = 0; dim_iter < 3; dim_iter++)
      cavityModel->setField(i * 3 + dim_iter, currentMyocardiumField[surfaceNodes[i] * 3 + dim_iter]);
  }

  // Set the z-field of the nodes along the mid-line
  // TODO: Hard-coded the values for the top and bottom of the cavities! Fix this ASAP!
  double cavityMinNode = 4148;
  double cavityMaxNode = 4068;
  double bottomMostNode = 1106;
  double cavityMinPos = currentMyocardiumField[cavityMinNode * 3 + 2];
  double cavityMaxPos = currentMyocardiumField[cavityMaxNode * 3 + 2];

  cout << "Min Pos:\t" << cavityMinPos << "\tMax Pos:\t" << cavityMaxPos << endl;

  // Linearly interpolate those nodes (assumes last ten nodes are the cavity!)
  int numberOfVerticalNodes = 10;
  for (int i = 0; i < numberOfVerticalNodes; i++)
  {
    double tempSlope = (cavityMaxPos - cavityMinPos)/(numberOfVerticalNodes-1);
    cavityModel->setField((bottomMostNode + i) * 3 + 2, tempSlope * (i) + cavityMinPos);
    // cout << "Node " << bottomMostNode + i << "\t Position: " << tempSlope * (i) + cavityMinPos << endl;
  }

  double EF = (referenceVolume - cavityModel->computeCurrentVolume())/referenceVolume;

  cout << "Reference Volume: \t" << referenceVolume << endl;
  cout << "Current Volume:   \t" << cavityModel->computeCurrentVolume() << endl;

  return EF;
}
