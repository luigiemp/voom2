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

double calculateCavityVolume(const vector<double> currentMyocardiumField, FEMesh* EndocardialSurfMesh, FEMesh* surfaceCapMesh, const vector<int> &endoBaseRingNodeSet);
void preConditionBoundaryCondition(MechanicsModel* myocardiumModel, int PbDoF, vector<MechanicsMaterial*>& plmaterials, EigenNRsolver* mySolver, double stiffness, int numSteps);
double Interpolate(double Field_i, double Time_i, double Field_ip1, double Time_ip1, double Time);

int main(int argc, char** argv)
{
  cout << string(50, '\n'); // Clear Screen

  // Timing
  time_t start, end;
  time(&start);

  // Activation Sequence? Set True to make activation a function of z
  bool activationSequence = false;

  // Constant Activation?
  bool constantActivation = false;
  double activationValue = 0.0;
  double minActivationFactor = 0.0;

  // Fibers across the wall thickness?
  bool fibersAcrossWall = true;

  // Calculate Ejection Fraction?
  bool calculateEjectionFractionFlag = true;

  // Use Conduction Velocity (if false - uses activation Time file)
  bool useConductionVelocity = true;

  // Conduction Velocity (cm/ms)
  double cv = 0.06;

  // Pressure File
  bool pressureFlag = false;
  string pressureFile = "InputFiles/Pressure_1ms_Inflation.dat";
  double pressureDivider = -1.0; // Divides the pressure by this value;
  
  // Simulation Time (in ms)
  double simTime = 2000;  

  // Time Step (in ms)
  double deltaT = 0.01;

  // OutputString
  string outputString = "/u/project/cardio/adityapo/ScratchResults/Scratch/Ellipsoid_";
  // string outputString = "/u/project/cardio/adityapo/ScratchResults/ContractionOnly/Ellipsoid";

  // Fiber Visualization String
  string fiberPlotString = outputString + "_Fiber.vtk";

  // Printing Volume to File
  string volumeFile = outputString + "_Volume.txt";

  // Spring BC:
  int SpringBCflag = 0;
  // Spring Stiffnes
  Real SpringK = 1.0e4; // 1.0e4;
  
  // Torsional Spring BC:
  int torsionalSpringBCflag = 1;
  // Torsional Spring Stiffness:
  int TorsionalSpringK = 1.0e2;

  
  // LJ Type Boundary Condition
  bool LJBoundaryConditionFlag = true;
  BCPotentialType potentialType = QUARTIC;		// QUADRATIC or QUARTIC
  FEMesh* LJBoundaryConditionMesh;
  if (LJBoundaryConditionFlag)
    LJBoundaryConditionMesh = new FEMesh("Mesh/EllipsoidMeshCoarse_Quadratic/EllipsoidCoarseQuadratic_LJBoundaryCondition0_5.node", "Mesh/EllipsoidMeshCoarse_Quadratic/EllipsoidCoarseQuadratic_LJBoundaryCondition0_5.ele");
  string LJBoundaryConditionEBCFile = "Mesh/EllipsoidMeshCoarse_Quadratic/EllipsoidCoarseQuadratic_LJBoundaryCondition0_5.EBCNodeSet";
  double LJsearchRadius = 1.50;
  double LJdepthPotentialWell = 1.0E3; // 4.0E2;
  double LJminDistance = 0.50;
  int numPreconditionSteps = 0;
  bool numNeighborFlag = true; int numNeighbors = 3;
  
  // Flexible:
  bool makeFlexible = false;
  double membraneStiffness = 1.0E2;
  double preTensionFactor = 0.98;                // Factor which controls pretension
  if (makeFlexible == true && LJBoundaryConditionFlag == false) {
    cout << "LJ Boundary Condition Flag needs to be on to make the BC Flexible" << endl;
    exit(EXIT_FAILURE);
  }

  // Lagrange Multiplier?
  bool LagrangeMultiplierApproach = true;
  double LagrangeMultiplier = 1.0;
  double vd = 0.5;

  // Windkessel Parameters
  bool windkesselFlag = true;
  double C = 5.16;
  double R1 = 0.030 * 1000;
  double R2 = 0.63 * 1000;
  double gamma = 0.15;
  double V_d = 85.0;
  double V_s = 60.0;
  double K = 1.0E4; // Bulk modulus of blood
  double P_LA = 15.0;	// mmHg (Left Atrial Pressure Constant)
  vector<double> P_LV; P_LV.push_back(15.0);	// mmHg (Initial left ventricular pressure)
  vector<double> V_LV; 


  // Initialize Mesh
  // Assumptions to use this main as is: strip has a face at z=0; tetrahedral mesh
  FEMesh Cube("Mesh/EllipsoidMeshCoarse_Quadratic/EllipsoidCoarseQuadratic.node", "Mesh/EllipsoidMeshCoarse_Quadratic/EllipsoidCoarseQuadratic.ele");
  FEMesh surfMesh("Mesh/EllipsoidMeshCoarse_Quadratic/EllipsoidCoarseQuadratic.node", "Mesh/EllipsoidMeshCoarse_Quadratic/EllipsoidCoarseQuadratic.EpicardiumElset");
  FEMesh innerSurfMesh("Mesh/EllipsoidMeshCoarse_Quadratic/EllipsoidCoarseQuadratic.node", "Mesh/EllipsoidMeshCoarse_Quadratic/EllipsoidCoarseQuadratic.EndocardiumElset");
  
  FEMesh surfaceCapMesh("Mesh/EllipsoidMeshCoarse_Quadratic/EllipsoidCoarseQuadratic_TopCapMesh.node", "Mesh/EllipsoidMeshCoarse_Quadratic/EllipsoidCoarseQuadratic_TopCapMesh.ele");
  string EndoBaseRingNodeSet = "Mesh/EllipsoidMeshCoarse_Quadratic/EllipsoidCoarseQuadratic_TopCapMesh.EndoBaseRingNodeset";

  string FiberFile = "Mesh/EllipsoidMeshCoarse_Quadratic/EllipsoidCoarseQuadratic.fiber";
  string BCfile = "Mesh/EllipsoidMeshCoarse_Quadratic/EllipsoidCoarseQuadratic.Null.bc";
  // string BCfile = "Mesh/EllipsoidMeshCoarse_Quadratic/EllipsoidCoarseQuadratic.Apex.bc";
  // string BCfile = "Mesh/EllipsoidMeshCoarse_Quadratic/EllipsoidCoarseQuadratic.BaseNodeset";
  string torsionalSpringBC = "Mesh/EllipsoidMeshCoarse_Quadratic/EllipsoidCoarseQuadratic.BaseNodeset";
  string ActivationTimeFile = "Mesh/EllipsoidMeshFiner/Small_B.activationTime";   // This is the Element Activation Time File
  string ActivationFile = "InputFiles/ActFunc_600ms_1msInterval.dat";  // This is the Calcium Transient
  ifstream FiberInp(FiberFile.c_str());
  ifstream ActTime(ActivationTimeFile.c_str());

  // Code for calculating Ejection Fraction
  vector <int> EndoBaseRingNodeSetVec;

  if (calculateEjectionFractionFlag)
  {
    ifstream EndoBaseRingNodeSetFileStream (EndoBaseRingNodeSet.c_str());

    if (EndoBaseRingNodeSetFileStream.is_open())
        cout << "** Opened EndoBaseRingNodeSet File Successfully." << endl;
    else
    {
        cout << "** Cannot Open EndoBaseRingNodeSet File." << endl;
        exit(1);
    }

    double numEndoBaseRingNodes;
    EndoBaseRingNodeSetFileStream >> numEndoBaseRingNodes;

    for (int i = 0; i < numEndoBaseRingNodes; i++)
    {
      double tempEndoBaseRingNode;
      EndoBaseRingNodeSetFileStream >> tempEndoBaseRingNode;
      EndoBaseRingNodeSetVec.push_back(tempEndoBaseRingNode);
    }
    EndoBaseRingNodeSetFileStream.close();
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
  uint NumEl =  Cube.getNumberOfElements();
  vector<GeomElement*> meshElements = Cube.getElements();
  int numQuadPoints = meshElements[0]->getNumberOfQuadPoints();
  uint NumMat = NumEl * numQuadPoints;
  
  vector<MechanicsMaterial * > PLmaterials;
  PLmaterials.reserve(NumEl * numQuadPoints);


  APForceVelPotential TestPotential(4.0, 1000.0 * 3.0, 3.0);	// 50.0 for 2nd parameter, force
  BlankViscousPotential ViscPotential;
  // NewtonianViscousPotential ViscPotential(0.005, 0.5);
  Vector3d HardParam(1.,1.,1.);

  // Visualize Fiber directions
  ofstream out;
  // Fibers VTK File Headers:
  out.open(fiberPlotString.c_str());
  out << "# vtk DataFile Version 3.1" << endl;
  out << "Fiber vector representation" << endl;
  out << "ASCII" << endl;
  out << "DATASET UNSTRUCTURED_GRID" << endl;
  out << "POINTS " << NumEl*numQuadPoints << " FLOAT" << endl;

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
  

  // Calculate the activation time for each quadrature point in each element
  vector <double> activationTimesQP;

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
    // Figure out Element Activation Time:
    double eleActivationTime = 0.0;
    if (ActTime.is_open())
      ActTime >> eleActivationTime;

    int el_numQuadPoints = meshElements[el_iter]->getNumberOfQuadPoints();
    int el_numNodes = meshElements[el_iter]->getNodesPerElement();
    vector<int> el_nodeIds = meshElements[el_iter]->getNodesID();
    for (int quadPt_iter = 0; quadPt_iter < el_numQuadPoints; quadPt_iter++) {
      // Compute the location of the quadrature point
      Real quadPointX = 0.0; Real quadPointY = 0.0; Real quadPointZ = 0.0;
      for (int el_node_iter = 0; el_node_iter < el_numNodes; el_node_iter++) {
	quadPointX += meshElements[el_iter]->getN(quadPt_iter, el_node_iter) * Cube.getX(el_nodeIds[el_node_iter], 0);
	quadPointY += meshElements[el_iter]->getN(quadPt_iter, el_node_iter) * Cube.getX(el_nodeIds[el_node_iter], 1);
	quadPointZ += meshElements[el_iter]->getN(quadPt_iter, el_node_iter) * Cube.getX(el_nodeIds[el_node_iter], 2);
      }
      
      // Activation Time Calculations:
      if(activationSequence) {
	if (useConductionVelocity)
	  activationTimesQP.push_back((quadPointZ - z_min)/cv);
	else
	  activationTimesQP.push_back(eleActivationTime);  // THIS IS WRONG BECAUSE ALL OF THE ELEMENT GETS ACTIVATED ONCE AND NOT BY EACH QUADRATURE POINT.
      }
      else
	activationTimesQP.push_back(0.0);  // Everything gets activated right away
      
      out << quadPointX << " " << quadPointY << " " << quadPointZ << endl;

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
	el_vectors[0] << 0., 0., 1.;
	el_vectors[1] << 1., 0., 0.;
	el_vectors[2] << 0., 1., 0.;
      }
      fiberVectors.push_back(el_vectors[0]);
      sheetVectors.push_back(el_vectors[1]);
      sheetNormalVectors.push_back(el_vectors[2]);

      // CompNeoHookean* PassiveMat = new CompNeoHookean(0, 1.0e4, 1.0e4);
      // CompNeoHookean* ActiveMat = new CompNeoHookean(0, 1.0e4, 1.0e4);
      // PLmaterials.push_back(PassiveMat);
      
      
      Humphrey_Compressible* PassiveMat = new Humphrey_Compressible(0, 15.98, 55.85, 0.0, -33.27, 30.21, 300.590, 6400.62, el_vectors);
      LinYinActive_Compressible* ActiveMat = new LinYinActive_Compressible(0, -38.70, 40.83, 25.12, 90.51, 1710.18, el_vectors);

      
      PlasticMaterial* PlMat = new PlasticMaterial(el_iter, ActiveMat, PassiveMat, &TestPotential, &ViscPotential);
      PlMat->setDirectionVectors(el_vectors);
      PlMat->setHardeningParameters(HardParam);
      PlMat->setActiveDeformationGradient(Matrix3d::Identity(3,3));
      PlMat->setTotalDeformationGradient(Matrix3d::Identity(3,3));
      
      PlMat->setTimestep(deltaT);
      PlMat->setActivationMultiplier(0.0);

      PLmaterials.push_back(PlMat);
      
    } // Quad point loop
  } // Element loop

  // Finish vtk file for plotting fiber data:
  out << "POINT_DATA " << meshElements.size() * numQuadPoints << endl;
  out << "VECTORS FiberDirection FLOAT" << endl;
  for (int el_iter = 0; el_iter < meshElements.size(); el_iter++)
    for (int quadPt_iter = 0; quadPt_iter < numQuadPoints; quadPt_iter++)
    out << fiberVectors[el_iter * numQuadPoints + quadPt_iter][0] << " " << fiberVectors[el_iter * numQuadPoints + quadPt_iter][1] << " " << fiberVectors[el_iter * numQuadPoints + quadPt_iter][2] << endl;
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
  
  cout << "** Applying Boundary Conditions" << endl;
  
  if (SpringBCflag)
    myModel.initSpringBC("Mesh/EllipsoidMeshCoarse_Quadratic/EllipsoidCoarseQuadratic.EpicardiumNodeset", &surfMesh, SpringK);

  if (torsionalSpringBCflag)
    myModel.initTorsionalSpringBC(torsionalSpringBC, TorsionalSpringK);
  
  if (LJBoundaryConditionFlag) {
    myModel.initializeLennardJonesBC(potentialType, "Mesh/EllipsoidMeshCoarse_Quadratic/EllipsoidCoarseQuadratic.EpicardiumNodeset", LJBoundaryConditionMesh, &surfMesh, LJsearchRadius, LJdepthPotentialWell, LJminDistance);
    myModel.toggleNumberNeighborFlag(numNeighborFlag, numNeighbors);
    if (makeFlexible) {
      myModel.initializeMembraneStiffness(membraneStiffness);
      myModel.initiatePreTensionInMembrane(preTensionFactor);
    }
  }

  // Impose a Lagrange Multiplier
  if (LagrangeMultiplierApproach) {
    cout << "Setting up Lagrange multiplier approach" << endl;
    myModel.initializeLagrangeMultiplierMethod(&innerSurfMesh, &surfaceCapMesh, EndoBaseRingNodeSetVec, vd);
    cout << "CAVITY VOLUME: " << myModel.computeCavityVolume() << endl;
    myModel.setLagrangeMultiplier(LagrangeMultiplier);
    myModel.setTargetVolume(myModel.computeCavityVolume());
  }

  cout << "Model Setup" << endl;

  // Initialize Result
  uint myRequest;
  uint PbDoF = (Cube.getNumberOfNodes())*myModel.getDoFperNode();
  uint mainDoF = PbDoF; // This is for adding EBC's in case it's flexible
  if (makeFlexible) PbDoF += LJBoundaryConditionMesh->getNumberOfNodes() * myModel.getDoFperNode();
  if (LagrangeMultiplierApproach) PbDoF++;
  EigenResult myResults(PbDoF, 0);

  // Run Consistency check

  Real perturbationFactor = 0.1;
  myRequest = 7; // Check both Forces and Stiffness
  Real myH = 1e-6;
  Real myTol = 1e-7;

  // Before checking consistency, the perturbed deformation state must be 
  // set to the current deformation state.

  // myModel.checkConsistency(myResults, perturbationFactor, myRequest, myH, myTol);

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
      // cout << BCid[ind] << " " <<  BCvalues[ind] << endl;
      ind++;
    }
  }
  BCinp.close();

  // Setting EBC for LJ Boundary Condition
  if (makeFlexible) {
    cout << "********" << " Setting Membrane EBCs " << "********" << endl;
    ifstream MembraneBCinp(LJBoundaryConditionEBCFile.c_str());

    if (MembraneBCinp.is_open())
      cout << "Membrane BC File Opened Successfully!" << endl;
    else
      cout << "ERROR: BC File Failed to Open!" << endl;

   MembraneBCinp >> NumBC;
    cout << "Number of Nodes with MembraneEBC: " << NumBC << endl;
    BCid.reserve(NumBC*3);
    BCvalues.reserve(NumBC*3);
    for(int i = 0; i < NumBC; i++) {
      MembraneBCinp >> node;
      BCnodes.push_back(node);
      for (int j = 0; j < 3; j++) {
        BCid.push_back(mainDoF + node*3 + j);
        BCvalues.push_back(LJBoundaryConditionMesh->getX(node, j));
        // cout << BCid[ind] << " " <<  BCvalues[ind] << endl;
        ind++;
      }
    }
    MembraneBCinp.close();
  }

  // Solver
  Real NRtol = 1.0e-4;
  uint NRmaxIter = 100;
  EigenNRsolver mySolver(&myModel, BCid, BCvalues, CHOL, NRtol, NRmaxIter);

  // ********************************* // 
  // SOLVE:

  ind = 0;
  myModel.finalizeCompute();

  // Precondition
  if (LJBoundaryConditionFlag && numPreconditionSteps > 0)
    preConditionBoundaryCondition(&myModel, PbDoF, PLmaterials, &mySolver, LJdepthPotentialWell, numPreconditionSteps);

  ofstream outVolume;
  outVolume.open(volumeFile.c_str());

  for (int s = 0; s < simTime/deltaT; s++)
  {
    cout << endl << "*** Step " << s << " ***" << endl;
    if (SpringBCflag) myModel.computeSpringNormals();

    // Update pressure:
    if (pressureFlag) {
	double appliedPressure = pressureData[s][1]/pressureDivider;
        myModel.updatePressure(appliedPressure);
	cout << "* Applying " << appliedPressure << " units of pressure to the endocardium." << endl;
    }

    for (int k = 0; k < meshElements.size(); k++)
    {
      for (int q = 0; q < numQuadPoints; q++) {
	(PLmaterials[k * numQuadPoints + q])->setTimestep(deltaT/1000);
	
	if (s * deltaT < activationTimesQP[k * numQuadPoints + q] || s * deltaT > activationTimesQP[k * numQuadPoints + q] + cycleLength)
	  (PLmaterials[k * numQuadPoints + q])->setActivationMultiplier(minActivationFactor);
	else
	{
	  // Figure out time in cycle
	  double tempNormalizedTime = s * deltaT - activationTimesQP[k];
          double tempActivationMultiplier = ActivationFactor[tempNormalizedTime/deltaT];
          if (tempActivationMultiplier < minActivationFactor)
            tempActivationMultiplier = minActivationFactor;
	  (PLmaterials[k * numQuadPoints + q])->setActivationMultiplier(tempActivationMultiplier);
	}
      }
    }
      
    ind++;
    mySolver.solve(DISP);
    // myModel.setPrevField();
    myModel.finalizeCompute();  // This sets the previous field for the spring normals.

    // Print out myocardium volume and ejection fraction (if requested)
    outVolume << s * deltaT << "\t" << myModel.computeCurrentVolume();

    if (calculateEjectionFractionFlag)
    {
      vector <double> currentMyocardiumField(PbDoF, 0.0);
      // vector <double> currentMyocardiumField(Cube.getNumberOfNodes() * Cube.getDimension(), 0.0);
      myModel.getField(currentMyocardiumField);
      double cavityVolume = calculateCavityVolume(currentMyocardiumField, &innerSurfMesh, &surfaceCapMesh, EndoBaseRingNodeSetVec);
      outVolume << "\t" << cavityVolume;
      cout << "Cavity Volume: " << cavityVolume << endl;
      if (LagrangeMultiplierApproach) outVolume << "\t" << myModel.getLagrangeMultiplier();
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
    cout << "Reference Volume: " << myModel.computeRefVolume() << "\t Current Volume: " << myModel.computeCurrentVolume() << endl;
    if (LagrangeMultiplierApproach) cout << "Internal Pressure (Lagrange Multiplier): " << myModel.getLagrangeMultiplier() << endl;
  }
  outVolume.close();

  // Timing
  time (&end);
  cout << endl << "BVP solved in " << difftime(end,start) << " s" << endl;

  return 0;
}



double calculateCavityVolume(const vector<double> currentMyocardiumField, FEMesh* EndocardialSurfMesh, FEMesh* surfaceCapMesh, const vector<int> &endoBaseRingNodeSet) {
  double volume = 0.0;
  
  // First deform the surfaceCapMesh
  vector <VectorXd> surfaceCapMeshNewX;
  VectorXd baseCenterNode(3);
  baseCenterNode << 0.0, 0.0, 0.0;

  // Preallocate the vector of vectorXd
  for (int node_iter = 0; node_iter < surfaceCapMesh->getNumberOfNodes(); node_iter++) {
    VectorXd tempVecXd(3);
    tempVecXd << 0.0, 0.0, 0.0;
    surfaceCapMeshNewX.push_back(tempVecXd);
  }

  for (int surfCapNode_Iter = 0; surfCapNode_Iter < endoBaseRingNodeSet.size(); surfCapNode_Iter++) {
    surfaceCapMeshNewX[surfCapNode_Iter](0) = currentMyocardiumField[endoBaseRingNodeSet[surfCapNode_Iter] * 3 + 0];
    surfaceCapMeshNewX[surfCapNode_Iter](1) = currentMyocardiumField[endoBaseRingNodeSet[surfCapNode_Iter] * 3 + 1];
    surfaceCapMeshNewX[surfCapNode_Iter](2) = currentMyocardiumField[endoBaseRingNodeSet[surfCapNode_Iter] * 3 + 2];
    baseCenterNode = baseCenterNode + surfaceCapMeshNewX[surfCapNode_Iter];
  }
  
  // Compute the center node of the base
  baseCenterNode = baseCenterNode/endoBaseRingNodeSet.size();
  surfaceCapMeshNewX[endoBaseRingNodeSet.size()] = baseCenterNode;
  int baseCenterNodeId = endoBaseRingNodeSet.size();
  
  vector<GeomElement*> surfaceCapElements = surfaceCapMesh->getElements();
  // Compute positions of midside nodes
  for (int el_iter = 0; el_iter < surfaceCapElements.size(); el_iter++) {
    vector<int> el_nodeIds = surfaceCapElements[el_iter]->getNodesID();

    if (el_nodeIds.size() > 3) {
      if (el_nodeIds.size() > 6) {
	cout << "ERROR: Can only compute midside node for a Quad triangle." << endl;
	return EXIT_FAILURE;
      }
      // Node 0 and Node 1
      if (el_nodeIds[0] == baseCenterNodeId || el_nodeIds[1] == baseCenterNodeId) {
        int edgeVertex = (el_nodeIds[0] == baseCenterNodeId) ? el_nodeIds[1] : el_nodeIds[0];
        VectorXd midsideNode = 0.5 * (surfaceCapMeshNewX[edgeVertex] + baseCenterNode);
	surfaceCapMeshNewX[el_nodeIds[3]] = midsideNode;
      }
      if (el_nodeIds[1] == baseCenterNodeId || el_nodeIds[2] == baseCenterNodeId) {
        int edgeVertex = (el_nodeIds[1] == baseCenterNodeId) ? el_nodeIds[2] : el_nodeIds[1];
        VectorXd midsideNode = 0.5 * (surfaceCapMeshNewX[edgeVertex] + baseCenterNode);
        surfaceCapMeshNewX[el_nodeIds[4]] = midsideNode;
      }
      if (el_nodeIds[2] == baseCenterNodeId || el_nodeIds[0] == baseCenterNodeId) {
        int edgeVertex = (el_nodeIds[2] == baseCenterNodeId) ? el_nodeIds[0] : el_nodeIds[2];
        VectorXd midsideNode = 0.5 * (surfaceCapMeshNewX[edgeVertex] + baseCenterNode);
        surfaceCapMeshNewX[el_nodeIds[5]] = midsideNode;
      }        
    }
  }
  surfaceCapMesh->setX(surfaceCapMeshNewX);

  // Compute first the volume associated with myocardiumModel
  vector<GeomElement*> endocardialElements = EndocardialSurfMesh->getElements();
  for (int el_iter = 0; el_iter < EndocardialSurfMesh->getNumberOfElements(); el_iter++) {
    vector <int> el_nodeIds = endocardialElements[el_iter]->getNodesID();

    int numQP = endocardialElements[el_iter]->getNumberOfQuadPoints();
    for (int quadPt_iter = 0; quadPt_iter < numQP; quadPt_iter++) {
      Vector3d a1 = Vector3d::Zero(), a2 = Vector3d::Zero(), a3 = Vector3d::Zero();
      Vector3d quadPtLocation = Vector3d::Zero();

      for (int a = 0; a < el_nodeIds.size(); a++) {
        int nodeID = el_nodeIds[a];
        Vector3d current_x;
	current_x << currentMyocardiumField[nodeID * 3], currentMyocardiumField[nodeID * 3 + 1], currentMyocardiumField[nodeID * 3 + 2];
	quadPtLocation += current_x * endocardialElements[el_iter]->getN(quadPt_iter, a);

        a1 += current_x * endocardialElements[el_iter]->getDN(quadPt_iter, a, 0);
        a2 += current_x * endocardialElements[el_iter]->getDN(quadPt_iter, a, 1);
      }
      a3 = a1.cross(a2);
      volume += 1./3. * quadPtLocation.dot(a3) * endocardialElements[el_iter]->getQPweights(quadPt_iter);
    }
  } // End of for loop endocardial elements


  // Compute the volume associated with surfaceCapMesh
  for (int el_iter = 0; el_iter < surfaceCapMesh->getNumberOfElements(); el_iter++) {
    vector <int> el_nodeIds = surfaceCapElements[el_iter]->getNodesID();

    int numQP = surfaceCapElements[el_iter]->getNumberOfQuadPoints();
    for (int quadPt_iter = 0; quadPt_iter < numQP; quadPt_iter++) {
      Vector3d a1 = Vector3d::Zero(), a2 = Vector3d::Zero(), a3 = Vector3d::Zero();
      Vector3d quadPtLocation = Vector3d::Zero();

      for (int a = 0; a < el_nodeIds.size(); a++) {
        int nodeID = el_nodeIds[a];
        Vector3d current_x;
        current_x << surfaceCapMeshNewX[nodeID](0), surfaceCapMeshNewX[nodeID](1), surfaceCapMeshNewX[nodeID](2);
        quadPtLocation += current_x * surfaceCapElements[el_iter]->getN(quadPt_iter, a);

        a1 += current_x * surfaceCapElements[el_iter]->getDN(quadPt_iter, a, 0);
        a2 += current_x * surfaceCapElements[el_iter]->getDN(quadPt_iter, a, 1);
      }
      a3 = a1.cross(a2);
      volume += 1./3. * quadPtLocation.dot(a3) * surfaceCapElements[el_iter]->getQPweights(quadPt_iter);
    }
  } // End of for loop endocardial elements
  
  return volume;
}

void preConditionBoundaryCondition(MechanicsModel* myocardiumModel, int PbDoF, vector<MechanicsMaterial*>& plmaterials, EigenNRsolver* mySolver, double stiffness, int numSteps) {
  // myocardiumModel->toggleConstantMinimumDistanceFlag(false);
  cout << " **Preconditioning the pericardium boundary condition..." << endl;
  for (int iter = 0; iter < numSteps; iter++) {
    // double currentStiffness = stiffness/numSteps * (iter + 1);
    double currentStiffness = pow(10, (log10(stiffness)/numSteps * (iter + 1))) - 1.0;
    // double currentStiffness = stiffness/numSteps * (1+1);
    cout << "** Preconditioning with stiffness = " << currentStiffness << endl;

    // Compute reference configuration
    // 1. Set stiffness
    myocardiumModel->setLJStiffness(currentStiffness);
    
    // 2. Solve with zero loads
    mySolver->solve(DISP);
    myocardiumModel->writeOutputVTK("PreconditionStep_", iter);
    
    // 3. Get current field and set it to X
    Mesh* myMesh = myocardiumModel->getMesh();
    int numNodes = myMesh->getNumberOfNodes();
    int dim = myMesh->getDimension();
    // vector<double> currentField(numNodes * dim);
    vector<double> currentField(PbDoF, 0.0);
    myocardiumModel->getField(currentField);
    
    vector<VectorXd> newX;
    for (int node_iter = 0; node_iter < numNodes; node_iter++) {
      VectorXd nodeX(3);
      nodeX << currentField[node_iter * dim], currentField[node_iter * dim + 1], currentField[node_iter * dim + 2];
      newX.push_back(nodeX);
    }

    // VectorXd tempX(3);
    // tempX << myMesh->getX(0,0), myMesh->getX(0,1), myMesh->getX(0,2);
    // tempX = tempX - newX[0];
    // cout << "TEMPX: " << tempX.norm() << endl;

    myMesh->setX(newX);
    myocardiumModel->writeOutputVTK("PreconditionStepPostSetX_", iter);
    
    // 4. Reset Q, field, and history variables
    for (int mat_iter = 0; mat_iter < plmaterials.size(); mat_iter++) {
      plmaterials[mat_iter]->resetState();
    }
    myocardiumModel->setPrevField();

   // 5. Call finalize model compute - computes the new nearest neighbors for pericardium boundary condition
   myocardiumModel->finalizeCompute(); 
   
   // 6. Recompute the distance to reset energy to 0. NOT SURE IF WE WANT TO DO THIS YET?
   // myocardiumModel->recomputeAverageMinDistance();
  } 
  cout << "** Finished Preconditioning the pericardium boundary condition." << endl;
}


double Interpolate(double Field_i, double Time_i, double Field_ip1, double Time_ip1, double Time) {
    double m = (Field_ip1 - Field_i)/(Time_ip1 - Time_i);
    double c = Field_i - m * Time_i;
    return m * Time + c;
}
