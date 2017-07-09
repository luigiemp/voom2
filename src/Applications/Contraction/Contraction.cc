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

// Boost Header Files for String manipulation
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

// VTK Header Files
#include <vtkVersion.h>
#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkTensor.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkCellType.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>

// Namespaces
using namespace voom;
using namespace boost;

// Struct to store secondary Meshes
struct secondaryMesh {
  string secondaryMeshName;
  string secondaryMeshNodeFileName;
  string secondaryMeshEleFileName;
  FEMesh* secondaryMeshObject;
};

// Struct to store nodesets
struct nodeSet {
  string nodeSetName;
  string nodeSetMesh;
  string nodeSetPath;
};

enum materialTypes {
  COMPNEOHOOKEAN, JACOBIAN, LINYINACTIVE_COMPRESSIBLE, HUMPHREY_COMPRESSIBLE,
  HILLFORCEVELPOTENTIAL, APFORCEVELPOTENTIAL,
  BLANKVISCOUSPOTENTIAL, NEWTONIANVISCOUSPOTENTIAL
};
// Struct for material
struct material {
  string materialResponseCaptured;
  materialTypes materialName;
  vector<double> materialProperties;
};

// enumerated Type for boundary condition type
enum boundaryConditionTypes {
  EBC, TORSIONALSPRING, SPRINGFOUNDATION, INTERACTION
};

// Struct for boundary conditions
struct boundaryCondition {
  boundaryConditionTypes BCType;
  nodeSet* BCNodeSet;
  Mesh* surfMesh; // Only used for SPRINGFOUNDATION and INTERACTION boundary conditions
  Mesh* BCMesh; // Only used for Interaction type boundary conditions
  vector<double> BCProperties;
};

// Helper Function Prototypes:
FEMesh* createFEMesh(string nodeFile, string eleFile, string meshName);
vector <vector <Vector3d> > readFiberData(string fiberFile, const FEMesh* mainMesh, bool writeFibers, string fiberPlotString);
materialTypes determineMaterialName(string materialName);
MechanicsMaterial* createMechanicsMaterial(materialTypes materialName, vector<double> materialProperties, vector<Vector3d> fibers);
Potential* createKineticPotentialMaterial(materialTypes materialName, vector<double> materialProperties);
ViscousPotential* createViscousPotentialMaterial(materialTypes materialName, vector<double> materialProperties);
boundaryConditionTypes determineBoundaryConditionType(string boundaryConditionName);
FEMesh* findMeshObject(string meshName, vector<secondaryMesh> &secondaryMeshObjects);
nodeSet* findNodeSet(string nodeSetName, vector<nodeSet> &nodeSetObjects);
void readTwoColumnDataFile(string twoColumnDataFilePath, vector <vector <double> > &twoColumnData);
double calculateEjectionFraction(MechanicsModel* cavityModel, const MechanicsModel* myocardiumModel, const vector<int> surfaceNodes, const vector<double> currentMyocardiumField);

int main(int argc, char** argv)
{
  // Timing
  time_t start, end;
  time(&start);
  

  cout << string(50, '\n'); // Clear Screen

  // Activation Sequence? Set True to make activation a function of z
  bool activationSequence = false;

  // Constant Activation?
  bool constantActivation = true;
  double activationValue = 0.0;
  double minActivationFactor = 0.0;

  // Fibers across the wall thickness?
  bool fibersAcrossWall = true;

  // Calculate Ejection Fraction?
  bool calculateEjectionFractionFlag = false;

  // Use Conduction Velocity (if false - uses activation Time file)
  bool useConductionVelocity = true;

  // Conduction Velocity (cm/ms)
  double cv = 0.06;

  // Read Input File
  // Check if input file has been supplied
  if (argc == 1) {
    cerr << "\033[41;37mERROR\033[0m\t: Input file not specified." << endl;
    cerr << "Usage\t: " << argv[0] << " InputFile" << endl;
    return EXIT_FAILURE;
  }

  // Check if Input File can be opened:
  ifstream inp;
  inp.open(argv[1]);
  if (!inp.is_open()) {
    cerr << "\033[41;37mERROR\033[0m\t: Cannot open file " << argv[1] << endl;
    return EXIT_FAILURE;
  }

  // ** Variables
  // * Mesh
  string nodeFile, eleFile, fiberFile;
  vector<secondaryMesh> secondaryMeshes;
  vector<nodeSet> nodeSets;

  // * Compute Options
  double simulationTime = 100;
  double deltaT = 0.01;
  double NRTOL = 1.0e-5;
  double NRIterations = 100;
  bool checkConsistencyFlag = false;

  // * Output
  bool writeFibers = false; string modelOutputName = "Contraction";
  string outputDirectory = "./";

  // * Material Properties
  vector<material> materials;

  // * Boundary Conditions
  vector <boundaryCondition> boundaryConditions;

  // * Pressure
  bool pressureFlag = false;
  FEMesh* pressureSurfMesh;
  vector <vector <double> > timePressureData;

  // * Activation Profile
  bool activationProfileProvided = false;
  vector <vector <double> > activationProfileData;

  // Read Input File:
  cout << "\033[1;33mReading Input File\033[0m\t: " << argv[1]  << endl << endl;
  string line;
  while (getline(inp, line)) {
    trim(line);
    if (find_first(line,"#") || find_first(line, "$") || find_first(line, "*"))
      continue; // comments
    vector<string> varArgs;
    split(varArgs, line, is_any_of("\t "), token_compress_on);
    to_upper(varArgs[0]);
    
    // Mesh Data
    if (varArgs[0] == "NODEFILE") {
      if (varArgs.size() >= 2) nodeFile = varArgs[1];
    }
    if (varArgs[0] == "ELEFILE") {
      if (varArgs.size() >= 2) eleFile = varArgs[1];
    }
    if (varArgs[0] == "FIBERFILE") {
      if (varArgs.size() >= 2) fiberFile = varArgs[1];
    }
    if (varArgs[0] == "SECONDARYMESH") {
      secondaryMeshes.push_back(secondaryMesh());
      secondaryMeshes[secondaryMeshes.size() - 1].secondaryMeshObject = (FEMesh *) malloc(sizeof(FEMesh));
      if (varArgs.size() >= 2) secondaryMeshes[secondaryMeshes.size() - 1].secondaryMeshName = varArgs[1];
      if (varArgs.size() >= 3) secondaryMeshes[secondaryMeshes.size() - 1].secondaryMeshNodeFileName = varArgs[2];
      if (varArgs.size() >= 4) secondaryMeshes[secondaryMeshes.size() - 1].secondaryMeshEleFileName = varArgs[3];

      secondaryMeshes[secondaryMeshes.size() - 1].secondaryMeshObject = createFEMesh(secondaryMeshes[secondaryMeshes.size() - 1].secondaryMeshNodeFileName, secondaryMeshes[secondaryMeshes.size() - 1].secondaryMeshEleFileName, secondaryMeshes[secondaryMeshes.size() - 1].secondaryMeshName);
    }
    if (varArgs[0] == "NODESET") {
      nodeSets.push_back(nodeSet());
      if (varArgs.size() >= 2) nodeSets[nodeSets.size() - 1].nodeSetName = varArgs[1];
      if (varArgs.size() >= 3) nodeSets[nodeSets.size() - 1].nodeSetMesh = varArgs[2];
      if (varArgs.size() >= 4) nodeSets[nodeSets.size() - 1].nodeSetPath = varArgs[3];
    }
    if (varArgs[0] == "WRITEFIBERS") {
      if (varArgs.size() >= 2) {
	to_upper(varArgs[1]);
	writeFibers = (varArgs[1] == "TRUE") ? true : false;
      }
    }
    if (varArgs[0] == "MODELOUTPUTNAME") {
      if (varArgs.size() >= 2) modelOutputName = varArgs[1];
    }
    if (varArgs[0] == "OUTPUTDIRECTORY") {
      if (varArgs.size() >= 2) outputDirectory = varArgs[1];
    }
    if (varArgs[0] == "SIMULATIONTIME") {
      if (varArgs.size() >= 2) simulationTime = atof(varArgs[1].c_str());
    }
    if (varArgs[0] == "DELTAT") {
      if (varArgs.size() >= 2) deltaT = atof(varArgs[1].c_str());
    }
    if (varArgs[0] == "MATERIAL") {
      materials.push_back(material());
      if (varArgs.size() >= 2) {to_upper(varArgs[1]); materials[materials.size() - 1].materialResponseCaptured = varArgs[1]; }
      if (varArgs.size() >= 3) materials[materials.size() - 1].materialName = determineMaterialName(varArgs[2]);
      if (varArgs.size() >= 4) {
	for (int matPropIter = 3; matPropIter < varArgs.size(); matPropIter++)
	  materials[materials.size() - 1].materialProperties.push_back(atof(varArgs[matPropIter].c_str()));
      }
    }
    if (varArgs[0] == "BOUNDARYCONDITION") {
      boundaryConditions.push_back(boundaryCondition());
      boundaryConditions[boundaryConditions.size() - 1].BCNodeSet = (nodeSet * ) malloc(sizeof(nodeSet));
      boundaryConditions[boundaryConditions.size() - 1].surfMesh = (FEMesh * ) malloc(sizeof(FEMesh));
      boundaryConditions[boundaryConditions.size() - 1].BCMesh = (FEMesh * ) malloc(sizeof(FEMesh));
      int startOfBCProps = 3;
      if (varArgs.size() >= 2) {boundaryConditions[boundaryConditions.size() - 1].BCType = determineBoundaryConditionType(varArgs[1]);}
      if (varArgs.size() >= 3) {boundaryConditions[boundaryConditions.size() - 1].BCNodeSet = findNodeSet(varArgs[2], nodeSets);}
      if (varArgs.size() >= 4) {
	if (boundaryConditions[boundaryConditions.size() - 1].BCType == INTERACTION) {
	  startOfBCProps = 5;
          boundaryConditions[boundaryConditions.size() - 1].surfMesh = findMeshObject(varArgs[3], secondaryMeshes);
	  boundaryConditions[boundaryConditions.size() - 1].BCMesh = findMeshObject(varArgs[4], secondaryMeshes);
        }
        else if (boundaryConditions[boundaryConditions.size() - 1].BCType == SPRINGFOUNDATION) {
          boundaryConditions[boundaryConditions.size() - 1].surfMesh = findMeshObject(varArgs[3], secondaryMeshes);
	  startOfBCProps = 4;
        }
	for (int BCPropIter = startOfBCProps; BCPropIter < varArgs.size(); BCPropIter++)
	  boundaryConditions[boundaryConditions.size() - 1].BCProperties.push_back(atof(varArgs[BCPropIter].c_str()));
      }
    } // End of BC
    if (varArgs[0] == "PRESSURE") {
      pressureFlag = true;
      cout << "\033[1;33mPressure flag is on.\033[0m\t " << endl;
      if (varArgs.size() >= 2) pressureSurfMesh = findMeshObject(varArgs[1], secondaryMeshes);
      if (varArgs.size() >= 3) {
        to_upper(varArgs[2]);
        if (varArgs[2] == "CONSTANT") {
	  vector<double> timePressureDataTemp(2, 0.0);
	  timePressureDataTemp[0] = 0.0; timePressureDataTemp[1] = atof(varArgs[3].c_str());
	  timePressureData.push_back(timePressureDataTemp);
	  cout << "\t\033[1;33mApplying a constant pressure of \033[0m"<< timePressureDataTemp[1] << " units." << endl;
        }
	else if (varArgs[2] == "FILE") {
	  readTwoColumnDataFile(varArgs[3], timePressureData);
        }
	else {
          cerr << "\033[41;37mERROR\033[0m\t: Unknown pressure type " << varArgs[1] << endl;
          exit(EXIT_FAILURE);
        }
        cout << endl;
      }
    }
    if (varArgs[0] == "ACTIVATIONPROFILE") {
      if (varArgs.size() >= 2) {
        to_upper(varArgs[1]);
	activationProfileProvided = true;
	if (varArgs[1] == "FILE") {
	  readTwoColumnDataFile(varArgs[2], activationProfileData);
        }
        else if (varArgs[1] == "CONSTANT") {
          vector<double> activationProfileDataTemp(2, 0.0);
          activationProfileDataTemp[0] = 0.0; activationProfileDataTemp[1] = atof(varArgs[3].c_str());
          activationProfileData.push_back(activationProfileDataTemp);
        }
	else {
	  cerr << "\033[41;37mERROR\033[0m\t: Unknown activation profile type " << varArgs[1] << endl;
          exit(EXIT_FAILURE);
        }
      }
    }
    if (varArgs[0] == "ACTIVATIONSEQUENCE") {
      if (varArgs.size() >= 2) {
        to_upper(varArgs[1]);
        if (varArgs[1] == "CONDUCTIONVELOCITY") {

        }
        else if (varArgs[1] == "INSTANTANEOUS") {

        }
        else if (varArgs[1] == "FILE") {

        }
        else {
	  cerr << "\033[41;37mERROR\033[0m\t: Unknown activation sequence type " << varArgs[1] << endl;
          exit(EXIT_FAILURE);
        }
      }
    }
    if (varArgs[0] == "NRTOL") {
      if (varArgs.size() >= 2) {NRTOL = atof(varArgs[1].c_str());}
    }
    if (varArgs[0] == "NRITERATIONS") {
      if (varArgs.size() >= 2) {NRIterations = atoi(varArgs[1].c_str());}
    }
    if (varArgs[0] == "CHECKCONSISTENCY") {
      if (varArgs.size() >= 2) {
        to_upper(varArgs[1]);
        checkConsistencyFlag = (varArgs[1] == "TRUE") ? true : false;
      }
    }

    varArgs.clear();
  }
  inp.close();
  // ****** // 
  // Setup Output
  string outputString = outputDirectory + "/" + modelOutputName;

  // Fiber Visualization String
  string fiberPlotString = outputString + "_Fiber.vtu";
  
  // Printing Volume to File
  string volumeFile = outputString + "_Volume.txt";

  // Initialize Mesh
  FEMesh* mainMesh = createFEMesh(nodeFile, eleFile, "Main-Mesh");
  vector <vector <Vector3d> > fiberData = readFiberData(fiberFile, mainMesh, writeFibers, fiberPlotString);

  // Print Compute Parameters to Screen


  // Initialize Mesh
  string ActivationTimeFile = "Mesh/EllipsoidMeshFiner/Small_B.activationTime";   // This is the Element Activation Time File
  string ActivationFile = "InputFiles/ActFunc_600ms_1msInterval.dat";  // This is the Calcium Transient
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
  for (int node_iter = 0; node_iter < mainMesh->getNumberOfNodes(); node_iter++)
  {
    VectorXd tempNodalPos = mainMesh->getX(node_iter);
    if(tempNodalPos[2] < z_min)
      z_min = tempNodalPos[2];
  }
  
  // ***************************************** //
  // **************Begin Computing************ //
  // ***************************************** //
  cout << endl;
  cout << "\033[1;32mNumber Of Nodes \t : \033[0m" << mainMesh->getNumberOfNodes() << endl;
  cout << "\033[1;32mNumber Of Element \t : \033[0m" << mainMesh->getNumberOfElements() << endl;
  cout << "\033[1;32mMesh Dimension \t\t : \033[0m" << mainMesh->getDimension() << endl << endl;

  // Initialize Material
  uint NumEl =  mainMesh->getNumberOfElements();
  vector<GeomElement*> meshElements = mainMesh->getElements();
  int numQuadPoints = meshElements[0]->getNumberOfQuadPoints();
  uint NumMat = NumEl * numQuadPoints;
  
  vector<MechanicsMaterial * > PLmaterials;
  PLmaterials.reserve(NumEl * numQuadPoints);



  // NewtonianViscousPotential ViscPotential(0.005, 0.5);
  Vector3d HardParam(1.,1.,1.);

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
      // Activation Time Calculations:
      if(activationSequence) {
	if (useConductionVelocity) {
	  // activationTimesQP.push_back((quadPointZ - z_min)/cv);
        }
	else
	  activationTimesQP.push_back(eleActivationTime);  // THIS IS WRONG BECAUSE ALL OF THE ELEMENT GETS ACTIVATED ONCE AND NOT BY EACH QUADRATURE POINT.
      }
      else
	activationTimesQP.push_back(0.0);  // Everything gets activated right away
      
      vector<Vector3d> el_vectors = fiberData[el_iter * el_numQuadPoints + quadPt_iter];


      // Create Passive, Active, Viscous, and Kinetic Potential Materials
      MechanicsMaterial* PassiveMat;
      MechanicsMaterial* ActiveMat;
      ViscousPotential* ViscPotential;
      Potential* TestPotential;

      for (int material_iter = 0; material_iter < materials.size(); material_iter++) {
	if (materials[material_iter].materialResponseCaptured == "PASSIVEMATERIAL") PassiveMat = createMechanicsMaterial(materials[material_iter].materialName, materials[material_iter].materialProperties, el_vectors);
        if (materials[material_iter].materialResponseCaptured == "ACTIVEMATERIAL") ActiveMat = createMechanicsMaterial(materials[material_iter].materialName, materials[material_iter].materialProperties, el_vectors);
        if (materials[material_iter].materialResponseCaptured == "KINETICPOTENTIAL") TestPotential = createKineticPotentialMaterial(materials[material_iter].materialName, materials[material_iter].materialProperties);
        if (materials[material_iter].materialResponseCaptured == "VISCOUSPOTENTIAL") ViscPotential = createViscousPotentialMaterial(materials[material_iter].materialName, materials[material_iter].materialProperties);
      }

      PLmaterials.push_back(PassiveMat);
      
      
      // PlasticMaterial* PlMat = new PlasticMaterial(el_iter, ActiveMat, PassiveMat, &TestPotential, &ViscPotential);
      // PlMat->setDirectionVectors(el_vectors);
      // PlMat->setHardeningParameters(HardParam);
      // PlMat->setActiveDeformationGradient(Matrix3d::Identity(3,3));
      // PlMat->setTotalDeformationGradient(Matrix3d::Identity(3,3));
      
      // PlMat->setTimestep(deltaT);
      // PlMat->setActivationMultiplier(0.0);

      // PLmaterials.push_back(PlMat);
    } // Quad point loop
  } // Element loop

  ActTime.close();

  // Initialize Model
  int NodeDoF = 3;
  int PressureFlag = 1;
  Real Pressure = 0.0;
  int NodalForcesFlag = 0;
  vector<int > ForcesID;
  vector<Real > Forces;
  MechanicsModel myModel(mainMesh, PLmaterials, NodeDoF, PressureFlag, pressureSurfMesh,
			 NodalForcesFlag, 0);
  myModel.updatePressure(Pressure);
  myModel.updateNodalForces(&ForcesID, &Forces);



  // Incorporate Boundary Conditions:
  int NumBC = 0, node = 0;
  vector<int > BCnodes;
  vector<int > BCid;
  vector<Real > BCvalues;

  cout << "\033[1;33mApplying boundary conditions to the model\033[0m\t"  << endl;
  for (int bc_iter = 0; bc_iter < boundaryConditions.size(); bc_iter++) {
    boundaryConditionTypes bc_type = boundaryConditions[bc_iter].BCType;
    switch(bc_type) {
      case EBC: {
	cout << "\t\033[1;33mApplying EBC.\033[0m" << endl;
	cout << "\t\033[1;33mEBC node set used\033[0m:\t" << boundaryConditions[bc_iter].BCNodeSet->nodeSetName << endl;
        
        ifstream BCinp(boundaryConditions[bc_iter].BCNodeSet->nodeSetPath.c_str());
        if (BCinp.is_open())
          cout << "\t\033[1;33mEBC nodeset file opened successfully\033[0m:\t" << boundaryConditions[bc_iter].BCNodeSet->nodeSetPath.c_str() << endl;
        else {
          cerr << "\t\033[41;37mERROR\033[0m\t: Unable to open NodeSet file " << boundaryConditions[bc_iter].BCNodeSet->nodeSetPath.c_str() << endl;
          exit(EXIT_FAILURE);
        }

        int tempNumBC;
        BCinp >> tempNumBC;
        NumBC += tempNumBC;
        BCid.reserve(NumBC*3);
        BCvalues.reserve(NumBC*3);
        for(int i = 0; i < tempNumBC; i++) {
          BCinp >> node;
          BCnodes.push_back(node);
          for (int j = 0; j < 3; j++) {
            BCid.push_back(node*3 + j);
            BCvalues.push_back(mainMesh->getX(node, j));
          }
        }
        BCinp.close();
        cout << "\t\033[1;33mEBC applied to \033[0m" << tempNumBC << " nodes." << endl;
        break;
      } // End of EBC
      case TORSIONALSPRING: {
	cout << "\t\033[1;33mApplying Torsional Spring BC.\033[0m" << endl;
        cout << "\t\033[1;33mTorsional node set used\033[0m:\t" << boundaryConditions[bc_iter].BCNodeSet->nodeSetName << endl;

        myModel.initTorsionalSpringBC(boundaryConditions[bc_iter].BCNodeSet->nodeSetPath.c_str(), boundaryConditions[bc_iter].BCProperties[0]);
        cout << "\t\033[1;33mTorsional node set path\033[0m:\t" << boundaryConditions[bc_iter].BCNodeSet->nodeSetPath << endl;
	cout << "\t\033[1;33mTorsional spring stiffness\033[0m:\t" << boundaryConditions[bc_iter].BCProperties[0] << endl;
        break;
      } // End of TORSIONALSPRING
      case SPRINGFOUNDATION: {
	cout << "\t\033[1;33mApplying Spring Foundation BC.\033[0m" << endl;
        cout << "\t\033[1;33mSpring foundation node set used\033[0m:\t" << boundaryConditions[bc_iter].BCNodeSet->nodeSetName << endl;
        // cout << boundaryConditions[bc_iter].surfMesh->getNumberOfNodes() << endl;
        myModel.initSpringBC(boundaryConditions[bc_iter].BCNodeSet->nodeSetPath.c_str(), boundaryConditions[bc_iter].surfMesh, boundaryConditions[bc_iter].BCProperties[0]);
        cout << "\t\033[1;33mSpring foundation node set path\033[0m:\t" << boundaryConditions[bc_iter].BCNodeSet->nodeSetPath << endl;
        // cout << "\t\033[1;33mSpring foundation surface mesh\033[0m:\t" << boundaryConditions[bc_iter].surfMesh->secondaryMeshName << endl;
        cout << "\t\033[1;33mSpring foundation stiffness\033[0m:\t" << boundaryConditions[bc_iter].BCProperties[0] << endl;
	break;
      } // End of SPRINGFOUNDATION
      case INTERACTION: {
        cout << "\t\033[1;33mApplying Interaction BC.\033[0m" << endl;
        cout << "\t\033[1;33mInteraction BC node set used\033[0m:\t" << boundaryConditions[bc_iter].BCNodeSet->nodeSetName << endl;

	myModel.initializeLennardJonesBC(boundaryConditions[bc_iter].BCNodeSet->nodeSetPath.c_str(), boundaryConditions[bc_iter].BCMesh, boundaryConditions[bc_iter].surfMesh, boundaryConditions[bc_iter].BCProperties[0], boundaryConditions[bc_iter].BCProperties[1], boundaryConditions[bc_iter].BCProperties[2]);
	cout << "\t\033[1;33mInteraction Search Radius\033[0m:\t" << boundaryConditions[bc_iter].BCProperties[0] << endl;
	cout << "\t\033[1;33mInteraction Energy\033[0m:\t" << boundaryConditions[bc_iter].BCProperties[1] << endl;
	cout << "\t\033[1;33mInteraction Minimum Distance\033[0m:\t" << boundaryConditions[bc_iter].BCProperties[2] << endl;
        break;
      }
      default:
	cerr << "\t\033[41;37mERROR\033[0m\t: Not sure how to implement boundary condition #" << bc_iter << "." << endl;
        exit(EXIT_FAILURE);
    } // end of switch statement
    cout << endl;
  } // end of bc_iter


  // Initialize Result
  uint myRequest;
  uint PbDoF = (mainMesh->getNumberOfNodes())*myModel.getDoFperNode();
  EigenResult myResults(PbDoF, 0);

  // Run Consistency check

  Real perturbationFactor = 0.1;
  myRequest = 7; // Check both Forces and Stiffness
  Real myH = 1e-6;
  Real myTol = 1e-7;
  // Before checking consistency, the perturbed deformation state must be 
  // set to the current deformation state.
  if (checkConsistencyFlag) {
    cout << "\033[1;33mChecking Model Consistency...\033[0m" << endl;
    myModel.checkConsistency(&myResults, perturbationFactor, myRequest, myH, myTol);
  }


  // Print initial configuration
  myModel.writeOutputVTK(outputString, 0);

  // Solver
  EigenNRsolver mySolver(&myModel, BCid, BCvalues, CHOL, NRTOL, NRIterations);

  /*
  // ********************************* // 
  // SOLVE:

  ind = 0;
  myModel.finalizeCompute();

  ofstream outVolume;
  outVolume.open(volumeFile.c_str());

  for (int s = 0; s < simulationTime/deltaT; s++)
  {
    cout << endl << "*** Step " << s << " ***" << endl;
    if (SpringBCflag) myModel.computeNormals();

    // Update pressure:
    if (pressureFlag) {
	double appliedPressure = timePressureData[s][1]/pressureDivider;
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
      vector <double> currentMyocardiumField(mainMesh->getNumberOfNodes() * mainMesh->getDimension(), 0.0);
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
    cout << "Reference Volume: " << myModel.computeRefVolume() << "\t Current Volume: " << myModel.computeCurrentVolume() << endl;
  }
  outVolume.close();

  */


  // Timing
  time (&end);
  cout << endl << "BVP solved in " << difftime(end,start) << " s" << endl;

  return 0;
} // End main()



// ******************** //
// Helper Functions
FEMesh* createFEMesh(string nodeFile, string eleFile, string meshName) {
  cout << "\033[1;33mCreating " << meshName << " FE-mesh object...\033[0m" << endl; 
  // Check to make sure Node and Element Files have been specified
  if (nodeFile.empty()) {
    cerr << "\033[41;37mERROR\033[0m\t: Node file not specified." << endl;
    exit(EXIT_FAILURE);
  }
  if (eleFile.empty()) {
    cerr << "\033[41;37mERROR\033[0m\t: Element connectivity file not specified." << endl;
    exit(EXIT_FAILURE);
  }
  FEMesh* mesh = new FEMesh(nodeFile, eleFile);
  cout << "\033[1;33mFE-mesh object Created successfully with\033[0m" << endl;
  cout << "\t\033[1;33mNode file\033[0m:\t\t" << nodeFile << endl;
  cout << "\t\033[1;33mEl-connectivity file\033[0m:\t" << eleFile << endl << endl;
  return mesh;
}


// Generating the Fiber Data
vector <vector <Vector3d> > readFiberData(string fiberFile, const FEMesh* mainMesh, bool writeFibers, string fiberPlotString) {
  cout << "\033[1;33mGenerating fiber data...\033[0m" << endl;

  // VTK File for Fibers
  vtkSmartPointer<vtkUnstructuredGrid> NodalPointGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkPoints> NodalPoints = vtkSmartPointer<vtkPoints>::New();

  vtkSmartPointer<vtkDoubleArray> fiberDirections = vtkSmartPointer<vtkDoubleArray>::New();
  fiberDirections->SetNumberOfComponents(3);
  fiberDirections->SetName("Fiber_Direction");
  
  vtkSmartPointer<vtkDoubleArray> sheetDirections = vtkSmartPointer<vtkDoubleArray>::New();
  sheetDirections->SetNumberOfComponents(3);
  sheetDirections->SetName("Sheet_Direction");

  vtkSmartPointer<vtkDoubleArray> sheetNormalDirections = vtkSmartPointer<vtkDoubleArray>::New();
  sheetNormalDirections->SetNumberOfComponents(3);
  sheetNormalDirections->SetName("Sheet_Normal_Direction");

  // Check to see if Fiber file has been specified
  bool fiberFileGiven = true;
  if (fiberFile.empty()) {
    cout << "\t\033[1;33mFiber file has not been specified.\033[0m" << endl;
    cout << "\t\033[1;33mFibers will be set to f = [1 0 0], s = [0 1 0], n = [0 0 1].\033[0m" << endl;
    fiberFileGiven = false;
  }

  // Open fiber file
  ifstream fiberFileStream;
  fiberFileStream.open(fiberFile.c_str());

  if (fiberFileStream.is_open())
    cout << "\t\033[1;33mFiber file\033[0m\t: " << fiberFile << endl;
  else {
    cout << "\t\033[41;37mERROR\033[0m\t: Cannot open fiber file " << fiberFile << endl;
    exit(EXIT_FAILURE);
  }

  vector <vector <Vector3d> > fiberData;
  
  for (int el_iter = 0; el_iter < mainMesh->getNumberOfElements(); el_iter++) {
    GeomElement* tempEle = mainMesh->getElements()[el_iter]; 
    vector<int> el_nodeIds = tempEle->getNodesID();
    for (int quad_iter = 0; quad_iter < tempEle->getNumberOfQuadPoints(); quad_iter++) {
      vector <Vector3d> el_vectors(3, Vector3d::Zero(3,1));
      if (fiberFileGiven)
      {
	string fiberLine;
	if (!getline(fiberFileStream, fiberLine)) {
	  cout << "\t\033[41;37mERROR\033[0m\t: Not enough lines in fiber file." << endl;
	  cout << "\t\t: Number of lines should equal NumElements * NumQuadPtsPerElement." << endl;
	  exit(EXIT_FAILURE);
	}

        trim(fiberLine);
        vector<string> fiberVals;
        split(fiberVals, fiberLine, is_any_of("\t "), token_compress_on);
        // Load Fiber Data:
        el_vectors[0][0] = atof(fiberVals[0].c_str()); el_vectors[0][1] = atof(fiberVals[1].c_str()); el_vectors[0][2] = atof(fiberVals[2].c_str());
	el_vectors[1][0] = atof(fiberVals[3].c_str()); el_vectors[1][1] = atof(fiberVals[4].c_str()); el_vectors[1][2] = atof(fiberVals[5].c_str());
        el_vectors[2][0] = atof(fiberVals[6].c_str()); el_vectors[2][1] = atof(fiberVals[7].c_str()); el_vectors[2][2] = atof(fiberVals[8].c_str());
      }
      else
      {
        el_vectors[0] << 1., 0., 0.;
        el_vectors[1] << 0., 1., 0.;
        el_vectors[2] << 0., 0., 1.;
      }
      fiberData.push_back(el_vectors);

      // Write Fibers to file?
      if(writeFibers) {
	// Compute the location of the quadrature point
	Real quadPointX = 0.0; Real quadPointY = 0.0; Real quadPointZ = 0.0;
 	float tempPoint[3] = {0.0};
	for (int el_node_iter = 0; el_node_iter < el_nodeIds.size(); el_node_iter++) {
	  tempPoint[0] += tempEle->getN(quad_iter, el_node_iter) * mainMesh->getX(el_nodeIds[el_node_iter], 0);
	  tempPoint[1] += tempEle->getN(quad_iter, el_node_iter) * mainMesh->getX(el_nodeIds[el_node_iter], 1);
	  tempPoint[2] += tempEle->getN(quad_iter, el_node_iter) * mainMesh->getX(el_nodeIds[el_node_iter], 2);
	}
	NodalPoints->InsertNextPoint(tempPoint);
	double tempFiberDirection[3] = {el_vectors[0][0], el_vectors[0][1], el_vectors[0][2]};
	double tempSheetDirection[3] = {el_vectors[1][0], el_vectors[1][1], el_vectors[1][2]};
	double tempSheetNormalDirection[3] = {el_vectors[2][0], el_vectors[2][1], el_vectors[2][2]};
	fiberDirections->InsertNextTuple(tempFiberDirection);
	sheetDirections->InsertNextTuple(tempSheetDirection);
	sheetNormalDirections->InsertNextTuple(tempSheetNormalDirection);
      } // writeFibers flag
    } // Loop over quad Points
  } // Loop over elements

  if (writeFibers) {
    cout << "\t\033[1;33mWriting fiber VTK file\033[0m\t: " << fiberPlotString << endl;
    NodalPointGrid->SetPoints(NodalPoints);
    NodalPointGrid->GetPointData()->AddArray(fiberDirections);
    NodalPointGrid->GetPointData()->AddArray(sheetDirections);
    NodalPointGrid->GetPointData()->AddArray(sheetNormalDirections);

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> NodalPointGridWriter = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    NodalPointGridWriter->SetFileName(fiberPlotString.c_str());
    NodalPointGridWriter->SetInput(NodalPointGrid);
    NodalPointGridWriter->Write();
  }

  cout << endl;
  
  return fiberData;  
}


materialTypes determineMaterialName(string materialName) {
  to_upper(materialName);
  if (materialName == "COMPNEOHOOKEAN") return COMPNEOHOOKEAN;
  else if (materialName == "JACOBIAN") return JACOBIAN;
  else if (materialName == "LINYINACTIVE_COMPRESSIBLE") return LINYINACTIVE_COMPRESSIBLE;
  else if (materialName == "HUMPHREY_COMPRESSIBLE") return HUMPHREY_COMPRESSIBLE;

  // ForceVel Potentials
  else if (materialName == "HILLFORCEVELPOTENTIAL") return HILLFORCEVELPOTENTIAL;
  else if (materialName == "APFORCEVELPOTENTIAL") return APFORCEVELPOTENTIAL;
  
  // Viscous Potentials
  else if (materialName == "BLANKVISCOUSPOTENTIAL") return BLANKVISCOUSPOTENTIAL;
  else if (materialName == "NEWTONIANVISCOUSPOTENTIAL") return NEWTONIANVISCOUSPOTENTIAL;
  else {
    cout << "\t\033[41;37mERROR\033[0m\t: Unknown material type " << materialName << endl;
    exit(EXIT_FAILURE);
  }
}


// Helper function to create the Material
MechanicsMaterial* createMechanicsMaterial(materialTypes materialName, vector<double> materialProperties, vector<Vector3d> fibers) {
  MechanicsMaterial* mechMaterial;
  switch(materialName) {
    case COMPNEOHOOKEAN :
      assert(materialProperties.size() == 2); // Check to make sure there are enough material properties specified
      mechMaterial = new CompNeoHookean(0, materialProperties[0], materialProperties[1]);
      break;
    case JACOBIAN :
      assert(materialProperties.size() == 0); // Check to make sure there are enough material properties specified
      mechMaterial = new Jacobian(0);
      break;
    case LINYINACTIVE_COMPRESSIBLE :
      assert(materialProperties.size() == 5); // Check to make sure there are enough material properties specified
      mechMaterial = new LinYinActive_Compressible(0, materialProperties[0], materialProperties[1], materialProperties[2], materialProperties[3], materialProperties[4], fibers);
      break;
    case HUMPHREY_COMPRESSIBLE :
      assert(materialProperties.size() == 7); // Check to make sure there are enough material properties specified
      mechMaterial = new Humphrey_Compressible(0, materialProperties[0], materialProperties[1], materialProperties[2], materialProperties[3], materialProperties[4], materialProperties[5], materialProperties[6], fibers);
      break;
    default:
      cout << "\t\033[41;37mERROR\033[0m\t: Unknown mechanics material type " << materialName << endl;
      exit(EXIT_FAILURE);
  }
  return mechMaterial;
}

Potential* createKineticPotentialMaterial(materialTypes materialName, vector<double> materialProperties) {
  Potential* kineticPotential;
  switch(materialName) {
    case HILLFORCEVELPOTENTIAL :
      assert(materialProperties.size() == 3); // Check to make sure there are enough material properties specified
      kineticPotential = new HillForceVelPotential(materialProperties[0], materialProperties[1], materialProperties[2]);
      break;
    case APFORCEVELPOTENTIAL :
      assert(materialProperties.size() == 3); // Check to make sure there are enough material properties specified
      kineticPotential = new APForceVelPotential(materialProperties[0], materialProperties[1], materialProperties[2]);
      break;
    default:
      cout << "\t\033[41;37mERROR\033[0m\t: Unknown kinetic potential material type " << materialName << endl;
      exit(EXIT_FAILURE);
  }
  return kineticPotential;
}


ViscousPotential* createViscousPotentialMaterial(materialTypes materialName, vector<double> materialProperties) {
  ViscousPotential* viscousPotential;
  switch(materialName) {
    case BLANKVISCOUSPOTENTIAL :
      assert(materialProperties.size() == 0); // Check to make sure there are enough material properties specified
      viscousPotential = new BlankViscousPotential;
      break;
    case NEWTONIANVISCOUSPOTENTIAL :
      assert(materialProperties.size() == 2); // Check to make sure there are enough material properties specified
      viscousPotential = new NewtonianViscousPotential(materialProperties[0], materialProperties[1]);
      break;
    default:
      cout << "\t\033[41;37mERROR\033[0m\t: Unknown viscous potential material type " << materialName << endl;
      exit(EXIT_FAILURE);
  }
  return viscousPotential;
}


boundaryConditionTypes determineBoundaryConditionType(string boundaryConditionName) {
  to_upper(boundaryConditionName);
  if (boundaryConditionName == "EBC") return EBC;
  else if (boundaryConditionName == "TORSIONALSPRING") return TORSIONALSPRING;
  else if (boundaryConditionName == "SPRINGFOUNDATION") return SPRINGFOUNDATION;
  else if (boundaryConditionName == "INTERACTION") return INTERACTION;
  else {
    cout << "\t\033[41;37mERROR\033[0m\t: Unknown boundary condition type " << boundaryConditionName << endl;
    exit(EXIT_FAILURE);
  }
}

FEMesh* findMeshObject(string meshName, vector<secondaryMesh> &secondaryMeshObjects) {
  for (int i = 0; i < secondaryMeshObjects.size(); i++) {
    to_upper(meshName);
    string tempMeshName = secondaryMeshObjects[i].secondaryMeshName;
    to_upper(tempMeshName);
    if(meshName == tempMeshName) {
      return secondaryMeshObjects[i].secondaryMeshObject;
    }
  }
  cout << "\t\033[41;37mERROR\033[0m\t: Unknown Mesh Name " << meshName << endl;
  exit(EXIT_FAILURE);
}

nodeSet* findNodeSet(string nodeSetName, vector<nodeSet> &nodeSetObjects) {
  for (int i = 0; i < nodeSetObjects.size(); i++) {
    to_upper(nodeSetName);
    string tempNodeSetName = nodeSetObjects[i].nodeSetName;
    to_upper(tempNodeSetName);
    if (nodeSetName == tempNodeSetName)
      return &nodeSetObjects[i];
  }
  cout << "\t\033[41;37mERROR\033[0m\t: Unknown NodeSet Name " << nodeSetName << endl;
  exit(EXIT_FAILURE);
}


void readTwoColumnDataFile(string twoColumnDataFilePath, vector <vector <double> > &twoColumnData) {
  ifstream pressureFileStream;
  pressureFileStream.open(twoColumnDataFilePath.c_str());
  
  double tempCol1, tempCol2;
  if (pressureFileStream.is_open())
  {
    cout << "\t\033[1;33mReading pressure file \033[0m" << twoColumnDataFilePath << endl;
    while (pressureFileStream >> tempCol1 >> tempCol2) {
      vector <double> tempTwoColData;
      tempTwoColData.push_back(tempCol1);
      tempTwoColData.push_back(tempCol2);
      twoColumnData.push_back(tempTwoColData);
    }
  }
  else {
    cout << "\t\033[41;37mERROR\033[0m\t: Unable to open file " << twoColumnDataFilePath << endl;
    exit(EXIT_FAILURE);
  }
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
