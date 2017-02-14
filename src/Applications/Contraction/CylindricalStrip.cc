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

int main(int argc, char** argv)
{
  cout << string(50, '\n'); // Clear Screen

  // Timing
  time_t start, end;
  time(&start);

  // Pressure File
  bool pressureFlag = true;
  string pressureFile = "Pressure_1ms.dat";
  
  // Simulation Time (in ms)
  double simTime = 200;
  
  // Time Step (in ms)
  double deltaT = 0.01;

  // OutputString
  string outputString = "/u/project/cardio/adityapo/ScratchResults/CylindricalStrip/Cylinder";

  // Spring BC:
  int SpringBCflag = 1;

  // Spring Stiffnes
  Real SpringK = 1.0e9; // 1.0e4;
  
  // Torsional Spring BC:
  int torsionalSpringBCflag = 0;

  // Torsional Spring Stiffness:
  int TorsionalSpringK = 0.0e6;

  // Initialize Mesh
  // Assumptions to use this main as is: strip has a face at z=0; tetrahedral mesh
  FEMesh myMesh("Mesh/CylindricalStripMesh/Cylindrical_Mesh.node", "Mesh/CylindricalStripMesh/Cylindrical_Mesh.ele");
  FEMesh endocardialMesh("Mesh/CylindricalStripMesh/Cylindrical_Mesh.node", "Mesh/CylindricalStripMesh/Cylindrical_Mesh.EndocardiumElset");
  FEMesh springBCSurfMesh("Mesh/CylindricalStripMesh/Cylindrical_Mesh.node", "Mesh/CylindricalStripMesh/Cylindrical_Mesh.EpiBotElset");
  string EBCfile = "Mesh/CylindricalStripMesh/Cylindrical_Mesh.Null.bc";
  // EBCfile = "Mesh/CylindricalStripMesh/Cylindrical_Mesh.TopSurfaceNodeset";
  string linearSpringBC = "Mesh/CylindricalStripMesh/Cylindrical_Mesh.EpiBotNodeset";
  string torsionalSpringBC = "Mesh/CylindricalStripMesh/Cylindrical_Mesh.BottomSurfaceNodeset"; // Mesh/EllipsoidMesh/Small_A.torsionalSpringNodes";

  cout << endl;
  cout << "Number Of Nodes   : " << myMesh.getNumberOfNodes() << endl;
  cout << "Number Of Element : " << myMesh.getNumberOfElements() << endl;
  cout << "Mesh Dimension    : " << myMesh.getDimension() << endl << endl;
  
  // Initialize Material
  vector <Vector3d> el_vectors(3, Vector3d::Zero(3,1));
  el_vectors[0] << 1., 0., 0.;
  el_vectors[1] << 0., 1., 0.;
  el_vectors[2] << 0., 0., 1.;

  uint NumMat =  myMesh.getNumberOfElements();
  vector<GeomElement*> meshElements = myMesh.getElements();
  
  vector<MechanicsMaterial * > PLmaterials;
  PLmaterials.reserve(NumMat);

  APForceVelPotential TestPotential(4.0, 0.0, 3.0);	// 50.0 for 2nd parameter, force
  BlankViscousPotential ViscPotential;
  Vector3d HardParam(0.,0.,0.);

  Humphrey_Compressible* PassiveMat = new Humphrey_Compressible(0, 15.98, 55.85, 0.0, -33.27, 30.21, 3.590, 64.62, el_vectors);
  LinYinActive_Compressible* ActiveMat = new LinYinActive_Compressible(0, -38.70, 40.83, 25.12, 9.51, 171.18, el_vectors);
  // CompNeoHookean* PassiveMat = new CompNeoHookean(0, 10.0, 10.0);  

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

  for (int el_iter = 0; el_iter < meshElements.size(); el_iter++) {
    PlasticMaterial* PlMat = new PlasticMaterial(el_iter, ActiveMat, PassiveMat, &TestPotential, &ViscPotential);
    PlMat->setDirectionVectors(el_vectors);
    PlMat->setHardeningParameters(HardParam);
    PlMat->setActiveDeformationGradient(Matrix3d::Identity(3,3));
    PlMat->setTotalDeformationGradient(Matrix3d::Identity(3,3));
    PlMat->setTimestep(deltaT);
    PlMat->setActivationMultiplier(0.0);

    PLmaterials.push_back(PlMat);
    // PLmaterials.push_back(PassiveMat);
  }

  // Initialize Model
  int NodeDoF = 3;
  int PressureFlag = 1;
  Real Pressure = 0.0;
  int NodalForcesFlag = 0;
  vector<int > ForcesID;
  vector<Real > Forces;
  MechanicsModel myModel(&myMesh, PLmaterials, NodeDoF, PressureFlag, &endocardialMesh,
			 NodalForcesFlag, SpringBCflag);
  myModel.updatePressure(Pressure);
  myModel.updateNodalForces(&ForcesID, &Forces);

  if (SpringBCflag)
    myModel.initSpringBC(linearSpringBC, &springBCSurfMesh, SpringK);

  if (torsionalSpringBCflag)
    myModel.initTorsionalSpringBC(torsionalSpringBC, TorsionalSpringK);



  // Initialize Result
  uint myRequest;
  uint PbDoF = (myMesh.getNumberOfNodes())*myModel.getDoFperNode();
  EigenResult myResults(PbDoF, NumMat*2);

  // Run Consistency check

  Real perturbationFactor = 0.1;
  myRequest = 7; // Check both Forces and Stiffness
  Real myH = 1e-6;
  Real myTol = 1e-7;

  // Before checking consistency, the perturbed deformation state must be 
  // set to the current deformation state.

  // myModel.checkConsistency(&myResults, perturbationFactor, myRequest, myH, myTol);

  // Print initial configuration
  myModel.writeOutputVTK(outputString, 0);

  // EBC
  cout << "********" << " Setting up EBCs " << "********" << endl;
  int NumBC = 0, node = 0, ind = 0;;
  vector<int > BCnodes; 
  vector<int > BCid;
  vector<Real > BCvalues;

  ifstream BCinp(EBCfile.c_str());
   
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
      BCvalues.push_back(myMesh.getX(node, j));
      cout << BCid[ind] << " " <<  BCvalues[ind] << endl;
      ind++;
    }
  }

  // Solver
  Real NRtol = 1.0e-4;
  uint NRmaxIter = 100;
  EigenNRsolver mySolver(&myModel, BCid, BCvalues, CHOL, NRtol, NRmaxIter);

  
  // SOLVE:
  ind = 0;

  for (int s = 0; s < 400; s++)
  {
    cout << "Step " << s << endl;
    if (SpringBCflag) myModel.computeNormals();

    // Update pressure	
    if (pressureFlag)
      myModel.updatePressure(1.0/20. * pressureData[s][1]);

    for (int k = 0; k < NumMat; k++)
	(PLmaterials[k])->setActivationMultiplier(0.0);
      
    ind++;
    mySolver.solve(DISP);
    myModel.setPrevField();

    // Update State Variables:     
    for (int k = 0; k < NumMat; k++)
      (PLmaterials[k])->updateStateVariables();

    cout << "State Variables Updated." << endl;
    
    // Write Output
    myModel.writeOutputVTK(outputString, ind);
    cout << "Output Written for step." << endl;
  }


  // Timing
  time (&end);
  cout << endl << "BVP solved in " << difftime(end,start) << " s" << endl;

  return 0;
}
