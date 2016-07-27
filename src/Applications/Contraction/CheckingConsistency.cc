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
#include "Humphrey.h"
#include "LinYinActive.h"

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
  bool fibersAcrossWall = true;
  double theta_min = 60.0;
  double theta_max = -60.0;

  // Conduction Velocity (cm/ms)
  double cv = 0.006;
  
  // Simulation Time (in ms)
  double simTime = 10;  
  
  // Time Step (in ms)
  double deltaT = 0.01;

  // Run Only Consistency Check?
  bool consistencyCheckOnly = true;

  // OutputString
  string outputString = "Strip_LV_CV_";

  // Fiber Visualization String
  string fiberPlotString = "Strip_LV_CV_Fiber.vtk";

  // Initialize Mesh
  // Assumptions to use this main as is: strip has a face at z=0; tetrahedral mesh
  FEMesh Cube("Mesh/Cube/Cube1.node", "Mesh/Cube/Cube1.ele");
  FEMesh surfMesh("Mesh/Cube/Cube1.node", "Mesh/Cube/Cube1Surf.ele");
  Real xmax = 1.0;
  // FEMesh Cube("Strip36.node", "Strip36.ele");
  // FEMesh surfMesh("Strip36.node", "Strip36Surf.ele");
  // FEMesh Cube("Strip144.node", "Strip144.ele");
  // FEMesh surfMesh("Strip144.node", "Strip144Surf.ele");
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
  // Humphrey PassiveMat(0, 1.0, 2.0, 3.0, 4.0, 5.0);
  // LinYinActive ActiveMat(0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
  CompNeoHookean ActiveMat(0, 4.0, 0.4);
  // APForceVelPotential TestPotential(1.0, 500.0);
  // HillForceVelPotential TestPotential(4.4*pow(10,-3), .01*0.59, 25);
  // BlankPotential TestPotential;

  // BlankViscousPotential ViscPotential;
  NewtonianViscousPotential ViscPotential(1000, 0.5);

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
    APForceVelPotential* TestPotential = new APForceVelPotential(4.0, 50000.0, 3.0);
    PlasticMaterial* PlMat = new PlasticMaterial(el_iter, &ActiveMat, &PassiveMat, TestPotential, &ViscPotential);
    PlMat->setDirectionVectors(el_vectors);
    PlMat->setHardeningParameters(HardParam);
    PlMat->setActiveDeformationGradient(Matrix3d::Identity(3,3));
    PlMat->setTotalDeformationGradient(Matrix3d::Identity(3,3));

    PlMat->setTimestep(deltaT);
    PlMat->setActivationMultiplier(0.0);

    PLmaterials.push_back(PlMat);
  }

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

 
  // Testing four different conditions for consistency:
  // 1. \mathbf{F} = \mathbf{I} and \mathbf{Q} = \mathbf{0}
  //    u_{ia} = 0
  cout << "********************************" << endl;
  cout << "Checking Consistency Test 1 of 4" << endl;
  cout << "u_{ia} = 0, mathbf{F} = mathbf{I}, and mathbf{Q} = mathbf{0}" << endl;
  cout << "********************************" << endl;
  for (int el_iter = 0; el_iter < meshElements.size(); el_iter++)
    {
      PlasticMaterial* tempPlasticMaterial = dynamic_cast<PlasticMaterial*>(PLmaterials[el_iter]);
      tempPlasticMaterial->setTotalDeformationGradient(Matrix3d::Identity(3,3));
    }
  myModel.checkConsistency(&myResults, 0.0, (FORCE|STIFFNESS), myH, myTol);

  // 2. \mathbf{F} = \mathbf{I} and \mathbf{Q} = \mathbf{0}
  //    u_{ia} = rand  - Perturb u but not \mathbf{F}
  cout << "********************************" << endl;
  cout << "Checking Consistency Test 2 of 4" << endl;
  cout << "u_{ia} = rand and mathbf{F} = mathbf{I} and mathbf{Q} = mathbf{0}" << endl;
  cout << "********************************" << endl;
  // Set random displacement:
  srand(time(NULL));
  vector<Real> perturb(Cube.getNumberOfNodes() * NodeDoF, 0.0);
  for (uint a = 0; a < Cube.getNumberOfNodes(); a++)
    {
      for (uint i = 0; i < NodeDoF; i++)
	{
	  Real randomNum = perturbationFactor*(Real(rand())/RAND_MAX - 0.5);
	  perturb[a*NodeDoF + i] = randomNum;
	  myModel.linearizedUpdate(a, i, randomNum);
	}
    }

  for (int el_iter = 0; el_iter < meshElements.size(); el_iter++)
    {
      PlasticMaterial* tempPlasticMaterial = dynamic_cast<PlasticMaterial*>(PLmaterials[el_iter]);
      tempPlasticMaterial->setTotalDeformationGradient(Matrix3d::Identity(3,3));
      tempPlasticMaterial->setHardeningParameters(HardParam);
      tempPlasticMaterial->setActiveDeformationGradient(Matrix3d::Identity(3,3));
      tempPlasticMaterial->setTotalDeformationGradient(Matrix3d::Identity(3,3));
    }
    
  myModel.checkConsistency(&myResults, 0.0, (FORCE|STIFFNESS), myH, myTol);
    
  // Reset field to initial values
  for(uint a = 0; a < Cube.getNumberOfNodes(); a++)
    for(uint i = 0; i < NodeDoF; i++)
      myModel.linearizedUpdate(a, i, -perturb[a*NodeDoF + i]);

    
  // 3. \mathbf{F} = \mathbf{I} and \mathbf{Q} = \mathbf{0}
  //    u_{ia} = rand -> \mathbf{F} and \mathbf{Q} match u_{ia}
  cout << "********************************" << endl;
  cout << "Checking Consistency Test 3 of 4" << endl;
  cout << "u_{ia} = rand and mathbf{F,Q} match u_{ia}" << endl;
  cout << "********************************" << endl;
  // Set random displacement:
  srand(time(NULL));
  // vector<Real> perturb(Cube.getNumberOfNodes() * NodeDoF, 0.0);
  std::fill(perturb.begin(), perturb.end(), 0.0);
  for (uint a = 0; a < Cube.getNumberOfNodes(); a++)
    {
      for (uint i = 0; i < NodeDoF; i++)
	{
	  Real randomNum = perturbationFactor*(Real(rand())/RAND_MAX - 0.5);
	  perturb[a*NodeDoF + i] = randomNum;
	  myModel.linearizedUpdate(a, i, randomNum);
	}
    }
    
  vector<double> tempField(Cube.getNumberOfNodes() * NodeDoF, 0.0);
  myModel.getField(tempField);

  for (int el_iter = 0; el_iter < meshElements.size(); el_iter++)
    {
      PlasticMaterial* tempPlasticMaterial = dynamic_cast<PlasticMaterial*>(PLmaterials[el_iter]);
      GeomElement* geomEl = meshElements[el_iter];
      

      
      // Set Flist to current deformation gradient Fn:
      for (int quadPoint_iter = 0; quadPoint_iter < geomEl->getNumberOfQuadPoints(); quadPoint_iter++)
	{
	  // Compute Deformation Gradient at Quad Point
	  Matrix3d Fcomp = Matrix3d::Zero();
	  // Compute F at all quadrature points
	  const uint numQP = geomEl->getNumberOfQuadPoints();
	  const uint dim   = Cube.getDimension();

	  const vector<int > & NodesID = geomEl->getNodesID();
	  const uint nodeNum = NodesID.size();

	  for(uint i = 0; i < dim; i++) 
	    for(uint J = 0; J < dim; J++)
	      for(uint a = 0; a < nodeNum; a++)
		Fcomp(i,J) += tempField[NodesID[a]*dim + i] * geomEl->getDN(quadPoint_iter, a, J);
	
	  MechanicsMaterial::FKresults FKres;
	  FKres.request = 0;
	  // TODO: Only works for one quad point because tempPlasticMaterial is not indexed by q
	  // TODO: doesn't use el_vectors[0]
	  tempPlasticMaterial->compute(FKres, Fcomp, NULL); 
	  tempPlasticMaterial->updateStateVariables();
	  cout << tempPlasticMaterial->getTotalDeformationGradient() << endl;
	  cout << tempPlasticMaterial->getActiveDeformationGradient() << endl;
	}
    }
    
  myModel.checkConsistency(&myResults, 0.0, (FORCE|STIFFNESS), myH, myTol);
  // Reset field to initial values
  for(uint a = 0; a < Cube.getNumberOfNodes(); a++)
    for(uint i = 0; i < NodeDoF; i++)
      myModel.linearizedUpdate(a, i, -perturb[a*NodeDoF + i]);

  // 4. \mathbf{F} = rand and \mathbf{Q} = rand
  //    u_{ia} = rand
  cout << "********************************" << endl;
  cout << "Checking Consistency Test 4 of 4" << endl;
  cout << "u_{ia} = rand, mathbf{F,Q} = rand" << endl;
  cout << "********************************" << endl;
  // Set random displacement:
  srand(time(NULL));
  // vector<double> perturb(Cube.getNumberOfNodes() * NodeDoF, 0.0);
  std::fill(perturb.begin(), perturb.end(), 0.0);
  for (uint a = 0; a < Cube.getNumberOfNodes(); a++)
    {
      for (uint i = 0; i < NodeDoF; i++)
	{
	  Real randomNum = perturbationFactor*(Real(rand())/RAND_MAX - 0.5);
	  perturb[a*NodeDoF + i] = randomNum;
	  myModel.linearizedUpdate(a, i, randomNum);
	}
    }
    
  // vector<double> tempField (Cube.getNumberOfNodes() * NodeDoF, 0.0);
  std::fill(tempField.begin(), tempField.end(), 0.0);
  myModel.getField(tempField);

  for (int el_iter = 0; el_iter < meshElements.size(); el_iter++)
    {
      PlasticMaterial* tempPlasticMaterial = dynamic_cast<PlasticMaterial*>(PLmaterials[el_iter]);
      GeomElement* geomEl = meshElements[el_iter];

      Matrix3d Fcomp = Matrix3d::Identity();
            
      srand(time(NULL));
      for (unsigned int i = 0; i<3; i++) 
	for (unsigned int J = 0; J<3; J++) 
	  Fcomp(i,J) += 0.1*(double(rand())/RAND_MAX); 

      Vector3d HardParamRandom(0.0, 0.0, 0.0);
      for (unsigned int i = 0; i < 3; i++)
	HardParamRandom[i] = 0.1*(double(rand())/RAND_MAX);
	
      tempPlasticMaterial->setTotalDeformationGradient(Fcomp);
      tempPlasticMaterial->setHardeningParameters(HardParamRandom);
    }
    
  myModel.checkConsistency(&myResults, 0.0, (FORCE|STIFFNESS), myH, myTol);
  // Reset field to initial values
  for(uint a = 0; a < Cube.getNumberOfNodes(); a++)
    for(uint i = 0; i < NodeDoF; i++)
      myModel.linearizedUpdate(a, i, -perturb[a*NodeDoF + i]);

  // 5. \mathbf{F} = \mathbf{I} and \mathbf{Q} = \mathbf{0}
  //    u_{ia} = rand -> \mathbf{F} and \mathbf{Q} match u_{ia}
  cout << "********************************" << endl;
  cout << "Checking Consistency Test 5 of 5" << endl;
  cout << "u_{ia} = rand, mathbf{F,Q} = rand" << endl;
  cout << "Model consistency followed by Material Consistency for same Fn+1" << endl;
  cout << "********************************" << endl;
  // Initialize field to be 0 displacement:
  cout << "** Setting Field to 0 " << endl;

  std::fill(tempField.begin(), tempField.end(), 0.0);
  myModel.getField(tempField);

  for (int el_iter = 0; el_iter < meshElements.size(); el_iter++)
  {
    PlasticMaterial* tempPlasticMaterial = dynamic_cast<PlasticMaterial*>(PLmaterials[el_iter]);
    GeomElement* geomEl = meshElements[el_iter];
    
    // Set Flist to current deformation gradient Fn:
    for (int quadPoint_iter = 0; quadPoint_iter < geomEl->getNumberOfQuadPoints(); quadPoint_iter++)
    {
      // Compute Deformation Gradient at Quad Point
      Matrix3d Fcomp = Matrix3d::Zero();

      // Compute F at all quadrature points
      const uint numQP = geomEl->getNumberOfQuadPoints();
      const uint dim   = Cube.getDimension();

      const vector<int > & NodesID = geomEl->getNodesID();
      const uint nodeNum = NodesID.size();
      
      for(uint i = 0; i < dim; i++) 
	for(uint J = 0; J < dim; J++)
	  for(uint a = 0; a < nodeNum; a++)
	    Fcomp(i,J) += tempField[NodesID[a]*dim + i] * geomEl->getDN(quadPoint_iter, a, J);
      
      MechanicsMaterial::FKresults FKres;
      FKres.request = 0;
      
      tempPlasticMaterial->setHardeningParameters(HardParam);
      tempPlasticMaterial->setTotalDeformationGradient(Matrix3d::Identity());
      tempPlasticMaterial->setActiveDeformationGradient(Matrix3d::Identity());
      
      // TODO: Only works for one quad point because tempPlasticMaterial is not indexed by q
      // TODO: doesn't use el_vectors[0]
      tempPlasticMaterial->compute(FKres, Fcomp, NULL); 
      tempPlasticMaterial->updateStateVariables();
      // cout << tempPlasticMaterial->getTotalDeformationGradient() << endl;
      // cout << tempPlasticMaterial->getActiveDeformationGradient() << endl;
      // cout << tempPlasticMaterial->getHardeningParameters() << endl;
    }
  }
  cout << "*** Initial State Set. F = I, Q = {0,0,0}, Fa = I." << endl;

  // Perturb the field

  srand(time(NULL));
  // vector<double> perturb(Cube.getNumberOfNodes() * NodeDoF, 0.0);
  std::fill(perturb.begin(), perturb.end(), 0.0);
  for (uint a = 0; a < Cube.getNumberOfNodes(); a++)
  {
    for (uint i = 0; i < NodeDoF; i++)
    {
      Real randomNum = perturbationFactor*(Real(rand())/RAND_MAX - 0.5);
      perturb[a*NodeDoF + i] = randomNum;
      myModel.linearizedUpdate(a, i, randomNum);
    }
  }

  std::fill(tempField.begin(), tempField.end(), 0.0);
  myModel.getField(tempField);
  
  for (int el_iter = 0; el_iter < meshElements.size(); el_iter++)
  {
    PlasticMaterial* tempPlasticMaterial = dynamic_cast<PlasticMaterial*>(PLmaterials[el_iter]);
    GeomElement* geomEl = meshElements[el_iter];
    
    // Set Flist to current deformation gradient Fn:
    for (int quadPoint_iter = 0; quadPoint_iter < geomEl->getNumberOfQuadPoints(); quadPoint_iter++)
    {
      // Compute Deformation Gradient at Quad Point
      Matrix3d Fcomp = Matrix3d::Zero();

      // Compute F at all quadrature points
      const uint numQP = geomEl->getNumberOfQuadPoints();
      const uint dim   = Cube.getDimension();

      const vector<int > & NodesID = geomEl->getNodesID();
      const uint nodeNum = NodesID.size();
      
      for(uint i = 0; i < dim; i++) 
	for(uint J = 0; J < dim; J++)
	  for(uint a = 0; a < nodeNum; a++)
	    Fcomp(i,J) += tempField[NodesID[a]*dim + i] * geomEl->getDN(quadPoint_iter, a, J);
      
      MechanicsMaterial::FKresults FKres;
      FKres.request = 7;

      cout << "Checking material consistency with perturbed deformation gradient. " << endl;
      tempPlasticMaterial->checkConsistency(FKres, Fcomp);
            
      // TODO: Only works for one quad point because tempPlasticMaterial is not indexed by q
      // TODO: doesn't use el_vectors[0]
      // tempPlasticMaterial->compute(FKres, Fcomp, NULL); 
      // tempPlasticMaterial->updateStateVariables();
      // cout << tempPlasticMaterial->getTotalDeformationGradient() << endl;
      // cout << tempPlasticMaterial->getActiveDeformationGradient() << endl;
      // cout << tempPlasticMaterial->getHardeningParameters() << endl;
    }
  }
  myResults.resetResidualToZero();
  myResults.resetStiffnessToZero();

  myModel.checkConsistency(&myResults, 0.0, (STIFFNESS), myH, myTol);
  // Reset field to initial values
  for(uint a = 0; a < Cube.getNumberOfNodes(); a++)
    for(uint i = 0; i < NodeDoF; i++)
      myModel.linearizedUpdate(a, i, -perturb[a*NodeDoF + i]);

  {
    // 6. \mathbf{F} = \mathbf{I} and \mathbf{Q} = \mathbf{0}
    //    u_{ia} = rand -> \mathbf{F} and \mathbf{Q} match u_{ia}
    cout << "********************************" << endl;
    cout << "Checking Consistency Test 6" << endl;
    cout << "u_{ia} = rand" << endl;
    cout << "Checking Model consistency for a normal mechanics material" << endl;
    cout << "********************************" << endl;

    // Test to ensure that Model Consistency passes for a normal material
    vector <MechanicsMaterial*> MechMaterialVector;
    MechMaterialVector.reserve(NumMat);

    CompNeoHookean* TestCompNeoHookeanMaterial = new CompNeoHookean(0, 4.0, 0.4);
    for (int el_iter = 0; el_iter < meshElements.size(); el_iter++) {
      // PLmaterials.push_back(&PassiveMat);
      MechMaterialVector.push_back(TestCompNeoHookeanMaterial);
    }
    MechanicsModel myModel2(&Cube, MechMaterialVector, NodeDoF, PressureFlag, &surfMesh,
			     NodalForcesFlag);
    myModel2.updatePressure(Pressure);
    myModel2.updateNodalForces(&ForcesID, &Forces);

    std::fill(tempField.begin(), tempField.end(), 0.0);
    myModel2.getField(tempField);

    // Perturb the field

    srand(time(NULL));
    // vector<double> perturb(Cube.getNumberOfNodes() * NodeDoF, 0.0);
    std::fill(perturb.begin(), perturb.end(), 0.0);
    for (uint a = 0; a < Cube.getNumberOfNodes(); a++)
      {
	for (uint i = 0; i < NodeDoF; i++)
	  {
	    Real randomNum = perturbationFactor*(Real(rand())/RAND_MAX - 0.5);
	    perturb[a*NodeDoF + i] = randomNum;
	    myModel2.linearizedUpdate(a, i, randomNum);
	  }
      }

    myModel2.checkConsistency(&myResults, 0.0, (FORCE|STIFFNESS), myH, myTol);
    // Reset field to initial values
    for(uint a = 0; a < Cube.getNumberOfNodes(); a++)
      for(uint i = 0; i < NodeDoF; i++)
	myModel2.linearizedUpdate(a, i, -perturb[a*NodeDoF + i]);
    
  }

  // Timing
  time (&end);
  cout << endl << "Consistency Check completed in " << difftime(end,start) << " s" << endl;



  return 0;
}
