#include "State.h"
#include "EigenNRsolver.h"
#include "LBFGSB.h"
#include "FEMesh.h"
#include "Model.h"
#include "MechanicsBody.h"
#include "PotentialBody.h"
#include "TorsionalBody.h"
#include "PreferredNormalBody.h"
#include "ImposedKinematics.h"
#include "Harmonic.h"


using namespace voom;

int main(int argc, char** argv)
{

  // Timing
  time_t start, end;
  time(&start);

  if(argc < 2){
    cout << "Input file missing." << endl;
    return(0);
  }

       



  // Problem parameters
     int Nsteps = 20;
     // Imposed Kinematics Material
     
     Real ERRa =  0.36, ERRb =  0.09;
     Real ECCa = -0.15, ECCb = -0.10;
     Real ELLa = -0.19, ELLb = -0.03;

     Real EFF = -0.15; 
     
     Real Beta  = 1.0;
     Real Gamma = 1.0;
     int Ndir = 4;
     vector<Real > Alphas(Ndir, 1.0), Stretches(Ndir, 1.0); 
     Real Repi = 4.5, Rendo = 3.5;
     Real thetaEpi = -M_PI*53.5/180, thetaEndo = 39.5*M_PI/180.0; 

     // BC
     Real kP  = 0.1, kB  = 0.1;
     Real r0P = 0.0, r0B = 0.0;
     Real SearchR_P = 1.0, SearchR_B = 0.5;

     Real TorK = 0.0;
     Vector3d Center = Vector3d::Zero();

     Real Knormal = 1.0;;



     string parameterFileName = argv[1];
     // Reading input from file passed as argument
     ifstream inp;
     inp.open(parameterFileName.c_str(), ios::in);
     if (!inp) {
       cout << "Cannot open input file: " << parameterFileName << endl;
       return(0);
     }
     string ModelName, OutFolder, temp;
  
     inp >> temp >> ModelName;
     inp >> temp >> OutFolder;
     inp >> temp >> Nsteps;
     inp >> temp >> Alphas[0];
     inp >> temp >> Alphas[1];
     inp >> temp >> Alphas[2];
     inp >> temp >> Alphas[3];
     inp >> temp >> Beta;
     inp >> temp >> Gamma;
     inp >> temp >> kP;
     inp >> temp >> kB;
     inp >> temp >> r0P;
     inp >> temp >> r0B;
     inp >> temp >> SearchR_P;
     inp >> temp >> SearchR_B;
     inp >> temp >> TorK;
     inp >> temp >> Knormal;
     
     inp.close();

     // List input parameters
     cout << " Model name    : " << ModelName << endl
	  << " OutputFolder  : " << OutFolder << endl
	  << " Nsteps        : " << Nsteps    << endl
	  << " Alphas[0]     : " << Alphas[0] << endl
	  << " Alphas[1]     : " << Alphas[1] << endl
	  << " Alphas[2]     : " << Alphas[2] << endl
	  << " Alphas[3]     : " << Alphas[3] << endl
	  << " Beta          : " << Beta      << endl
	  << " Gamma         : " << Gamma     << endl
	  << " kP            : " << kP        << endl
	  << " kB            : " << kB        << endl
	  << " r0P           : " << r0P       << endl
	  << " r0B           : " << r0B       << endl
	  << " SearchR_P     : " << SearchR_P << endl
	  << " SearchR_B     : " << SearchR_B << endl
	  << " TorK          : " << TorK      << endl
	  << " Knormal       : " << Knormal   << endl;
     
   





  // Initialize State and Mesh
  State myState, AuxState;
  bool CheckOverlap = false;
  int dofPerNode = 3;
  ostringstream inpA, inpB;
  inpA << ModelName << ".node";   inpB << ModelName << ".ele";
  FEMesh SliceMesh(inpA.str(), inpB.str(), &myState, dofPerNode, CheckOverlap);
  // Membrane meshes
  inpB.str(""); inpB.clear(); inpB << ModelName << ".EpiEleSet";
  FEMesh PericardiumMesh(inpA.str(), inpB.str(), &AuxState, dofPerNode, CheckOverlap);
  CheckOverlap = true;
  inpB.str(""); inpB.clear(); inpB << ModelName << ".BottomEleSet";
  FEMesh BottomMesh(inpA.str(), inpB.str(), &AuxState, dofPerNode, CheckOverlap);
  inpB.str(""); inpB.clear(); inpB << ModelName << ".TopEleSet";
  FEMesh TopMesh(inpA.str(), inpB.str(), &AuxState, dofPerNode, CheckOverlap);

 



  // Initialize bodies
  // Initialize Materials   
  // Cube 
  vector<MechanicsMaterial * > materials;
  vector<GeomElement* > Elements = SliceMesh.getElements();
  int indMat = 0;
  vector<Vector3d > Directions(Ndir, Vector3d::Zero());
  Vector3d er, ec, ez; 
  ez << 0.0, 0.0, 1.0;
  vector<vector<Real > > StretchGrad;
  
  for (int e = 0; e < Elements.size(); e++) {
    const vector<int > & NodesID = Elements[e]->getNodesID();
    for (int q = 0; q < Elements[e]->getNumberOfQuadPoints(); q++) { // One material object per QP
      Vector3d Xq = Vector3d::Zero();
      for (int n = 0; n < NodesID.size(); n++) {
	Xq += Elements[e]->getN(q, n)*myState.getX(NodesID[n]);
      }
      
      // Normalized QP position
      Real r = Xq.norm();
      Real alphaNorm = (Repi-r)/(Repi-Rendo);
      // Compute angle to rotate fiber 
      Real theta = thetaEpi + alphaNorm*(thetaEndo - thetaEpi);
    
      // Tangent fiber in the xy plane
      Vector3d f0; f0 << Xq(1), -Xq(0), 0.0;
      f0 = f0 / f0.norm();
    
      // Rodrigues' rotation formula 
      // Rotation axis
      Vector3d k0; k0 << Xq(0), Xq(1), 0.0;
      k0 = k0 / k0.norm();
    
      // Compute fiber angle at Xq
      Vector3d fiber = Vector3d::Zero();
      fiber = f0*cos(theta) + k0.cross(f0)*sin(theta);
      fiber = fiber / fiber.norm();

      // cout << fiber << endl;

      // Directions[0] = fiber;
      // Directions[1] = ez;
      er << Xq(0), Xq(1), 0.0; er = er/er.norm();
      ec << Xq(1),-Xq(0), 0.0; ec = ec/ec.norm();
      Directions[0] = er;
      Directions[1] = ec;
      Directions[2] = ez;
      Directions[3] = fiber;

      vector<Real > QPStretchGrad(Ndir, 0.0);
      QPStretchGrad[0] = 2.0*(ERRa + alphaNorm*ERRb) / double(Nsteps);
      QPStretchGrad[1] = 2.0*(ECCa + alphaNorm*ECCb) / double(Nsteps);
      QPStretchGrad[2] = 2.0*(ELLa + alphaNorm*ELLb) / double(Nsteps);
      QPStretchGrad[3] = 2.0*EFF / double(Nsteps);
      StretchGrad.push_back(QPStretchGrad);

      materials.push_back(new ImposedKinematics(indMat, Alphas, Stretches, Directions, Beta, Gamma));
      indMat++;
    }
  }
    
  // Initialize Mech Body One
  MechanicsBody Slice(&SliceMesh, &myState, materials);

       int PbDoF = myState.getDOFcount();
       EigenResult TestResult(PbDoF, 0);
       TestResult.initializeResults(PbDoF*4);
 
       TestResult.setRequest(ENERGY);
       TestResult.resetResults(ENERGY);
       Slice.compute(&TestResult);
       cout << endl << "Slice energy is  = " << TestResult.getEnergy() << endl << endl;
  
  ifstream EpiNodeFile, BottomNodeFile, TopNodeFile;
  int NumNodes = 0;
  inpB.str(""); inpB.clear(); inpB << ModelName << ".EpiNodeSet";
  EpiNodeFile.open(inpB.str().c_str());
  EpiNodeFile >> NumNodes;
  cout << "Epi Num Nodes : " << NumNodes << endl;
  vector<int > EpiNodes(NumNodes, 0);
  for(int i = 0; i < NumNodes; i++) {
    EpiNodeFile >> EpiNodes[i];
    // cout << EpiNodes[i] << endl;
  }

  inpB.str(""); inpB.clear(); inpB << ModelName << ".BottomNodeSet";
  BottomNodeFile.open(inpB.str().c_str());
  BottomNodeFile >> NumNodes;
  cout << "Bottom Num Nodes : " << NumNodes << endl;
  vector<int > BottomNodes(NumNodes, 0);
  for(int i = 0; i < NumNodes; i++) {
    BottomNodeFile >> BottomNodes[i];
    // cout << BottomNodes[i] << endl;
  }

  inpB.str(""); inpB.clear(); inpB << ModelName << ".TopNodeSet";
  TopNodeFile.open(inpB.str().c_str());
  TopNodeFile >> NumNodes;
  cout << "Top Num Nodes : " << NumNodes << endl;
  vector<int > TopNodes(NumNodes, 0);
  for(int i = 0; i < NumNodes; i++) {
    TopNodeFile >> TopNodes[i];
    // cout << TopNodes[i] << endl;
  }

  Harmonic PairMatP( kP, r0P );	
  PotentialBody Pericardium(&PericardiumMesh, &myState, &PairMatP, EpiNodes, SearchR_P);
  TestResult.resetResults(ENERGY);
  Pericardium.compute(&TestResult);
  cout << endl << "Pericardium energy is  = " << TestResult.getEnergy() << endl;
  Harmonic PairMatB( kB, r0B );
  PotentialBody Bottom(&BottomMesh, &myState, &PairMatB, BottomNodes, SearchR_B);
  TestResult.resetResults(ENERGY);
  Bottom.compute(&TestResult);
  cout << endl << "Bottom energy is  = " << TestResult.getEnergy() << endl;
  TorsionalBody TorBody(&SliceMesh, &myState, TopNodes, TorK, Center);
  TestResult.resetResults(ENERGY);
  TorBody.compute(&TestResult);
  cout << endl << "Torsional energy is  = " << TestResult.getEnergy() << endl;
  
  PreferredNormalBody PnormBody(&TopMesh, &myState, Knormal);



  
  // Initialize Model
  vector<Body *> Bodies;
  Bodies.push_back(&Slice);
  Bodies.push_back(&Pericardium);
  Bodies.push_back(&Bottom);
  // Bodies.push_back(&TorBody);
  Bodies.push_back(&PnormBody);
  Model myModel(Bodies, &myState);
  




  // Initialize Solver
  vector<int > DoFid;
  vector<Real > DoFvalues;

  EigenNRsolver mySolver(&myState, &myModel, 
			 DoFid, DoFvalues,
			 CHOL, 1.0e-7, 50);

  
  Slice.writeOutputVTK(OutFolder, 0);
  inpB.str(""); inpB.clear(); inpB << OutFolder << "_QP";
  Slice.writeQPdataVTK(inpB.str(), 0);
  inpB.str(""); inpB.clear(); inpB << OutFolder << "_Pericardium";
  Pericardium.writeOutputVTK(inpB.str(), 0);
  inpB.str(""); inpB.clear(); inpB << OutFolder << "_Bottom";
  Bottom.writeOutputVTK(inpB.str(), 0);

  // Incrementally increase contraction
  for (int i = 1; i <= Nsteps; i++) {
    cout << endl << "Step = " << i << endl;    
    
    for (int m = 0; m < materials.size(); m++) {
      for (int j = 0; j < Ndir; j++) {
	Stretches[j] = 1.0 + StretchGrad[m][j] * double(i);
	// cout << Stretches[j] << endl;
      }
      materials[m]->setInternalParameters(Stretches);
    }

    // Update BC bodies
    // Pericardium.setSearchR(SearchR_P*(1.0 + double(i)/Nsteps));
    // Pericardium.setInteractions();
    // TorBody.updateXprev();
    // Bottom.setSearchR(SearchR_B*(1.0 + double(i)/Nsteps));
    // Bottom.setInteractions();

    mySolver.solve(DISP);
    Slice.writeOutputVTK(OutFolder, i);
    inpB.str(""); inpB.clear(); inpB << OutFolder << "_QP";
    Slice.writeQPdataVTK(inpB.str(), i);
    inpB.str(""); inpB.clear(); inpB << OutFolder << "_Pericardium";
    Pericardium.writeOutputVTK(inpB.str(), i);
    inpB.str(""); inpB.clear(); inpB << OutFolder << "_Bottom";
    Bottom.writeOutputVTK(inpB.str(), i);
  }

  
  EigenResult myR(myState.getDOFcount(), 0);
  Real ResX = 0.0, ResY = 0.0, ResZ = 0.0;
  myR.setRequest(FORCE);
  myModel.compute(&myR);
  for (int n = 0; n < myState.getXsize(); n++) {
    ResX +=  myR.getResidual(n*3);
    ResY +=  myR.getResidual(n*3+1);
    ResZ +=  myR.getResidual(n*3+2);
  }
  cout << endl << "ResX = " << ResX << " - ResY = " << ResY << " - ResZ = " << ResZ << endl;

  

  // Timing
  time (&end);
  cout << endl << "BVP solved in " << difftime(end,start) << " s" << endl;



  return 0;
}
