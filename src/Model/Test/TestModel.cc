// Test for Model

// #include "PoissonModel.h"
// #include "IsotropicDiffusion.h"
#include "MechanicsModel.h"
#include "FEMesh.h"
#include "PassMyoA.h"
#include "CompNeoHookean.h"

#include "EigenEllipticResult.h"

#include "LoopShellMesh.h"
#include "SCElastic.h"
#include "LoopShellModel.h"


using namespace voom;

int main(int argc, char** argv) {

  
  
  // //Test Mechanics model
  // {
  //   cout << " ------------------------------- " << endl;
  //   cout << " TEST OF MECHANICS MODEL " << endl << endl;
    
  //   // FEMesh myFEmesh("../../Solver/Test/QuadTet.node", "../../Solver/Test/QuadTet.ele");
  //   // FEMesh surfMesh("../../Solver/Test/QuadTet.node", "../../Solver/Test/QuadTet.surf");
  //   // FEMesh myFEmesh("../../Mesh/Test/CoarseLV.node", "../../Mesh/Test/CoarseLV.ele");
  //   // FEMesh surfMesh("../../Mesh/Test/CoarseLV.node", "../../Mesh/Test/CoarseLV.surf");
  //   FEMesh myFEmesh("../../Mesh/Test/CubeQuad.node", "../../Mesh/Test/CubeQuad.ele");
  //   FEMesh surfMesh("../../Mesh/Test/CubeQuad.node", "../../Mesh/Test/SurfCubeQuad.ele");
  //   // FEMesh myFEmesh("../../Mesh/Test/Cube.node", "../../Mesh/Test/Cube.ele");
  //   // FEMesh surfMesh("../../Mesh/Test/Cube.node", "../../Mesh/Test/SurfCube.ele");
  //   // FEMesh myFEmesh("../../Mesh/Test/NodeFile.dat", "../../Mesh/Test/ElFile.dat");
    
  //   // Initialize Model
  //   uint NodeDoF = 3;

  //   uint NumMat = myFEmesh.getNumberOfElements();
  //   vector<MechanicsMaterial * > materials;
  //   materials.reserve(NumMat);
  //   for (int k = 0; k < NumMat; k++) {
  //     // // PassMyoA* Mat = new PassMyoA(1.0+double(rand())/RAND_MAX, 3.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX, 2.0+double(rand())/RAND_MAX, 2.0+double(rand())/RAND_MAX);
  //     // PassMyoA* Mat = new PassMyoA(k, 1.0+double(rand())/RAND_MAX, 3.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX,  1.0+double(rand())/RAND_MAX, 2.0+double(rand())/RAND_MAX);
  //     // Vector3d N; N << 1.0, 0.0, 0.0;
  //     // Mat->setN(N);
  //     // materials.push_back(Mat);
  //    materials.push_back(new CompNeoHookean(k, 10.0, 3.0) );
  //   }
  //   // CompNeoHookean *Mat = new CompNeoHookean(0, 10.0, 3.0);
  //   // for (int k = 0; k < NumMat; k++) {
  //   //   materials.push_back(Mat);
  //   // }

  
  //   int PressureFlag = 1;
  //   Real Pressure = 1.0;
  //   MechanicsModel myModel(&myFEmesh, materials, NodeDoF, PressureFlag, Pressure, &surfMesh);
    
  //   // Run consistency test
  //   uint PbDoF = (myFEmesh.getNumberOfNodes())*myModel.getDoFperNode();
  //   int TotNumMatProp = NumMat*2;
  //   EigenEllipticResult myResults(PbDoF, TotNumMatProp);

  //   Real perturbationFactor = 0.1;
  //   uint myRequest = 6; // Check both Forces and Stiffness
  //   Real myH = 1e-6;
  //   Real myTol = 1e-7;

  //   myModel.checkConsistency(myResults, perturbationFactor, myRequest, myH, myTol);
  
  //   myModel.checkDmat(myResults, perturbationFactor, myH, myTol);
    
  //   cout << endl << " END OF TEST OF MECHANICS MODEL " << endl;
  //   cout << " ------------------------------ " << endl << endl;
  // }

  {
    cout << " ------------------------------- " << endl;
    cout << " Testing Loop Shell Model " << endl << endl;
    //LoopShellMesh icosa_mesh("sphere_nodes_1SD.dat","sphere_conn_1SD.dat");
    //LoopShellMesh icosa_mesh("sphere_loop_nodes.dat","sphere_loop_conn.dat");
    LoopShellMesh icosa_mesh("T5sphere_nodes.dat","T5sphere_conn.dat");
    //LoopShellMesh icosa_mesh("nonicosa_sphere_nodes.dat","nonicosa_sphere_conn.dat");
    uint NumMat = icosa_mesh.getNumberOfElements();
    uint NodeDoF = 3;
    vector<SCElastic *> materials;
    materials.reserve(NumMat);
    for(int k = 0; k < NumMat; k++)
      materials.push_back(new SCElastic(1,0,0));

    LoopShellModel model( &icosa_mesh, materials, NodeDoF);

    // Run consistency test
    uint PbDoF = (icosa_mesh.getNumberOfNodes())*model.getDoFperNode();
    int TotNumMatProp = NumMat*2;
    EigenEllipticResult myResults(PbDoF, TotNumMatProp);
    ComputeRequest myRequest = FORCE;
    myResults.setRequest(myRequest);
    //model.printField();
    model.compute( myResults);
    //cout << "Energy : " << myResults.getEnergy()<<endl;
    //model.checkConsistency(myResults, perturbationFactor, myRequest, myH, myTol);
    cout << std::scientific;
    //cout << (*(myResults._residual)).cwiseAbs().colwise().maxCoeff() << endl;
    cout << (*(myResults._residual)) << endl;
    int counter = 0;
    for (int i=0; i< PbDoF; i++ ){
      if (abs((*(myResults._residual))(i)) > .0001){
	counter++;
	//cout << "Detected at " << i << endl;
      }
    }
    //cout << "counter = " << counter << endl;
    cout << endl << " END OF TEST OF MECHANICS MODEL " << endl;
    cout << " ------------------------------ " << endl << endl;


  }

  // //Test Poisson model
  // {
  //   cout << " ----------------------------- " << endl;
  //   cout << " TEST OF POISSON MODEL " << endl << endl;

  //   // Finish to initialize Mesh
  //   vector<int > LocalDoF(12, 0); // DoF - not nodal - mapping
  //   vector<int > GhostDoF;
  //   for (uint i = 0; i < LocalDoF.size(); i++)
  //     LocalDoF[i] = i;

  //   FEMesh myFEmesh(Positions, Connectivity, LocalDoF, GhostDoF, ElementType, QuadOrder);

  //   // Initialize Model
  //   vector<string > ElMatType(2, "IsoTropic");

  //   map<string, DiffusionMaterial* > ElMaterials;
  //   IsotropicDiffusion IsoDiffMat(3.24);
  //   ElMaterials.insert(make_pair("IsoTropic", &IsoDiffMat));

  //   PoissonModel myPoissonModel(&myFEmesh, 
  // 				ElMatType, 
  // 				ElMaterials);

    
    
  //   // Run consistency test
  //   vector<int> myLocalDoF = myFEmesh.getLocalDoF();
  //   vector<int> myGhostDoF = myFEmesh.getGhostDoF();
  //   EpetraEllipticResult myResults(mpicomm, myLocalDoF, myGhostDoF);

  //   Real perturbationFactor = 0.1;
  //   uint myRequest = 4;
  //   Real myH = 1e-6;
  //   Real myTol = 1e-8;

  //   myPoissonModel.checkConsistency(myResults, perturbationFactor, myRequest, myH, myTol);

  //   cout << endl << " END OF TEST OF POISSON MODEL " << endl;
  //   cout << " ---------------------------- " << endl << endl;
  // }
  


}
