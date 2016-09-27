// Test for Model

#include "PoissonModel.h"
#include "IsotropicDiffusion.h"
#include "MechanicsModel.h"
#include "FEMesh.h"
#include "PassMyoA.h"
#include "CompNeoHookean.h"
#include "EigenResult.h"

#include "LoopShellMesh.h"
#include "SCElastic.h"
#include "LoopShellModel.h"
#include "Jacobian.h"
#include "EigenEllipticResult.h"


using namespace voom;

int main(int argc, char** argv) {

  
  
  //Test Mechanics model
  {
    cout << " ------------------------------- " << endl;
    cout << " TEST OF MECHANICS MODEL " << endl << endl;
  
    cout << " TEST OF 1st MECHANICS MODEL " << endl << endl;
    
    // FEMesh myFEmesh("../../Solver/Test/QuadTet.node", "../../Solver/Test/QuadTet.ele");
    // FEMesh surfMesh("../../Solver/Test/QuadTet.node", "../../Solver/Test/QuadTet.surf");
    // FEMesh myFEmesh("../../Mesh/Test/CoarseLV.node", "../../Mesh/Test/CoarseLV.ele");
    // FEMesh surfMesh("../../Mesh/Test/CoarseLV.node", "../../Mesh/Test/CoarseLV.surf");
    FEMesh myFEmesh("../../Mesh/Test/CubeQuad.node", "../../Mesh/Test/CubeQuad.ele");
    FEMesh surfMesh("../../Mesh/Test/CubeQuad.node", "../../Mesh/Test/SurfCubeQuad.ele");
    // FEMesh myFEmesh("../../Mesh/Test/Cube.node", "../../Mesh/Test/Cube.ele");
    // FEMesh surfMesh("../../Mesh/Test/Cube.node", "../../Mesh/Test/SurfCube.ele");
    // FEMesh myFEmesh("../../Mesh/Test/NodeFile.dat", "../../Mesh/Test/ElFile.dat");

    // FEMesh myFEmesh("/u/home/l/luigiemp/project-cardio/CardiacMesh/Contraction/SmallHeart/Small_A.node", 
    // 		    "/u/home/l/luigiemp/project-cardio/CardiacMesh/Contraction/SmallHeart/Small_A.ele");
    // FEMesh surfMesh("/u/home/l/luigiemp/project-cardio/CardiacMesh/Contraction/SmallHeart/Small_A.node", 
    // 		    "/u/home/l/luigiemp/project-cardio/CardiacMesh/Contraction/SmallHeart/Small_A.surfEle");

    
    // Initialize Model
    uint NodeDoF = 3;
    uint NumQP = 4;

    uint NumMat = myFEmesh.getNumberOfElements()*NumQP;
    vector<MechanicsMaterial * > materials;
    materials.reserve(NumMat);
    
    for (int k = 0; k < NumMat; k++) {
      // vector<Vector3d > Fibers;
      // Vector3d N; N << 1.0, 0.0, 0.0;
      // Fibers.push_back(N);
      // PassMyoA* Mat = new PassMyoA(k, 1.0+double(rand())/RAND_MAX, 3.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX,  1.0+double(rand())/RAND_MAX, 2.0+double(rand())/RAND_MAX, Fibers);

      CompNeoHookean* Mat = new CompNeoHookean(k, 1.0+double(rand())/RAND_MAX, 3.0+double(rand())/RAND_MAX);
      materials.push_back(Mat);
    }
   

    // Apply pressure
    int PressureFlag = 1;
    MechanicsModel myModel(&myFEmesh, materials, NodeDoF, PressureFlag, &surfMesh);

    // Apply spring BC
    Real SpringK = 1.0;
    myModel.initSpringBC("../../Mesh/Test/SurfCubeQuad.nodes", &surfMesh, SpringK);



<<<<<<< HEAD
    int PressureFlag = 1;
    Real Pressure = 1.0;
    MechanicsModel myModel(&myFEmesh, materials, NodeDoF, PressureFlag, Pressure, &surfMesh);
  
=======
    
>>>>>>> 29309bbfa22d9ab99c53776e10f3040558987a65
    // Run consistency test
    uint PbDoF = (myFEmesh.getNumberOfNodes())*myModel.getDoFperNode();
    int TotNumMatProp = NumMat*2;
    EigenResult myResults(PbDoF, TotNumMatProp);

    Real perturbationFactor = 0.1;
    uint myRequest = 6; // Check both Forces and Stiffness
    Real myH = 1e-6;
    Real myTol = 1e-7;

    myModel.checkConsistency(myResults, perturbationFactor, myRequest, myH, myTol);
  
    myModel.checkDmat(myResults, perturbationFactor, myH, myTol);

    // Check model volume
    cout << endl << "Model reference volume is = " << myModel.computeRefVolume() << endl;
    cout << endl << "Model current volume is   = " << myModel.computeCurrentVolume() << endl;

    // Change field and recompute volume
    vector<Real > x(PbDoF, 0.0);
    myModel.getField(x);
    for (int i=0; i<PbDoF; i++) {
      myModel.setField(i, x[i]*0.5);
    }
    cout << endl << "Model current volume is   = " << myModel.computeCurrentVolume() << endl;
    
    
    cout << endl << " END OF TEST OF 1st MECHANICS MODEL " << endl;
    cout << " ------------------------------ " << endl << endl;
  }






  ////////////////////////////////////////////////







  //Test Mechanics model
  {
    cout << " ------------------------------- " << endl;
    cout << " TEST OF 2nd MECHANICS MODEL " << endl << endl;

    FEMesh myFEmesh("/u/home/l/luigiemp/project-cardio/CardiacModels/Contraction/SmallHeart/Small_A_Cavity.node", 
    		    "/u/home/l/luigiemp/project-cardio/CardiacModels/Contraction/SmallHeart/Small_A_Cavity.ele");
    
    // Initialize Model
    uint NodeDoF = 3;
    uint NumQP = 4;

    uint NumMat = myFEmesh.getNumberOfElements()*NumQP;
    vector<MechanicsMaterial * > materials;
    materials.reserve(NumMat);
    
    for (int k = 0; k < NumMat; k++) {
      Jacobian* Mat = new Jacobian(k);
      materials.push_back(Mat);
    }
   
    MechanicsModel myModel(&myFEmesh, materials, NodeDoF);

 
    // Run consistency test
    uint PbDoF = (myFEmesh.getNumberOfNodes())*myModel.getDoFperNode();
    int TotNumMatProp = NumMat;
    EigenEllipticResult myResults(PbDoF, TotNumMatProp);
 
    Real perturbationFactor = 0.1;
    uint myRequest = 2; // Check Forces
    Real myH = 1e-6;
    Real myTol = 1e-7;

    myModel.checkConsistency(myResults, perturbationFactor, myRequest, myH, myTol);

    // Check model volume
    cout << endl << "Model reference volume is       = " << myModel.computeRefVolume() << endl;
    cout << endl << "Model current volume is         = " << myModel.computeCurrentVolume() << endl;
    myResults.setRequest(1);
    myModel.compute(myResults);
    cout << endl << "Model energy = model volume is  = " << myResults.getEnergy() << endl;

    // Change field and recompute volume
    vector<Real > x(PbDoF, 0.0);
    myModel.getField(x);
    for (int i=0; i<PbDoF; i++) {
      myModel.setField(i, x[i]*0.5);
    }
    cout << endl << "Model current volume is   = " << myModel.computeCurrentVolume() << endl;
    myModel.compute(myResults);
    cout << endl << "Model energy = model volume is  = " << myResults.getEnergy() << endl;
    
    
<<<<<<< HEAD
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
    EigenResult myResults(PbDoF, TotNumMatProp);
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
=======

    cout << endl << " END OF TEST OF 2nd MECHANICS MODEL " << endl;
>>>>>>> 29309bbfa22d9ab99c53776e10f3040558987a65
    cout << " ------------------------------ " << endl << endl;


  }




  ///////////////////////////////////////////////////////////







  //Test Poisson model
  {
    cout << " ----------------------------- " << endl;
    cout << " TEST OF POISSON MODEL " << endl << endl;

    // Initialize Model
    FEMesh myFEmesh("../../Mesh/Test/Cube.node", "../../Mesh/Test/Cube.ele");
    
    uint NodeDoF = 1;
    uint NumQP = 1;

    uint NumMat = myFEmesh.getNumberOfElements()*NumQP;
    vector<DiffusionMaterial * > materials;
    materials.reserve(NumMat);
    
    for (int k = 0; k < NumMat; k++) {
      IsotropicDiffusion* Mat = new IsotropicDiffusion( double(rand())/RAND_MAX );
      materials.push_back(Mat);
    }

    PoissonModel myModel(&myFEmesh, materials, NodeDoF);



    // Run consistency test
    uint PbDoF = (myFEmesh.getNumberOfNodes())*myModel.getDoFperNode();
    int TotNumMatProp = NumMat*2;
    EigenEllipticResult myResults(PbDoF, TotNumMatProp);

    Real perturbationFactor = 0.1;
    uint myRequest = 4; // Check stiffness only - There is no energy in the current Poisson model
    Real myH = 1e-6;
    Real myTol = 1e-7;

    myModel.checkConsistency(myResults, perturbationFactor, myRequest, myH, myTol);

    // Print results to Paraview
    myModel.writeOutputVTK("PoissonTest", 0);

    cout << endl << " END OF TEST OF POISSON MODEL " << endl;
    cout << " ---------------------------- " << endl << endl;
  }
  


}
