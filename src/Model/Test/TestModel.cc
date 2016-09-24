// Test for Model

#include "PoissonModel.h"
#include "IsotropicDiffusion.h"
#include "MechanicsModel.h"
#include "FEMesh.h"
#include "PassMyoA.h"
#include "CompNeoHookean.h"
#include "Jacobian.h"
#include "EigenEllipticResult.h"


using namespace voom;

int main(int argc, char** argv) {

  
  
  //Test Mechanics model
  {
    cout << " ------------------------------- " << endl;
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



    
    // Run consistency test
    uint PbDoF = (myFEmesh.getNumberOfNodes())*myModel.getDoFperNode();
    int TotNumMatProp = NumMat*2;
    EigenEllipticResult myResults(PbDoF, TotNumMatProp);

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
    
    

    cout << endl << " END OF TEST OF 2nd MECHANICS MODEL " << endl;
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
