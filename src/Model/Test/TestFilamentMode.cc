// Test for Model

// #include "PoissonModel.h"
// #include "IsotropicDiffusion.h"
#include "FilamentsModel.h"
#include "FEMesh.h"
#include "EigenResult.h"


using namespace voom;

int main(int argc, char** argv) {

  
  
  //Test Mechanics model
  {
    cout << " ------------------------------- " << endl;
    cout << " TEST OF GEL MODEL " << endl << endl;
    
    
    FEMesh TestFEmesh("../../Mesh/Test/Filament.node", "../../Mesh/Test/Filament.ele");
    
    // Initialize Model
    uint NodeDoF = 3;
    uint NumMat = 1;
    uint NumEle = TestFEmesh.getNumberOfElements();
    Real k = 1.0;
    Vector3d d0;

    d0 << 1.0, 0.0, 0.0;


    std::cout << "Number of elements is " << NumEle << std::endl;
    
    vector<FilamentsMaterial * > materials;
    materials.reserve(NumMat);
    
    Spring* MatSpring = new Spring(0, k, d0);

    materials.push_back(MatSpring);
    
    /*
    for (int k = 0; k < NumMat; k++) {
      PassMyoA* Mat = new PassMyoA(k, 1.0+double(rand())/RAND_MAX, 3.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX,  1.0+double(rand())/RAND_MAX, 2.0+double(rand())/RAND_MAX);

      materials.push_back(Mat);

      Vector3d N; N << 1.0, 0.0, 0.0;
      Fibers.push_back(N);
      // materials.push_back(new CompNeoHookean(k, 10.0, 3.0) );
    }
    
    x
    myFEmesh.setFibers(Fibers);


    int PressureFlag = 1;
    Real Pressure = 1.0;
    MechanicsModel myModel(&myFEmesh, materials, NodeDoF, PressureFlag, Pressure, &surfMesh);
    
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
    */

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
