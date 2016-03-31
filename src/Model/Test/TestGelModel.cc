// Test for Model

// #include "PoissonModel.h"
// #include "IsotropicDiffusion.h"
#include "GelModel.h"
#include "GelMesh.h"
#include "EigenResult.h"
#include "Spring.h"

using namespace voom;

int main(int argc, char** argv) {

  
  
  //Test Mechanics model
  {
    cout << " ------------------------------- " << endl;
    cout << " TEST OF GEL MODEL " << endl << endl;
    
    
    GelMesh TestFEmesh("../../Mesh/Test/Filament2.node", "../../Mesh/Test/Filament2.ele");
    
    // Initialize Model
    uint NodeDoF = 3;

    uint NumEle = TestFEmesh.getNumberOfFilaments();
    
    std::cout << "Number of elements is " << NumEle << std::endl;
    
    vector<FilamentMaterial * > materials;
    materials.reserve(NumEle);
    
    Real k = 1.0;
    Vector3d d0;
    d0 << 1.0, 0.0, 0.0;
    
    for (int k = 0; k < NumEle; k++) {        
      Spring* Mat = new Spring(0,k,d0);
      materials.push_back(Mat);
    }
    
    GelModel myModel(&TestFEmesh, materials, NodeDoF,0,1);
    
     
    myModel.initializeField(2.0);    
    const vector<GeomFilament* > elements = TestFEmesh.getFilaments();    
    GeomFilament* geomFil = elements[0];  
    
    
    vector<Vector3d > dlist(1,Vector3d::Zero());
     
    myModel.computeDeformation(dlist,geomFil);

    //cout << dlist[0] << endl;
    


    /* 
    for (int k = 0; k < NumMat; k++) {
      // PassMyoA* Mat = new PassMyoA(1.0+double(rand())/RAND_MAX, 3.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX, 2.0+double(rand())/RAND_MAX, 2.0+double(rand())/RAND_MAX);
      PassMyoA* Mat = new PassMyoA(k, 1.0+double(rand())/RAND_MAX, 3.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX, 1.0+double(rand())/RAND_MAX,  1.0+double(rand())/RAND_MAX, 2.0+double(rand())/RAND_MAX);

      materials.push_back(Mat);

      Vector3d N; N << 1.0, 0.0, 0.0;
      Fibers.push_back(N);
      // materials.push_back(new CompNeoHookean(k, 10.0, 3.0) );
    }
    // CompNeoHookean *Mat = new CompNeoHookean(0, 10.0, 3.0);
    // for (int k = 0; k < NumMat; k++) {
    //   materials.push_back(Mat);
    // }

  
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
