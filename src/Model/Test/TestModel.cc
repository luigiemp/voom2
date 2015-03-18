// Test for Model

// #include "PoissonModel.h"
// #include "IsotropicDiffusion.h"
#include "MechanicsModel.h"
#include "FEMesh.h"
#include "PassMyoA.h"
#include "EigenEllipticResult.h"


using namespace voom;

int main(int argc, char** argv) {

  
  
  //Test Mechanics model
  {
    cout << " ------------------------------- " << endl;
    cout << " TEST OF MECHANICS MODEL " << endl << endl;
    
    FEMesh myFEmesh("../../Mesh/Test/NodeFile.dat", "../../Mesh/Test/ElFile.dat");
    
    // Initialize Model
    uint NodeDoF = 3;

    uint NumMat = 2;
    vector<MechanicsMaterial * > materials;
    materials.reserve(NumMat);
    for (uint k = 0; k < NumMat; k++) {
      materials.push_back(new PassMyoA(1.0+double(rand())/RAND_MAX, 3.0+double(rand())/RAND_MAX, 2.0+double(rand())/RAND_MAX, 2.0+double(rand())/RAND_MAX));
      Vector3d N; N << 1.0, 0.0, 0.0;
      materials[k]->setN(N);
    }

    MechanicsModel myModel(&myFEmesh, materials, NodeDoF);
    
    // Run consistency test
    uint PbDoF = (myFEmesh.getNumberOfNodes())*myModel.getDoFperNode();
    EigenEllipticResult myResults(PbDoF, NumMat);

    Real perturbationFactor = 0.1;
    uint myRequest = 6; // Check both Forces and Stiffness
    Real myH = 1e-6;
    Real myTol = 1e-7;

    myModel.compute(myResults);
    cout << "Energy = " << myResults.getEnergy() << endl;

    myModel.checkConsistency(myResults, perturbationFactor, myRequest, myH, myTol);
    
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
