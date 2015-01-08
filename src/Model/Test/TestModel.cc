// Test for Model
#include "FEMesh.h"
#include "PoissonModel.h"
#include "IsotropicDiffusion.h"
#include "MechanicsModel.h"
#include "CompNeoHookean.h"


using namespace voom;

int main(int argc, char** argv) {

  MPI_Init(&argc, &argv);

  // Initialize Epetra
  Epetra_MpiComm* mpicomm = new Epetra_MpiComm(MPI_COMM_WORLD);

  
  // Create manually an FEmesh
  vector<VectorXd > Positions(12, Vector3d::Zero());
  
    Positions[0] <<-1.0,-1.0,-1.0; 
    Positions[1] << 1.0,-1.0,-1.0; 
    Positions[2] << 1.0, 1.0,-1.0;   
    Positions[3] <<-1.0, 1.0,-1.0;
    Positions[4] <<-1.0,-1.0, 1.0;   
    Positions[5] << 1.0,-1.0, 1.0;
    Positions[6] << 1.0, 1.0, 1.0;
    Positions[7] <<-1.0, 1.0, 1.0;
    Positions[8] <<-1.0,-1.0, 2.0;   
    Positions[9] << 1.0,-1.0, 2.0;   
    Positions[10]<< 1.0, 1.0, 2.0;
    Positions[11]<<-1.0, 1.0, 2.0;
    for (uint i = 0; i < Positions.size(); i++)
      Positions[i] += (Vector3d::Random()*0.1);
    
    vector<vector<int > > Connectivity(2, vector<int >(8,0));
    Connectivity[0][0] = 0;  Connectivity[0][1] = 1;  Connectivity[0][2] = 2;  Connectivity[0][3] = 3;
    Connectivity[0][4] = 4;  Connectivity[0][5] = 5;  Connectivity[0][6] = 6;  Connectivity[0][7] = 7;
    
    Connectivity[1][0] = 4;  Connectivity[1][1] = 5;  Connectivity[1][2] = 6;  Connectivity[1][3] = 7;
    Connectivity[1][4] = 8;  Connectivity[1][5] = 9;  Connectivity[1][6] = 10; Connectivity[1][7] = 11;
    
    string ElementType = "Hexa";
    uint QuadOrder = 2;
    
    
    
  //Test Poisson model
  {
    cout << " ----------------------------- " << endl;
    cout << " TEST OF POISSON MODEL " << endl << endl;

    // Finish to initialize Mesh
    vector<int > LocalDoF(12, 0); // DoF - not nodal - mapping
    vector<int > GhostDoF;
    for (uint i = 0; i < LocalDoF.size(); i++)
      LocalDoF[i] = i;

    FEMesh myFEmesh(Positions, Connectivity, LocalDoF, GhostDoF, ElementType, QuadOrder);

    // Initialize Model
    vector<string > ElMatType(2, "IsoTropic");

    map<string, DiffusionMaterial* > ElMaterials;
    IsotropicDiffusion IsoDiffMat(3.24);
    ElMaterials.insert(make_pair("IsoTropic", &IsoDiffMat));

    PoissonModel myPoissonModel(&myFEmesh, 
				ElMatType, 
				ElMaterials);

    
    
    // Run consistency test
    vector<int> myLocalDoF = myFEmesh.getLocalDoF();
    vector<int> myGhostDoF = myFEmesh.getGhostDoF();
    EpetraEllipticResult myResults(mpicomm, myLocalDoF, myGhostDoF);

    Real perturbationFactor = 0.1;
    uint myRequest = 4;
    Real myH = 1e-6;
    Real myTol = 1e-8;

    myPoissonModel.checkConsistency(myResults, perturbationFactor, myRequest, myH, myTol);

    cout << endl << " END OF TEST OF POISSON MODEL " << endl;
    cout << " ---------------------------- " << endl << endl;
  }



  //Test Mechanics model
  {
    cout << " ------------------------------- " << endl;
    cout << " TEST OF MECHANICS MODEL " << endl << endl;
    
    // Finish to initialize Mesh
    vector<int > LocalDoF(36, 0); // DoF - not nodal - mapping
    vector<int > GhostDoF;
    for (uint i = 0; i < LocalDoF.size(); i++)
      LocalDoF[i] = i;
    
    FEMesh myFEmesh(Positions, Connectivity, LocalDoF, GhostDoF, ElementType, QuadOrder);
    
    // Initialize Model
    vector<string > ElMatType(2, "CompNeoHookean");

    map<string, MechanicsMaterial* > ElMaterials;
    CompNeoHookean NeoHookMat(1.23, 4.02);
    ElMaterials.insert(make_pair("CompNeoHookean", &NeoHookMat));

    MechanicsModel myMechanicsModel(&myFEmesh, 
				    ElMatType, 
				    ElMaterials);

    
    
    // Run consistency test
    vector<int> myLocalDoF = myFEmesh.getLocalDoF();
    vector<int> myGhostDoF = myFEmesh.getGhostDoF();
    EpetraEllipticResult myResults(mpicomm, myLocalDoF, myGhostDoF);

    Real perturbationFactor = 0.1;
    uint myRequest = 6; // Check both Forces and Stiffness
    Real myH = 1e-6;
    Real myTol = 1e-7;

    myMechanicsModel.checkConsistency(myResults, perturbationFactor, myRequest, myH, myTol);
    
    cout << endl << " END OF TEST OF MECHANICS MODEL " << endl;
    cout << " ------------------------------ " << endl << endl;
  }
  
  // Finilize MPI
  MPI_Finalize();
}
