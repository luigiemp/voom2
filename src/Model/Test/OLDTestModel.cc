// Test for Model
//#include "MechanicsModel.h"
#include "PoissonModel.h"
#include "FEMesh.h"

//#include "CompNeoHookean.h"

using namespace voom;

int main(int argc, char** argv) {

  MPI_Init(&argc, &argv);

  // Initialize Epetra
  Epetra_MpiComm* mpicomm = new Epetra_MpiComm(MPI_COMM_WORLD);

  
  // //Test mechanics model
  // {
  //   // Initialize material
  //   Real A[3] = {1.0,0.3,1.0}; // [E, nu. rho]
  //   CompNeoHookean Mat(A);
  //   std::vector<MechanicMaterial*> matData(1, &Mat);

  //   // Initialize model
  //   MechanicsModel<3>* myMechModel = new MechanicsModel<3>( std::string(argv[1]),
  // 							    matData);
  //   // Run consistency test
  //   myMechModel->checkConsistency(mpicomm);

  //   // Delete models
  //   delete myMechModel;
  // }
  
  cout << " ----- *** ----- " << endl;
  
  //Test Poisson model
  {
    Mesh* myMesh = Mesh::New( argv[1] );
    Model* myPoissonModel = Model::New( myMesh, argv[1] );
    // // Initialize material
    // vector<double > A(1, 1.0); // [\alpha]
   
    // // Initialize model
    // PoissonModel<3>* myPoissonModel = new PoissonModel<3>( std::string(argv[1]),
							   // A );
    // Run consistency test
    myPoissonModel->checkConsistency(mpicomm);

    // Delete model
    delete myMesh;
    delete myPoissonModel;
  }
   
  // Finilize MPI
  MPI_Finalize();
}
