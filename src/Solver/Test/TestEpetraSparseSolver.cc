#include "FEMesh.h"
#include "PoissonModel.h"
#include "IsotropicDiffusion.h"
#include "EpetraEllipticResult.h"
#include "EpetraSparseSolver.h"


#include "EpetraSparseSolver.h"

using namespace voom;

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  // Initialize Epetra
  Epetra_MpiComm* mpicomm = new Epetra_MpiComm(MPI_COMM_WORLD);
  int myRank = 0;
  MPI_Comm_rank( MPI_COMM_WORLD, &myRank );
  cout<<"Processor "<< myRank << " started"<<endl;

  // Initialize Mesh
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

     vector<vector<int > > Connectivity(2, vector<int >(8,0));
     Connectivity[0][0] = 0;  Connectivity[0][1] = 1;  Connectivity[0][2] = 2;  Connectivity[0][3] = 3;
     Connectivity[0][4] = 4;  Connectivity[0][5] = 5;  Connectivity[0][6] = 6;  Connectivity[0][7] = 7;
     
     Connectivity[1][0] = 4;  Connectivity[1][1] = 5;  Connectivity[1][2] = 6;  Connectivity[1][3] = 7;
     Connectivity[1][4] = 8;  Connectivity[1][5] = 9;  Connectivity[1][6] = 10; Connectivity[1][7] = 11;
     
     vector<int > LocalDoF(12, 0); // DoF - not nodal - mapping
     vector<int > GhostDoF;
     for (uint i = 0; i < LocalDoF.size(); i++)
       LocalDoF[i] = i;
     
     string ElementType = "Hexa";
     
     uint QuadOrder = 2;
     
     FEMesh myFEmesh(Positions, Connectivity, LocalDoF, GhostDoF, ElementType, QuadOrder);


    
  // Initialize Poisson Model
     vector<string > ElMatType(2, "IsoTropic");

     map<string, DiffusionMaterial* > ElMaterials;
     IsotropicDiffusion IsoDiffMat(2.0);
     ElMaterials.insert(make_pair("IsoTropic", &IsoDiffMat));
     
     PoissonModel myPoissonModel(&myFEmesh, 
				 ElMatType, 
				 ElMaterials);
     


  // Initialize EpetraElliptic Results
     EpetraEllipticResult myResults(mpicomm,
				    myFEmesh.getLocalDoF(),
				    myFEmesh.getGhostDoF()); 

     

  // Initialize EpetraSparseSolver
     EpetraSparseSolver mySolver(mpicomm,
				 1.0e-12, 1000000,
				 myResults.getResidual(), myResults.getGhostResidual(),
				 myResults.getStiffness(),
				 myResults.getTargetMap(),
				 myResults.getSourceMap());


     
  // Apply boundary conditions
     vector<uint > DoFids(8,0);
     DoFids[0] = 0;  DoFids[1] = 1;  DoFids[2] = 2;  DoFids[3] = 3;
     DoFids[4] = 8;  DoFids[5] = 9;  DoFids[6] = 10; DoFids[7] = 11;
     vector<Real > DoFvalues(8,0.0);
     DoFvalues[0] = 0.0;  DoFvalues[1] = 0.0;  DoFvalues[2] = 0.0;  DoFvalues[3] = 0.0;
     DoFvalues[4] = 10.0; DoFvalues[5] = 10.0; DoFvalues[6] = 10.0; DoFvalues[7] = 10.0;
     
     for (uint i = 0; i < DoFids.size(); i++) {
       myPoissonModel.setField(DoFids[i], DoFvalues[i]);
     }
     


  // Model compute - Fill in RHS and 'Stiffness' matrix
     myResults.setRequest(6);
     myPoissonModel.compute(myResults);
     //cout << *(myResults.getResidual()) << endl;
     
     
  // Solve
     mySolver.applyBC(&myFEmesh, DoFids, DoFvalues);
     mySolver.solve(-1.0);
     vector<Real > Zeros(DoFids.size(), 0.0);
     mySolver.setSolution(DoFids, Zeros);
     mySolver.updateModel(&myPoissonModel);
     mySolver.PrintSolution();


  // Write solution
     myPoissonModel.writeOutput("TestSparseSolver.out","ASCII");

}
