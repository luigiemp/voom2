#include "FEMesh.h"
#include "PoissonModel.h"
#include "IsotropicDiffusion.h"
#include "MechanicsModel.h"
#include "CompNeoHookean.h"
#include "EpetraEllipticResult.h"
#include "EpetraNRsolver.h"



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
     
     string ElementType = "Hexa";
     uint QuadOrder = 2;



     // Poisson model example
     {
       // Finish to initialize Mesh - DoF are problem dependent
       vector<int > LocalDoF(12, 0); // DoF - not nodal - mapping
       vector<int > GhostDoF;
       for (uint i = 0; i < LocalDoF.size(); i++)
	 LocalDoF[i] = i;
     
       FEMesh myFEmesh(Positions, Connectivity, LocalDoF, GhostDoF, ElementType, QuadOrder);
       
       // Define boundary conditions
       vector<uint > DoFids(8,0);
       DoFids[0] = 0;  DoFids[1] = 1;  DoFids[2] = 2;  DoFids[3] = 3;
       DoFids[4] = 8;  DoFids[5] = 9;  DoFids[6] = 10; DoFids[7] = 11;
       vector<Real > DoFvalues(8,0.0);
       DoFvalues[0] = 0.0;  DoFvalues[1] = 0.0;  DoFvalues[2] = 0.0;  DoFvalues[3] = 0.0;
       DoFvalues[4] = 10.0; DoFvalues[5] = 10.0; DoFvalues[6] = 10.0; DoFvalues[7] = 10.0;
       
       // Initialize Poisson Model
       vector<string > ElMatType(2, "IsoTropic");
       
       map<string, DiffusionMaterial* > ElMaterials;
       IsotropicDiffusion IsoDiffMat(2.0);
       ElMaterials.insert(make_pair("IsoTropic", &IsoDiffMat));
       
       PoissonModel myPoissonModel(&myFEmesh, 
				   ElMatType, 
				   ElMaterials);
       
       // Initialize Epetra NR solver
       EpetraNRsolver myNRsolver(mpicomm,
				 &myFEmesh,
				 &myPoissonModel,
				 DoFids, 
				 DoFvalues,
				 1.0e-8, 10);
       
       // Solve
       myNRsolver.solve(); 
        
       // Write solution
       myPoissonModel.writeOutput("TestNRsolver.out","ASCII");
       myPoissonModel.PrintField();

     } // End of Poisson model example



     // Mechanics model example
     {
       // Finish to initialize Mesh - DoF are problem dependent
       vector<int > LocalDoF(36, 0); // DoF - not nodal - mapping
       vector<int > GhostDoF;
       for (uint i = 0; i < LocalDoF.size(); i++)
	 LocalDoF[i] = i;
     
       FEMesh myFEmesh(Positions, Connectivity, LocalDoF, GhostDoF, ElementType, QuadOrder);
       
       // Define boundary conditions
       vector<uint > DoFids;
       uint EasyAssign[] = {0, 1, 2, 4, 5, 8, 9, 11, 26, 29, 32, 35};
       DoFids.assign(EasyAssign, EasyAssign+12);
       vector<Real > DoFvalues(12,0.0);
       DoFvalues[0] = -1.0; DoFvalues[1] = -1.0; DoFvalues[2] = -1.0; DoFvalues[3] = -1.0;
       DoFvalues[4] = -1.0; DoFvalues[5] = -1.0; DoFvalues[6] = -1.0; DoFvalues[7] = -1.0;
       DoFvalues[8] =  2.1; DoFvalues[9] =  2.1; DoFvalues[10] = 2.1; DoFvalues[11] = 2.1;

       // Initialize Mechanics Model
       vector<string > ElMatType(2, "CompNeoHookean");
        
       map<string, MechanicsMaterial* > ElMaterials; 
       CompNeoHookean NeoHookMat(1.0, 1.0);
       ElMaterials.insert(make_pair("CompNeoHookean", &NeoHookMat));
       
       MechanicsModel myMechanicsModel(&myFEmesh, 
				       ElMatType, 
				       ElMaterials);
        
       // Initialize Epetra NR solver
       EpetraNRsolver myNRsolver(mpicomm,
				 &myFEmesh,
				 &myMechanicsModel,
				 DoFids, 
				 DoFvalues,
				 1.0e-8, 10);
        
       // Solve 
       myNRsolver.solve(); 
       myMechanicsModel.PrintField();
        
       // Write solution
       myMechanicsModel.writeOutput("TestNRsolver.out","ASCII");
       
     } // End of Mechanics model example

}
