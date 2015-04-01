//-*-C++-*-
/*!
  \file EpetraSparseSolver.h
  \brief A sparse solver implemented using AztecOO from Trilinos. Solves the
  system \f$ Ax = b \$ where \f$A\f$ is a Epetra_FECrsMatrix, \f$x,b\f$ are
  Epetra_Vectors. The type of Sparse solver used can be changed by the user.
*/

#ifndef __EpetraSparseSolver_h__
#define __EpetraSparseSolver_h__

#include "voom.h"
#include "Model.h"
#include "Mesh.h"

// Epetra Headers
#include "Epetra_ConfigDefs.h"
// #ifdef HAVE_MPI
// #include "mpi.h"
// #include "Epetra_MpiComm.h"
// #else
#include "Epetra_SerialComm.h"
// #endif

#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsMatrix.h"

// Aztec gives us Iterative Solvers
#include "AztecOO_config.h"
#include "AztecOO.h"
#include "AztecOO_Operator.h"

#include "Amesos.h"

namespace voom{

  //! Solver options from AztecOO
  enum AZ_SolverType { CG, GMRES, CGS, TFQMR, BICGSTAB, LU };
  //! Precondition options from AztecOO
  enum AZ_PreConditionerType { NONE, JACOBI, NEUMANN, LS, SYM_GS, DOM_DECOMP, ILU};
  //! Convergence options
  enum AZ_ConvergenceType { R0, RHS, ANORM, NOSCALED, SOL, WEIGHTED};
			       
 //  class EpetraSparseSolver
//   {
//   public:
//     //! Constructor
//     EpetraSparseSolver(Epetra_MpiComm* comm,
// 		       Real linearSolverTolerance, int linearSolverMaxIterations,
// 		       Epetra_Vector* RHS, Epetra_Vector* ghostRHS,
// 		       Epetra_FECrsMatrix* SystemMatrix,
// 		       Epetra_Map * TargetMap,
// 		       Epetra_Map * SourceMap);

//     //! Destructor
//     ~EpetraSparseSolver() {
//       delete _X; 
//       delete _ghostX;
//       delete _exporter; 
//       delete _problem; 
//       delete _solver;
//     }

//     //! Apply Natural and essential BC
//     void applyBC(Mesh * myMesh,
// 		 const vector<uint > & EBC_Dof,
// 		 const vector<Real> & EBC_Values);
    
//     // Set manually value of _X and _ghostX
//     void setSolution(vector<uint > DoFid, vector<Real > DoFvalues);



//     //! Solve function
//     void solve(Real gamma = 1.0) {
//       // LEP: We should only need to zero _X, not _ghostX (not 100% sure when we solve in parallel)
//       this->ZeroX();
//       _B->Scale(gamma);
//       _ghostB->Scale(gamma);
//       _solver->Iterate(_linearSolverMaxIterations, _linearSolverTolerance);
//     }

//     void ZeroX() {
//       for (uint i = 0; i < _X->MyLength(); i++)
// 	_extractedX[i] = 0.0;
//     }
    
//     //! Update solution values through Model-specific linear update function
//     void updateModel(Model * myModel) {
//       myModel->linearizedUpdate( _extractedX, _extractedGhostX );
//     }

//     Real SolutionNorm() {
//       Real Xnorm;
//       _X->Norm2(&Xnorm);
//       return Xnorm;
//     }

//     Real DeltaXnorm(vector<uint > ExcludedDoFid) {
//       Real FreeSolutionNorm = 0.0;
//       set<uint > ExcludedSet;
//       for (uint i = 0; i < ExcludedDoFid.size(); i++)
// 	ExcludedSet.insert(ExcludedDoFid[i]);
	
//       for (uint i = 0; i < _X->MyLength(); i++) {
// 	if (ExcludedSet.find(i) == ExcludedSet.end()) {
// 	  FreeSolutionNorm += pow(_extractedX[i], 2.0);
// 	}
//       }
      
//       return pow(FreeSolutionNorm, 0.5);
//     }

//     Real DeltaRHSnorm(vector<uint > ExcludedDoFid) {
//       Real ForcingTermNorm = 0.0;
//       set<uint > ExcludedSet;
//       for (uint i = 0; i < ExcludedDoFid.size(); i++)
// 	ExcludedSet.insert(ExcludedDoFid[i]);
      
//       Real *extractedB;
//       _B->ExtractView(&extractedB);
//       for (uint i = 0; i < _B->MyLength(); i++) {
// 	if (ExcludedSet.find(i) == ExcludedSet.end()) {
// 	  ForcingTermNorm += pow(extractedB[i], 2.0);
// 	}
//       }
      
//       return pow(ForcingTermNorm, 0.5);
//     }

//     void PrintSolution() {
//       cout << *_X << endl;
//     }

    

//     //! Set Solver Type
//     void setSolverType(AZ_SolverType type) {
//       switch(type) {
//       case 0: _solver->SetAztecOption(AZ_solver, AZ_cg);          break;
//       case 1: _solver->SetAztecOption(AZ_solver, AZ_cgs);	  break;
//       case 2: _solver->SetAztecOption(AZ_solver, AZ_gmres);	  break;
//       case 3: _solver->SetAztecOption(AZ_solver, AZ_tfqmr);	  break;
//       case 4: _solver->SetAztecOption(AZ_solver, AZ_bicgstab);	  break;
//       case 5: _solver->SetAztecOption(AZ_solver, AZ_ilu);  	  break;
//       default:							 
// 	{
// 	  std::cout << "Unknown Solver Type: " << type << "\n"
// 		    << "Exiting...\n";
// 	}
//       }
//     }	

//     //! Set Preconditioner type
//     void setPreconditionerType(AZ_PreConditionerType type) {
//       switch(type) {
//       case 0: _solver->SetAztecOption(AZ_precond, AZ_none);        break;
//       case 1: _solver->SetAztecOption(AZ_precond, AZ_Jacobi);      break;
//       case 2: _solver->SetAztecOption(AZ_precond, AZ_Neumann);     break;
//       case 3: _solver->SetAztecOption(AZ_precond, AZ_ls);	   break;
//       case 4: _solver->SetAztecOption(AZ_precond, AZ_sym_GS);      break;
//       case 5: _solver->SetAztecOption(AZ_precond, AZ_dom_decomp);  break;
//       case 6: _solver->SetAztecOption(AZ_precond, AZ_lu);          break;
//       default:
// 	{
// 	  cout << "Unknown Preconditioner Type: " << type << "\n"
// 	       << "Exiting...\n";
// 	}
//       }
//     }

//     //! Set convergence criterion
//     void setConvergenceCriterion(AZ_ConvergenceType type) {
//       switch(type) {
//       case 0: _solver->SetAztecOption(AZ_conv, AZ_r0);         break;
//       case 1: _solver->SetAztecOption(AZ_conv, AZ_rhs);	       break;
//       case 2: _solver->SetAztecOption(AZ_conv, AZ_Anorm);      break;
//       case 3: _solver->SetAztecOption(AZ_conv, AZ_noscaled);   break;
//       case 4: _solver->SetAztecOption(AZ_conv, AZ_sol);	       break;
//       case 5: _solver->SetAztecOption(AZ_conv, AZ_weighted);   break;
//       default:						      
// 	{
// 	  cout << "Unknown Convergence Type: " << type << "\n"
// 	       << "Exiting...\n";
// 	}
//       }
//     }



// protected:
//     //! Mpicomm
//     Epetra_MpiComm* _comm;

//     //! Linear System Solver tolerance
//     Real            _linearSolverTolerance;
    
//     //! Linear Solver max Iterations
//     int             _linearSolverMaxIterations;

//     //! Solution vector
//     Epetra_Vector*        _X;
//     Epetra_Vector*        _ghostX;
//     Real*                 _extractedX;
//     Real*                 _extractedGhostX;

//     //! RHS 
//     Epetra_Vector*        _B;
//     Epetra_Vector*        _ghostB;

//     //! System matrix
//     Epetra_FECrsMatrix*   _A;
    
//     //! Epetra maps
//     Epetra_Map *          _targetMap;
//     Epetra_Map *          _sourceMap;

//     //! Epetra Exporter
//     Epetra_Export*        _exporter;

//     //! Linear Problem
//     Epetra_LinearProblem* _problem;
//     AztecOO*              _solver;
//   };
}

#endif // __SparseSolver_h__

