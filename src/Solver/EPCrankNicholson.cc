#include "EP_CrankNicholson.h"

namespace voom{

  EP_CrankNicholson(Epetra_MpiComm *comm, const Model& model, 
		    EPetraEPResult& result)
  {
    EPSolver(Epetra_MpiComm *comm, model, result);
    
    // RESULT STRUCT TYPE TBD
    _result.request(0);

    // Set alpha, beta, gamma
    _result.alpha = 0;
    _result.beta = 0;
    _result.gamma = 0;

    _myModel->Compute(_result);

    

    // Get Epetra Matrices and Vectors

    // Set up Diffusion Problem Solver
    _diffusionProblem = new Epetra_LinearProblem(_LHS, _Voltage, 
					_RHS);
    
    _diffusionSolver = new AztecOO(*_diffusionProblem);
    _diffusionSolver->SetAztecOption(AZ_precond, AZ_Jacobi);
    _diffusionSolver->SetAztecOption(AZ_solver, AZ_cg);
    _diffusionSolver->SetAztecOption(AZ_diagnostics, AZ_none);
    _diffusionSolver->SetAztecOption(AZ_output, 0);

    // Ionic Problem Solver
    _ionicProblem = new Epetra_LinearProblem(_CL, _Voltage, 
					_diffusiveCurrent);
    
    _ionicSolver = new AztecOO(*_ionicProblem);
    _ionicSolver->SetAztecOption(AZ_precond, AZ_Jacobi);
    _ionicSolver->SetAztecOption(AZ_solver, AZ_cg);
    _ionicSolver->SetAztecOption(AZ_diagnostics, AZ_none);
    _ionicSolver->SetAztecOption(AZ_output, 0);

  }
  
  void EP_CrankNicholson::Solve() {

    // Step 1: Compute the diffusive Current from body
    
    // Step 2: Get Voltage from Ionic current only first body only
    
    // Step 3: Same as step 1
    
  } // End of Solve routine

} // namespace
