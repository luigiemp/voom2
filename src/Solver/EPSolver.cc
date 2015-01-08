//-*-C++-*-
#include "EPSolver.h"
// #include "EP_ExplicitSolver.h"
// #include "EP_BackwardEuler.h"
#include "EP_CrankNicholson.h"

namespace voom{
  EPSolver::EPSolver(Epetra_MpiComm *comm, Model* model,
		     EpetraEPResult* result, LumpingScheme lumpingScheme) {
    _comm             = comm;
    _myModel          = model;
    _result           = result;
    _lumpingScheme    = lumpingScheme;
    _maxIter          = 20;
    

    
    // Define Epetra vectors
    _exporter = new Epetra_Export(result->getSourceMap(), result->getTargetMap());
    
    // Local Values
    // _diffusiveCurrent = new Epetra_Vector(*_targetMap);

    // Ghost Values
    // _ghostDiffusiveCurrent = new Epetra_Vector(*_sourceMap);

    // Linear Solver Tolerance
    _linearSolverTolerance = 1E-5;

  }

  //! Write Voltage Output
  void EPSolver::WriteOutput(const Real time, std::string path)
  {

  }

  //! Gateway to Solver Class
  EPSolver* EPSolver::New(Epetra_MpiComm *comm, Model& model, 
			 EpetraEPResult& EPResult,
			 const LumpingScheme LScheme, const Method SolverType)
  {
    EPSolver* mySolver;
    switch(SolverType)
    {
      case EXPLICIT:
	// mySolver = new EP_ExplicitSolver(comm, model, EPResult, LScheme);
	break;
      case BACKWARDEULER:
	// mySolver = new EP_BackwardEuler(comm, model, EPResult, LScheme);
	break;
      case CRANKNICHOLSON:
	mySolver = new EP_CrankNicholson(comm, model, EPResult, LScheme);
	break;
      default:
	std::cout << "Unknown Solver specified\n";
	exit(0);
     }
    return mySolver;
  }

}// End of namespace
