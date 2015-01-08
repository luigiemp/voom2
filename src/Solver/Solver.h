//-*-C++-*-
/*!
  \file Solver.h

  \brief A base solver class.
*/
#ifndef __Solver_h__
#define __Solver_h__

#include "voom.h"
#include "Model.h"
#include "EpetraEllipticResult.h"

// Epetra Headers
#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif


namespace voom
{
  class Solver
  {
  

  public:
    //! Constructor
    Solver(Epetra_MpiComm* comm, Model* myModel);

    //! Set linear solver maximum iterations
    void setLinearSolverMaxIterations(const int iter) {
      _linearSolverMaxIterations = iter;
    }

    //! Set linear solver tolerance
    void setLinearSolverTolerance(const Real toler) {
      _linearSolverTolerance = toler;
    }

    //! Compute function
    virtual void compute(const int loadStepId = 0) = 0;

protected:
    //! Mpicomm
    Epetra_MpiComm* _comm;

    //! Pointer to Model being used
    Model*          _myModel;

    //! Linear System Solver tolerance
    Real            _linearSolverTolerance;

    //! Linear Solver max Iterations
    int                   _linearSolverMaxIterations;

    //! Model Results Struct
    Model::ModelResults   _results;

    //! Solve the system
    virtual void _Solve() = 0;

    //! Update solution
    virtual void _Update() = 0;

  };
}

#endif
