//-*-C++-*-
#ifndef _EP_CrankNicholson_h_
#define _EP_CrankNicholson_h_
#include "EPSolver.h"

namespace voom{
  
  class EP_CrankNicholson : public EPSolver {
  private:
    //! Diffusion Problem variables
    Epetra_LinearProblem* _diffusionProblem;
    AztecOO*              _diffusionSolver;

    //! Ionic Problem Variables
    Epetra_LinearProblem* _ionicProblem;
    AztecOO*              _ionicSolver;

  public:
    //! Constructor
    EP_CrankNicholson(Epetra_MpiComm *comm, const Model& model, 
		      EPetraEPResult& Result, const Mesh* mesh);

    //! Destructor
    ~EP_CrankNicholson()
    {
      delete _diffusionProblem; delete _diffusionSolver; 
      delete _ionicProblem; delete _ionicSolver;
    }

    //! Crank Nicholson solve function
    void Solve();
  };
}

#endif
