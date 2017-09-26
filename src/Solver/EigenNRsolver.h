//-*-C++-*-
/*!
  \file EigenNRsolver.h
  \brief Newton Raphson solver based on Eigen sparse solver
*/

#ifndef __EigenNRsolver_h__
#define __EigenNRsolver_h__

#include "voom.h"
#include "State.h"
#include "EigenResult.h"
#include "Mesh.h"
#include "Model.h"
#include "MechanicsMaterial.h"
#include "MechanicsBody.h"


namespace voom{
			    
  //! Enumerator for requesting computed results 
  enum SolverType {
    CHOL = 0,
    CG   = 1,
    LU   = 2
  };

  //! Enumerator for requesting computed results 
  enum SolveFor {
    DISP = 0,
    MAT  = 1
  };
   
  class EigenNRsolver
  {
  public:
    //! Constructor
    EigenNRsolver(State *myState,
		  Model *myModel, 
		  vector<int > & DoFid,
		  vector<Real > & DoFvalues,
		  SolverType LinSolType = CHOL,
		  Real NRtol = 1.0e-8, int NRmaxIter = 100):
      _myState(myState), _myModel(myModel),
      _DoFid(DoFid), _DoFvalues(DoFvalues),
      _linSolType(LinSolType),
      _NRtol(NRtol), _NRmaxIter(NRmaxIter) {};

    //! Destructor
    ~EigenNRsolver() {};
      
    //! Solve function
    void solve(SolveFor = DISP);

    //! Apply essential BC
    void applyEBC(EigenResult & R);

  protected:
    State*          _myState;
    Model*          _myModel;
    vector<int  > & _DoFid;
    vector<Real > & _DoFvalues;
    SolverType      _linSolType;
    Real            _NRtol; 
    int             _NRmaxIter;
    
  };

}

#endif // __EigenNRsolver_h__

