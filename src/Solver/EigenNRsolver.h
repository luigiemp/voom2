//-*-C++-*-
/*!
  \file EigenNRsolver.h
  \brief Newton Raphson solver based on Eigen sparse solver
*/

#ifndef __EigenNRsolver_h__
#define __EigenNRsolver_h__

#include "voom.h"
#include "Mesh.h"
#include "EllipticModel.h"
#include "EigenEllipticResult.h"

namespace voom{
			       
  class EigenNRsolver
  {
  public:
    //! Constructor
    EigenNRsolver( EllipticModel *myModel, 
		   vector<int > & DoFid,
		   vector<Real > & DoFvalues,
		   string LinSolType,
		   Real NRtol = 1.0e-8, uint NRmaxIter = 100):
      _myModel(myModel),
      _DoFid(DoFid), _DoFvalues(DoFvalues),
      _linSolType(LinSolType),
      _NRtol(NRtol), _NRmaxIter(NRmaxIter) {};

    //! Destructor
    ~EigenNRsolver() {};
      
    //! Solve function
    void solve();

    //! Apply essential BC
    void applyEBC(EigenEllipticResult & myResults);

  protected:
    EllipticModel*  _myModel;
    vector<int > &  _DoFid;
    vector<Real > & _DoFvalues;
    string          _linSolType;
    Real            _NRtol; 
    uint            _NRmaxIter;
    
  };

}

#endif // __EigenNRsolver_h__

