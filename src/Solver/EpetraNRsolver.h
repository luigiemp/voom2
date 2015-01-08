//-*-C++-*-
/*!
  \file EpetraNRsolver.h
  \brief Newton Raphson solver based on Epetra sparse solver
*/

#ifndef __EpetraNRsolver_h__
#define __EpetraNRsolver_h__

#include "voom.h"
#include "Mesh.h"
#include "EllipticModel.h"
#include "EpetraEllipticResult.h"
#include "EpetraSparseSolver.h"

namespace voom{
			       
  class EpetraNRsolver
  {
  public:
    //! Constructor
    EpetraNRsolver(Epetra_MpiComm* mpicomm,
		   Mesh *myMesh,
		   EllipticModel *myModel,
		   const vector<uint >& DoFids, 
		   const vector<Real >& DoFvalues,
		   Real NRtol = 1.0e-8, uint NRmaxIter = 100):
      _mpicomm(mpicomm),
      _myMesh(myMesh),
      _myModel(myModel),
      _EBCdof(DoFids),
      _EBCvalues(DoFvalues),
      _NRtol(NRtol), _NRmaxIter(NRmaxIter) {};

    //! Destructor
    ~EpetraNRsolver() {};
      
    //! Solve function
    void solve();

  protected:
    Epetra_MpiComm* _mpicomm;
    Mesh* _myMesh;
    EllipticModel* _myModel;
    const vector<uint >& _EBCdof;
    const vector<Real >& _EBCvalues;
    Real _NRtol; 
    uint _NRmaxIter;
  };
}

#endif // __EpetraNRsolver_h__

