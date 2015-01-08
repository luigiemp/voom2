//-*-C++-*-
/*!
  \file MechancsSolver.h

  \brief Derived class for solving mechanics problem
*/

#ifndef __Solver_h__
#define __Solver_h__
#include "Solver.h"

namespace voom{
  template<class nodeType, unsigned int dim>
  class MechanicsSolver:public Solver {
  private:
    //! Displacement
    Teuchos::RCP<Epetra_Vector>     _U;
    Teuchos::RCP<Epetra_Vector>     _ghostU;

    //! Internal Force
    Teuchos::RCP<Epetra_Vector>     _FInternal;
    Teuchos::RCP<Epetra_Vector>     _ghostFInternal
;
    //! External Force
    Teuchos::RCP<Epetra_Vector>     _FExternal;
    Teuchos::RCP<Epetra_Vector>     _ghostFExternal

    //! Stiffness matrix
    Teuchos::RCP<Epetra_CrsMatrix>  _K;

    //! Epetra Maps
    Epetra_Map*                     _sourceMap;
    Epetra_Map*                     _targetMap;

    //! Epetra Exporter
    Epetra_Export*                  _exporter;

  public:
    //! Constructor
    MechanicsSolver(Epetra_Mpicomm* comm, Model<nodeType,dim>* myModel);

    //! Destructor
    ~MechanicsSolver() {
      delete _sourceMap; delete _targetMap; delete _exporter;
    }

    //! Solver the model using Newton Raphson Solver
    void Solve();
  };
}
#endif
