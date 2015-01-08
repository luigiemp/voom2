#include "EpetraNRsolver.h"

namespace voom
{
  void EpetraNRsolver::solve()
  {
    Real error = 1.0;
    uint iter = 0;

    // Set field in model based on EBC - needed in non-linear problems
    for (uint i = 0; i < _EBCdof.size(); i++) {
      _myModel->setField(_EBCdof[i], _EBCvalues[i]);
    } 

    // NR loop
    while (iter < _NRmaxIter && error > _NRtol)
    {
      // Create Epetra Elliptic results
      EpetraEllipticResult myResults(_mpicomm,
      				     _myMesh->getLocalDoF(),
       				     _myMesh->getGhostDoF()); 
      
      // Create Epetra Elliptic Sparse Solver
      // LEP: it would be nice to create this and Elliptic results once and for 
      // all outside the while loop but this would produce Epetra internal errors. 
      // I do not know why yet.
      EpetraSparseSolver mySparseSolver(_mpicomm,
					1.0e-12, 1000000,
					myResults.getResidual(), myResults.getGhostResidual(),
					myResults.getStiffness(),
					myResults.getTargetMap(),
					myResults.getSourceMap());

      // Compute stiffness and residual
      myResults.setRequest(6); 
      _myModel->compute(myResults);

      // Apply BC
      mySparseSolver.applyBC(_myMesh, _EBCdof, _EBCvalues);

      // Solve
      mySparseSolver.solve(-1.0);

      // Change solver solution so that EBC are not added to field in Model multiple times
      vector<Real > Zeros(_EBCdof.size(), 0.0);
      mySparseSolver.setSolution(_EBCdof, Zeros);
      // Update field
      mySparseSolver.updateModel(_myModel);
    
      // Update iter and error
      iter++;
      // Compute residual
      myResults.setRequest(2);
      _myModel->compute(myResults);

      error = mySparseSolver.DeltaRHSnorm(_EBCdof);
      cout << "NR iter = " << iter << "   - NR error = " << error << endl;

    } // while loop

  } // solve function

} // namespace voom
