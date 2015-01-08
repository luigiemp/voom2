#include "EpetraSparseSolver.h"

namespace voom{
  // Constructor
  EpetraSparseSolver::EpetraSparseSolver(Epetra_MpiComm* Comm,
					 Real LinearSolverTolerance, int LinearSolverMaxIterations,
					 Epetra_Vector* RHS, Epetra_Vector* ghostRHS,
					 Epetra_FECrsMatrix* SystemMatrix,
					 Epetra_Map * TargetMap,
					 Epetra_Map * SourceMap):
    _comm(Comm), 
    _linearSolverTolerance(LinearSolverTolerance),
    _linearSolverMaxIterations(LinearSolverMaxIterations),
    _B(RHS), _ghostB(ghostRHS), _A(SystemMatrix),
    _targetMap(TargetMap), _sourceMap(SourceMap)
  {
    // create Epetra_Exporter
    _exporter  = new Epetra_Export(*_sourceMap, *_targetMap);

    // Create Epetra solution vectors
    _X       = new Epetra_Vector(*_targetMap);
    _ghostX  = new Epetra_Vector(*_sourceMap);
    _X->ExtractView( &_extractedX );
    _ghostX->ExtractView(&_extractedGhostX);

    // Create Solver instances
    _problem = new Epetra_LinearProblem(_A, _X, _B);
    _solver  = new AztecOO( *_problem );
    // Solver settings can be modified after creating a Solver object
    _solver->SetAztecOption(AZ_solver, AZ_gmres);
    _solver->SetAztecOption(AZ_precond, AZ_Jacobi);
    _solver->SetAztecOption(AZ_diagnostics, AZ_none);
    _solver->SetAztecOption(AZ_output, AZ_warnings);
    _solver->SetAztecOption(AZ_conv, AZ_noscaled);
  }
 


  // Apply boundary condition
  void EpetraSparseSolver::applyBC(Mesh * myMesh,
				   const vector<uint > & EBC_DoF,
				   const vector<Real> & EBC_Values)
  {
    // Maps to be computed to modify system matrix A and apply BC
    map<int, int > globalToLocalID;      // In model we already have local to global mapping
    map<int, vector<int > > zeroColumn;  // Column to be zeroed
    map<int, vector<Real > > zeroValues; // Column values before are zeroed

    // Construct globalToLocalID
    vector<int > localDoF = myMesh->getLocalDoF();
    vector<int > ghostDoF = myMesh->getGhostDoF();
    const int nlocalDoF = localDoF.size();

    for(int i = 0; i < nlocalDoF; i++) 
      globalToLocalID[ localDoF[i] ] = i;
    for(int i = 0; i < ghostDoF.size(); i++)
      globalToLocalID[ ghostDoF[i] ] = i + nlocalDoF;
    
    // All EBC are ONLY for local nodes.
    Real *extracted_B, *extracted_ghostB;
    _B->ExtractView(&extracted_B);
    _ghostB->ExtractView(&extracted_ghostB);


    // int myRank = 0;
    // MPI_Comm_rank( MPI_COMM_WORLD, &myRank );
    
    for(int i = 0; i < EBC_DoF.size(); i++) {
      int index = localDoF[ EBC_DoF[i] ];     // Global Row ID
      int numvals = _A->NumGlobalEntries( index );
      vector<Real > rowVals(numvals);
      vector<int >  rowInds(numvals);

      // Get non zero entries in a row
      _A->ExtractGlobalRowCopy( index, numvals, numvals, 
				&rowVals[0], &rowInds[0]); 
      // Loop over the global column Ids
      vector<Real > newVals(numvals, 0.0);
      for(int colId = 0; colId < numvals; colId++)
      {
    	const int columnId = rowInds[colId];
    	// If not diagonal make it 0 or else 1
    	zeroColumn[ columnId ].push_back( index );
    	if ( columnId != index) {
	  // const int localColId = globalToLocalID[ columnId ];
    	  // if ( localColId  < nlocalDoF)
    	  //   extracted_B[localColId] -= rowVals[ colId ]*EBC_Values[i];
    	  // else
    	  //   extracted_ghostB[ localColId - nlocalDoF] 
    	  //     -= rowVals[ colId ]*EBC_Values[i];
    	  zeroValues[ columnId ].push_back( 0.0 );
    	} else {
    	  zeroValues[ columnId ].push_back( 1.0 );
    	  newVals[colId] = 1.;
    	}
      } // Loop over column
      // Replace row with new values
      _A->ReplaceGlobalValues( index, numvals,
			       &(newVals[0]),&(rowInds[0]));
    }

   

    for(int i = 0; i < EBC_DoF.size(); i++) {
      int index = EBC_DoF[i];  // local Row ID
      if ( index < nlocalDoF )
    	extracted_B[index] = EBC_Values[i];
      else
    	extracted_ghostB[index - nlocalDoF] = EBC_Values[i];
    } // Loop over EBC_DoF

    // Loop over zerocolumnID vector
    //       for(int rowId = 0; rowId < zeroColumn.size(); rowId++) {
    for(map<int, vector<int> >::iterator it = zeroColumn.begin();
    	it != zeroColumn.end(); it++)
    {
      // Zero out entires
      const int rowId = it->first;
      _A->ReplaceGlobalValues(rowId, zeroColumn[rowId].size(), 
			      &(zeroValues[rowId][0]), 
			      &(zeroColumn[rowId][0]) );
    } // Loop over zerocolumnID vector

  } // End of applyEBC



  void EpetraSparseSolver::setSolution(vector<uint > DoFids, vector<Real > DoFvalues)
  { 
    int nLocalDoF = _X->MyLength();
    for(uint i = 0; i < DoFids.size(); i++) {
      int index = DoFids[i];
      if ( index < nLocalDoF )
    	_extractedX[index] = DoFvalues[i];
      else
    	_extractedGhostX[index - nLocalDoF] = DoFvalues[i];
    }
  }

} // namespace voom
