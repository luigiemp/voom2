#include "EigenNRsolver.h"

namespace voom
{
  void EigenNRsolver::solve()
  { 
    // Identify number of DoF and unique material models
    uint PbDoF = ( (_myModel->getMesh())->getNumberOfNodes() )*( _myModel->getDoFperNode() );
    uint NumMat =  _myModel->getNumMat();

    // Initialize field vector
    for (int i = 0; i < _DoFid.size(); i++) {
      _myModel->setField(_DoFid[i], _DoFvalues[i]); // Set known displacements
    }
    
    // Initialize solver parameters
    Real error = 1.0;
    uint iter = 0;
    
    // NR loop
    while (iter < _NRmaxIter && error > _NRtol)
    {
      // Create Eigen Elliptic results
      EigenEllipticResult myResults( PbDoF, NumMat );

      // Compute stiffness and residual
      myResults.setRequest(6); 
      _myModel->compute(myResults);

      // Apply essential boundary conditions
      this->applyEBC(myResults);

      // Solve
      SimplicialCholesky<SparseMatrix<Real > > chol(*(myResults._stiffness));  // performs a Cholesky factorization of _stiffness
      VectorXd Deltax = chol.solve(*(myResults._residual));         // use the factorization to solve for the given right hand side
 

      // Change solver solution so that EBC are not added to field in Model multiple times
      for (int i = 0; i < _DoFid.size(); i++) {
	Deltax(_DoFid[i]) = 0.0; // So we do not add the known EBC mutiple times
      }
      _myModel->linearizedUpdate(Deltax.data());
    
      // Update iter and error
      iter++;
      // Compute NR error
         // myResults.setRequest(2);
         // _myModel->compute(myResults);
         // error = (myResults._residual).norm();
      error = Deltax.norm();
      cout << "NR iter = " << iter << "   - NR error = " << error << endl;

    } // while loop

  } // solve function



  void EigenNRsolver::applyEBC(EigenEllipticResult & myResults) {
    
// | A    B |  | x       |   | f |
// |        |  |         | = |   |
// | B^T  C |  | \bar{x} |   | r |

    // Build auxiliary set with EBC dof id
    int PbDoF = (myResults._stiffness)->rows();
    set<int > EBC;
    VectorXd X = VectorXd::Zero(PbDoF);
    for (int i = 0; i < _DoFid.size(); i++) {
      EBC.insert(_DoFid[i]);
      X(_DoFid[i]) = _DoFvalues[i];
    }
    
    // Initialize sparse matrix B
    SparseMatrix<Real > B(PbDoF, PbDoF);
    vector<Triplet<Real > > BtripletList;
    BtripletList.reserve(PbDoF*3);

    // Fill in matrix B and change matrix A
    for (int k = 0; k < (myResults._stiffness)->outerSize(); k++) {
      for (SparseMatrix<Real>::InnerIterator it(*(myResults._stiffness), k); it; ++it) {	 
	  // int col = it.col(); == k here because we are using the default column major storing
	  int row = it.row();
	  if ( EBC.find(k) != EBC.end() ) {
	    if ( EBC.find(row) == EBC.end() ) { // extracting B and setting it = 0 afterward
	      BtripletList.push_back( Triplet<Real >(row, k, it.value() ) );
	      it.valueRef() = 0;
	    }
	    else if ( row != k ) { // setting C = I
	      it.valueRef() = 0;
	    }
	    else {
	      it.valueRef() = 1;
	    }
	  }

	  if ( EBC.find(k) == EBC.end() ) { // setting B^T = 0
	    it.valueRef() = 0;
	  }
	    
	} // end of loop over rows
      } // end of loop over columns

    // Clean A
    (myResults._stiffness)->prune(0.0);
      
    // Construct B
    B.setFromTriplets(BtripletList.begin(), BtripletList.end());
    B.makeCompressed();

    // Add -(B * \bar{x}) to RHS
    *(myResults._residual) -= B*X;
    
  } // applyEBC function



} // namespace voom
