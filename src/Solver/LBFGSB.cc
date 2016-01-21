// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          Adapted from voom
//
// Reference:
//  Nonlinear quasi-Newton solver using L-BFGS-B code of Jorge Nocedal.
//  http://www.ece.northwestern.edu/~nocedal/lbfgsb.html
//
/////////////////////////////////////////////////////////////////////////

#include "LBFGSB.h"

extern "C" void setulb_(int * n, int *m, 
			double * x, double * l, double * u, 
			int * nbd, 
			double * f, double * g, 
			double * factr, double * pgtol,
			double * wa, 
			int * iwa,
			char * task, 
			int * iprint,  
			char * csave, 
			int * lsave, 
			int * isave, 
			double * dsave );

namespace voom
{

  //! Constructor
  LBFGSB::LBFGSB( EllipticModel *myModel,
		  EigenEllipticResult* myResults, 
		  int m,
		  double factr, double pgtol,
		  int iprint, int maxIterations):
    _myModel(myModel), _myResults(myResults),
    _m(m),
    _factr(factr), _pgtol(pgtol),
    _iprint(iprint), _maxIterations(maxIterations), _iterNo(0) {
    
    _counter = 0; //********************** Temporary *********************
    _n = ( (_myModel->getMesh())->getNumberOfNodes() )*( _myModel->getDoFperNode() );
    this->resize(_n);
    _projg = 0.0;

    // Initialize x vector
    _myModel->getField(_x);
  }; // Constructor



  void LBFGSB::solve() {
    
    // set up misc. arrays and data
    char task[60], csave[60];
    for(int i=0; i<60; i++) task[i] = csave[i] = '\0';
    
    int lsave[4];
    for(int i=0; i<4; i++) lsave[i]=0;
    
    double dsave[29];
    for(int i=0; i<29; i++) dsave[i]=0.0;
    
    int isave[44];
    for(int i=0; i<44; i++) isave[i]=0;
    
    cout << "================================================================================"
	 << endl << endl
	 << "Starting BFGS iterations."
	 << endl << endl;

    cout << setw(14) << scientific << right << "|proj g|" 
	 << setw(14) << scientific << right << "|g|" 
	 << setw(14) << scientific  << right << "f" 
	 << setw(14) << right << "iterations" 
	 << setw(14) << right << "evaluations" 
	 << endl
	 << "--------------------------------------------------------------------------------"
	 << endl;

    // We start the iteration by initializing task.
    sprintf(task,"START");

    this->_setFieldAndCompute();

    // ------- the beginning of the loop ----------
    while(true) {
      
      // This is the call to the L-BFGS-B code.
      setulb_(&_n, &_m, _x.data(), _l.data(), _u.data(), _nbd.data(), 
	      &_f, _g.data(), &_factr, &_pgtol, _wa.data(),_iwa.data(), 
	      &(task[0]), &_iprint, &(csave[0]),
	      &(lsave[0]),&(isave[0]),&(dsave[0]));

      if(strncmp(task,"FG",2) == 0) {
	// The minimization routine has returned to request the
	// function f and gradient g values at the current x.

	this->_setFieldAndCompute();

	// Go back to the minimization routine.
	continue;
      }
      else if( strncmp(task,"NEW_X",5) == 0 ) {
	// stop if maximum number of iterations has been reached
	if ( _maxIterations > 0 && isave[29] > _maxIterations ) {
	  break;
	}

	// The minimization routine has returned with a new iterate,
	// and we have opted to continue the iteration.
	if ( _iprint > 0 && isave[29]%_iprint == 0 ){
	  cout << setw(14) << scientific << right << dsave[12]
	       << setw(14) << scientific << right << _g.lpNorm<Infinity>()
	       << setw(14) << scientific  << right << _f
	       << setw(14) << right << isave[29] 
	       << setw(14) << right << isave[33] 
	       << endl;
	}

	continue;
      }
      else if( strncmp(task,"CONV",4) == 0 ) {
	this->_setFieldAndCompute();
	cout << task << endl;
	break;
      }
      else if( strncmp(task,"ABNORM",6) == 0 ) {
	this->_setFieldAndCompute();
	cout << task << endl;
	break;
      }

      // If task is neither FG nor NEW_X we terminate execution.
      else {
	cout << task << endl;
	break;
      }

      // ---------- the end of the loop -------------
    }
    cout << setw(14) << scientific << right << dsave[12]
	 << setw(14) << scientific << right << _g.lpNorm<Infinity>()
	 << setw(14) << scientific << right << _f
	 << setw(14) << right << isave[29] 
	 << setw(14) << right << isave[33] 
	 << endl << endl
	 << "================================================================================"
	 << endl << endl;
    
    cout.unsetf(ios_base::scientific);


    // ---------- the end of solve() -------------
    _projg = dsave[12];
    _iterNo = isave[29];

    this->_setFieldAndCompute();
    
  } // solve 

}; // end namespace
