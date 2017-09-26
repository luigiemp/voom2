//-*-C++-*-
#include "Model.h"

namespace voom {

  // Compute //
  void Model::compute(Result* R) {

    // Empty results before new calculations
    R->resetResults(R->getRequest());
    // Sum contributions from all bodies
    for (int i = 0; i < _bodies.size(); i++) {
      _bodies[i]->compute(R);
    }
    // Finalize Model stiffness only after summing contributions from all bodies
    if ( R->getRequest() & STIFFNESS ) {
      R->FinalizeGlobalStiffnessAssembly(); 
    }

  }



  // Consistency Checks //
  void Model::checkConsistency(Result* R, Real perturbationFactor, int request, 
			       Real h, Real tol)
  {
    // Perturb field randomly to change from reference configuration
    // Save perturbed field to set the configuration back to reference at the end of the test
    int DoFnum = _myState->getDOFcount();
    vector<Real > perturb(DoFnum, 0.0);
    srand( time(NULL) );
    for(int i = 0; i < DoFnum; i++) {
      Real randomNum =  perturbationFactor*(Real(rand())/RAND_MAX - 0.5);
      perturb[i] = randomNum;
      _myState->linearizedUpdate(i, randomNum);
    }      

    // Force Check
    if ( request & FORCE ) {
      Real error = 0.0, norm = 0.0;

      R->setRequest(FORCE); // First compute forces numerically
      R->resetResults(FORCE);
      this->compute(R);

      R->setRequest(ENERGY); // Reset result request so that only energy is computed 
      R->resetResults(ENERGY);
      this->compute(R);

      cout << "Model energy at test start = " <<  R->getEnergy() << endl;

      for(int i = 0; i < DoFnum; i++) {
	// Perturb +h
	_myState->linearizedUpdate(i, h);
	R->resetResults(ENERGY);
	this->compute(R);
	Real Wplus = R->getEnergy();
	
	// Perturb -2h
	_myState->linearizedUpdate(i, -2.0*h);
	R->resetResults(ENERGY);
	this->compute(R);
	Real Wminus = R->getEnergy();
	
	// Bring back to original position
	_myState->linearizedUpdate(i, h);
	
	error += pow( (Wplus-Wminus)/(2.0*h) - 
		      R->getResidual(i), 2);
	norm += pow(R->getResidual(i), 2);
      } // Loop over dimension
      error = sqrt(error);
      norm  = sqrt(norm);
      
      if ( abs(error) < norm * tol) {
	cout << "** Elliptic Model Force consistency check PASSED" << endl;
	cout << "** Error: " << error << " Norm: " << norm  << " Norm*tol: " << norm*tol << endl;
      }
      else {
	cout << "** Elliptic Model Force consistency check FAILED" << endl;
	cout << "** Error: " << error << " Norm: " << norm << " Norm*tol: " << norm*tol << endl;
      }
    } // Check Forces loop

    // Stiffness check
    if ( request & STIFFNESS ) {
      Real error = 0.0, norm = 0.0;

      R->setRequest(STIFFNESS); // First compute stiffness numerically
      R->resetResults(STIFFNESS);
      this->compute(R);

      R->setRequest(FORCE); // Reset result request so that only forces are computed 
      for(int i = 0; i < DoFnum; i++) {
	for(int j = 0; j < DoFnum; j++) {
	  // Perturb +h
	  _myState->linearizedUpdate(j, h);
	  R->resetResults(FORCE);
	  this->compute(R);
	  Real Fplus = R->getResidual(i);

	  // Perturb -2h
	  _myState->linearizedUpdate(j, -2*h);
	  R->resetResults(FORCE);
	  this->compute(R);
	  Real Fminus = R->getResidual(i);

	  // Bring back to original position
	  _myState->linearizedUpdate(j, h);

	  // Computing Error and Norm;
	  error += pow((Fplus - Fminus)/(2.*h) - R->getStiffness(i, j), 2.0);
	  norm += pow( R->getStiffness(i, j), 2.0); 
	  // cout << R->getStiffness(a*nodeDoF+i, b*nodeDoF+j) << "\t" << (Fplus - Fminus)/(2.*h) << endl;
	} // j loop
      } // i loop
      
      error = sqrt(error);
      norm  = sqrt(norm);
      if ( abs(error) < norm * tol) {
	cout << "** Elliptic Model Hessian consistency check PASSED" << endl;
	cout << "** Error: " << error << " Norm: " << norm << " Norm*tol: " << norm*tol << endl << endl;
      }
      else {
	cout << "** Elliptic Model Hessian consistency check FAILED" << endl;
	cout << "** Error: " << error << " Norm: " << norm << " Norm*tol: " << norm*tol << endl << endl;
      }
    }// Check Stiffness    

    // Reset field to initial values
    for(int i = 0; i < DoFnum; i++) {
      _myState->linearizedUpdate(i, -perturb[i]);
    } 
    
  } // Check consistency

} // Namespace voom
