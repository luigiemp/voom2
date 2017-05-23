//-*-C++-*-
#include "Model.h"

namespace voom {

  // Compute //
  void Model::compute(Result* R) {
    for (int i = 0; i < _bodies.size(); i++) {
      _bodies[i]->compute(R);
    }
  }



  // Consistency Checks //
  void Model::checkConsistency(Result* R,  int nodeNum, int nodeDoF, 
			       Real perturbationFactor, 
			       int request = 6, Real h = 1e-6, Real tol = 1e-6)
  {
    // Perturb field randomly to change from reference configuration
    // Save perturbed field to set the configuration back to reference 
    // at the end of the test
    vector<Real > perturb(nodeNum * nodeDoF, 0.0);
    srand( time(NULL) );
    for(int a = 0; a < nodeNum; a++) {
      for(int i = 0; i < nodeDoF; i++) {
	Real randomNum =  perturbationFactor*(Real(rand())/RAND_MAX - 0.5);
	perturb[a*nodeDoF + i] = randomNum;
	R->linearizedUpdate(a*nodeDoF + i, randomNum);
      }
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

      for(int a = 0; a < nodeNum; a++) {
	for(int i = 0; i < nodeDoF; i++) {
	  // Perturb +h
	  R->linearizedUpdate(a*nodeDoF + i, h);
	  R->resetResults(ENERGY);
	  this->compute(R);
	  Real Wplus = R->getEnergy();
	  
	  // Perturb -2h
	  R->linearizedUpdate(a*nodeDoF + i, -2.0*h);
	  R->resetResults(ENERGY);
	  this->compute(R);
	  Real Wminus = R->getEnergy();
	  
	  // Bring back to original position
	  R->linearizedUpdate(a*nodeDoF + i, h);
	  
	  error += pow( (Wplus-Wminus)/(2.0*h) - 
			R->getResidual(a*nodeDoF + i), 2);
	  norm += pow(R->getResidual(a*nodeDoF + i), 2);
	} // Loop over dimension
      } // Loop over nodes
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
      R->FinalizeGlobalStiffnessAssembly();

      R->setRequest(FORCE); // Reset result request so that only forces are computed 
      for(int a = 0; a < nodeNum; a++) {
	for(int i = 0; i < nodeDoF; i++) {
	  for(int b = 0; b < nodeNum; b++) {
	    for(int j = 0; j < nodeDoF; j++) {
	      // Perturb +h
	      R->linearizedUpdate(b*nodeDoF + j, h);
	      R->resetResults(FORCE);
	      this->compute(R);
	      Real Fplus = R->getResidual(a*nodeDoF + i);

	      // Perturb -2h
	      R->linearizedUpdate(b*nodeDoF + j, -2*h);
	      R->resetResults(FORCE);
	      this->compute(R);
	      Real Fminus = R->getResidual(a*nodeDoF + i);

	      // Bring back to original position
	      R->linearizedUpdate(b*nodeDoF + j, h);

	      // Computing Error and Norm;
	      error += pow((Fplus - Fminus)/(2.*h) - R->getStiffness(a*nodeDoF+i, b*nodeDoF+j), 2.0);
	      norm += pow( R->getStiffness(a*nodeDoF+i, b*nodeDoF+j), 2.0); 
	      // cout << R->getStiffness(a*nodeDoF+i, b*nodeDoF+j) << "\t" << (Fplus - Fminus)/(2.*h) << endl;
	    } // j loop
	  } // b loop
	} // i loop
      } // a loop

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
    } // Check Stiffness    

    // Reset field to initial values
    for(int a = 0; a < nodeNum; a++) {
      for(int i = 0; i < nodeDoF; i++) {
	R->linearizedUpdate(a*nodeDoF + i, -perturb[a*nodeDoF + i]);
      }
    } 
    
  } // Check consistency



} // Namespace voom
