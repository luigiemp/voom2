//-*-C++-*-
#include "Model.h"

namespace voom {

  // Consistency Checks //
  void Model::checkConsistency(Result & R, Real perturbationFactor, int request,
				       Real h, Real tol)
  {
    // Check only for local nodes
    const uint nodeNum   = this->getNumberOfNodes();
    const uint nLocalDoF = nodeNum*_nodeDoF;

    // Perturb field randomly to change from reference configuration
    // Save perturbed field to set the configuration back to reference 
    // at the end of the test
    vector<Real > perturb(nLocalDoF, 0.0);
    srand( time(NULL) );
    for(uint a = 0; a < nodeNum; a++) {
      for(uint i = 0; i < _nodeDoF; i++) {
	Real randomNum =  perturbationFactor*(Real(rand())/RAND_MAX - 0.5);
	perturb[a*_nodeDoF + i] = randomNum;
	this->linearizedUpdate(a, i, randomNum);
      }
    }      
    


    // Force Check
    if ( request & FORCE ) {
      Real error = 0.0, norm = 0.0;

      R.setRequest(3); // First compute forces numerically
      this->compute(R);

      R.setRequest(1); // Reset result request so that only energy is computed 
      // this->compute(R);

      cout << "Model energy at test start = " <<  R.getEnergy() << endl;

      for(int a = 0; a < nodeNum; a++) {
	for(int i = 0; i < _nodeDoF; i++) {
	  // Perturb +h
	  this->linearizedUpdate(a, i, h);
	  this->compute(R);
	  Real Wplus = R.getEnergy();
	  
	  // Perturb -2h
	  this->linearizedUpdate(a, i, -2*h);
	  this->compute(R);
	  Real Wminus = R.getEnergy();
	  
	  // Bring back to original position
	  this->linearizedUpdate(a, i, h);
	  
	  error += pow( (Wplus-Wminus)/(2.*h) - 
			R.getResidual(a*_nodeDoF + i), 2);
	  norm += pow(R.getResidual(a*_nodeDoF + i), 2);
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

      R.setRequest(4); // First compute stiffness numerically
      this->compute(R);

      R.setRequest(2); // Reset result request so that only forces are computed 
      for(int a = 0; a < nodeNum; a++) {
	for(int i = 0; i < _nodeDoF; i++) {
	  for(int b = 0; b < nodeNum; b++) {
	    for(int j = 0; j < _nodeDoF; j++) {
	      // Perturb +h
	      this->linearizedUpdate(b, j, h);
	      this->compute(R);
	      Real Fplus = R.getResidual(a*_nodeDoF + i);

	      // Perturb -2h
	      this->linearizedUpdate(b, j, -2*h);
	      this->compute(R);
	      Real Fminus = R.getResidual(a*_nodeDoF + i);

	      // Bring back to original position
	      this->linearizedUpdate(b, j, h);

	      // Computing Error and Norm;
	      error += pow((Fplus - Fminus)/(2.*h) - R.getStiffness(a*_nodeDoF+i, b*_nodeDoF+j), 2.0);
	      norm += pow( R.getStiffness(a*_nodeDoF+i, b*_nodeDoF+j), 2.0); 
	      cout << R.getStiffness(a*_nodeDoF+i, b*_nodeDoF+j) << "\t" << (Fplus - Fminus)/(2.*h) << endl;
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
    for(uint a = 0; a < nodeNum; a++) {
      for(uint i = 0; i < _nodeDoF; i++) {
	this->linearizedUpdate(a, i, -perturb[a*_nodeDoF + i]);
      }
    } 
    
  } // Check consistency


} // Namespace voom
