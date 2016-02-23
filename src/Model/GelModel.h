//-*-C++-*-
/*!
  \file EPModel.h

  \brief Implementation of Cardiac Electrophysiology. This is the parent
  class for cardiac EP. Monodomain and Bidomain classes will derive from 
  this class
*/

#ifndef __GelModel_h__
#define __GelModel_h__
#include "Model.h"
#include "GelResult.h"

#include "Spring.h"
#include "EntropicSpring.h"
#include "AngleSpring.h"
#include "BrownianRod.h"

#include "Constraint.h"
#include "voom.h"
#include "Node.h"
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include "Crosslink.h"
#include "Motor.h"
#include "Grid.h"
#include "TwoBodyPotential.h"
#include "IntersectionFinder.h"
#include "PeriodicTie.h"
#include "PeriodicBox.h"
#include "LeesEdwards.h"
#include "PinchForce.h"
#include "NematicProbTable.h"
#include "GelInput.h"
#include "Model.h"

namespace voom{
  class GelModel : public Model {
  private:
    //! Periodic box for Lees-Edward BC
    PeriodicBox * _box;

  public:
    //! Constructor
    GelModel(Mesh* Filaments, const GelInput inputFile, const uint NodeDof);

    //! Destructor
    virtual ~GelModel() {;}
    
    
    //! Initialize field
    // From constant value
    void initializeField(const Real value = 1.0) {
      const uint numNodes = _myMesh->getNumberOfNodes();
      const uint dim = _myMesh->getDimension();

      for (uint i = 0; i < numNodes; i++) 
	for (uint j = 0; j < dim; j++)
	  _field[i*dim+j] = _myMesh->getX(i,j)*value; // value = isotropic expansion/shrinking
    };

    //! From array
    void initializeField(const Real* value) {
      _field.assign(value, value+_field.size());
    };

    //! Linearized update
    void linearizedUpdate(const Real* localValues, Real fact) {
      const int nLocalDof = (_myMesh->getNumberOfNodes())*_nodeDoF;
      for(uint i = 0; i < nLocalDof; i++)
	_field[i] += fact*localValues[i];
    };

    // One value at the time (Node ID, dof index, value)
    void linearizedUpdate(const int id, const int dof, const Real value) {
      // const uint dim = _myMesh->getDimension();
      // assert( id < _field.size() && dof < dim );
      _field[id*_nodeDoF + dof] += value;
    }

    // One value at the time (Node ID, dof index, value)
    void linearizedUpdate(const int dof, const Real value) {
      _field[dof] += value;
    }
    
    void setField(uint dof, Real value) {
      _field[dof] = value;
    }
    void setField(const Real* value) {
      _field.assign(value, value+_field.size());
    };
    void getField(vector<double> & x) {
      assert(x.size() == _field.size());
      x = _field;
    }

    void setPrevField() {
      _prevField = _field;
    };

    void printField() {
      int i = 0;
      while (i < _field.size()) {
	for (uint j = 0; j < _nodeDoF; j++) {
	  cout << _field[i] << " ";
	  i++;
	}
	cout << endl;
      }
    }

    void writeField(string OutputFile, int step) {
      // Create outputFile name
      stringstream FileNameStream;
      FileNameStream << OutputFile << step << ".dat";
      ofstream out;
      out.open( (FileNameStream.str()).c_str() );
  
      out << _field.size() << endl;
      for (uint i = 0; i < _field.size(); i++) {
	out << setprecision(15) << _field[i] << endl;
      }
      out.close();
    }



    //! compute function call
    void compute(GelResult& R) {;}
  };
}

#endif
