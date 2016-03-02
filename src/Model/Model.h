//-*-C++-*-
/*!
  \file Model.h
  \brief A base model class creates basic entities needs for solving in parallel
  and stores generic boundary conditions. 
 */
#ifndef __Model_h__
#define __Model_h__

#include "voom.h"
#include "Mesh.h"
#include "Result.h"

namespace voom {
  class Model {
  private :
    Model() {}; // Default constructor is private because it should never be used by a derived class.
  
  public:
    //! Constructor using an input file ... to be changed after implementation of input class - no longer current
    // Model(Mesh* myMesh, const string inputFile, const uint NodeDoF);

    //! Constructor to be used
    Model(Mesh* aMesh, const uint NodeDoF)
      : _myMesh(aMesh), _nodeDoF(NodeDoF) {
      // WSK: evetually we want to copy the mesh object that is
      // pointed to.  To do this we need to implement a virtual copy
      // method to be used like this:
      //
    }

    //! Destructor
    virtual ~Model() {};

    //! Model compute function
    virtual void compute(voom::Result * R) = 0;

    //! Check consistency - common to all models
    void checkConsistency(Result* R, Real perturbationFactor, int request = 6,
			  Real h = 1e-6, Real tol = 1e-6);

    //! Return number of unique material pointers
    virtual uint getNumMat() = 0;

    //! Get DoF per Node
    uint getDoFperNode() {
      return _nodeDoF;
    }
    
    //! GetMesh
    Mesh * getMesh() {
      return _myMesh;
    }

    //! Initialize field
    // From constant value
    virtual void initializeField(const Real value = 0.0) = 0;
 
    //! From array
    virtual void initializeField(const Real* value) = 0;

    //! Linearized update (local and ghost solution)
    // From solution array (Used by Solver)
    // virtual void linearizedUpdate(const Real* locaValues, 
    // 				  const Real* ghostValues) = 0;
    virtual void linearizedUpdate(const Real* locaValues, Real fact) = 0;

    // One value at the time (Node ID, DoF index, value)
    virtual void linearizedUpdate(const int id, const int dof, const Real value) = 0;

    // One value at the time (DoF index, value)
    virtual void linearizedUpdate(const int dof, const Real value) = 0;

    // Set solution field
    virtual void setField(uint dof, Real value) = 0;
    virtual void setField(const Real* value) = 0;

    // Get solution field
    virtual void getField(vector<double> & x) = 0;

    //! Set Previous Field
    virtual void setPrevField() = 0;

    // Print field
    virtual void printField() = 0;

    //! Output generation
    virtual void writeOutputVTK(const string OutputFile, int step) = 0;



  protected:
    Mesh*        _myMesh;

    //! DoF per node
    uint         _nodeDoF;
  
  }; // Model

} // namespace voom

#endif
