//-*-C++-*-

#ifndef __GelMesh_h__
#define __GelMesh_h__
#include "voom.h"
#include "GelElement.h"

namespace voom{
  
  class GelMesh
  {
  public:
    //! Constructor from input file
    GelMesh(const string Nodes, const string ConnTable);

    //! Destructor
    virtual ~GelMesh() {
      for(uint i = 0; i < _elements.size(); i++)
        delete _elements[i];
    }
    
    //! Get mesh dimension   
    uint getDimension() { return _X[0].size(); }

    //! Get position data 
    VectorXd getX() {
      uint Xsize = _X.size();
      uint dim = _X[0].size();
      VectorXd X0 = VectorXd::Zero(Xsize*dim);
      for (int i = 0; i < Xsize; i++) {
        for (int j = 0; j < dim; j++) {
          X0(i*dim + j) = _X[i](j);
        }
      }
      return X0;
    }

    const VectorXd & getX(const int nodeId) {
      return _X[nodeId];
    }

    //! Get position component
    Real getX(const int nodeId, const uint dof) {
      return getX(nodeId)(dof);
    }

    //! Get number of Nodes
    int getNumberOfNodes() { return _X.size(); }

    //! Get number of elements
    int getNumberOfElements() { return _elements.size(); }

    //! Get list of elements
    const vector<GelElement* > & getElements() {
      return _elements;
    }

  protected:
    
    //! Nodal positions
    vector<VectorXd >     _X;

    //! List of Elements
    vector<GelElement* > _elements;
    

  };
}

#endif // __FEMesh_h__
