// -*- C++ -*-
/*!
  \file Mesh.h
  \brief Class for a mesh object.
*/

#if !defined(__Mesh_h__)
#define __Mesh_h__

#include "voom.h"
#include "GeomElement.h"
#include "FEgeomElement.h"
#include "Quadrature.h"
#include "Shape.h"

namespace voom
{
  //! A mesh base class, defining node positions and element connectivities
  /*!  This class defines an interface for accessing the nodal
    position vectors and element connectivities of a generic
    finite-element-style mesh.  The implementation details are left
    for derived classes to define, but the basic idea is that some
    number of nodes with positions in 3-dimensional space are
    connected by "elements" which are generalized to mean any set of
    nodes.  Meshes in 1-D and 2-D, are embedded in 3-D space.

    In the case of meshfree, the "elements" will correspond to
    material points where fields and their gradients are computed, so
    the element connectivity will then be a list of neighbor nodes,
    i.e., nodes within the "search radius" or basis function support
    size.

    Mesh classes are also responsible for i/o.  They provide methods to
    read and write input files, including the case of restarts.
  */
  // enum MESHTYPE { FINITEELEMENT, MESHFREERKPM, MESHFREELME };

  class Mesh
  {
  public:
    //! Input-file based Constructor
    Mesh(const string Nodes, const string ConnTable);

    //! Position only based constructor
    Mesh(const vector<VectorXd > &  X): _X(X) {};

    //! Destructor
    // Destructor made virtual due to inheritance
    virtual ~Mesh() {
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
    const vector<GeomElement* > & getElements() {
      return _elements;
    }

    //! Get element type (assumed one element type per mesh - no mixed meshes for now)
    string getElementType() { return _elementType; }

    // //! Return mapping between local and global DoF and ghost DoF
    // const vector<int > & getLocalDoF() { return _localDoF; };
    // const vector<int > & getGhostDoF() { return _ghostDoF; };

    //! Gateway to mesh class
    // static Mesh* New(const string inputFile);

  protected:
    //! Default constructor is protected because it should be called only by derived classes, not from outside
    Mesh() {};
    // //! Protected because elements are not initialized here and should be used only from a derived class
    // Mesh(const vector<VectorXd > &  Positions,
    // 	 const vector<int > & LocalDoF,
    // 	 const vector<int > & GhostDoF):
    //   _positions(Positions), _localDoF(LocalDoF), _ghostDoF(GhostDoF) {};

    //! Nodal positions
    vector<VectorXd >     _X;

    //! List of Elements
    vector<GeomElement* > _elements;

    //! Element type
    string _elementType;

    // //! Map of Degrees of freedom local ID to global ID
    // vector<int >          _localDoF;
    // vector<int >          _ghostDoF;

  }; // Class Mesh

} //namespace voom

#endif // __Mesh_h__
