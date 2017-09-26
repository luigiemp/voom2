// -*- C++ -*-
/*!
  \file Mesh.h
  \brief Class for a mesh object.
*/

#if !defined(__Mesh_h__)
#define __Mesh_h__

#include "voom.h"
#include "State.h"
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
  private:
    Mesh();

  public:
    //! Input-file based Constructor
    Mesh(string Nodes, string ConnTable, State* myState, int dofPerNode, bool CheckOverlap);
    Mesh(const vector<VectorXd > &  Positions, State* myState, int dofPerNode, bool CheckOverlap);

    //! Destructor
    // Destructor made virtual due to inheritance
    virtual ~Mesh() {
      for(uint i = 0; i < _elements.size(); i++)
	delete _elements[i];
    }

    //! Get number of elements
    int getNumberOfElements() { return _elements.size(); }

    //! Get list of elements
    const vector<GeomElement* > & getElements() { return _elements; }

    //! Get element type (assumed one element type per mesh - no mixed meshes for now)
    string getElementType() { return _elementType; }

    const vector<int > & getLocalToGlobal() { return _LocaltoGlobal; }
   
    //! Get X - useful if _gState is not the same as the one in model (e.g., see Potential Body)
    Vector3d getX(int nodeID) {return _gState->getX(nodeID); }

    // //! Return mapping between local and global DoF and ghost DoF
    // const vector<int > & getLocalDoF() { return _localDoF; };
    // const vector<int > & getGhostDoF() { return _ghostDoF; };

  protected:
    //! State where nodal positions and dof are stored
    State * _gState;

    //! List of Elements
    vector<GeomElement* > _elements;

    //! Element type
    string _elementType;

    //! Map between local nodes and global nodes: needed to fill in State and elements from connTable
    vector<int > _LocaltoGlobal; // Local nodal i to global node k - k = _LocaltoGlobal[i];

    // //! Map of Degrees of freedom local ID to global ID
    // vector<int >          _localDoF;
    // vector<int >          _ghostDoF;

  }; // Class Mesh

} //namespace voom

#endif // __Mesh_h__
