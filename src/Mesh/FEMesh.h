//-*-C++-*-

#ifndef __FEMesh_h__
#define __FEMesh_h__

#include "Mesh.h"
#include "HexQuadrature.h"
#include "TetQuadrature.h"
#include "LineQuadrature.h"
#include "QuadQuadrature.h"
#include "TriQuadrature.h"

#include "HexShape.h"
#include "BarShape.h"
#include "LinTriShape.h"
#include "QuadTriShape.h"
#include "LinTetShape.h"
#include "QuadTetShape.h"

namespace voom{

  class FEMesh: public Mesh {

  public:
    //! Constructor from input file
    FEMesh(const string Nodes, const string ConnTable);

    //! Constructor from nodes and connectivities
    //! (assume only one type of elements and quadrature in one mesh)
    // FEMesh(const vector<VectorXd > &  Positions,
    // 	   const vector<vector<int > > & Connectivity,
    // 	   string ElementType,
    // 	   uint QuadOrder);

    // Build a mesh using the same nodes that are already used by another mesh
    // FEMesh(FEMesh*, const string ConnTable);

    //! Destructor
    // Destructor made virtual due to inheritance
    virtual ~FEMesh() {
      for(uint i = 0; i < _shapes.size(); i++)
	delete _shapes[i];

      delete _quadrature;
    }

  protected:
    //! List of shape objects
    // One vector of shape functions per each element type
    vector<Shape*> _shapes;

    //! One Quadrature rule per each element type
    Quadrature*    _quadrature;

    //! Helper function to determine type of element and fills in
    //! \param _shapes and \param _quadrature and returns \return NumNodesEl
    int createElementShapeAndQuadrature(const string ElType);
  };
}

#endif // __FEMesh_h__
