//-*-C++-*-

#ifndef __FEMesh_h__
#define __FEMesh_h__

#include "Mesh.h"
#include "HexQuadrature.h"
#include "LineQuadrature.h"
#include "QuadQuadrature.h"

#include "HexShape.h"
#include "BarShape.h"

namespace voom{
  
  class FEMesh: public Mesh {
  
  public:
    //! Constructor from input file
    FEMesh(const string inputFile);

    //! Constructor from nodes and connectivities
    //! (assume only one type of elements and quadrature in the mesh)
    FEMesh(const vector<VectorXd > &  Positions,
	   const vector<vector<int > > & Connectivity, 
	   const vector<int > & LocalDoF,
	   const vector<int > & GhostDoF,
	   string ElementType,
	   uint QuadOrder);

    //! Destructor
    // Destructor made virtual due to inheritance
    virtual ~FEMesh() {
      for(uint i = 0; i < _shapes.size(); i++)
	for(uint j = 0; j < _shapes[i].size(); j++)
	  delete _shapes[i][j];

      for(uint i = 0; i < _quadrature.size(); i++) 
	delete _quadrature[i];
    }
  
  protected:
    //! List of shape objects 
    // One vector of shape functions per each element type
    vector<vector<Shape*> > _shapes;
    
    //! One Quadrature rule per each element type
    vector<Quadrature*>     _quadrature;
  };
}

#endif // __FEMesh_h__
