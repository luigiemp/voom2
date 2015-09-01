//-*-C++-*-

#ifndef __LoopShellMesh_h__
#define __LoopShellMesh_h__

#include "Mesh.h"
#include "TriQuadrature.h"
#include "LoopShellShape.h"

namespace voom{
  
  class LoopShellMesh: public Mesh {
  
  public:
    
    typedef Eigen::Vector3i TriangleConnectivity;
    typedef std::vector<int> CornerValences; //Valences of nodes at three corners

    //! Constructor from input file
    LoopShellMesh(const string Nodes, const string ConnTable);
    //! Constructor from nodes and connectivities
    //! (assume only one type of elements and quadrature in one mesh)
    // LoopShellMesh(const vector<VectorXd > &  Positions,
    // 	   const vector<vector<int > > & Connectivity, 
    // 	   string ElementType,
    // 	   uint QuadOrder);

    // Build a mesh using the same nodes that are already used by another mesh
    // LoopShellMesh(LoopShellMesh*, const string ConnTable);

    //! Destructor
    // Destructor made virtual due to inheritance
    virtual ~LoopShellMesh() {
      for(map < CornerValences, vector<LoopShellShape*> >::iterator it = _shapes.begin();
	  it != _shapes.end(); ++it){
	//loop over map
	for(uint j=0; j < (it->second).size(); j++){
	  //loop over vector
	  delete (it->second)[j];
	}
      }

      delete _quadrature;
    }
    void initializeMesh(vector<TriangleConnectivity> & connectivities);
    
  protected:
    //! Map of shape objects The integer key is the Valence
    //! corresponding to a regular or an irregular patch. The value
    //! stored is a vector of type Shape*. Each entry of the vector is
    //! the shape functions at a quadrature point
    map< CornerValences ,vector<LoopShellShape*> > _shapes;
    
    //! One Quadrature rule per each element type
    Quadrature*    _quadrature;
  };
}

#endif // __LoopShellMesh_h__
