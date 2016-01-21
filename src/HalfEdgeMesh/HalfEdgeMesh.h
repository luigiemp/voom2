// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                    (C) 2007 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__HEMesh_h__)
#define __HEMesh_h__

#include <vector>
#include "voom.h"

namespace voom {

  struct HalfEdge;
  struct Vertex;
  struct Face;

  //! A HalfEdge is the Left half of an edge, pointing CCW around a Face.
  struct HalfEdge {
    //! Default constructor
    HalfEdge() {
      opposite = next = prev = 0;
      face = 0; 
      vertex = 0;
      id = -1;
    }

    //! Adjacent HalfEdge
    HalfEdge * opposite;
    //! HalfEdge Counter-Clockwise around Face on left
    HalfEdge * next;
    //! HalfEdge Clockwise around Face on left (stored for convenience)
    HalfEdge * prev;
    //! Face on left
    Face * face;
    //! vertex at the end of the HalfEdge
    Vertex * vertex;

    //! HalfEdge id
    int id;
  };
  
  struct Face {
    //! Default constructor
    Face() { halfEdges[0] = halfEdges[1] = halfEdges[2] = 0; id=-1; }

    //! HalfEdges looping Counter-Clockwise around Face
    HalfEdge * halfEdges[3];

    //! Face id
    int id;
  };

  struct Vertex {
    //! Default constructor
    Vertex() { id = -1; boundary = false; }

    //! vertex id
    int id;

    //! is this a boundary vertex?
    bool boundary;

    //! list of incident HalfEdges (pointing towards the vertex)
    std::vector< HalfEdge* > halfEdges;
  };


  //! HalfEdge data structure to store a triangle mesh
  /*!  This class is built to convert a connectivity array into a
    HalfEdge data structure.
  */
  struct HalfEdgeMesh {

    std::vector< HalfEdge* > halfEdges;
    std::vector< Face* > faces;
    std::vector< Vertex* > vertices;

    //! Connectivity of a single triangle
    //typedef tvmet::Vector<int,3> TriangleConnectivity;
    typedef Eigen::Vector3i TriangleConnectivity;

    //! Container of connectivities for a mesh of triangles
    typedef std::vector<TriangleConnectivity> ConnectivityContainer;

    // construct HalfEdge data structure from a connectivity array
    HalfEdgeMesh(const ConnectivityContainer & connectivities, int nVertices);

    // Delete structs
    virtual ~HalfEdgeMesh() {
      for(int f=0; f<faces.size(); f++) delete faces[f];
      for(int h=0; h<halfEdges.size(); h++) delete halfEdges[h];
      for(int v=0; v<vertices.size(); v++) delete vertices[v];
    }

  };

  
}
#endif // __HEMesh_h__
