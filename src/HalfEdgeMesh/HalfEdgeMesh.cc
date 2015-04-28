// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                    (C) 2007 All Rights Reserved
//
//----------------------------------------------------------------------

#include "HalfEdgeMesh.h"

namespace voom {

HalfEdgeMesh::HalfEdgeMesh(const ConnectivityContainer & connectivities, 
			   int nVertices)
{

  std::cout << "Building HalfEdgeMesh..." << std::endl;

  // size containers
  faces.resize( connectivities.size() );
  for(int f=0; f<faces.size(); f++) faces[f] = new Face();
  halfEdges.resize( 3*connectivities.size() );
  for(int h=0; h<halfEdges.size(); h++) halfEdges[h] = new HalfEdge();
  vertices.resize(nVertices);
  for(int v=0; v<vertices.size(); v++) vertices[v] = new Vertex();

  // set vertex ids      
  for(int i=0; i<nVertices; i++) vertices[i]->id = i;

  // Link HalfEdges, Vertices, and Faces for each element; 
  // opposites remain undefined
  for(int f=0; f<connectivities.size(); f++) {
    Face * F = faces[f];
    F->id = f;
    for(int i=0; i<3; i++) {
      HalfEdge * HE = halfEdges[3*f+i];
      HE->id = 3*f+i;
      int id=connectivities[f](i);
      Vertex * V = vertices[id];

      // set vertex from connectivities
      HE->vertex = V;
      V->halfEdges.push_back( HE );

      // link HalfEdges to Face
      HE->face = F;
      F->halfEdges[i] = HE;

      // identify next HalfEdge (i+1 mod 3)
      int j = (i==2) ? 0 : i+1;
      HE->next = halfEdges[3*f+j];

      // identify prev HalfEdge (i-1 mod 3)
      j = (i==0) ? 2 : i-1;
      HE->prev = halfEdges[3*f+j];
    }
  }
  std::cout << "Linked " 
	    << halfEdges.size() << " half-edges, "
	    << faces.size() << " faces, and "
	    << vertices.size() << " vertices. "
	    << std::endl;

  std::cout << "Determining half-edge opposites." << std::endl;
  // For each  HalfEdge, find its opposite
  for(int h=0; h<halfEdges.size(); h++) {
    HalfEdge * H = halfEdges[h];
    // check if next is already set
    if( H->opposite != 0 ) continue;

    // if not, then look at all HalfEdges incident to my vertex and
    // check to see if they duplicate me (i.e., orientations are not
    // consistent) or if their next is my opposite

    for(int hv=0; hv<H->vertex->halfEdges.size(); hv++) {
      HalfEdge * Hv = H->vertex->halfEdges[hv];

      if( H == Hv ) continue;
    
      // Hv is a duplicate of H if both previous HalfEdges point to
      // the same vertex
      if( H->prev->vertex->id == Hv->prev->vertex->id ) {
	std::cout << "Error:  half-edges " << H->id << " and " << Hv->id
		  << " are duplicates.  " << std::endl
		  << "Adjacent faces " << H->face->id
		  << " and " << Hv->face->id 
		  << " must have opposite orientations." << std::endl;
	exit(0);
      }

      // Hv's next is my opposite if it points to the vertex of
      // my previous HalfEdge
      if( H->prev->vertex->id == Hv->next->vertex->id ) {
	H->opposite = Hv->next;
	break;
      }
    }

    if(H->opposite == 0 ) {
      std::cout << "  Warning."
		<< "  Couldn't find HalfEdge opposite to HalfEdge " << H->id 
		<< ".  This HalfEdge is on the boundary."
		<< std::endl;

      H->vertex->boundary = true;

//       for(int hv=0; hv<H->vertex->halfEdges.size(); hv++) {
// 	HalfEdge * Hv = H->vertex->halfEdges[hv];
// 	std::cout << "h: "<< Hv->id << std::endl
// 		  << "v: "<< Hv->vertex->id << std::endl
// 		  << "n: "<< Hv->next->id << std::endl
// 		  << "nv: "<< Hv->next->vertex->id << std::endl
// 		  << "p: "<< Hv->prev->id << std::endl
// 		  << "pv: "<< Hv->prev->vertex->id << std::endl;
// 	std::cout << std::endl;
//       }
//      exit(0);
    }
	
  }

  std::cout << "HalfEdgeMesh built." << std::endl;
	  
}

}
