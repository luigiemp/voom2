#include "HalfEdgeMesh.h"

using namespace voom;

int main(){

  HalfEdgeMesh::ConnectivityContainer connect(2);
  //
  //
  // 0---3
  // | \ |
  // 1---2
  connect[0] << 0, 1, 2;
  connect[1] << 2, 3, 0;
  int nverts=4;

  HalfEdgeMesh mesh(connect,nverts);

  std::cout << "vertices:"<< std::endl;
  for(int a=0; a<mesh.vertices.size(); a++) {
    std::cout << mesh.vertices[a]->id << std::endl;
  }

  std::cout << "HalfEdges:"<< std::endl;
  for(int h=0; h<mesh.halfEdges.size(); h++) {
    std::cout << "h: "<< mesh.halfEdges[h]->id << std::endl
  	      << "v: "<< (mesh.halfEdges[h]->vertex)->id << std::endl
  	      << "n: "<< (mesh.halfEdges[h]->next)->id << std::endl;
    if(mesh.halfEdges[h]->opposite) 
      std::cout << "o: "<< (mesh.halfEdges[h]->opposite)->id << std::endl;
    else
      std::cout << "o: "<< -1 << std::endl;
    std::cout << "f: "<< (mesh.halfEdges[h]->face)->id << std::endl
  	      << std::endl;
      
  }

}
