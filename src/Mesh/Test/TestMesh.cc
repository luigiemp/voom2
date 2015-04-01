#include "FEMesh.h"

using namespace voom;

int main(int argc, char** argv)
{
  cout << endl << "Testing voom mesh class ... " << endl;
  cout << "................................... " << endl << endl;
  
  // Test FEmesh constructor from nodes and conn table
  {
    cout << endl << "Test FEmesh constructor from nodes and conn table files " << endl;

    FEMesh TestFEmesh("NodeFile.dat", "ElFile.dat");
    cout << "Number Of Nodes   : " << TestFEmesh.getNumberOfNodes() << endl;
    cout << "Number Of Element : " << TestFEmesh.getNumberOfElements() << endl;
    cout << "Mesh Dimension    : " << TestFEmesh.getDimension() << endl;
    cout << endl << "Nodal position" << endl;
    for(uint i = 0; i < TestFEmesh.getNumberOfNodes(); i++) {
      for(uint m = 0 ; m < TestFEmesh.getDimension(); m++)
	cout << setw(10) << TestFEmesh.getX(i, m) << " " ;
      cout << endl;
    }   
    cout << endl << "Conn table" << endl;
    vector<GeomElement* > Els = TestFEmesh.getElements();
    for(uint e = 0; e < TestFEmesh.getNumberOfElements(); e++) {
      vector<int > NodesID = Els[e]->getNodesID();
      for(uint n = 0 ; n < NodesID.size(); n++)
	cout << NodesID[n] << " " ;
      cout << endl;
    }   
    
    cout << endl << "END of Test FEmesh constructor from nodes and conn table files " << endl; 
  }



  // Test FEmesh constructor from nodes and conn table
  {
    cout << endl << "Test FEmesh constructor from nodes and conn table files " << endl;

    FEMesh TestFEmesh("Cube.node", "Cube.ele");
    cout << "Number Of Nodes   : " << TestFEmesh.getNumberOfNodes() << endl;
    cout << "Number Of Element : " << TestFEmesh.getNumberOfElements() << endl;
    cout << "Mesh Dimension    : " << TestFEmesh.getDimension() << endl;
    cout << endl << "Nodal position" << endl;
    for(uint i = 0; i < TestFEmesh.getNumberOfNodes(); i++) {
      for(uint m = 0 ; m < TestFEmesh.getDimension(); m++)
	cout << setw(10) << TestFEmesh.getX(i, m) << " " ;
      cout << endl;
    }   
    cout << endl << "Conn table" << endl;
    vector<GeomElement* > Els = TestFEmesh.getElements();
    for(uint e = 0; e < TestFEmesh.getNumberOfElements(); e++) {
      vector<int > NodesID = Els[e]->getNodesID();
      for(uint n = 0 ; n < NodesID.size(); n++)
	cout << NodesID[n] << " " ;
      cout << endl;
    }   
    
    cout << endl << "END of Test FEmesh constructor from nodes and conn table files " << endl; 
  }



  // Test FEmesh constructor from nodes and conn table - Lin tet cube
  {
    cout << endl << "Test FEmesh constructor from nodes and conn table files " << endl;

    FEMesh TestFEmesh("Cube.node", "SurfCube.ele");
    cout << "Number Of Nodes   : " << TestFEmesh.getNumberOfNodes() << endl;
    cout << "Number Of Element : " << TestFEmesh.getNumberOfElements() << endl;
    cout << "Mesh Dimension    : " << TestFEmesh.getDimension() << endl;
    cout << endl << "Nodal position" << endl;
    for(uint i = 0; i < TestFEmesh.getNumberOfNodes(); i++) {
      for(uint m = 0 ; m < TestFEmesh.getDimension(); m++)
	cout << setw(10) << TestFEmesh.getX(i, m) << " " ;
      cout << endl;
    }   
    cout << endl << "Conn table" << endl;
    vector<GeomElement* > Els = TestFEmesh.getElements();
    for(uint e = 0; e < TestFEmesh.getNumberOfElements(); e++) {
      vector<int > NodesID = Els[e]->getNodesID();
      for(uint n = 0 ; n < NodesID.size(); n++)
	cout << NodesID[n] << " " ;
      cout << endl;
    }   
    
    cout << endl << "END of Test FEmesh constructor from nodes and conn table files " << endl; 
  }



  // Test FEmesh constructor from nodes and conn table
  {
    cout << endl << "Test FEmesh constructor from nodes and conn table files " << endl;

    FEMesh TestFEmesh("CubeQuad.node", "CubeQuad.ele");
    cout << "Number Of Nodes   : " << TestFEmesh.getNumberOfNodes() << endl;
    cout << "Number Of Element : " << TestFEmesh.getNumberOfElements() << endl;
    cout << "Mesh Dimension    : " << TestFEmesh.getDimension() << endl;
    cout << endl << "Nodal position" << endl;
    for(uint i = 0; i < TestFEmesh.getNumberOfNodes(); i++) {
      for(uint m = 0 ; m < TestFEmesh.getDimension(); m++)
	cout << setw(10) << TestFEmesh.getX(i, m) << " " ;
      cout << endl;
    }   
    cout << endl << "Conn table" << endl;
    vector<GeomElement* > Els = TestFEmesh.getElements();
    for(uint e = 0; e < TestFEmesh.getNumberOfElements(); e++) {
      vector<int > NodesID = Els[e]->getNodesID();
      for(uint n = 0 ; n < NodesID.size(); n++)
	cout << NodesID[n] << " " ;
      cout << endl;
    }   
    
    cout << endl << "END of Test FEmesh constructor from nodes and conn table files " << endl; 
  }



  // Test FEmesh constructor from nodes and conn table - Quad cube
  {
    cout << endl << "Test FEmesh constructor from nodes and conn table files " << endl;

    FEMesh TestFEmesh("CubeQuad.node", "SurfCubeQuad.ele");
    cout << "Number Of Nodes   : " << TestFEmesh.getNumberOfNodes() << endl;
    cout << "Number Of Element : " << TestFEmesh.getNumberOfElements() << endl;
    cout << "Mesh Dimension    : " << TestFEmesh.getDimension() << endl;
    cout << endl << "Nodal position" << endl;
    for(uint i = 0; i < TestFEmesh.getNumberOfNodes(); i++) {
      for(uint m = 0 ; m < TestFEmesh.getDimension(); m++)
	cout << setw(10) << TestFEmesh.getX(i, m) << " " ;
      cout << endl;
    }   
    cout << endl << "Conn table" << endl;
    vector<GeomElement* > Els = TestFEmesh.getElements();
    for(uint e = 0; e < TestFEmesh.getNumberOfElements(); e++) {
      vector<int > NodesID = Els[e]->getNodesID();
      for(uint n = 0 ; n < NodesID.size(); n++)
	cout << NodesID[n] << " " ;
      cout << endl;
    }   
    
    cout << endl << "END of Test FEmesh constructor from nodes and conn table files " << endl; 
  }



  // Test FEmesh constructor from nodes and conn table - Coarse LV mesh
  {
    cout << endl << "CoarseLV " << endl;

    FEMesh TestFEmesh("CoarseLV.node", "CoarseLV.ele");
    cout << "Number Of Nodes   : " << TestFEmesh.getNumberOfNodes() << endl;
    cout << "Number Of Element : " << TestFEmesh.getNumberOfElements() << endl;
    cout << "Mesh Dimension    : " << TestFEmesh.getDimension() << endl;
    // cout << endl << "Nodal position" << endl;
    // for(uint i = 0; i < TestFEmesh.getNumberOfNodes(); i++) {
    //   for(uint m = 0 ; m < TestFEmesh.getDimension(); m++)
    // 	cout << setw(10) << TestFEmesh.getX(i, m) << " " ;
    //   cout << endl;
    // }   
    // cout << endl << "Conn table" << endl;
    // vector<GeomElement* > Els = TestFEmesh.getElements();
    // for(uint e = 0; e < TestFEmesh.getNumberOfElements(); e++) {
    //   vector<int > NodesID = Els[e]->getNodesID();
    //   for(uint n = 0 ; n < NodesID.size(); n++)
    // 	cout << NodesID[n] << " " ;
    //   cout << endl;
    // }   
    
    cout << endl << "END of CoarseLV test " << endl; 
  }



 
  cout << endl << "........................ " << endl;
  cout << "Test of voom mesh class completed" << endl;
  return 0;
}
