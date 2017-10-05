#include "FEMesh.h"
// #include "LoopShellMesh.h"

using namespace voom;

int main()
{
  
  cout << endl << "Testing voom mesh class ... " << endl;
  cout << "................................... " << endl << endl;
  
  // Test FEmesh constructor from nodes and conn table
  {
    cout << endl << "Test FEmesh constructor from nodes and conn table files " << endl;

    State myState;
    bool CheckOverlap = false;
    int dofPerNode = 3;
    FEMesh TestFEmesh("../../../../voom2/src/Mesh/Test/NodeFile.dat", "../../../../voom2/src/Mesh/Test/ElFile.dat", &myState, dofPerNode, CheckOverlap);
    cout << "Number Of Element : " << TestFEmesh.getNumberOfElements() << endl;
    cout << endl << "Nodal position" << endl;
    for(int i = 0; i < myState.getXsize(); i++) {
      cout << myState.getX(i)(0) << " " << myState.getX(i)(1) << " " << myState.getX(i)(2) << " " ;
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

    State myState;
    bool CheckOverlap = false;
    int dofPerNode = 3;
      FEMesh TestFEmesh("../../../../voom2/src/Mesh/Test/Cube.node", "../../../../voom2/src/Mesh/Test/Cube.ele", &myState, dofPerNode, CheckOverlap);
    cout << "Number Of Element : " << TestFEmesh.getNumberOfElements() << endl;
    cout << endl << "Nodal position" << endl;
    for(int i = 0; i < myState.getXsize(); i++) {
      cout << myState.getX(i)(0) << " " << myState.getX(i)(1) << " " << myState.getX(i)(2) << " " ;
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

    State myState;
    bool CheckOverlap = false;
    int dofPerNode = 3;
    FEMesh TestFEmesh("../../../../voom2/src/Mesh/Test/Cube.node", "../../../../voom2/src/Mesh/Test/SurfCube.ele", &myState, dofPerNode, CheckOverlap);
    cout << "Number Of Element : " << TestFEmesh.getNumberOfElements() << endl;
    cout << endl << "Nodal position" << endl;
    for(int i = 0; i < myState.getXsize(); i++) {
      cout << myState.getX(i)(0) << " " << myState.getX(i)(1) << " " << myState.getX(i)(2) << " " ;
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

    State myState;
    bool CheckOverlap = false;
    int dofPerNode = 3;
    FEMesh TestFEmesh("../../../../voom2/src/Mesh/Test/CubeQuad.node", "../../../../voom2/src/Mesh/Test/CubeQuad.ele", &myState, dofPerNode, CheckOverlap);
    cout << "Number Of Element : " << TestFEmesh.getNumberOfElements() << endl;
    cout << endl << "Nodal position" << endl;
    for(int i = 0; i < myState.getXsize(); i++) {
      cout << myState.getX(i)(0) << " " << myState.getX(i)(1) << " " << myState.getX(i)(2) << " " ;
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

    State myState;
    bool CheckOverlap = false;
    int dofPerNode = 3;
    FEMesh TestFEmesh("../../../../voom2/src/Mesh/Test/CubeQuad.node", "../../../../voom2/src/Mesh/Test/SurfCubeQuad.ele", &myState, dofPerNode, CheckOverlap);
    cout << "Number Of Element : " << TestFEmesh.getNumberOfElements() << endl;
    cout << endl << "Nodal position" << endl;
    for(int i = 0; i < myState.getXsize(); i++) {
      cout << myState.getX(i)(0) << " " << myState.getX(i)(1) << " " << myState.getX(i)(2) << " " ;
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

    State myState;
    bool CheckOverlap = false;
    int dofPerNode = 3;
    FEMesh TestFEmesh("../../../../voom2/src/Mesh/Test/CoarseLV.node", "../../../../voom2/src/Mesh/Test/CoarseLV.ele", &myState, dofPerNode, CheckOverlap);
    cout << "Number Of Element : " << TestFEmesh.getNumberOfElements() << endl;
    cout << endl << "Nodal position" << endl;
    for(int i = 0; i < myState.getXsize(); i++) {
      cout << myState.getX(i)(0) << " " << myState.getX(i)(1) << " " << myState.getX(i)(2) << " " ;
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
    
    cout << endl << "END of CoarseLV test " << endl; 
  }



  // Test FEmesh constructor with shared nodes
  {
    cout << endl << "Double cube " << endl;

    State myState;
    bool CheckOverlap = true;
    int dofPerNode = 3;
    FEMesh CubeMeshOne("../../../../voom2/src/Mesh/Test/CubeQuad.node", "../../../../voom2/src/Mesh/Test/CubeQuad.ele", &myState, dofPerNode, CheckOverlap);
    CheckOverlap = true;
    FEMesh CubeMeshTwo("../../../../voom2/src/Mesh/Test/CubeQuadTZ.node", "../../../../voom2/src/Mesh/Test/CubeQuad.ele", &myState, dofPerNode, CheckOverlap);

    cout << "Number Of Element : " << CubeMeshTwo.getNumberOfElements() << endl;
    cout << endl << "Nodal position" << endl;
    for(int i = 0; i < myState.getXsize(); i++) {
      cout << myState.getX(i)(0) << " " << myState.getX(i)(1) << " " << myState.getX(i)(2) << " " ;
      cout << endl;
    } 
    cout << endl << "Conn table" << endl;
    vector<GeomElement* > Els = CubeMeshTwo.getElements();
    for(uint e = 0; e < CubeMeshTwo.getNumberOfElements(); e++) {
      vector<int > NodesID = Els[e]->getNodesID();
      for(uint n = 0 ; n < NodesID.size(); n++)
  	cout << NodesID[n] << " " ;
      cout << endl;
    }     
    
    cout << endl << "END of Double Cube Test " << endl; 
  }



  // Test FEmesh constructor with linear hexa mesh plus shared nodes
  {
    cout << endl << "Double hexa cube " << endl;

    State myState;
    bool CheckOverlap = true;
    int dofPerNode = 3;
    FEMesh CubeMeshOne("../../../../voom2/src/Mesh/Test/CubeHexa.node", "../../../../voom2/src/Mesh/Test/CubeHexa.ele", &myState, dofPerNode, CheckOverlap);
    CheckOverlap = true;
    FEMesh CubeMeshTwo("../../../../voom2/src/Mesh/Test/CubeHexa.node", "../../../../voom2/src/Mesh/Test/CubeHexa.ele", &myState, dofPerNode, CheckOverlap);

    cout << "Number Of Element : " << CubeMeshTwo.getNumberOfElements() << endl;
    cout << endl << "Nodal position" << endl;
    for(int i = 0; i < myState.getXsize(); i++) {
      cout << myState.getX(i)(0) << " " << myState.getX(i)(1) << " " << myState.getX(i)(2) << " " ;
      cout << endl;
    } 
    cout << endl << "Conn table" << endl;
    vector<GeomElement* > Els = CubeMeshTwo.getElements();
    for(uint e = 0; e < CubeMeshTwo.getNumberOfElements(); e++) {
      vector<int > NodesID = Els[e]->getNodesID();
      for(uint n = 0 ; n < NodesID.size(); n++)
  	cout << NodesID[n] << " " ;
      cout << endl;
    }     
    
    cout << endl << "END of Double Cube Test " << endl; 
  }



  // LoopShellMesh icosa("T7nodes.dat","T7connectivity.dat");
  


  cout << endl << "........................ " << endl;
  cout << "Test of voom mesh class completed" << endl;
  
}
