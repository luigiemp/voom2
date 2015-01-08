#include "FEMesh.h"

using namespace voom;

int main(int argc, char** argv)
{
  cout << endl << "Testing voom mesh class ... " << endl;
  cout << "................................... " << endl << endl;
  
  // Test FEmesh constructor from nodes and conn table
  {
    cout << endl << "Test FEmesh constructor from nodes and conn table " << endl;

    vector<VectorXd > Positions(12, Vector3d::Zero());
  
    Positions[0] <<-1.0,-1.0,-1.0; 
    Positions[1] << 1.0,-1.0,-1.0; 
    Positions[2] << 1.0, 1.0,-1.0;   
    Positions[3] <<-1.0, 1.0,-1.0;
    Positions[4] <<-1.0,-1.0, 1.0;   
    Positions[5] << 1.0,-1.0, 1.0;
    Positions[6] << 1.0, 1.0, 1.0;
    Positions[7] <<-1.0, 1.0, 1.0;
    Positions[8] <<-1.0,-1.0, 2.0;   
    Positions[9] << 1.0,-1.0, 2.0;   
    Positions[10]<< 1.0, 1.0, 2.0;
    Positions[11]<<-1.0, 1.0, 2.0;
    for (uint i = 0; i < Positions.size(); i++)
      Positions[i] += (Vector3d::Random()*0.1);

    vector<vector<int > > Connectivity(2, vector<int >(8,0));
    Connectivity[0][0] = 0;  Connectivity[0][1] = 1;  Connectivity[0][2] = 2;  Connectivity[0][3] = 3;
    Connectivity[0][4] = 4;  Connectivity[0][5] = 5;  Connectivity[0][6] = 6;  Connectivity[0][7] = 7;
    
    Connectivity[1][0] = 4;  Connectivity[1][1] = 5;  Connectivity[1][2] = 6;  Connectivity[1][3] = 7;
    Connectivity[1][4] = 8;  Connectivity[1][5] = 9;  Connectivity[1][6] = 10; Connectivity[1][7] = 11;
    
    vector<int > LocalDoF(36, 0); // DoF - not nodal - mapping
    vector<int > GhostDoF;
    for (uint i = 0; i < LocalDoF.size(); i++)
      LocalDoF[i] = i;
    
    string ElementType = "Hexa";
    
    uint QuadOrder = 2;
    
    FEMesh TestFEmesh(Positions, Connectivity, LocalDoF, GhostDoF, ElementType, QuadOrder);
    cout << "Number Of Nodes   : " << TestFEmesh.getNumberOfNodes() << endl;
    cout << "Number Of Element : " << TestFEmesh.getNumberOfElements() << endl;
    cout << "Mesh Dimension    : " << TestFEmesh.getDimension() << endl;
    cout << endl << "Nodal position" << endl;
    for(uint i = 0; i < TestFEmesh.getNumberOfNodes(); i++) {
      for(uint m = 0 ; m < TestFEmesh.getDimension(); m++)
	cout << setw(10) << TestFEmesh.getX(i, m) << " " ;
      cout << endl;
    }   
    
    cout << endl << "END of Test FEmesh constructor from nodes and conn table " << endl; 
  }
  /*
  // Test Mesh constructor from input file
  {
    cout << endl << "Test Mesh constructor from input file " << endl;
    Mesh* myMesh = Mesh::New( string( argv[1] ) );
    cout << "Number Of Nodes   : " << myMesh->getNumberOfNodes() << endl;
    cout << "Number Of Element : " << myMesh->getNumberOfElements() << endl;
    cout << "Mesh Dimension    : " << myMesh->getDimension() << endl;
    cout << endl << "Nodal position" << endl;
    for(uint i = 0; i < myMesh->getNumberOfNodes(); i++) {
      for(uint m = 0 ; m < myMesh->getDimension(); m++)
	cout << setw(10) << myMesh->getX(i, m) << " " ;
      cout << endl;
    }    
    delete myMesh; 
    cout << endl << "END of Test Mesh constructor from input file " << endl;
  }
  */
  cout << endl << "........................ " << endl;
  cout << "Test of voom mesh class completed" << endl;
  return 0;
}
