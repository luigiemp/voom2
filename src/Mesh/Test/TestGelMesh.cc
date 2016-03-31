#include "GelMesh.h"
        
using namespace voom;
 
int main(int argc, char** argv)
{
  cout << endl << "Testing voom mesh class ... " << endl;
  cout << "................................... " << endl << endl;
  {
    cout << endl << "Test Gelmesh constructor for filament from nodes and conn table files " << endl;

    GelMesh TestGelmesh("Filament2.node", "Filament2.ele");
    cout << "Number Of Nodes   : " << TestGelmesh.getNumberOfNodes() << endl;
    cout << "Number Of Filament : " << TestGelmesh.getNumberOfFilaments() << endl;
    cout << "Mesh Dimension    : " << TestGelmesh.getDimension() << endl;
    cout << endl << "Nodal position" << endl;
    for(uint i = 0; i < TestGelmesh.getNumberOfNodes(); i++) {
      for(uint m = 0 ; m < TestGelmesh.getDimension(); m++)
        cout << setw(10) << TestGelmesh.getX(i, m) << " " ;
      cout << endl;
    }
    cout << endl << "Conn table" << endl;
    vector<GeomFilament* > Els = TestGelmesh.getFilaments();
    for(uint e = 0; e < TestGelmesh.getNumberOfFilaments(); e++) {
      vector<int > NodesID = Els[e]->getNodesID();
      for(uint n = 0 ; n < NodesID.size(); n++)
        cout << NodesID[n] << " " ;
      cout << endl;
    }

    cout << endl << "END of Test GelMesh constructor for filament from nodes and conn table files " << endl;
  }


 
  cout << endl << "........................ " << endl;
  cout << "Test of voom mesh class completed" << endl;
  return 0;
}
