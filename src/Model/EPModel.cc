#include "EPModel.h"

namespace voom {
  EPModel::EPModel(Mesh* myMesh, const string inputFile, const uint NodeDof)
    : Model(myMesh, inputFile, NodeDof) {
    // Read diffusion tensor data
    ifstream inp(inputFile.c_str() );
    string line;
    _diffusion.resize( _myMesh->getNumberOfElements() );
    while( getline(inp, line) ) {
      // Read diffusion tensor data
      if (line.find("*DIFFUSION") == 0) 
	for(uint i = 0; i < _diffusion.size(); i++) {
	  getline(inp, line);
	  vector<string> strs = splitString(line, " \t");
	  for(uint m = 0; m < 3; m++)
	    for(uint n = 0; n < 3; n++) 
	      _diffusion[i](m,n) = atof( strs[3*m + n].c_str() );
	} // Loop over number of diffusion values

      // Reading in of load steps
      
    } // While getline loop
    inp.close();
  } // Constructor

}
