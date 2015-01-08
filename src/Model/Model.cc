//-*-C++-*-
#include "Model.h"
// #include "PoissonModel.h"
// #include "MechanicsModel.h"

namespace voom {
  // Constructor
  // Model::Model(Mesh* myMesh, const string inputFile, const uint NodeDoF):
  //   _myMesh(myMesh), _inputFile(inputFile), _nodeDoF(NodeDoF) {
  // };



  // Static function which returns derived class pointer
  // Model* Model::New(Mesh* myMesh, const string inputFile, const uint NodeDoF)
  // { 
  //   // The static constructor only initialize the analysis type
  //   Model* myModel;
    // // Reading control deck data
    // ifstream inp(inputFile.c_str());
    // vector<string > analysisTypes;
    // string line;

    // while ( inp.good() ) {
    //   getline(inp, line);
    //   if ( line.find("*STEP") == 0 ) {
    // 	while( getline(inp, line) ) { // Reading step info
    // 	  if ( line.find( "**" ) == 0) continue;
    // 	  else {
    // 	    analysisTypes.push_back(line); break;
    // 	  }
    // 	}
    //   } // End of find *STEP
    // } // In ctrl deck file is found
    // inp.close();

    // // Check that analysis type has been declared
    // if ( analysisTypes.size() == 0 ) {
    //   cerr << "**ERROR: No analysis defined in model" << endl;
    //   cerr << "Exiting....\n";
    //   exit(EXIT_FAILURE);
    // }
    
    // // Define type of model to be created
    // if ( analysisTypes[0].find("*HEAT") == 0) {
    //   vector<string > strs = splitString( analysisTypes[0], ",");
    //   if (strs[1].find("STEADY STATE") != string::npos ) {
    // 	myModel = new PoissonModel(myMesh, inputFile, NodeDoF);
    //   }
    // } else if ( analysisTypes[0].find("STATIC")== 0) {
    //   myModel = new MechanicsModel(myMesh, inputFile, NodeDoF);
    //   /*vector<string > strs = splitString( analysisTypes[0], ",");
    //   if (strs[1].find("LINEARIZED KINEMATICS") != string::npos ) {
    // 	myModel->setKinematics(0);
    //   }
    //   else if (strs[1].find("FINITE KINEMATICS") != string::npos ) {
    // 	myModel->setKinematics(1);
    // 	}*/
    // } else {
    //   cerr << "**ERROR:  Analysis type " << analysisTypes[0]
    // 	   << " not implemented"<< endl;
    //   cerr << "Exiting...\n";
    //   exit(EXIT_FAILURE);
  //   // }
  //   return myModel;
  // } // Static Function class



  // // Parse input file to read list of material per each element
  // void Model::_getMaterialNamesList(vector<string > & names)
  // {
  //   ifstream inp( _inputFile.c_str() );
  //   string line;
  //   names.resize( _myMesh->getNumberOfElements() );
  //   // Read file until find the element data block
  //   while( getline(inp, line) ) 
  //     if ( line.find("*ELEMENT") == 0) break;
  //   // Fill in a list with one material per each element
  //   for(uint i = 0; i < names.size(); i++) {
  //     getline(inp, line);
  //     vector<string > strs = splitString(line, " \t");
  //     names[i] = strs.back();
  //   }
  //   inp.close();
  // } // end of _getMaterialNamesList


  /*
  // Parse input deck to get local to global and ghost to global ID map
  void Model::_getNodalMapping(vector<int > & localToGlobalID,
			       vector<int > & ghostToGlobalID) 
  {
    ifstream inp( _inputFile.c_str() );
    string line;

    // Read number of local nodes data
    while( getline(inp, line) )
      if ( line.find("*NUMLOCALNODES") == 0 ) break;

    getline(inp, line);
    const int numLocalNodes = atoi(line.c_str());
    
    // Read local to global mapping data
    while( getline(inp, line) )
      if ( line.find("*GLOBALID") == 0 ) break;
    
    localToGlobalID.resize(numLocalNodes);
    ghostToGlobalID.resize(_myMesh->getNumberOfNodes() - numLocalNodes);
    // First read Local to Global nodal mapping
    for(uint i = 0; i < numLocalNodes; i++) {
      getline(inp, line);
      localToGlobalID[i] = atoi(line.c_str());
    }
    // second read Ghost to Global nodal mapping
    for(uint i = 0; i < ghostToGlobalID.size(); i++) {
      getline(inp, line);
      ghostToGlobalID[i] = atoi(line.c_str());
    }

    inp.close();
  } // end of _getNodalMapping
  */


  // Initialize output file
  void Model::initializeOutput(const int DOF_PER_NODE   , const char* DOF_NAMES,
			       const int OUTPUT_PER_ELEM, const char* OUTPUT_NAMES,
			       const int NSTEPS,
			       const string OutputFile, const string format)
  {
    // Write Output Header Information
    string fileName(OutputFile);

    ofstream out;
    char header[500];
    const int nNodes = _myMesh->getNumberOfNodes();
    const int nElems = _myMesh->getNumberOfElements();

    strcpy(header, "NLOCALNODES DOF_PER_NODE DOF_NAMES\nNELEMENTS OUTPUT_PER_ELEM OUTPUT_NAMES\nNSTEPS\n");
    // Writing header information
    if (format == "BINARY") {
      out.open(fileName.c_str(), ios::out|ios::binary);
      out.write(header, 500);
      out.write((char*)(&nNodes), sizeof(int));
      out.write((char*)(&DOF_PER_NODE), sizeof(int));
      out.write((char*)(DOF_NAMES), sizeof(int)); // Each out name should not be more than 30 chars
      out.write((char*)(&nElems), sizeof(int));
      out.write((char*)(&OUTPUT_PER_ELEM), sizeof(int));
      out.write((char*)(OUTPUT_NAMES), sizeof(int));
      out.write((char*)(&NSTEPS), sizeof(int));
    } 
    else {
      out.open( fileName.c_str() );
      out << header << endl;
      out << _myMesh->getNumberOfNodes() << " " << DOF_PER_NODE << " " << DOF_NAMES << endl;
      out << _myMesh->getNumberOfElements() << " " << OUTPUT_PER_ELEM << " " << OUTPUT_NAMES << endl;
      out << NSTEPS << endl;
    }

    // Writing results is done in the derived classes
    out.close();

  } // end of initializeOutput



} // Namespace voom
