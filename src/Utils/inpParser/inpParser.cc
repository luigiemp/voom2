//-*-C++-*-
/*
  Parses input mesh data and creates a .v2i input to be used with voom2
  for analysis.

  Input
  ******
  A text file which defines the following
  MODEL  - Input FE model Abaqus or vtk format
  ELEMENTMAP - What elemental mapping has to be used. 

  Output
  ******
  The output v2i files writes information in pseudo abaqus format. We use
  most of the cards from Abaqus. However node and element numbers start 
  from 0. The deck contains
  1. Nodal ID/position information
        We first write the number of nodes (for that cpu) in the model. 
	This include local and global nodes. Example
	*NUMBEROFNODES
	4
	*NODE
	0.0 0.0 0.0
	0.0 0.0 1.0
	0.0 1.0 1.0
	0.0 1.0 0.0

	We cannot mix 2D and 3D. All nodes should be same dimensions
  2. Element definition
        Similar to nodes we first write the number of nodes (for that cpu)
	data. Elements are only local. Then when defining an element we 
	first write a string which defines the element type (Abaqus format + 
	our own format), connectivity table, Material data for this 
	element. For FE we use Abaqus naming convention. For Meshfree we will 
	use our own naming convention. Example
	*NUMBEROFELEMENTS
	2
	*ELEMENT
	S3 0 1 2 STEEL
	S3 1 2 3 WOOD

	Here we create two tria elements with material name STEEL and WOOD.
  3. Node Sets
        Needs to be implemented
  4. Element Sets
        Needs to be implemented
  5. Nodal Data. Scalars/vectors/tensors to be stored at nodes
        Needs to be implemented
  6. Elemental Data. Scalars/vector/tensors to be stored at elements
        Needs to be implemented
  7. Material
        We use *MATERIAL card similar to Abaqus. The name of the material
	will be same as what appears in the last column of the element.
	NOTE: To achieve this when using hypermesh name the component and 
	material the same. Then nothing more needs to be done.
  8. Load Steps
        We use *STEP and *ENDSTEP to denote a loading. Current *BOUNDARY is 
	used to contrain displacement and temperature. Abaqus uses DOF 
	starting from 1. Voom2 uses it from 0. Hence inpParser will decrement 
	values to start it from 0. If the input format is vtk then this 
	decrement is NOT done. The vtkbc file should be written by the 
	use assuming DOF starts at 0.


 Refer to user manual for more information (To be written)
 Added documentation

*/
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <dirent.h>
#include <vector>
#include <set>
#include <map>

// Vtk Headers
#include <vtkSmartPointer.h>
#include <vtkIdList.h>
#include <vtkDataSet.h>
#include <vtkDataSetReader.h>
#include <vtkPoints.h>

using namespace std;
using namespace boost;

//***************************************************************************
// Data structure to read in EBC/NBC 
//***************************************************************************
struct NodalData {
  // Node ID at which constraint is applied
  vector<int> nodalID;
  // DOF data, constraint values are stored as a string. Easier handling
  vector<string> data;
};

struct ElemData {
  // Element ID at which constraint is applied
  vector<int> elementID;
  // DOF data, constraint values are stored as a string. Easier handling
  vector<string> data;
};

struct BC {
  // Analysis definition
  string    analysisType;
  // EBC on nodes
  NodalData boundary;
  // NBC on nodes
  NodalData cload;
  // NBC on elements
  ElemData  dload;
};

struct MaterialData {
  vector<string> data;
};

//***************************************************************************
//                            Main Program
//***************************************************************************
int main(int argc, char** argv){
  cout << "File     : " << __FILE__ << endl;
  cout << "Compiled : " << __DATE__ << " " << __TIME__ << endl << endl;
  if (argc != 2) {
    cout << "Usage : " << argv[0] << " Inputfile\n";
    exit(1);
  }

  ifstream inp( argv[1] );
  if (!inp.is_open() ) {
    cerr << "ERROR opening input file " << argv[1] << endl;
    cerr << "Exiting....\n";
    exit(1);
  }

  //***************************************************************************
  // Read the input file for inpParser
  //***************************************************************************
  string file, fileName, extension, line, metisFile, ctrlFile;
  const int width = 15;
  map<string, string> elementMap;     vector< string > strs;
  while( inp.good() ) {
    getline(inp, line); trim(line); strs.clear();
    split( strs, line, is_any_of("\t "));
    to_upper(strs[0]);
    // Input File Name
    if ( strs[0] == "MODEL" ) { 
      fileName = strs[1];  file = fileName.substr(0, fileName.length()-4);
      extension = fileName.substr( fileName.length()-3, fileName.length());
      metisFile = fileName.substr(0, fileName.length()-4) + ".metis.";
    }
    // Elemental Map information
    if ( strs[0] == "ELEMENTMAP")  elementMap[ strs[1] ] = strs[2];
  }
  inp.close();

  if ( elementMap.size() == 0) {
    cerr << "** ERROR: No Elemental Mapping card found."
	 << " Please provide this data\n"
	 << "Exiting...\n";
    exit(1);
  }

  // Data structure to hold nodes and element numbers - For each CPU
  vector< set<int> > localNodes, localElements, ghostNodes;
  // For each CPU map of Global to Local Node ID
  vector< map<int, int> > nodalRenum;
  // For each global Node ID w know which cpu it lives on
  vector<int> globalNodeIdToCpuId;
  // Node/Element set names
  vector< string > nSets, eSets;
  // Global Load Steps
  vector< BC > loadSteps; 
  // Material data
  vector< MaterialData> matData;
  // For a global element ID we get which CPU it lives on and its
  // local element number. Use by Distributed load boundary condition.
  map<int, vector<int> > globalToLocalElemId;

  //******************************************************************
  // Parse Control deck if present
  //******************************************************************
  if (extension == "vtk") ctrlFile = file + ".vtkbc"; else ctrlFile = fileName;
  cout << "** Reading Control deck data\n";
  inp.open( ctrlFile.c_str() );
  while (inp.good() ) {
    getline(inp, line);
    to_upper(line);
    // Read Material Data
    if (find_first(line, "*MATERIAL")) {
      MaterialData material; trim(line);
      material.data.push_back(line);
      while( getline(inp, line) ) { // Reading material info
	to_upper(line);
	if ( find_first(line, "**" )) continue;
	if ( find_first(line, "*MATERIAL" ) || 
	     find_first(line, "*STEP") ) { 
	  inp.seekg( int(inp.tellg()) - line.length() - 1); line = "";
	  break;
	}
	trim(line);
	material.data.push_back( line );
      } // Read material data
      matData.push_back( material );
    } // If * material

    // If * Boundary
    if ( find_first(line, "*STEP") ) {
      BC loadStep;
      // Read next line which says analysis type	
      while( getline(inp, line) ) { // Reading step info
	to_upper(line);
	if ( find_first(line, "**" )) continue;
	else {
	  loadStep.analysisType = line; break;
	}
      };
      // Read step information
      while( getline(inp, line) ) { // Reading step info
	to_upper(line);
	if ( find_first(line, "*END") ) break;
	if ( find_first(line, "*BOUNDARY") ) {
	  while( getline(inp, line) ) {
	    if ( find_first(line, "**" )) continue;
	    if ( find_first(line, "*" )) { 
	      inp.seekg( int(inp.tellg()) - line.length() - 1);
	      break;
	    }
	    strs.clear(); trim(line); split( strs, line, is_any_of(","));
	    // Decrement node number by 1 if input is abaqus
	    if ( extension == "inp")
	      loadStep.boundary.nodalID.push_back( atoi(strs[0].c_str())-1 );
	    else
	      loadStep.boundary.nodalID.push_back( atoi(strs[0].c_str()) );
	      
	    string a;
	    for(vector<string>::iterator it = strs.begin() + 1;
		it != strs.end() - 1; it++) {
	      trim(*it); 
	      int dof = atoi( it->c_str() ); 
	      // Incase 2nd dof value is not written use the first dof value
	      if (it->length() == 0) dof = atoi( (it-1)->c_str() );
	      if (extension == "inp" ) dof--; char buf[20];
	      sprintf(buf, "%10d", dof);
	      a += " " +  string(buf);
	    }
	    char buf[20]; sprintf(buf, "%10s", strs.back() );
	    a += " " + string( buf );
	    loadStep.boundary.data.push_back( a );
	  }
	} // End of boundary check
	if ( find_first(line, "*CLOAD") ) {
	  while( getline(inp, line) ) {
	    if ( find_first(line, "**" )) continue;
	    if ( find_first(line, "*" )) { 
	      inp.seekg( int(inp.tellg()) - line.length() - 1);
	      break;
	    }
	    strs.clear(); trim(line); split( strs, line, is_any_of(","));
	    // Decrement node number by 1 if input is abaqus
	    if (extension == "inp")
	      loadStep.cload.nodalID.push_back( atoi(strs[0].c_str())-1 );
	    else
	      loadStep.cload.nodalID.push_back( atoi(strs[0].c_str()) );

	    int dof = atoi(strs[1].c_str());
	    if ( extension == "inp") dof --;
	    char buf[20]; sprintf( buf, "%10d", dof);
	    string a = string(buf);
	    sprintf(buf, "%15s", strs.back() );
	    a += string(buf);
	    loadStep.cload.data.push_back( a );
	  }
	} // End of Cload check
	if ( find_first(line, "*DLOAD") ) {
	  while( getline(inp, line) ) {
	    if ( find_first(line, "**" )) continue;
	    if ( find_first(line, "*" )) { 
	      inp.seekg( int(inp.tellg()) - line.length() - 1);
	      break;
	    }
	    strs.clear(); trim(line);
	    split( strs, line, is_any_of(","));
	    if (extension == "inp")
	      loadStep.dload.elementID.push_back( atoi(strs[0].c_str()) - 1);
	    else
	      loadStep.dload.elementID.push_back( atoi(strs[0].c_str()) );

	    string a;
	    for(vector<string>::iterator it = strs.begin() + 1;
		it != strs.end(); it++) {
	      char buf[20]; sprintf(buf, "%15s", *it );
	      a += string( buf );
	    }
	    loadStep.dload.data.push_back( a );
	  }
	} // End of Dload check
      } // Read data per step
      loadSteps.push_back( loadStep );
    } // End of find *STEP
  } // End of reading input file
  inp.close();
  cout << "** Number of material data in Ctrl deck : " << matData.size() 
       << endl;
  cout << "** Number of load steps in Ctrl deck    : " << loadSteps.size() 
       << endl;

  //**************************************************************************
  // Read Model Information, Nodal Position, element connectivity
  //**************************************************************************
  vector< vector<double> > position; // Store position for all nodes
  vector< vector<int> > connectivity; // Store connectivity for all elements
  vector< string > elementType; // Element types for all elements
  if ( extension == "vtk" ) {
    cout << "** Parsing input VTK File\n";
    vtkSmartPointer<vtkDataSetReader> reader = 
      vtkSmartPointer<vtkDataSetReader>::New();
    reader->SetFileName( fileName.c_str() );
    reader->Update();
    reader->GetOutput()->Register(reader);
    vtkDataSet* dataSet = vtkDataSet::SafeDownCast( reader->GetOutput() );
    position.resize( dataSet->GetNumberOfPoints() );
    for(uint i = 0; i < dataSet->GetNumberOfPoints(); i++) {
      vector<double> Position(3, 0.);
      dataSet->GetPoint(i, &Position[0]);
      position[i] = Position;
    }
    connectivity.resize( dataSet->GetNumberOfCells() );
    elementType.resize( dataSet->GetNumberOfCells() );
    // In Vtk all elements belong to only one element Set
    eSets.resize( dataSet->GetNumberOfCells(), "ALL" );
    for(int i = 0; i < dataSet->GetNumberOfCells(); i++) {
      vtkSmartPointer<vtkIdList> ptsIds = vtkSmartPointer<vtkIdList>::New();
      dataSet->GetCellPoints(i, ptsIds);
      vector<int> conn( ptsIds->GetNumberOfIds());
      for(uint m = 0; m < conn.size(); m++) conn[m] = ptsIds->GetId( m );
      connectivity[i] = conn;
      elementType[i] = lexical_cast<string>( dataSet->GetCellType(i) );
    }
  } else if ( extension == "inp" ) {
    inp.open( ctrlFile.c_str() );
    while( getline(inp, line) ) {
      // Read nodal info
      trim(line), to_upper(line);
      if (find_first(line, "*NODE")) {
	while( getline(inp, line) ) {
	  if ( find_first(line, "**" )) continue;
	  if ( find_first(line, "*" )) { 
	    inp.seekg( int(inp.tellg()) - line.length() - 1);
	    break;
	  }
	  strs.clear(); trim(line);
	  split( strs, line, is_any_of(","));
	  vector<double> pos;
	  for(uint m = 1; m < strs.size(); m++) 
	    pos.push_back( atof( strs[m].c_str() ) );
	  position.push_back( pos );
	}
      } // Nodal reading
      // Read Element Info
      if (find_first(line, "*ELEMENT")) {
	// Get element Type first
	strs.clear(); split( strs, line, is_any_of(",") );
	vector<string> data; split( data, strs[1], is_any_of("="));
	string eType = data[1]; data.clear();
	split( data, strs[2], is_any_of("="));
	string setName = data[1];
	while( getline(inp, line) ) {
	  if ( find_first(line, "**" )) continue;
	  if ( find_first(line, "*" )) { 
	    inp.seekg( int(inp.tellg()) - line.length() - 1);
	    break;
	  }
	  if ( eType == "C3D8" || eType == "C3D8R" ) {
	    string nLine; getline( inp, nLine); line += nLine;
	  }
	  strs.clear(); trim(line);
	  split( strs, line, is_any_of(","));
	  vector<int> conn;
	  elementType.push_back( eType );
	  for(uint m = 1; m < strs.size(); m++) 
	    conn.push_back( atoi( strs[m].c_str() ) - 1 );
	  eSets.push_back( setName );
	  connectivity.push_back( conn );
	} // Inner getline loop
      } // Elemental reading
    } // Reading the file
    inp.close();
  } else {
    cerr << "** ERROR: Unknown extension " << extension << endl;
    cerr << "** Exiting...\n";
    exit(1);
  }

  //**************************************************************************
  // See if metis file exists
  //**************************************************************************
  char buff[200];  getcwd(buff, 200);
  DIR* cwd = opendir( buff );  struct dirent *DirEntry;
  int nCpu = 1;
  globalNodeIdToCpuId.resize( position.size(), 0 ); // Resize array
  while( DirEntry =readdir(cwd) ) {
    string str(DirEntry->d_name);
    if ( str.find( metisFile) != string::npos ){
      split( strs, str, is_any_of(".") );
      nCpu = atoi( strs[ strs.size() - 1].c_str() ); // Find the number of cpu
      break;
    }
  }
  localNodes.resize( nCpu );    localElements.resize( nCpu );
  nodalRenum.resize( nCpu );    ghostNodes.resize( nCpu );  
  // Read metis file if nCpu is not 1.
  if ( nCpu != 1 ) {
    string nMetisFile = file + ".metis.npart." + lexical_cast<string>(nCpu);
    string eMetisFile = file + ".metis.epart." + lexical_cast<string>(nCpu);
    cout << "** Reading domain decomposition data\n";
    cout << "** Number of Cpus: " << nCpu << endl;
    ifstream nfile(nMetisFile.c_str()), efile(eMetisFile.c_str());
    int id = 0, cpuNumber;
    // Read local nodes
    while( nfile >> cpuNumber ) { 
      // IN a metis file line number corresponds to the node number and 
      // the cpu entry in that line corresponds to the cpu where this 
      // node lives.
      localNodes[cpuNumber].insert( id ); 
      globalNodeIdToCpuId[id] = cpuNumber;
      id++; // Line Number counter
    }
    // Read local elements
    id = 0;
    while( efile >> cpuNumber ) { 
      localElements[cpuNumber].insert( id );
      vector<int> data(2); 
      data[0] = cpuNumber; data[1] = localElements[cpuNumber].size();
      globalToLocalElemId[id] = data;
      id++;
    }
    nfile.close(); efile.close();
    // For a given local element find nodes which make up this element. Then
    // look for these nodes in the localNodes[cpuId] set. If not present it 
    // is a ghost node.
    for(int cpuId = 0; cpuId < nCpu; cpuId++) {
      for( set<int>::iterator it = localElements[cpuId].begin(); 
	   it != localElements[cpuId].end(); it++){
	// Parse over connectivity to find ghost nodes
	for(unsigned int nodeId = 0; nodeId < connectivity[*it].size(); 
	    nodeId++) {
	  int node = connectivity[*it][nodeId];
	  if( localNodes[cpuId].find(node) == localNodes[cpuId].end() )
	    ghostNodes[cpuId].insert(node);
	}
      }
    }
    // Renumber local and ghost Nodes for each CPU.
    for(int cpuId = 0; cpuId < nCpu; cpuId++) {
      int localNumber = 0;
      for(set<int>::iterator it = localNodes[cpuId].begin();
	  it != localNodes[cpuId].end(); it++)
	nodalRenum[cpuId][*it] = localNumber++;
      for(set<int>::iterator it = ghostNodes[cpuId].begin();
	  it != ghostNodes[cpuId].end(); it++)
	nodalRenum[cpuId][*it] = localNumber++;
    }
  } // Loop to read metis file
  else {
    // Job is on one cpu so all data is trivial.
    // Fill Local Node Info
    for(int i = 0; i < position.size(); i++) {
      localNodes[0].insert(i);
      nodalRenum[0][i] = i;
    }
    // Fill local elements info
    for(int i = 0; i < connectivity.size(); i++) {
      localElements[0].insert(i);
      vector<int> data(2); data[0] = 0; data[1] = i;
      globalToLocalElemId[i] = data;
    }    
  } // End of domain decomposition data parsing

  //******************************************************************
  // Writing model and domain decomposition data
  //******************************************************************    
  // Writing input for voom2
  for(unsigned int cpuId = 0; cpuId < nCpu; cpuId++) {
    char outfile[200];
    sprintf(outfile, "%s.%03d.v2i", file.c_str(), cpuId);
    ofstream out( outfile );
    out<< "**\n** Voom2 input for Processor = " << cpuId << endl;

    // Nodal Data      
    out << "** Nodes are assumed to be sequentially numbered from "
        << "0\n";
    out << "**\n*NUMBEROFNODES\n"  
        << localNodes[cpuId].size() + ghostNodes[cpuId].size() << endl;
    out << "*NODE\n";
    // Write Local Nodes first
    for(set<int>::iterator it = localNodes[cpuId].begin();
	it != localNodes[cpuId].end(); ++it) 
      out << setw(width) << position[*it][0] << setw(width)
	  << position[*it][1] << setw(width) 
	  << position[*it][2] << endl;
    // Write ghost nodes next
    for(set<int>::iterator it = ghostNodes[cpuId].begin();
	it != ghostNodes[cpuId].end(); ++it) 
      out << setw(width) << position[*it][0] << setw(width)
	  << position[*it][1] << setw(width) 
	  << position[*it][2] << endl;
      
    out << "** Elements are assumed to be sequentially numbered from "
        << "0\n";    
    out << "**\n*NUMBEROFELEMENTS\n" << localElements[cpuId].size() 
        << endl;
    out << "*ELEMENT\n";
    // Filling Cell Data
    for( set<int>::iterator it = localElements[cpuId].begin(); 
	 it != localElements[cpuId].end(); it++){
      string eType = elementType[*it];
      if ( elementMap.find(eType) != elementMap.end() ) 
	eType = elementMap[eType];
      out << setw(width) << eType;
      for(unsigned int nodeId = 0; nodeId < connectivity[*it].size(); 
	  nodeId++) 
	out << setw(width) << nodalRenum[cpuId][connectivity[*it][nodeId]];
      out << setw(width) << eSets[*it];
      out << endl;
    }
    // Number of Local Nodes in the cpu
    out << "**\n** Number of Local Nodes\n**\n";
    out << "*NUMLOCALNODES\n" << localNodes[cpuId].size() << endl;

    //  A map of local node number to global node number
    out << "**\n** MPI Info LocalNodeID - GlobalNodeID\n**\n";
    out << "*GLOBALID\n";
    for(set<int>::iterator it = localNodes[cpuId].begin();
	it != localNodes[cpuId].end(); ++it) 
      out << setw(width) << *it << endl;
    
    for(set<int>::iterator it = ghostNodes[cpuId].begin();
	it != ghostNodes[cpuId].end(); ++it) 
      out << setw(width) << *it << endl;

    // Writing Material data
    out << "**\n** Material data\n**\n";
    for(uint m = 0; m < matData.size(); m++)
      for(vector<string>::iterator it = matData[m].data.begin();
	  it != matData[m].data.end(); it++)
	out << *it << endl;
    out.close();
  }// Cpu Loop

  //**************************************************************************
  // Writing/Appending Control Deck Data
  //**************************************************************************
  ofstream out[nCpu];
  for(unsigned int cpuId = 0; cpuId < nCpu; cpuId++) {
    char outfile[200];
    sprintf(outfile, "%s.%03d.v2i", file.c_str(), cpuId);
    out[cpuId].open( outfile, ios::app );
    out[cpuId] << "**\n** Analysis information\n**\n";
  }
  for(unsigned int stepId = 0; stepId < loadSteps.size(); stepId++) {
    for(unsigned int cpuId = 0; cpuId < nCpu; cpuId++) 
      out[cpuId] << "**\n** Load step: " << stepId << "\n*STEP\n"
		 << loadSteps[stepId].analysisType << "\n"
		 << "*BOUNDARY\n";
    map<int, int>::iterator it;
    // Write *BOUNDARY cards
    if (loadSteps[stepId].boundary.data.size() != 0 ) 
      for(int m = 0; m < loadSteps[stepId].boundary.data.size(); m++) {
	int nodeID = loadSteps[stepId].boundary.nodalID[m];
	int cpuId =  globalNodeIdToCpuId[ nodeID ];
	int localId = nodalRenum[cpuId][ nodeID ];
	out[cpuId] << setw(width) << localId << setw(width) 
		   << loadSteps[stepId].boundary.data[m]
		   << endl;
      } 

    // Write *CLOAD cards
    for(unsigned int cpuId = 0; cpuId < nCpu; cpuId++) 
      out[cpuId] << "*CLOAD\n";
    if (loadSteps[stepId].cload.data.size() != 0 ) 
      for(int m = 0; m < loadSteps[stepId].cload.data.size(); m++) {
	int nodeID = loadSteps[stepId].cload.nodalID[m];
	int cpuId =  globalNodeIdToCpuId[ nodeID ];
	int localId = nodalRenum[cpuId][ nodeID ];
	out[cpuId] << setw(width) << localId << setw(width) 
		   << loadSteps[stepId].cload.data[m]
		   << endl;
      } 

    // Write *DLOAD cards
    for(unsigned int cpuId = 0; cpuId < nCpu; cpuId++) 
      out[cpuId] << "*DLOAD\n";
    if (loadSteps[stepId].dload.data.size() != 0 ) 
      for(int m = 0; m < loadSteps[stepId].dload.data.size(); m++) {
	int elemID = loadSteps[stepId].dload.elementID[m];
	int cpuId = globalToLocalElemId[ elemID ][0];
	int localId = globalToLocalElemId[ elemID ][1];
	out[cpuId] << setw(width) << localId << setw(width)
		   << loadSteps[stepId].dload.data[m]
		   << endl;
      } 
    
    for(unsigned int cpuId = 0; cpuId < nCpu; cpuId++) 
      out[cpuId] << "*END STEP\n";
  } // Loop over steps

  // Closing file handles
  for(unsigned int cpuId = 0; cpuId < nCpu; cpuId++) out[cpuId].close();
}
