#include "PoissonModel.h"

namespace voom {

  // Constructor
  // PoissonModel::PoissonModel(Mesh* myMesh, const string inputFile, const uint NodeDoF):
  //   EllipticModel(myMesh, inputFile, NodeDoF) {
    
  //   // Resize and initialize (default function) _field vector
  //   _field.resize(  _localDoF.size() + _ghostDoF.size());
  //   this->initializeField(0.0);

  //   // Read only relevant information from the input file. Will ignore
  //   // unrelated cards
  //   ifstream inp(_inputFile.c_str());
  //   string line;

  //   map<string, Real> matData;   // String to Material Property
  //   while (getline(inp, line) ) {
  //     // Read Material Data
  //     if ( line.find("*MATERIAL") == 0) {
  // 	vector<string> strs = splitString(line, "=");
  // 	// Save part of line which contains material information
  // 	string matName = strs.back();
  // 	// prepare string for next use
  // 	strs.clear(); 
  // 	// Read all conductivity data
  // 	while( getline(inp, line) ) {
  // 	  if ( line.find("*CONDUCTIVITY") == 0) {
  // 	    getline(inp, line);
  // 	    strs = splitString(line, ",");
  // 	    // Store material conductivity in the matData vector
  // 	    // One conductivity per each material type
  // 	    matData[matName] = atof( strs[0].c_str() );
  // 	    break;
  // 	  } // if ( line.find("*CONDUCTIVITY") == 0) loop
  // 	} // while loop
  //     } // if ( line.find("*MATERIAL") == 0) loop

  //     // If STEP is found all materials have been read
  //     if ( line.find("*STEP") == 0 ) {
  // 	// Reposition file reading to beginning of STEP line
  // 	inp.seekg( int(inp.tellg()) - line.length() - 1);
  // 	break;
  //     }

  //   } // while loop - End of reading material data

  //   // Create one material per each element based on:
  //   // - the list of material in matData
  //   // - the list of materials per each element
  //   vector<string> names;
  //   this->_getMaterialNamesList(names);
  //   // Fill the list of material objects
  //   _material.reserve(names.size() ); // one material object per each element
  //   for(uint i = 0; i < names.size(); i++) 
  //   {
  //     // For now we assume isotropic material only
  //     // - TO BE CHANGED to include anysotropic diffusion
  //     // (Model constructor and inputparser need also to be changed)
  //     _material.push_back(new IsotropicDiffusion(matData[ names[i] ]) );
  //   }


  //   // Read step data
  //   while (getline(inp, line) ) {
  //     // Find line which starts STEP type data
  //     if ( line.find("*STEP") == 0 ) {
  // 	// Read step information
  // 	while( getline(inp, line) )
  // 	{ 
  // 	  if ( line.find("*END") == 0 ) break;
  // 	  // Scalar BC defined using BOUNDARY card
  // 	  if ( line.find("*BOUNDARY") == 0 ) {
  // 	    while( getline(inp, line) ) {
  // 	      if ( line.find("**" ) == 0) continue;
  // 	      if ( line.find("*") == 0) { 
  // 		inp.seekg( int(inp.tellg()) - line.length() - 1);
  // 		break;
  // 	      }
  // 	      // Line in input file contains 
  // 	      vector<string> strs = splitString( line, " ");
  // 	      int nodeID = atoi(strs[0].c_str()), 
  // 		start = atoi(strs[1].c_str()), stop = atoi(strs[2].c_str());
  // 	      Real value = atof(strs[3].c_str());

  // 	      // Scalar field denoted as 10. Model has only one dof at node
  // 	      // hence subtract d - start
  // 	      for(unsigned int d = start; d <= stop; d++) {
  // 		_EBC.DoF.push_back(nodeID  + d - start);
  // 		_EBC.Value.push_back(value);
  // 	      } // d loop
  // 	    } // while getline loop
  // 	  } // if boundary check loop
  // 	} // while 'Read step information' loop
  //     } // End of find *STEP
  //   } // While getline loop
  //   inp.close();

  // }



  PoissonModel::PoissonModel(Mesh* aMesh, 
			     const vector<string > & ElMatType, 
			     const map<string, DiffusionMaterial* > & ElMaterials):
    EllipticModel(aMesh, 1)
  {
    // Assume NodeDoF = 1 //

    // Initialize field
    _field.resize(  (_myMesh->getLocalDoF()).size() + (_myMesh->getGhostDoF()).size() );
    this->initializeField(0.0);

    // Fill in material vectors - ONE MATERIAL PER EACH QUADRATURE POINT
    vector<GeomElement* > Elements = _myMesh->getElements() ;
    uint NumEl = Elements.size();
    assert( ElMatType.size() == NumEl );
    
    // Create new material pointers so that _materials does not store  
    // pointers created outside Model and on which it has no control
    map<string, DiffusionMaterial* >::const_iterator it;
    map<string, DiffusionMaterial* > ElMaterialsCopy;
    for (it = ElMaterials.begin(); it != ElMaterials.end(); it++)
    {
      // If material has no history variables, then store the altready created pointer
      // It will be cloned anyway in the next loop over elements and quad point
      if ( (it->second)->HasHistoryVariables() ) {
	ElMaterialsCopy.insert(make_pair(it->first, it->second));
      }
      // Clone material so that in the next loop we store in _materials a Material* 
      // created here
      else {
	ElMaterialsCopy.insert( make_pair( it->first, (it->second)->clone() ) );
      }
    }
    
    for (uint e = 0; e <  NumEl; e++) // Loop over all the elements in the mesh
    { 
      DiffusionMaterial* ElMat = ElMaterialsCopy[ElMatType[e]];
      if ( ElMat->HasHistoryVariables() ) {
	for (uint j = 0; j < Elements[e]->getNumberOfQuadPoints(); j++) {
	  _materials.push_back( ElMat->clone() );
	}
      } else {
	for (uint j = 0; j < Elements[e]->getNumberOfQuadPoints(); j++) {
	  _materials.push_back( ElMat );
	}
      }
    }
    
  }
  


  // Compute Function - Compute Energy, Force, Stiffness
  void PoissonModel::compute(EllipticResult & R) 
  {
    // -----------------------------------
    // LEP: BC TO BE IMPLEMENTED !!!
    // -----------------------------------
    
    const vector<GeomElement* > elements = _myMesh->getElements();
    const int dim = _myMesh->getDimension();
    // Reset values in result struct
    if ( R.getRequest() & ENERGY ) {
      cout << "Error: Trying to compute energy in Poisson model" << endl;
      exit (EXIT_FAILURE);
    }
    if ( R.getRequest() & FORCE )  {
      R.resetResidualToZero();
    }
    if ( R.getRequest() & STIFFNESS ) { 
      R.resetStiffnessToZero();
    }

    // Loop through elements, also through material points array, which is unrolled
    uint index = 0;
    for(int e = 0; e < elements.size(); e++) {
      const vector<int  >& NodesID = elements[e]->getNodesID();
      const int numQP    = elements[e]->getNumberOfQuadPoints();
      const int numNodes = NodesID.size();
      
      // One material per quadrature point
      // Assume A is independent of T
      // Compute Residual
      
      for(int q = 0; q < numQP; q++) {
	DiffusionMaterial::DiffusionResults DiffRes;
	_materials[index]->compute(DiffRes, 0.0);
	index++;
	if (R.getRequest() & FORCE) {
    	  for(int a = 0; a < numNodes; a++) {
    	    Real tempResidual = 0.0;
    	    for(uint b = 0; b < numNodes; b++) {
    	      Real ScalarField = _field[NodesID[b]]; //nodes[b]->getV
    	      for (int I = 0; I < dim; I++) {
    		for (int J = 0; J < dim; J++) { 
    		  tempResidual += DiffRes.A(I,J) * elements[e]->getDN(q, a, I) * 
    		    elements[e]->getDN(q, b, J) * ScalarField;
    		} // J loop
    	      } // I loop
    	    } // b loop
    	    tempResidual *= elements[e]->getQPweights(q);
    	    R.addResidual(NodesID[a], tempResidual);
    	  } // a loop
	} // Internal force loop
       
	// Compute Hessian
	if (R.getRequest() & STIFFNESS) {
    	  for(int a = 0; a < numNodes; a++) {
    	    for(int b = 0; b < numNodes; b++) {
    	      Real sum = 0.0; 
    	      for (int I = 0; I < dim; I++) {
    		for(int J = 0; J < dim; J++)  {
    		  sum += DiffRes.A(I,J)*elements[e]->getDN(q, a, I)* 
    		    elements[e]->getDN(q, b, J);
    		} // J loop
    	      } // I loop
    	      sum *= elements[e]->getQPweights(q); 
    	      R.addStiffness(NodesID[a], NodesID[b], sum);
    	    } // b loop
    	  } //a loop
	} // Compute hessian matrix
      } // q loop
  
    } // element loop
    
    // Sum up all stiffness entries with the same indices
    // LEP: not sure if should be moved to Solver ...
    if ( R.getRequest() & STIFFNESS ) 
      R.FinalizeGlobalStiffnessAssembly(); 
  
  } // Poisson Model compute
  


  // Writing output
  void PoissonModel::writeOutput(const string OutputFile, const string format) 
  {
    // Create outputFile name
    string fileName(OutputFile);
    ofstream out;

    // Writing field values to output file
    if (format == "BINARY") {
      out.open( fileName.c_str(), ios::out|ios::binary|ios::app);
      for(uint i = 0; i < (_myMesh->getLocalDoF()).size(); i++)
	out.write((char*)(&_field[i]), sizeof(Real));
    } else {
      out.open( fileName.c_str(), ios::out|ios::app );
      for(uint i = 0; i < (_myMesh->getLocalDoF()).size(); i++)
	out << _field[i] << endl;
    }
    // Close file
    out.close();
  } // writeOutput

} // namespace voom




  
