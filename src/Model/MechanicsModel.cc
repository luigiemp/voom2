#include "MechanicsModel.h"

namespace voom {

  // Constructor
  MechanicsModel::MechanicsModel(Mesh* aMesh, 
				 const vector<string > & ElMatType, 
				 const map<string, MechanicsMaterial* > & ElMaterials): 
    EllipticModel(aMesh, aMesh->getDimension())
  {
    // Resize and initialize (default function) _field vector
    _field.resize(  (_myMesh->getLocalDoF()).size() + (_myMesh->getGhostDoF()).size() );
    this->initializeField();

    // Fill in material vectors - ONE MATERIAL PER EACH QUADRATURE POINT
    vector<GeomElement* > Elements = _myMesh->getElements() ;
    uint NumEl = Elements.size();
    assert( ElMatType.size() == NumEl );
    
    // Create new material pointers so that _materials does not store  
    // pointers created outside Model and on which Model has no control
    map<string, MechanicsMaterial* >::const_iterator it;
    map<string, MechanicsMaterial* > ElMaterialsCopy;
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
      MechanicsMaterial* ElMat = ElMaterialsCopy[ElMatType[e]];
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
  


    

    // // Resize and initialize (default function) _field vector
    // _field.resize(_localDoF.size() + _ghostDoF.size());
    // // If finite kinematics, _field needs to be initialized in the main function with the initial nodal positions
    // this->initializeField(1.0); 
    
    // // Read only relevant information from the input file. Will ignore
    // // unrelated cards
    // ifstream inp(_inputFile.c_str());
    // string line;

    // map<string, Real> matData;   // String to Material Property
    // while (getline(inp, line) ) {
    //   // Read Material Data
    //   if ( line.find("*MATERIAL") == 0) {
    // 	vector<string> strs = splitString(line, "=");
    // 	// Save part of line which contains material information
    // 	string matName = strs.back();
    // 	// prepare string for next use
    // 	strs.clear(); 
    // 	// Read all material mechanics coefficients
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
    //   } // if ( line.find("*MATERIAL") == 0) loop

    //   // If STEP is found all materials have been read
    //   if ( line.find("*STEP") == 0 ) {
    // 	// Reposition file reading to beginning of STEP line
    // 	inp.seekg( int(inp.tellg()) - line.length() - 1);
    // 	break;
    //   }

    // } // while loop - End of reading material data

    // // Create one material per each element based on:
    // // - the list of material in matData
    // // - the list of materials per each element
    // vector<string> names;
    // this->_getMaterialNamesList(names); // names has the size of element number
    // // Fill the list of material objects
    // _material.reserve(names.size() ); // one material object per each element
    // for(uint i = 0; i < names.size(); i++) 
    // {
    //   // For now we assume isotropic material only
    //   // - TO BE CHANGED to include anysotropic diffusion
    //   // (Model constructor and inputparser need also to be changed)
    //   _material.push_back(new IsotropicDiffusion(matData[ names[i] ]) );
    // }
    // /*

    


    // // Read step data
    // while (getline(inp, line) ) {
    //   // Find line which starts STEP type data
    //   if ( line.find("*STEP") == 0 ) {
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
/*! 
      Dirchilet boundary condition. While parsing input loadstep for Essential
      boundary conditions we have node_id, start_dof, end_dof, value. The 
      data from start_dof to end_dof is rolled into a single array. So the 
      entries in EBCDof array will be start_dof + nodal_id*dim to 
      end_dof + nodal_id*dim.
    */
   
    // 	      // hence subtract d - start
    // 	      for(unsigned int d = start; d <= stop; d++) {
    // 		_EBC.DoF.push_back(nodeID  + d - start);
    // 		_EBC.Value.push_back(value);
    // 	      } // d loop
    // 	    } // while getline loop
    // 	  } // if boundary check loop
    // 	} // while 'Read step information' loop
    //   } // End of find *STEP
    // } // While getline loop
    // inp.close();
    // */



  // Compute deformation gradient
  void MechanicsModel::computeDeformationGradient(vector<Matrix3d > & Flist, 
						  GeomElement* geomEl)
  {
    // Compute F at all quadrature points
    const uint numQP = geomEl->getNumberOfQuadPoints();
    const uint dim   = _myMesh->getDimension();

    const vector<int > & NodesID = geomEl->getNodesID();
    const uint nodeNum = NodesID.size();

    for(uint q = 0; q < numQP; q++) {
      // Initialize F to zero
      Flist[q] = Matrix3d::Zero();
	for(uint I = 0; I < dim; I++) 
	  for(uint j = 0; j < dim; j++)
	    for(uint a = 0; a < nodeNum; a++)
	    Flist[q](I,j) += 
	      _field[NodesID[a]*dim + I] * geomEl->getDN(q, a, j);
    } // loop over quadrature points

  }



  // Compute Function - Compute Energy, Force, Stiffness
  void MechanicsModel::compute(EllipticResult & R)
 {

   // -----------------------------------
   // LEP: BC TO BE IMPLEMENTED !!!
   // ----------------------------------- 

    const vector<GeomElement* > elements = _myMesh->getElements();
    const int dim = _myMesh->getDimension();

    // Reset values in result struct
    if ( R.getRequest() & ENERGY ) {
      R.setEnergy(0.0);
    }
    if ( R.getRequest() & FORCE )  {
      R.resetResidualToZero();
    }
    if ( R.getRequest() & STIFFNESS ) { 
      R.resetStiffnessToZero();
    }

    // Loop through elements, also through material points array, which is unrolled
    uint index = 0;
    MechanicsMaterial::FKresults FKres;
    FKres.request = R.getRequest();
    for(int e = 0; e < elements.size(); e++)
    {
      GeomElement* geomEl = elements[e];
      const vector<int  >& NodesID = geomEl->getNodesID();
      const int numQP    = geomEl->getNumberOfQuadPoints();
      const int numNodes = NodesID.size();
      
      // F at each quadrature point are computed at the same time in one element
      vector<Matrix3d > Flist(numQP, Matrix3d::Zero());
      // Compute deformation gradients for current element
      this->computeDeformationGradient(Flist, geomEl);
      
      // Loop over quadrature points
      for(int q = 0; q < numQP; q++) {
	_materials[index]->compute(FKres, Flist[q]);
	index++;

	// Volume associated with QP q
	Real Vol = geomEl->getQPweights(q);

	// Compute energy
	if (R.getRequest() & ENERGY) {
	  R.addEnergy(FKres.W*Vol);
	}

	// Compute Residual
	if (R.getRequest() & FORCE) {	
	  for(uint a = 0; a < numNodes; a++) {
	    for(uint i = 0; i < dim; i++) {
	      Real tempResidual = 0.0;
	      for (uint J = 0; J < dim; J++) { 
		tempResidual += FKres.P(i,J) * geomEl->getDN(q, a, J);
	      } // J loop
	      tempResidual *= Vol;
	      R.addResidual(NodesID[a]*dim+i, tempResidual);
	    } // i loop
	  } // a loop
	} // Internal force loop

	// Compute Stiffness
	if (R.getRequest() & STIFFNESS) {
	  for(uint a = 0; a < numNodes; a++) {
	    for(uint i = 0; i < dim; i++) {
	      for(uint b = 0; b < numNodes; b++) {
		for(uint j = 0; j < dim; j++) {
		  Real tempStiffness = 0.0; 
		  for(uint M = 0; M < dim; M++) {
		    for(uint N = 0; N < dim; N++) {
		      tempStiffness += FKres.K.get(i, M, j, N)*elements[e]->getDN(q, a, M)*
			                     elements[e]->getDN(q, b, N);
		    } // N loop
		  } // M loop
		  tempStiffness *= Vol;
		  R.addStiffness(NodesID[a]*dim + i, NodesID[b]*dim + j, tempStiffness);
		} // j loop
	      } // b loop
	    } // i loop
	  } // a loop
	} // Compute stiffness matrix

      } // QP loop
    } // Element loop

    // Sum up all stiffness entries with the same indices
      // LEP: should/can we move it to Solver?
    if ( R.getRequest() & STIFFNESS ) {
      R.FinalizeGlobalStiffnessAssembly(); 
    }

  } // Compute Mechanics Model



  // Writing output
  void MechanicsModel::writeOutput(const string OutputFile, const string format) 
  {
    // Create outputFile name
    string fileName(OutputFile);
    ofstream out;
    
    // Writing field values to output file
    if (format == "BINARY") {
      out.open( fileName.c_str(), ios::out|ios::binary|ios::app);
      for(uint i = 0; i < (_myMesh->getLocalDoF()).size(); i++)
	out.write((char*)(&_field[i]), sizeof(Real));
    } 
    else {
      out.open( fileName.c_str(), ios::out|ios::app );
      uint i = 0;
      const uint dim   = _myMesh->getDimension();
      while ( i < (_myMesh->getLocalDoF()).size() ) {
	for ( uint j = 0; j < dim; j++ ) {      
	  out << _field[i] << " ";
	  i++;
	}
	out << endl; // 1 nodal position per line
      }
    }
    // Close file
    out.close();
  } // writeOutput



} // namespace voom
