#include "PoissonModel.h"

namespace voom {

  PoissonModel::PoissonModel(Mesh* myMesh, 
			     vector<DiffusionMaterial * > Materials, 
			     const uint NodeDoF):
    EllipticModel(myMesh, NodeDoF), _materials(Materials)
  {
    // Assume NodeDoF = 1 //

    // Initialize field
    _field.resize(  (_myMesh->getNumberOfNodes() )*_nodeDoF );
    this->initializeField();

    this->setPrevField();
  }
  


  // Compute Function - Compute Energy, Force, Stiffness
  void PoissonModel::compute(EllipticResult & R) 
  {
    const vector<GeomElement* > elements = _myMesh->getElements();
    const int AvgNodePerEl = ((elements[0])->getNodesID()).size();
    const int NumEl = elements.size();
    const int dim = _myMesh->getDimension();
    vector<Triplet<Real > > KtripletList;
    
    int PbDoF = R.getPbDoF();

    // Reset values in result struct
    if ( R.getRequest() & ENERGY ) {
      R.setEnergy(0.0);
    }
    if ( R.getRequest() & FORCE )  {
      R.resetResidualToZero();
    }
    if ( R.getRequest() & STIFFNESS ) { 
      R.resetStiffnessToZero();
      KtripletList.reserve(NumEl*AvgNodePerEl*AvgNodePerEl);
    }
    
    

    // Loop through elements, also through material points array, which is unrolled
    DiffusionMaterial::DiffusionResults DiffRes;
    for(int e = 0; e < NumEl; e++)
    {
      GeomElement* geomEl = elements[e];
      const vector<int  >& NodesID = geomEl->getNodesID();
      const int numQP    = geomEl->getNumberOfQuadPoints();
      const int numNodes = NodesID.size();
      MatrixXd Kele = MatrixXd::Zero(numNodes, numNodes);

      // Compute Residual
      
      for(int q = 0; q < numQP; q++) {
	_materials[e*numQP + q]->compute(DiffRes, 0.0);
	Real Vol = elements[e]->getQPweights(q);
    
	if (R.getRequest() & FORCE) {
    	  for(int a = 0; a < numNodes; a++) {
    	    Real tempResidual = 0.0;
    	    for(uint b = 0; b < numNodes; b++) {
    	      Real ScalarField = _field[NodesID[b]];
    	      for (int I = 0; I < dim; I++) {
    		for (int J = 0; J < dim; J++) { 
    		  tempResidual += DiffRes.A(I,J) * elements[e]->getDN(q, a, I) * 
    		    elements[e]->getDN(q, b, J) * ScalarField;
    		} // J loop
    	      } // I loop
    	    } // b loop
    	    tempResidual *= Vol;
    	    R.addResidual(NodesID[a], tempResidual);
    	  } // a loop
	} // Internal force loop
       
	// Compute Hessian
	if (R.getRequest() & STIFFNESS) {
    	  for(int a = 0; a < numNodes; a++) {
    	    for(int b = 0; b < numNodes; b++) {
    	      Real tempStiffness = 0.0; 
    	      for (int I = 0; I < dim; I++) {
    		for(int J = 0; J < dim; J++)  {
    		  tempStiffness += DiffRes.A(I,J)*elements[e]->getDN(q, a, I)* 
    		    elements[e]->getDN(q, b, J);
    		} // J loop
    	      } // I loop
    	      tempStiffness *= Vol;
	      KtripletList.push_back( Triplet<Real >( NodesID[a], NodesID[b], tempStiffness ) );
    	    } // b loop
    	  } //a loop
	} // Compute hessian matrix
      } // q loop
  
    } // element loop
    
    // Sum up all stiffness entries with the same indices
    if ( R.getRequest() & STIFFNESS ) {
      R.setStiffnessFromTriplets(KtripletList);
      R.FinalizeGlobalStiffnessAssembly(); 
      // cout << "Stiffness assembled" << endl;
    }
  
  } // Poisson Model compute
  


  // Writing output
// Writing output
  void PoissonModel::writeOutputVTK(const string OutputFile, int step) 
  {
    ///// 
    // NEED TO BE REWRITTEN TAKING INTO ACCOUNT MULTIPLE QUADRATURE POINTS PER ELEMENT !!!
    /////
    // Create outputFile name
    stringstream FileNameStream;
    FileNameStream << OutputFile << step << ".vtk";
    ofstream out;
    out.open( (FileNameStream.str()).c_str() );
    int NumNodes = _myMesh->getNumberOfNodes();
  
    // Header
    char spchar = '#';
    out << spchar << " vtk DataFile Version 3.1" << endl;
    out << "vtk output" << endl;
    out << "ASCII" << endl;
    out << "DATASET UNSTRUCTURED_GRID" << endl;
    out << "POINTS " << NumNodes << " FLOAT" << endl;

    // For now we assumed dim == 3 !
    for (int i = 0; i < NumNodes; i++ ) {
      out << _myMesh->getX(i,0) << " " << _myMesh->getX(i,1) << " " << _myMesh->getX(i,2) << endl;
    }

    vector<GeomElement* > elements = _myMesh->getElements();
    int NumEl = elements.size();
    int NodePerEl = (elements[0])->getNodesPerElement();
    // To be adjusted - not general at all!
    int CellType = 0;
    switch (NodePerEl) {
    case 4:
      CellType = 10;
      break;
    case 10:
      CellType = 24;
      break;
    default:
      cout << "Error cell type not implemented in PoissonModel writeOutput. " << endl;
    }

    out << endl << "CELLS " << NumEl << " " << NumEl*(NodePerEl+1) << endl;

    for (int e = 0; e < NumEl; e++) {
      out << NodePerEl << " ";
      const vector<int > & NodesID = (elements[e])->getNodesID();
      for (int n = 0; n < NodePerEl; n++) {
	out << NodesID[n] << " ";
      }
      out << endl;
    }
    
    out << endl << "CELL_TYPES " << NumEl << endl;
    for (int e = 0; e < NumEl; e++) {
      out << CellType << " " << endl;
    }  // end of printing conntable

    // Point data section
    // Displacements
    out << endl << "POINT_DATA " << NumNodes << endl
	<< "SCALARS field double" << endl
	<< "LOOKUP_TABLE default" << endl;
  
    int dim = _myMesh->getDimension();
    for (int i = 0; i < NumNodes; i++ ) {
      out << _field[i] << endl;
    }    

    // Close file
    out.close();
  } // writeOutputVTK

  void PoissonModel::writeOutput(const string OutputFile, const string format) 
  {
    // Create outputFile name
    string fileName(OutputFile);
    ofstream out;

    // Writing field values to output file
    if (format == "BINARY") {
      out.open( fileName.c_str(), ios::out|ios::binary|ios::app);
      for(uint i = 0; i < _field.size(); i++)
	out.write((char*)(&_field[i]), sizeof(Real));
    } else {
      out.open( fileName.c_str(), ios::out|ios::app );
      for(uint i = 0; i < _field.size(); i++)
	out << _field[i] << endl;
    }
    // Close file
    out.close();
  } // writeOutput

} // namespace voom




  
