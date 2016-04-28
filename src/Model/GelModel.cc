#include "GelModel.h"
#include "Filament.h"
namespace voom {

  // Constructor
  GelModel::GelModel(GelMesh* aGelMesh, vector<FilamentMaterial * > springs,
		     vector<FilamentMaterial * > angleSprings,
		     const uint NodeDoF,
		     int NodalForcesFlag,
		     int ResetFlag):
    Model( NodeDoF),_myGelMesh(aGelMesh) , _springs(springs),_angleSprings(angleSprings), 
    _nodalForcesFlag(NodalForcesFlag), _resetFlag(ResetFlag)
  {
    // THERE IS ONE MATERIAL PER ELEMENT - CAN BE CHANGED - DIFFERENT THAN BEFORE
    // Resize and initialize (default function) _field vector
    _field.resize(  (_myGelMesh->getNumberOfNodes() )*_nodeDoF );
    this->initializeField();

    VectorXd _X0 = _myGelMesh->getX();

    const vector<GeomFilament* > filElement = _myGelMesh->getFilaments();
    uint NumFil = _myGelMesh->getNumberOfFilaments();
    
    for (int nFil = 0 ; nFil < NumFil ; nFil++){

      const vector<int > & NodesID = filElement[nFil]->getNodesID();
      const uint nodeNum = NodesID.size();
      vector<Vector3d > xlist(nodeNum,Vector3d::Zero());

      getFilamentx(xlist,filElement[nFil]);
      
      Filament* Fil =  new Filament(filElement[nFil], _springs[nFil],_angleSprings[nFil]);

      _gel.push_back(Fil);
    }
    
  }
  

  void GelModel::getFilamentx(vector<Vector3d > & xlist, GeomFilament* geomEl)
  {
    // Compute all segment extensions 
    const uint dim   = _myGelMesh->getDimension();

    const vector<int > & NodesID = geomEl->getNodesID();
    const uint nodeNum = NodesID.size();
    
    for(uint n = 0; n < nodeNum; n++)
      {
	for(uint i = 0; i < dim; i++)
	  {
	    xlist[n](i) = _field[NodesID[n]*dim+i];
	  }
      }
  }
  
  void GelModel::compute(Result & R)
  {
    const vector<GeomFilament* > filElement = _myGelMesh->getFilaments();
    uint NumFil = _myGelMesh->getNumberOfFilaments();

    if ( R.getRequest() & ENERGY ) {
      R.setEnergy(0.0);
    }
    if ( R.getRequest() & FORCE )  {
      R.resetResidualToZero();
    } 

    FilamentMaterial::Filresults Rf;
    Rf.request = R.getRequest();

    for (int nFil = 0 ; nFil < NumFil ; nFil++){

      int nodeNum =  _gel[nFil]->getNumberOfNodes();
      vector<Vector3d > xlist(nodeNum,Vector3d::Zero());
      
      getFilamentx(xlist,filElement[nFil]);
      
      _gel[nFil]->compute(R, xlist);
      
    }

  }

  void GelModel::writeOutput(const string OutputFile, int step) 
  {
    // Create outputFile name
    stringstream FileNameStream;
    FileNameStream << OutputFile << step << ".vtk";
    ofstream out;
    out.open( (FileNameStream.str()).c_str() );
    int NumNodes = _myGelMesh->getNumberOfNodes();
    int dim = _myGelMesh->getDimension();

    for (int i = 0; i < NumNodes; i++ ) {
      for (int j = 0; j < dim; j++) {
	out << _myGelMesh->getX(i,j) << " ";
      }
      out << endl;
    }
    out.close();
  }
  
  
  // Writing output
  void GelModel::writeOutputVTK(const string OutputFile, int step) 
  {
    /*
    // Create outputFile name
    stringstream FileNameStream;
    FileNameStream << OutputFile << step << ".vtk";
    ofstream out;
    out.open( (FileNameStream.str()).c_str() );
    int NumNodes = _myGelMesh->getNumberOfNodes();
  
    // Header
    char spchar = '#';
    out << spchar << " vtk DataFile Version 3.1" << endl;
    out << "vtk output" << endl;
    out << "ASCII" << endl;
    out << "DATASET UNSTRUCTURED_GRID" << endl;
    out << "POINTS " << NumNodes << " FLOAT" << endl;

    // For now we assumed dim == 3 !
    for (int i = 0; i < NumNodes; i++ ) {
      out << _myGelMesh->getX(i) << endl;
    }

    vector<GeomElement* > elements = _myGelMesh->getElements();
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
      cout << "Error cell type not implemented in MechanicsModel writeOutput. " << endl;
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
	<< "VECTORS displacements double" << endl;
  
    int dim = _myGelMesh->getDimension();
    for (int i = 0; i < NumNodes; i++ ) {
      VectorXd x = VectorXd::Zero(dim);
      for (int j = 0; j < dim; j++) {
	x(j) = _field[i*dim + j];
      }
      out << x - _myGelMesh->getX(i) << endl;
    }    

    // Residuals
    out << "VECTORS residual double" << endl;
  
    // Compute Residual
    uint PbDoF = ( _myGelMesh->getNumberOfNodes())*this->getDoFperNode();
    EigenResult myResults(PbDoF, 2);
    int myRequest = 2;
    myResults.setRequest(myRequest);
    this->compute(myResults);
    VectorXd R = *(myResults._residual);
    
    for (int i = 0; i < NumNodes; i++ ) {
      for (int j = 0; j < dim; j++) {
	out << R(i*dim + j) << endl;
      }
    }    

  

    // Cell data section
    // Alpha_1 material property
    out << endl << "CELL_DATA " << NumEl << endl
	<< "SCALARS Alpha1 double" << endl
	<< "LOOKUP_TABLE default" << endl;
    for ( int e = 0; e < NumEl; e++ ) {
      vector<Real > MatProp = _materials[e]->getMaterialParameters();
      out << MatProp[0] << endl;
    }    

    // Alpha_2 material property
    out << "SCALARS Alpha2 double" << endl
	<< "LOOKUP_TABLE default" << endl;
    for ( int e = 0; e < NumEl; e++ ) {
      vector<Real > MatProp = _materials[e]->getMaterialParameters();
      out << MatProp[1] << endl;
    }

    // Material internal variable
    out << "SCALARS InternalVariable double" << endl
	<< "LOOKUP_TABLE default" << endl;
    for ( int e = 0; e < NumEl; e++ ) {
      vector<Real > IntProp = _materials[e]->getInternalParameters();
      // for ( int i = 0; i < IntProp.size(); i++ ) {
 
//	     WARNING 
//	     THIS ONLY PRINTS THE FIRST INTERNAL VARIABLE
//	     BAD - NEED TO BE CHANGED!!
 
      if ( IntProp.size() > 0 ) {
	out << IntProp[0] << endl;
      } else {
	out << 0.0 << endl;
      }

	//}
    }
    
    // Material stress tensor
    out << "TENSORS P double" << endl;
    
      // Loop through elements, also through material points array, which is unrolled
      // uint index = 0;
      FilamentMaterial::FKresults FKres;
      FKres.request = 2;
      for(int e = 0; e < NumEl; e++)
      {
	GeomElement* geomEl = elements[e];
	const int numQP = geomEl->getNumberOfQuadPoints();
	Vector3d Fiber = geomEl->getFiber();
	
	// F at each quadrature point are computed at the same time in one element
	vector<Matrix3d > Flist(numQP, Matrix3d::Zero());
	// Compute deformation gradients for current element
	this->computeDeformationGradient(Flist, geomEl);
     
	// Loop over quadrature points
	// for(int q = 0; q < numQP; q++) {
	// _materials[e]->compute(FKres, Flist[q], &Fiber);

//	     WARNING 
//	     THIS ONLY WORKS WITH 1QP PER ELEMENT
//	     BAD - NEED TO BE CHANGED!!

	  // }
	_materials[e]->compute(FKres, Flist[0], &Fiber);
	out << FKres.P(0,0) << " " <<  FKres.P(0,1) << " " << FKres.P(0,2) << endl << 
	       FKres.P(1,0) << " " <<  FKres.P(1,1) << " " << FKres.P(1,2) << endl << 
	       FKres.P(2,0) << " " <<  FKres.P(2,1) << " " << FKres.P(2,2) << endl;
      }
      
    // Close file
    out.close();
    */
  } // writeOutput
    






} // namespace voom
