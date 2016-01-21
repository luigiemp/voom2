#include "PhaseModel.h"
#include "LoopShellElement.h"

namespace voom {

  // Constructor
  PhaseModel::PhaseModel(Mesh* aMesh, vector<CahnHilliard * > materials, 
				 const uint NodeDoF):
    EllipticModel(aMesh, NodeDoF), _materials(materials)
  {
    // THERE IS ONE MATERIAL PER ELEMENT - CAN BE CHANGEd - DIFFERENT THAN BEFORE
    // Resize and initialize (default function) _field vector
    _field.resize(  _myMesh->getNumberOfNodes()  );
    _lagMult = 0.0;
    this->initializeField(0.5-double(rand())/double(RAND_MAX));
    //this->initializeField();
  }
  

  // Compute Function - Compute Energy, Force, Stiffness
  void PhaseModel::compute(EllipticResult & R)
  {
    const vector<GeomElement* > elements = _myMesh->getElements();
    const int AvgNodePerEl = ((elements[0])->getNodesID()).size();
    const int NumEl = elements.size();
    const int dim = _myMesh->getDimension();
    vector<Triplet<Real > > KtripletList, HgtripletList;
    
    int PbDoF = R.getPbDoF();
    int TotNumMatProp = R.getNumMatProp();
    vector<VectorXd > dRdalpha;
    VectorXd Residual;
    //int NumPropPerMat = (_materials[0]->getMaterialParameters()).size(); // Assume all materials have the same number of material properties



    // Reset values in result struct
    if ( R.getRequest() & ENERGY ) {
      R.setEnergy(0.0);
    }
    if ( R.getRequest() & FORCE || R.getRequest() & DMATPROP )  {
      R.resetResidualToZero();
      R.resetConstraintToZero();
    }
    if ( R.getRequest() & STIFFNESS ) { 
      R.resetStiffnessToZero();
      KtripletList.reserve(dim*dim*NumEl*AvgNodePerEl*AvgNodePerEl);
    }
   
    // Loop through elements, also through material points array, which is unrolled
    // uint index = 0;
    CahnHilliard::Scalarresults ShRes;
    ShRes.request = R.getRequest();
    for(int e = 0; e < NumEl; e++)
      {
	LoopShellElement* element = dynamic_cast<LoopShellElement*>(elements[e]);
	const vector<int  >& NodesID = element->getNodesID();
	const int numQP    = element->getNumberOfQuadPoints();
	const int numNodes = NodesID.size();
      
	// Loop over quadrature points
	for(int q = 0; q < numQP; q++) {
	  // Initialize a and aPartials
	  vector<Vector3d> a(2, Vector3d::Zero());
	  vector<Vector3d> aPartials(3, Vector3d::Zero());
	  Vector3d temp;
	  Real phi;
	  Vector2d phiPartials = Vector2d::Zero();
	  for(uint J = 0; J < 2; J++) { //J =0,1 (represents two surface coordinates
	    temp = Vector3d::Zero();
	    phi = 0;
	    for(uint n = 0; n < numNodes; n++) {
	      for(uint i = 0; i < 3; i++) { //i =0,1,2 represents x,y,z
		temp(i) += _myMesh->getX(NodesID[n],i) * element->getN(q,n);
		a[J](i) += _myMesh->getX(NodesID[n],i) * element->getDN(q, n, J);
		aPartials[J](i) += _myMesh->getX(NodesID[n],i) * element->getDDN(q,n,J,J); //a_tt, a_pp
	      }
	      
	      phi += _field[NodesID[n]] * element->getN(q,n);
	      
	      phiPartials(J) += _field[NodesID[n]] * element->getDN(q,n,J);
	    }
	  }
	  
	  for(uint i = 0; i < dim; i++) { //i =0,1,2 represents three displacement indices
	    for(uint n = 0; n < numNodes; n++) {
	      aPartials[2](i) += _myMesh->getX(NodesID[n],i) *element->getDDN(q,n,0,1); //a_tt, a_pp
	    }
	  }
	  
	  ShellGeometry geometry = ShellGeometry(a,aPartials);
	  Real metric = geometry.metric();
	  
	  _materials[e]->compute(ShRes, phi, phiPartials, _lagMult, geometry);
	  
	  // Volume associated with QP q
	  Real Vol = element->getQPweights(q) * metric;
	  
	  // Compute energy
	  if (R.getRequest() & ENERGY) {
	    R.addEnergy(ShRes.W * Vol);
	  }
	  
	  // Compute Residual
	    if ( (R.getRequest() & FORCE) ) {	
	    
	    for(uint n = 0; n < numNodes; n++) {
	      Real tempResidual = 0.0;
	      tempResidual = ShRes.n * element->getN(q,n) +
		ShRes.m[0]*element->getDN(q,n,0) + ShRes.m[1]*element->getDN(q,n,1);
	      R.addResidual(NodesID[n], tempResidual*Vol);
	    } // n loop
	    R.addConstraint( ShRes.constr * Vol );
	    } // Internal force loop
	}
  

      }// Element loop
 
  } // Compute Mechanics Model




  // Writing output
  void PhaseModel::writeOutputVTK(const string OutputFile, int step) 
  {
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
      out << _myMesh->getX(i) << endl;
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
      cout << "Error cell type not implemented in PhaseModel writeOutput. " << endl;
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
    out << endl << "POINT_DATA " << NumNodes << endl
	<< "VECTORS displacements double" << endl;
  
    int dim = _myMesh->getDimension();
    for (int i = 0; i < NumNodes; i++ ) {
      VectorXd x = VectorXd::Zero(dim);
      for (int j = 0; j < dim; j++) {
	x(j) = _field[i*dim + j];
      }
      out << x - _myMesh->getX(i) << endl;
    }    

    // Close file
    out.close();
  } // writeOutput







} // namespace voom
