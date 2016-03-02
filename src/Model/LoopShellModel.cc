#include "LoopShellModel.h"
#include "LoopShellElement.h"
namespace voom {

  // Constructor
  LoopShellModel::LoopShellModel(Mesh* aMesh, vector<SCElastic * > materials, 
				 const uint NodeDoF):
    Model(aMesh, NodeDoF), _materials(materials)
  {
    // THERE IS ONE MATERIAL PER ELEMENT - CAN BE CHANGEd - DIFFERENT THAN BEFORE
    // Resize and initialize (default function) _field vector
    _field.resize(  (_myMesh->getNumberOfNodes() )*_nodeDoF );
    this->initializeField();

  }
  

  // Compute Function - Compute Energy, Force, Stiffness
  void LoopShellModel::compute(Result & R)
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
    }
    if ( R.getRequest() & STIFFNESS ) { 
      R.resetStiffnessToZero();
      KtripletList.reserve(dim*dim*NumEl*AvgNodePerEl*AvgNodePerEl);
    }
   
    // Loop through elements, also through material points array, which is unrolled
    // uint index = 0;
    SCElastic::Shellresults ShRes;
    ShRes.request = R.getRequest();
    for(int e = 0; e < NumEl; e++)
      {
	LoopShellElement* element = dynamic_cast<LoopShellElement*>(elements[e]);
	const vector<int  >& NodesID = element->getNodesID();
	const int numQP    = element->getNumberOfQuadPoints();
	const int numNodes = NodesID.size();
	//MatrixXd Kele = MatrixXd::Zero(numNodes*dim, numNodes*dim);
      
	// Loop over quadrature points
	for(int q = 0; q < numQP; q++) {
	  // Initialize a and aPartials
	  vector<Vector3d> a(2, Vector3d::Zero());
	  vector<Vector3d> aPartials(3, Vector3d::Zero());
	  Vector3d temp;
	  Real theta, phi;
	  for(uint J = 0; J < 2; J++) { //J =0,1 (represents two surface coordinates
	    temp = Vector3d::Zero();
	    for(uint i = 0; i < dim; i++) { //i =0,1,2 represents three displacement indices
	      for(uint n = 0; n < numNodes; n++) {
		temp(i) += _field[NodesID[n]*dim + i] * element->getN(q,n);
		a[J](i) += _field[NodesID[n]*dim + i] * element->getDN(q, n, J);
		// aPartials[J](i) += _displacements[n](i) *
		// aPartials[J](i) += _displacements(i,n) * 
		aPartials[J](i) += _field[NodesID[n]*dim + i] * element->getDDN(q,n,J,J); //a_tt, a_pp
	      }
	    }
	    
	  }
	  
	  //cout << (normal/normal.norm()).dot( temp/temp.norm()) << " works" <<endl;
	  //cout << temp.norm() << " ... " << numNodes<< endl;
	  //cout << aPartials[1].dot(aPartials[2]) <<endl;
	  // a_tp component
	  for(uint i = 0; i < dim; i++) { //i =0,1,2 represents three displacement indices
	    for(uint n = 0; n < numNodes; n++) {
	      aPartials[2](i) += _field[NodesID[n]*dim + i]*element->getDDN(q,n,0,1); //a_tt, a_pp
	    }
	  }

	  // phi = atan2(temp(1),temp(2)); theta = atan2(sqrt(temp(0)*temp(0)+temp(1)*temp(1)),temp(2));
	  // //cout << "theta = " << theta * 180/M_PI << " phi = "<< phi*180/M_PI <<endl;

	  // a[0] << cos(phi)*cos(theta), sin(phi)*cos(theta), -sin(theta);  //a_theta
	  // a[1] << -sin(phi)*sin(theta), cos(phi)*sin(theta), 0;       //a_phi
	  
	  // aPartials[0] << -cos(phi)*sin(theta), -sin(phi)*sin(theta), -cos(theta); // a_tt
	  // aPartials[1] << -cos(phi)*sin(theta), -sin(phi)*sin(theta), 0;      // a_pp
	  // aPartials[2] << -sin(phi)*cos(theta), cos(phi)*cos(theta), 0;        // a_tp
    	  

	  ShellGeometry geometry = ShellGeometry(a,aPartials);
	  Real metric = geometry.metric();
	  _materials[e]->compute(ShRes, geometry);

	  // Volume associated with QP q
	  Real Vol = element->getQPweights(q) * metric;
	
	  // Compute energy
	  if (R.getRequest() & ENERGY) {
	    R.addEnergy(ShRes.W*Vol);
	  }
	  
	  Real MC = - 0.5 * ( geometry.aDual()[0].dot( geometry.dPartials()[0] ) + \
			      geometry.aDual()[1].dot( geometry.dPartials()[1]) );
	  
	  //if (numNodes == 11 && q==0) cout << MC << endl;
	  //cout << MC <<endl;
	  //cout << geometry.metric()/sin(theta) << endl;
	  //if ( abs(temp.norm()-1)*100>.1) cout << "element = " << e << "Nodes : "<<NodesID.size() << endl;
	  //if (numNodes == 11 ) cout << MC << " ..." <<endl;
	  // Compute Residual
	    if ( (R.getRequest() & FORCE) ) {	
	    for(uint n = 0; n < numNodes; n++) {
	      for(uint i = 0; i < dim; i++) {
		Real tempResidual = 0.0;
		for ( int alpha = 0; alpha < 2; alpha++){
		  // Stress Resultant part
		  tempResidual += ShRes.n[alpha](i) *  element->getDN(q,n,alpha);
		  // Moment Resultant part
		  for(int beta=0; beta<2; beta++) {
		    uint index = -1;
		    if (alpha==0 && beta == 0) index = 0;
		    else if (alpha==1 && beta ==1) index = 1;
		    else if (alpha==1 && beta ==0) index =2;
		    else if (alpha==0 && beta ==1) index =3;
		    else exit(-1);
		    tempResidual +=  -  ShRes.m[alpha].dot(geometry.aDual()[beta]) *
		      ( element->getDDN(q,n,alpha,beta)*geometry.d()(i) + element->getDN(q,n,beta)*geometry.dPartials()[alpha](i) )
		      -  ShRes .m[alpha].dot(geometry.aDualPartials()[index])*element->getDN(q,n,beta)*geometry.d()(i);
		  }//beta loop
		}// alpha loop
		R.addResidual(NodesID[n]*dim+i, tempResidual*Vol);
		
	      } //i loop
	    } // n loop
	  } // Internal force loop
	}
  

      }// Element loop
 
  } // Compute Mechanics Model




  // Writing output
  void LoopShellModel::writeOutputVTK(const string OutputFile, int step) 
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
      cout << "Error cell type not implemented in LoopShellModel writeOutput. " << endl;
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
