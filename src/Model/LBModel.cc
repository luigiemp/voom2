#include "LBModel.h"
#include "LoopShellElement.h"

namespace voom {

  // Constructor
  LBModel::LBModel(Mesh* aMesh, vector<LandauBrazovskii * > materials, 
				 const uint NodeDoF):
    EllipticModel(aMesh, NodeDoF), _materials(materials)
  {
    // THERE IS ONE MATERIAL PER ELEMENT - CAN BE CHANGEd - DIFFERENT THAN BEFORE
    // Resize and initialize (default function) _field vector
    _field.resize(  _myMesh->getNumberOfNodes()  );
    //this->initializeField( (-double(rand())/double(RAND_MAX)) );
    this->initializeField();
    _implicitDynamics = false; //default
  }
  

  // Compute Function - Compute Energy, Force, Stiffness
  void LBModel::compute(EllipticResult & R)
  {
    const vector<GeomElement* > elements = _myMesh->getElements();
    const int AvgNodePerEl = ((elements[0])->getNodesID()).size();
    const int NumEl = elements.size();
    const int dim = _myMesh->getDimension();
    int PbDoF = R.getPbDoF();
    int TotNumMatProp = R.getNumMatProp();
    
    // Reset values in result struct
    if ( R.getRequest() & ENERGY ) {
      R.setEnergy(0.0);
    }
    if ( R.getRequest() & FORCE || R.getRequest() & DMATPROP )  {
      R.resetResidualToZero();
    }
    if ( R.getRequest() & STIFFNESS ) { 
      
    }
   
    // Loop through elements, also through material points array, which is unrolled
    LandauBrazovskii::Scalarresults ShRes;
    ShRes.request = R.getRequest();
    
#ifdef _OPENMP  
#pragma omp parallel for schedule(static) firstprivate(ShRes)
#endif  
    for(int e = 0; e < NumEl; e++)
      {
	//cout << e << endl;
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
	  Matrix2d phi2Partials = Matrix2d::Zero();
	  Real phiPrev;
	  	  	    
	  for(uint J = 0; J < 2; J++) { //J =0,1 (represents two surface coordinates
	    temp = Vector3d::Zero();
	    phi = 0;
	    phiPrev = 0;
	    for(uint n = 0; n < numNodes; n++) {
	      for(uint i = 0; i < 3; i++) { //i =0,1,2 represents x,y,z
		temp(i) += _myMesh->getX(NodesID[n],i) * element->getN(q,n);
		a[J](i) += _myMesh->getX(NodesID[n],i) * element->getDN(q, n, J);
		aPartials[J](i) += _myMesh->getX(NodesID[n],i) * element->getDDN(q,n,J,J); //a_tt, a_pp
	      }
	      
	      phi += _field[NodesID[n]] * element->getN(q,n);
	      if (_implicitDynamics == true){
		// Implict time stepping ***************************************************** XXXXXXXXXXXXXXXXXXXXXXX
		phiPrev += _prevField[NodesID[n]] * element->getN(q,n); 
	      }
	      phiPartials(J) += _field[NodesID[n]] * element->getDN(q,n,J);
	      phi2Partials(J,J) += _field[NodesID[n]] * element->getDDN(q,n,J,J);
	    }
	  }
	  
	  for(uint n = 0; n < numNodes; n++) {
	    for(uint i = 0; i < 3; i++) { //i =0,1,2 represents three displacement indices
	      aPartials[2](i) += _myMesh->getX(NodesID[n],i) *element->getDDN(q,n,0,1); //a_tp
	    }
	    phi2Partials(0,1) += _field[NodesID[n]] * element->getDDN(q,n,0,1);
	  }
	  phi2Partials(1,0) = phi2Partials(0,1);
	  ShellGeometry geometry = ShellGeometry(a,aPartials);
	  Real metric = geometry.metric();
	  Matrix2d gab = geometry.metricTensorInverse();
	  Vector2d GammaI = geometry.gklGammaI_kl();
	  /*
	  //**************** Check *********************************************
	  Real Phi = atan2(temp(1),temp(0)); Real theta = atan2(sqrt(temp(0)*temp(0)+temp(1)*temp(1)),temp(2));
	  Real X = temp(0); Real Y=temp(1); Real Z=temp(2);
	  // Real theta_u = (1/(X*X + Y*Y + Z*Z )) *
	  //                (1/sqrt(X*X+Y*Y))*
	  //                (Z*X*a[0](0) + Z*Y*a[0](1) - (X*X+Y*Y)*a[0](2));
	  // Real theta_v = (1/(X*X + Y*Y + Z*Z )) *
	  //                (1/sqrt(X*X+Y*Y))*
	  //                (Z*X*a[1](0) + Z*Y*a[1](1) - (X*X+Y*Y)*a[1](2));
	  
	  Real theta_u = -a[0](2)/sin(theta);
	  Real theta_v = -a[1](2)/sin(theta);
	  
	  Real theta_uu = -(aPartials[0](2) + cos(theta)*theta_u*theta_u)/sin(theta);
	  Real theta_vv = -(aPartials[1](2) + cos(theta)*theta_v*theta_v)/sin(theta);
	  Real theta_uv = -(aPartials[2](2) + cos(theta)*theta_u*theta_v)/sin(theta);
	  
	  Matrix2d D2theta;
	  D2theta << theta_uu, theta_uv, theta_uv, theta_vv;
	  //cout << "error = " << theta_u-Theta_u << " and " <<theta_v-Theta_v << endl;
	  Real phi_u = 1/(X*X+Y*Y)*(X*a[0](1)-Y*a[0](0));
	  Real phi_v = 1/(X*X+Y*Y)*(X*a[1](1)-Y*a[1](0));
	  
	  Real phi_uu = ((X*X+Y*Y)*(X*aPartials[0](1) - Y*aPartials[0](0))-2*(X*a[0](1)-Y*a[0](0))*(X*a[0](0)+Y*a[0](1)))/pow(X*X+Y*Y,2);
	  Real phi_vv = ((X*X+Y*Y)*(X*aPartials[1](1) - Y*aPartials[1](0))-2*(X*a[1](1)-Y*a[1](0))*(X*a[0](0)+Y*a[1](1)))/pow(X*X+Y*Y,2);
	  Real phi_uv = ((X*X+Y*Y)*(X*aPartials[2](1) + a[1](0)*a[0](1) - Y*aPartials[2](0) - a[1](1)*a[0](0)) -
			 2*(X*a[0](1)-Y*a[0](0))*(X*a[1](0)+Y*a[1](1)))/(pow(X*X+Y*Y,2));
	  Matrix2d D2phi;
	  D2phi << phi_uu, phi_uv, phi_uv, phi_vv;

	  Matrix2d Jac;
	  Jac << theta_u, phi_u,  theta_v, phi_v;
	  Vector2d phi_tp;
	  phi_tp = (Jac.inverse())*phiPartials;
	  
	  Matrix2d phi_2;
	  phi_2 = Jac.inverse() *(phi2Partials -phi_tp(0)*D2theta - phi_tp(1)*D2phi) *(Jac.transpose()).inverse();
	  
	  Real Lappsi = phi_2(0,0)+cos(theta)/sin(theta)*phi_tp(0)+(1/pow(sin(theta),2))*phi_2(1,1);
	 
	  // if (abs(Lappsi/phi+20) > 3.0){
	  //   cout << "phi_tp: " << phi_tp.transpose() << endl;
	  //   cout << "phi_2: " << endl << phi_2 << endl << endl;
	  //   cout <<"Lappsi/psi : " << Lappsi/phi << endl;
	  //   cout << "cot term : " << cos(theta)/sin(theta)*phi_tp(0) << " theta: " << theta*180/M_PI<< " phi: " << Phi*180/M_PI << endl;
	  //   cout << numNodes << endl;
	  // }
	  */
	  _materials[e]->compute(ShRes, phi, phiPartials, phi2Partials, geometry);
	  // Volume associated with QP q
	  Real Vol = element->getQPweights(q) * metric;
	  
	  Real dt = 0.01;
	  // Compute energy

	  if (R.getRequest() & ENERGY) {
	    if (_implicitDynamics == true){
	      // Implict time stepping ***************************************************** XXXXXXXXXXXXXXXXXXXXXXX
	      ShRes.W = ShRes.W + dt/2.0*pow( (phi-phiPrev)/dt, 2.0 );
	    }

	    R.addEnergy(ShRes.W * Vol);
	  }
	  
	  // Compute Residual
	  if ( (R.getRequest() & FORCE) ) {	
	      Real tempResidual = 0.0;
	      for(uint n = 0; n < numNodes; n++) {
		
		tempResidual = ShRes.n1 * element->getN(q,n) +
						     ShRes.n2(0)*element->getDN(q,n,0) + ShRes.n2(1)*element->getDN(q,n,1) + 
						     ShRes.n3 * ( gab(0,0) * element->getDDN(q,n,0,0)  + 
								  gab(0,1) * element->getDDN(q,n,0,1) +
								  gab(1,0) * element->getDDN(q,n,1,0) +
								  gab(1,1) * element->getDDN(q,n,1,1) +
								  - GammaI(0)*element->getDN(q,n,0) 
								  - GammaI(1)*element->getDN(q,n,1) );
		if (_implicitDynamics == true){
		  //Implicit time stepping ************************************ XXXXXXXXXXXXXXXX
		  tempResidual += (phi-phiPrev)/(dt)*element->getN(q,n);
		}
		R.addResidual(NodesID[n], tempResidual*Vol);
	    } // n loop
	      //cout << endl << e << " " << numNodes << " " << output << endl;	      
	    } // Internal force loop
	}
  

      }// Element loop
  } // Compute Mechanics Model




  // Writing output
  void LBModel::writeOutputVTK(const string OutputFile, int step) 
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
    out << "DATASET POLYDATA" << endl;
    out << "POINTS " << NumNodes << " FLOAT" << endl;

    // For now we assumed dim == 3 !
    for (int i = 0; i < NumNodes; i++ ) {
      out << (_myMesh->getX(i)).transpose() << endl;
    }

    int noOfElements = _myMesh->getNumberOfElements();
    out << "POLYGONS " << noOfElements << " " 
	<< noOfElements * 4 << endl;
    vector<GeomElement*> elements = _myMesh->getElements();
    for(int e=0; e< noOfElements; e++){
      vector<int> nodesID = elements[e]->getNodesID();
      out << "3 "<<nodesID[0]<< " " << nodesID[1] <<" "<< nodesID[2] << endl;
    }

    out << endl << "POINT_DATA " << NumNodes << endl
	<< "SCALARS Field double" << endl
        << "LOOKUP_TABLE default" << endl;
  
    int dim = 1;
    for (int i = 0; i < NumNodes; i++ ) {
      VectorXd x = VectorXd::Zero(dim);
      for (int j = 0; j < dim; j++) {
	x(j) = _field[i*dim + j];
      }
      out << x  << endl;
    }    

    // Close file
    out.close();
  } // writeOutput
} // namespace voom
