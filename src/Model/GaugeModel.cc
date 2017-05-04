#include "GaugeModel.h"
#include "LoopShellElement.h"
namespace voom {
  // Constructor
  GaugeModel::GaugeModel(Mesh* aMesh, vector<GaugeLipid * > materials, 
				 const uint NodeDoF):
    EllipticModel(aMesh, NodeDoF), _materials(materials)
  {
    // THERE IS ONE MATERIAL PER ELEMENT - CAN BE CHANGEd - DIFFERENT THAN BEFORE
    // Resize and initialize (default function) _field vector
    _field.resize(  (_myMesh->getNumberOfNodes() )*_nodeDoF );
    this->initializeField(0.001);

  }
  

  // Compute Function - Compute Energy, Force, Stiffness
  void GaugeModel::compute(EllipticResult &Res)
  {
    const vector<GeomElement* > elements = _myMesh->getElements();
    const int NumEl = elements.size();
    const int dim = _myMesh->getDimension();
    
    int PbDoF =Res.getPbDoF();
    VectorXd Residual;

    // Reset values in result struct
    if (Res.getRequest() & ENERGY ) {
     Res.setEnergy(0.0);
    }
    if (Res.getRequest() & FORCE ||Res.getRequest() & DMATPROP )  {
     Res.resetResidualToZero();
    }
    if (Res.getRequest() & STIFFNESS ) { 
     Res.resetStiffnessToZero();
    }
    // Loop through elements to compute area
    double area = 0;
    double RefArea = 8.06;
    double volume = 0;
    double gaugeConstX = 0;
    double gaugeConstY = 0;
    double gaugeConstZ = 0;
    double press = 0;
    double epsilonA = 1e-3;
    //double epsilonX = 1e-3;
    double epsilonY = 1e3;;
    //double epsilonZ = 1e-3;
    Real trX = 0.0;
    Real trY = 0.1;
    Real trZ = 0.0;
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
	  Vector3d R;
	  vector<Vector3d> Refa(2, Vector3d::Zero());
	  vector<Vector3d> RefaPartials(3, Vector3d::Zero());
	  Vector3d RefR;
	  
	  for(uint J = 0; J < 2; J++) { //J =0,1 (represents two surface coordinates
	    R = Vector3d::Zero();
	    RefR = Vector3d::Zero();
	    for(uint i = 0; i < dim; i++) { //i =0,1,2 represents three displacement indices
	      for(uint n = 0; n < numNodes; n++) {
		R(i) += _field[NodesID[n]*dim + i] * element->getN(q,n);
		a[J](i) += _field[NodesID[n]*dim + i] * element->getDN(q, n, J);
		aPartials[J](i) += _field[NodesID[n]*dim + i] * element->getDDN(q,n,J,J); //a_tt, a_pp

		RefR(i) += _myMesh->getX(NodesID[n],i) * element->getN(q,n);
		Refa[J](i) += _myMesh->getX(NodesID[n],i) * element->getDN(q, n, J);
		RefaPartials[J](i) += _myMesh->getX(NodesID[n],i) * element->getDDN(q,n,J,J); //a_tt, a_pp
	      }
	    }
	    
	  }	 
	  // a_tp component
	  for(uint i = 0; i < dim; i++) { //i =0,1,2 represents three displacement indices
	    for(uint n = 0; n < numNodes; n++) {
	      aPartials[2](i) += _field[NodesID[n]*dim + i]*element->getDDN(q,n,0,1); //a_tp
	      RefaPartials[2](i) += _myMesh->getX(NodesID[n],i) *element->getDDN(q,n,0,1); //a_tp
	    }
	  }
	  ShellGeometry geometry = ShellGeometry(a,aPartials);
	  ShellGeometry Refgeometry = ShellGeometry(Refa,RefaPartials);
	  Real metric = geometry.metric();
	  // Volume associated with QP q
	  Real da = element->getQPweights(q) * metric;
	  Real dA = element->getQPweights(q) * Refgeometry.metric();
	  area += da;
	  //RefArea += dA;
	  volume += da*(geometry.d()).dot(R)/3.;
	  //gaugeConstX += (RefR(0)-trX)*da;
	  gaugeConstY += (RefR(1)-trY)*da;
	  //gaugeConstZ += (RefR(2)-trZ)*da;
	}
      }
    // End: Total area computation
    // Loop through elements, also through material points array, which is unrolled
    GaugeLipid::Shellresults ShRes;
    ShRes.request =Res.getRequest();
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
	  Vector3d R;
	  vector<Vector3d> Refa(2, Vector3d::Zero());
	  vector<Vector3d> RefaPartials(3, Vector3d::Zero());
	  Vector3d RefR;
	  
	  for(uint J = 0; J < 2; J++) { //J =0,1 (represents two surface coordinates
	    R = Vector3d::Zero();
	    RefR = Vector3d::Zero();
	    for(uint i = 0; i < dim; i++) { //i =0,1,2 represents three displacement indices
	      for(uint n = 0; n < numNodes; n++) {
		R(i) += _field[NodesID[n]*dim + i] * element->getN(q,n);
		a[J](i) += _field[NodesID[n]*dim + i] * element->getDN(q, n, J);
		aPartials[J](i) += _field[NodesID[n]*dim + i] * element->getDDN(q,n,J,J); //a_tt, a_pp

		RefR(i) += _myMesh->getX(NodesID[n],i) * element->getN(q,n);
		Refa[J](i) += _myMesh->getX(NodesID[n],i) * element->getDN(q, n, J);
		RefaPartials[J](i) += _myMesh->getX(NodesID[n],i) * element->getDDN(q,n,J,J); //a_tt, a_pp
	      }
	    }
	    
	  }	 
	  // a_tp component
	  for(uint i = 0; i < dim; i++) { //i =0,1,2 represents three displacement indices
	    for(uint n = 0; n < numNodes; n++) {
	      aPartials[2](i) += _field[NodesID[n]*dim + i]*element->getDDN(q,n,0,1); //a_tp
	      RefaPartials[2](i) += _myMesh->getX(NodesID[n],i) *element->getDDN(q,n,0,1); //a_tp
	    }
	  }
	  ShellGeometry geometry = ShellGeometry(a,aPartials);
	  ShellGeometry Refgeometry = ShellGeometry(Refa,RefaPartials);
	  _materials[e]->setRefGeometry(Refgeometry);
	  Real metric = geometry.metric();
	  _materials[e]->compute(ShRes, geometry);
	  Matrix2d gab = geometry.metricTensorInverse();
	  // Volume associated with QP q
	  Real Vol = element->getQPweights(q) * metric;
	  Real RefVol = element->getQPweights(q) * Refgeometry.metric();
	  double gamma = (area-RefArea)/epsilonA;
	  
	  //double gammaGaugeX = gaugeConstX/epsilonX;
	  double gammaGaugeY = gaugeConstY*epsilonY;
	  //double gammaGaugeZ = gaugeConstZ/epsilonX;

	  // Compute energy
	  if (Res.getRequest() & ENERGY) {
	   Res.addEnergy(ShRes.W*Vol);
	    //Add work done due to pressure
	   Res.addEnergy(-press/3*Vol*(geometry.d()).dot(R));
	   
	  }
	  Vector3d tempResidual;
	  tempResidual = Vector3d::Zero();
	  // Compute Residual
	    if ( (Res.getRequest() & FORCE) ) {
	      for(uint n = 0; n < numNodes; n++) {
	      for(uint i = 0; i < dim; i++) {
		Real tempResidual = 0.0;
		for ( int alpha = 0; alpha < 2; alpha++){
		  // Stress Resultant part
		  tempResidual += (ShRes.n[alpha](i)
				   + (gamma+(RefR(1)-trY)*gammaGaugeY)*(geometry.aDual())[alpha](i)) *  element->getDN(q,n,alpha);
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
		      -  ShRes.m[alpha].dot(geometry.aDualPartials()[index])*element->getDN(q,n,beta)*geometry.d()(i);
		  }//beta loop
		  // ------------------ volume variation the alpha contribution
		  tempResidual += -press/3.0* (R.dot(geometry.d())*geometry.aDual()[alpha] - R.dot(geometry.aDual()[alpha])*geometry.d())(i)*element->getDN(q,n,alpha);
		}// alpha loop
		// Volume variation
		// ------------------ volume variation 
		tempResidual += -press/3.0*(geometry.d())(i) * element->getN(q,n);
		Res.addResidual(NodesID[n]*dim+i, tempResidual*Vol);
	      } //i loop
	    } // n loop
	    } // Internal force loop
	}
      }// Element loop
    // Add global energy integral at the end
    Res.addEnergy(pow((area-RefArea),2)/(2*epsilonA));
    //Res.addEnergy(pow((volume-Vbar),2)/(2*epsilonV)); 
    //Res.addEnergy(pow((gaugeConstX),2)/(2*epsilonX));
    Res.addEnergy(pow((gaugeConstY),2)/(2)*epsilonY);
    //Res.addEnergy(pow((gaugeConstZ),2)/(2*epsilonZ));
    
    
  } // Compute Mechanics Model  

  // Writing output
  void GaugeModel::writeOutputVTK(const string OutputFile, int step) 
  {
    // Create outputFile name
    stringstream FileNameStream;
    FileNameStream << OutputFile << step << ".vtk";
    ofstream out;
    out.open( (FileNameStream.str()).c_str() );
    int NumNodes = _myMesh->getNumberOfNodes(); //_myMesh->getNumberOfNodes()-8*2; //56*2
    // Header
    char spchar = '#';
    out << spchar << " vtk DataFile Version 3.1" << endl;
    out << "vtk output" << endl;
    out << "ASCII" << endl;
    out << "DATASET POLYDATA" << endl;
    out << "POINTS " << NumNodes << " FLOAT" << endl;

    int dim = 3;
    for (int i = 0; i < NumNodes; i++ ) {
      VectorXd x = VectorXd::Zero(dim);
      for (int j = 0; j < dim; j++) {
        x(j) = _field[i*dim + j];
      }
      out << (x).transpose()  << endl;
    }    

    int noOfElements =  _myMesh->getNumberOfElements(); 
    out << "POLYGONS " << noOfElements << " " 
        << noOfElements * 4 << endl;
    vector<GeomElement*> elements = _myMesh->getElements();
    for(int e=0; e< noOfElements; e++){
      vector<int> nodesID = elements[e]->getNodesID();
      out << "3 "<<nodesID[0]<< " " << nodesID[1] <<" "<< nodesID[2] << endl;
    }

    // Close file
    out.close();
  } // writeOutput
} // namespace voom
