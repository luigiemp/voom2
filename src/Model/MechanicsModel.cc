#include "MechanicsModel.h"

namespace voom {

  // Constructor
  MechanicsModel::MechanicsModel(Mesh* aMesh, vector<MechanicsMaterial * > materials, 
				 const uint NodeDoF,
				 int PressureFlag, Mesh* SurfaceMesh,
				 int NodalForcesFlag,
				 int ResetFlag):
    EllipticModel(aMesh, NodeDoF), _materials(materials), 
    _pressureFlag(PressureFlag), _pressure(0.0), _surfaceMesh(SurfaceMesh),
    _nodalForcesFlag(NodalForcesFlag), _forcesID(NULL), _forces(NULL), _resetFlag(ResetFlag)
  {
    // THERE IS ONE MATERIAL PER ELEMENT - CAN BE CHANGED - DIFFERENT THAN BEFORE
    // Resize and initialize (default function) _field vector
    _field.resize(  (_myMesh->getNumberOfNodes() )*_nodeDoF );
    this->initializeField();

    if (_pressureFlag == 1) {
      _prevField.resize( _field.size() );
      this->setPrevField();
    }

  }
  






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
	for(uint i = 0; i < dim; i++) 
	  for(uint J = 0; J < dim; J++)
	    for(uint a = 0; a < nodeNum; a++)
	      Flist[q](i,J) += 
		_field[NodesID[a]*dim + i] * geomEl->getDN(q, a, J);

    } // loop over quadrature points

  }







  // Compute Function - Compute Energy, Force, Stiffness
  void MechanicsModel::compute(EllipticResult & R)
 {
    const vector<GeomElement* > elements = _myMesh->getElements();
    const int AvgNodePerEl = ((elements[0])->getNodesID()).size();
    const int NumEl = elements.size();
    const int dim = _myMesh->getDimension();
    vector<Triplet<Real > > KtripletList, HgtripletList;
    
    int PbDoF = R.getPbDoF();
    int TotNumMatProp = R.getNumMatProp();
    vector<VectorXd > dRdalpha;
    int NumPropPerMat = (_materials[0]->getMaterialParameters()).size(); // Assume all materials have the same number of material properties


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
    
    if ( R.getRequest() & DMATPROP ) {
      dRdalpha.assign( TotNumMatProp, VectorXd::Zero(PbDoF) );
      if ( _resetFlag == 1 ) {
	R.resetGradgToZero();
	R.resetHgToZero();
	HgtripletList.reserve(TotNumMatProp*TotNumMatProp);
      }
    }
    
    

    // Loop through elements, also through material points array, which is unrolled
    // uint index = 0;
    MechanicsMaterial::FKresults FKres;
    FKres.request = R.getRequest();
    for(int e = 0; e < NumEl; e++)
    {
      GeomElement* geomEl = elements[e];
      const vector<int  >& NodesID = geomEl->getNodesID();
      const int numQP    = geomEl->getNumberOfQuadPoints();
      Vector3d Fiber = geomEl->getFiber();
      const int numNodes = NodesID.size();
      MatrixXd Kele = MatrixXd::Zero(numNodes*dim, numNodes*dim);
      
      // F at each quadrature point are computed at the same time in one element
      vector<Matrix3d > Flist(numQP, Matrix3d::Zero());
      // Compute deformation gradients for current element
      this->computeDeformationGradient(Flist, geomEl);
     
      // Loop over quadrature points
      for(int q = 0; q < numQP; q++) {
	_materials[e]->compute(FKres, Flist[q], &Fiber);

	// Volume associated with QP q
	Real Vol = geomEl->getQPweights(q);

	// Compute energy
	if (R.getRequest() & ENERGY) {
	  R.addEnergy(FKres.W*Vol);
	}
 
	// Compute Residual
	if ( (R.getRequest() & FORCE) || (R.getRequest() & DMATPROP) ) {
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
		  // R.addStiffness(NodesID[a]*dim + i, NodesID[b]*dim + j, tempStiffness);
		  Kele(a*dim + i, b*dim + j) += tempStiffness;
		} // j loop
	      } // b loop
	    } // i loop
	  } // a loop
	} // Compute stiffness matrix

	if ( R.getRequest() & DMATPROP ) { 
	  // Assemble dRdalpha
	  for (uint alpha = 0; alpha < NumPropPerMat; alpha++) {
	    for(uint a = 0; a < numNodes; a++) {
	      for(uint i = 0; i < dim; i++) {
		Real tempdRdalpha = 0.0;
		for (uint J = 0; J < dim; J++) { 
		  tempdRdalpha += FKres.Dmat.get(alpha,i,J) * geomEl->getDN(q, a, J);
		} // J loop
		tempdRdalpha *= Vol;
		(dRdalpha[ (_materials[e]->getMatID())*NumPropPerMat + alpha])( NodesID[a]*dim + i) += tempdRdalpha;
	      } // i loop
	    } // a loop
	  } // alpha loop
	  
	} // Compute DMATPROP
	
      } // QP loop

      if ( R.getRequest() & STIFFNESS ) { 
	// Transform in triplets Kele
	for(uint a = 0; a < numNodes; a++) {
	  for(uint i = 0; i < dim; i++) {
	    for(uint b = 0; b < numNodes; b++) {
	      for(uint j = 0; j < dim; j++) {
		KtripletList.push_back( Triplet<Real >( NodesID[a]*dim + i, NodesID[b]*dim + j, Kele(a*dim + i, b*dim + j) ) );
	      }
	    }
	  }
	}
      }

    } // Element loop

    // Sum up all stiffness entries with the same indices
    if ( R.getRequest() & STIFFNESS ) {
      R.setStiffnessFromTriplets(KtripletList);
      R.FinalizeGlobalStiffnessAssembly(); 
      // cout << "Stiffness assembled" << endl;
    }
 
    // Add pressure/external load if any
    if (_pressureFlag == 1) {
      this->applyPressure(R);
    }
 
    // Add external load if any
    if (_nodalForcesFlag == 1 && (R.getRequest() & FORCE || R.getRequest() & DMATPROP) ) {
      for (int i = 0; i < _forcesID->size(); i++) {
	R.addResidual( (*_forcesID)[i], (*_forces)[i] );
      }
    }
 
    // Compute Gradg and Hg
    if ( R.getRequest() & DMATPROP ) { 

      // First extract residual
      VectorXd Residual = VectorXd::Zero(PbDoF);
      for (int i = 0; i < PbDoF; i++) {
	Residual(i) = R.getResidual(i); // local copy
      }

      if (_resetFlag == 1) {
	for (uint alpha = 0; alpha < TotNumMatProp; alpha++) {
	  R.addGradg(alpha, 2.0*dRdalpha[alpha].dot(Residual) );
	  for (uint beta = 0; beta < TotNumMatProp; beta++) {
	    // !!!
	    // WARNING : assume W is linear in alpha!!!!!
	    // !!!
	    HgtripletList.push_back( Triplet<Real >( alpha, beta,  2.0*dRdalpha[alpha].dot(dRdalpha[beta]) ) );
	  } // alpha loop
	} // beta loop
	R.setHgFromTriplets(HgtripletList);
      }
      else {
	for (uint alpha = 0; alpha < TotNumMatProp; alpha++) {
	  R.addGradg(alpha, 2.0*dRdalpha[alpha].dot(Residual) );
	  for (uint beta = 0; beta < TotNumMatProp; beta++) {
	    // !!!
	    // WARNING : assume W is linear in alpha!!!!!
	    // !!!
	    R.addHg(alpha, beta,  2.0*dRdalpha[alpha].dot(dRdalpha[beta]) );
	  } // alpha loop
	} // beta loop
      }
      



    } // Compute Gradg and Hg
    
 } // Compute Mechanics Model
  






  void MechanicsModel::applyPressure(EllipticResult & R) {

    const vector<GeomElement* > elements = _surfaceMesh->getElements();

    // Loop through elements
    for(int e = 0; e < elements.size(); e++)
    {
      GeomElement* geomEl = elements[e];
      const vector<int  >& NodesID = geomEl->getNodesID();
      const int numQP    = geomEl->getNumberOfQuadPoints();
      const int numNodes = NodesID.size();
      
      // Loop over quadrature points
      for(int q = 0; q < numQP; q++) {

	// Compute normal based on _prevField and displacement
	Vector3d a1 = Vector3d::Zero(), a2 = Vector3d::Zero(), a3 = Vector3d::Zero(), u = Vector3d::Zero();
	for (int a = 0; a < NodesID.size(); a++) {
	  int nodeID = NodesID[a];
	  Vector3d xa_prev, xa_curr;
	  xa_prev << _prevField[nodeID*3], _prevField[nodeID*3+1], _prevField[nodeID*3+2];
	  xa_curr << _field[nodeID*3], _field[nodeID*3+1], _field[nodeID*3+2];
	  a1 += xa_prev*geomEl->getDN(q, a, 0);
	  a2 += xa_prev*geomEl->getDN(q, a, 1);

	  u += (xa_curr - _surfaceMesh->getX(nodeID))*geomEl->getN(q, a);
	}
	a3 = a1.cross(a2);

	// Surface associated with QP q
	Real Area = a3.norm();
	a3 /= Area;
	Area *= geomEl->getQPweights(q);
	// cout << "Area = " << Area << endl;

	// Compute energy
	if (R.getRequest() & ENERGY) { 
	  R.addEnergy( _pressure*Area*a3.dot(u) );  
	}
       
	// Compute Residual
	if ( (R.getRequest() & FORCE) || (R.getRequest() & DMATPROP) ) {	
	  for(uint a = 0; a < NodesID.size(); a++) {
	    for(uint i = 0; i < 3; i++) {
	      R.addResidual(NodesID[a]*3+i, _pressure * Area * a3(i) * geomEl->getN(q, a));
	      // cout << NodesID[a] << " " << geomEl->getN(q, a) << endl;
	    } // i loop
	  } // a loop
	} // Internal force loop
    
      } // loop over QP
    } // loop over elements

  } // apply pressure







  void MechanicsModel::checkDmat(EigenEllipticResult & R, Real perturbationFactor, Real hM, Real tol) 
  {

    // Perturb initial config - gradg and Hg are zero at F = I
    const uint nodeNum   = _myMesh->getNumberOfNodes();
    const uint nLocalDoF = nodeNum*_nodeDoF;
    
    // Perturb field randomly to change from reference configuration
    // Save perturbed field to set the configuration back to reference 
    // at the end of the test
    vector<Real > perturb(nLocalDoF, 0.0);
    srand( time(NULL) );
    for(uint a = 0; a < nodeNum; a++) {
      for(uint i = 0; i < _nodeDoF; i++) {
    	Real randomNum =  perturbationFactor*(Real(rand())/RAND_MAX - 0.5);
    	perturb[a*_nodeDoF + i] = randomNum;
    	this->linearizedUpdate(a, i, randomNum);
      }
    }  
    
    Real error = 0.0, norm = 0.0;

    set<MechanicsMaterial *> UNIQUEmaterials;
    for (uint i = 0; i < _materials.size(); i++) 
      UNIQUEmaterials.insert(_materials[i]);
    
    cout << "Number of unique materials = " <<  UNIQUEmaterials.size() << endl;
    R.setRequest(8); // First compute gradg and Hg numerically
    this->compute(R);

    // Test gradg //
    R.setRequest(2); // Reset result request so that only forces are computed
    for ( set<MechanicsMaterial *>::iterator itMat =  UNIQUEmaterials.begin();
	  itMat != UNIQUEmaterials.end(); itMat++ )
    {
      vector<Real > MatProp = (*itMat)->getMaterialParameters();
      int MatID = (*itMat)->getMatID();
      for(int m = 0; m < MatProp.size(); m++) {
	// Perturb +hM the material property alpha
	MatProp[m] += hM;
	// Reset matProp in the materials with MatProp[m]
	(*itMat)->setMaterialParameters(MatProp);	  
	// Compute R
	this->compute(R);
	Real RTRplus = (R._residual)->dot(*R._residual);
	// Perturb -2hM the material property alpha
	MatProp[m] -= 2.0*hM;
	// Reset matProp in the materials with MatProp[m]
	(*itMat)->setMaterialParameters(MatProp);
	// Compute R
	this->compute(R);
	Real RTRminus = (R._residual)->dot(*R._residual);
	
	// Bring back to original value of alpha
	MatProp[m] += hM;
	// Reset matProp in all the materials with MatProp[m]
	(*itMat)->setMaterialParameters(MatProp);
	
	error += pow( (RTRplus-RTRminus)/(2.0*hM) - 
		      (*R._Gradg)( MatID*MatProp.size() + m), 2.0 );
	norm += pow( (*R._Gradg)(  MatID*MatProp.size() + m), 2.0 );
	
	// cout << (*R._Gradg)( MatID*MatProp.size() + m) << " " << (RTRplus-RTRminus)/(2.0*hM) << endl; 
      } // Loop over m
    } // Loop over unique materials
    error = sqrt(error);
    norm  = sqrt(norm);
    
    if ( abs(error) < norm * tol) {
      cout << "** Gradg consistency check PASSED" << endl;
      cout << "** Error: " << error << " Norm: " << norm  << " Norm*tol: " << norm*tol << endl;
    }
    else {
      cout << "** Gradg consistency check FAILED" << endl;
      cout << "** Error: " << error << " Norm: " << norm << " Norm*tol: " << norm*tol << endl;
    }


    
    // Test Hg //
    error = 0.0; norm = 0.0;
    R.setRequest(8);
    this->compute(R);
    SparseMatrix<Real > HgAn = *R._Hg;
    for ( set<MechanicsMaterial *>::iterator itMatA =  UNIQUEmaterials.begin(); itMatA != UNIQUEmaterials.end(); itMatA++ )
    {
      vector<Real > MatPropA = (*itMatA)->getMaterialParameters();
      int MatIDA = (*itMatA)->getMatID();

      for(int mA = 0; mA < MatPropA.size(); mA++) {
	
	for ( set<MechanicsMaterial *>::iterator itMatB =  UNIQUEmaterials.begin(); itMatB != UNIQUEmaterials.end(); itMatB++ ) 
	{
	    vector<Real > MatPropB = (*itMatB)->getMaterialParameters();
	    int MatIDB = (*itMatB)->getMatID();
	    
	    for(int mB = 0; mB < MatPropB.size(); mB++) {

	      // Perturb +hM the material property alpha
	      MatPropB[mB] += hM;
	      // Reset matProp in the materials with MatProp[m]
	      (*itMatB)->setMaterialParameters(MatPropB);	  
	      // Compute R
	      this->compute(R);
	      Real GradPlus = R.getGradg( MatIDA*MatPropA.size() + mA);
	
	      // Perturb -2hM the material property alpha
	      MatPropB[mB] -= 2.0*hM;
	      // Reset matProp in the materials with MatProp[m]
	      (*itMatB)->setMaterialParameters(MatPropB);
	      // Compute R
	      this->compute(R);
	      Real GradMinus = R.getGradg( MatIDA*MatPropA.size() + mA);
	
	      // Bring back to original value of alpha
	      MatPropB[mB] += hM;
	      // Reset matProp in all the materials with MatPropB[m]
	      (*itMatB)->setMaterialParameters(MatPropB);
	
	      error += pow( (GradPlus - GradMinus)/(2.0*hM) - 
			    HgAn.coeff( MatIDA*MatPropA.size() + mA, MatIDB*MatPropB.size() + mB ), 2.0);
	      norm += pow( HgAn.coeff( MatIDA*MatPropA.size() + mA, MatIDB*MatPropB.size() + mB ) , 2.0 );
	
	      // cout << (*R._Gradg)( MatID*MatProp.size() + m) << " " << (RTRplus-RTRminus)/(2.0*hM) << endl;

	    } // Loop over mB
	} // Loop over unique materials B
      } // Loop over mA
    } // Loop over unique materials A
    error = sqrt(error);
    norm  = sqrt(norm);
    
    if ( abs(error) < norm * tol) {
      cout << "** Hg consistency check PASSED" << endl;
      cout << "** Error: " << error << " Norm: " << norm  << " Norm*tol: " << norm*tol << endl;
    }
    else {
      cout << "** Hg consistency check FAILED" << endl;
      cout << "** Error: " << error << " Norm: " << norm << " Norm*tol: " << norm*tol << endl;
    }
    

    
    // Reset field to initial values
    for(uint a = 0; a < nodeNum; a++) {
      for(uint i = 0; i < _nodeDoF; i++) {
	this->linearizedUpdate(a, i, -perturb[a*_nodeDoF + i]);
      }
    } 
  
  } // Check consistency of gradg and Hg - checkDmat







  // Writing output
  void MechanicsModel::writeOutputVTK(const string OutputFile, int step) 
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
  
    int dim = _myMesh->getDimension();
    for (int i = 0; i < NumNodes; i++ ) {
      VectorXd x = VectorXd::Zero(dim);
      for (int j = 0; j < dim; j++) {
	x(j) = _field[i*dim + j];
      }
      out << x - _myMesh->getX(i) << endl;
    }    

    // Residuals
    out << "VECTORS residual double" << endl;
  
    // Compute Residual
    uint PbDoF = ( _myMesh->getNumberOfNodes())*this->getDoFperNode();
    EigenEllipticResult myResults(PbDoF, 2);
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

    // Close file
    out.close();
  } // writeOutput







} // namespace voom
