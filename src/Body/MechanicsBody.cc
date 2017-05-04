#include "MechanicsBody.h"

namespace voom {

  // Constructor
  MechanicsBody::MechanicsBody(Mesh* myMesh, const int NodeDoF,
			       vector<MechanicsMaterial * > Materials,
			       Result* R):
    Body(myMesh, NodeDoF), _materials(Materials)
    {
      // Initialize _field vector for the part corresponding to this Body
      this->initializeField(R);
    } // Constructor



  // Compute deformation gradient
  void MechanicsBody::computeDeformationGradient(vector<Matrix3d > & Flist, GeomElement* geomEl, Result* R)
  {
    // Compute F at all quadrature points
    const int numQP = geomEl->getNumberOfQuadPoints();
    const int dim   = _myMesh->getDimension();

    const vector<int > & NodesID = geomEl->getNodesID();
    const int nodeNum = NodesID.size();

    for(int q = 0; q < numQP; q++) {
      // Initialize F to zero
      Flist[q] = Matrix3d::Zero();
      for(int i = 0; i < dim; i++) {
	for(int J = 0; J < dim; J++) {
	  for(int a = 0; a < nodeNum; a++) {
	    Flist[q](i,J) += R->getField(NodesID[a]*dim + i) * geomEl->getDN(q, a, J);
	  }
	}
      }
    } // loop over quadrature points
  };  // computeDeformationGradient





  void MechanicsBody::computeGreenLagrangianStrainTensor(vector<Matrix3d> & Elist, GeomElement* geomEl, Result* R) 
  {
    const int numQP = geomEl->getNumberOfQuadPoints();
    vector<Matrix3d > Flist(numQP, Matrix3d::Zero());
    this->computeDeformationGradient(Flist, geomEl, R);

    for (int q = 0; q < numQP; q++) {
      Elist[q] = 0.5 * (Flist[q].transpose() * Flist[q] - Matrix3d::Identity());
    }
  };  // computeGreenLagrangianStrainTensor





  // Compute Function - Compute Energy, Force, Stiffness
  void MechanicsBody::compute(Result* R)
  {
    const vector<GeomElement* > elements = _myMesh->getElements();
    const int AvgNodePerEl = ((elements[0])->getNodesID()).size(); // It is just to pre-allocate memory space
    const int NumEl = elements.size();
    const int dim = _myMesh->getDimension();
    

    int PbDoF = R->getPbDoF();
    int TotNumMatProp = R->getNumMatProp();
    vector<VectorXd > dRdalpha;
    int NumPropPerMat = (_materials[0]->getMaterialParameters()).size(); // Assume all materials have the same number of material properties


    // Loop through elements, also through material points array, which is unrolled
    // int index = 0;
    MechanicsMaterial::FKresults FKres;
    FKres.request = R->getRequest();
    for(int e = 0; e < NumEl; e++)
    {
      GeomElement* geomEl = elements[e];
      const vector<int  >& NodesID = geomEl->getNodesID();
      const int numQP    = geomEl->getNumberOfQuadPoints();
      const int numNodes = NodesID.size();
      MatrixXd Kele = MatrixXd::Zero(numNodes*dim, numNodes*dim);

      // F at each quadrature point are computed at the same time in one element
      vector<Matrix3d > Flist(numQP, Matrix3d::Zero());
      // Compute deformation gradients for current element
      this->computeDeformationGradient(Flist, geomEl, R);

      // Loop over quadrature points
      for(int q = 0; q < numQP; q++) {
        _materials[e*numQP + q]->compute(FKres, Flist[q]);

        // Volume associated with QP q
        Real Vol = geomEl->getQPweights(q);

        // Compute energy
        if (R->getRequest() & ENERGY) {
	  R->addEnergy(FKres.W*Vol);
        }

        // Compute Residual
        if ( (R->getRequest() & FORCE) || (R->getRequest() & DMATPROP) ) {
          for(int a = 0; a < numNodes; a++) {
            for(int i = 0; i < dim; i++) {
              Real tempResidual = 0.0;
              for (int J = 0; J < dim; J++) {
                tempResidual += FKres.P(i,J) * geomEl->getDN(q, a, J);
              } // J loop
              tempResidual *= Vol;
              R->addResidual(NodesID[a]*dim+i, tempResidual);
            } // i loop
          } // a loop
        } // Internal force loop

        // Compute Stiffness
        if ( R->getRequest() & STIFFNESS ) {
          for(int a = 0; a < numNodes; a++) {
            for(int i = 0; i < dim; i++) {
              for(int b = 0; b < numNodes; b++) {
                for(int j = 0; j < dim; j++) {
                  Real tempStiffness = 0.0;
                  for(int M = 0; M < dim; M++) {
                    for(int N = 0; N < dim; N++) {
                      tempStiffness += FKres.K.get(i, M, j, N)*elements[e]->getDN(q, a, M)*
                      elements[e]->getDN(q, b, N);
                    } // N loop
                  } // M loop
                  tempStiffness *= Vol;
                  Kele(a*dim + i, b*dim + j) += tempStiffness;
                } // j loop
              } // b loop
            } // i loop
          } // a loop
        } // Compute stiffness matrix

        if ( R->getRequest() & DMATPROP ) {
          // Assemble dRdalpha
          for (int alpha = 0; alpha < NumPropPerMat; alpha++) {
            for(int a = 0; a < numNodes; a++) {
              for(int i = 0; i < dim; i++) {
                Real tempdRdalpha = 0.0;
                for (int J = 0; J < dim; J++) {
                  tempdRdalpha += FKres.Dmat.get(alpha,i,J) * geomEl->getDN(q, a, J);
                } // J loop
                tempdRdalpha *= Vol;
                (dRdalpha[ (_materials[e*numQP + q]->getMatID())*NumPropPerMat + alpha])( NodesID[a]*dim + i) += tempdRdalpha;
              } // i loop
            } // a loop
          } // alpha loop

        } // Compute DMATPROP

      } // QP loop

      if ( R->getRequest() & STIFFNESS ) {
        // Transform in triplets Kele
        for(int a = 0; a < numNodes; a++) {
          for(int i = 0; i < dim; i++) {
            for(int b = 0; b < numNodes; b++) {
              for(int j = 0; j < dim; j++) {
		R->addStiffness(NodesID[a]*dim + i, NodesID[b]*dim + j, Kele(a*dim + i, b*dim + j));
                // KtripletList.push_back( Triplet<Real >( NodesID[a]*dim + i, NodesID[b]*dim + j, Kele(a*dim + i, b*dim + j) ) );
              }
            }
          }
        }
      }

    } // Element loop



    // Compute Gradg and Hg
    if ( R->getRequest() & DMATPROP ) {

      // First extract residual
      VectorXd Residual = VectorXd::Zero(PbDoF);
      for (int i = 0; i < PbDoF; i++) {
        Residual(i) = R->getResidual(i); // local copy
      }

      for (int alpha = 0; alpha < TotNumMatProp; alpha++) {
	R->addGradg(alpha, 2.0*dRdalpha[alpha].dot(Residual) );
	for (int beta = 0; beta < TotNumMatProp; beta++) {
	  // !!!
	  // WARNING : assume W is linear in alpha!!!!!
	  // !!!
	  R->addHg(alpha, beta,  2.0*dRdalpha[alpha].dot(dRdalpha[beta]));
	  // HgtripletList.push_back( Triplet<Real >( alpha, beta,  2.0*dRdalpha[alpha].dot(dRdalpha[beta]) ) );
	} // alpha loop
      } // beta loop
    } // Compute Gradg and Hg



  }; // Compute Mechanics Body







  // Check consistency on Hg needed for material parameter identification
  void MechanicsBody::checkDmat(Result* R, Real perturbationFactor, Real hM, Real tol)
  {

    // Perturb initial config - gradg and Hg are zero at F = I
    const int nodeNum   = _myMesh->getNumberOfNodes();
    const int nLocalDoF = nodeNum*_nodeDoF;
    const int PbDoF = R->getPbDoF();

    // Perturb field randomly to change from reference configuration
    // Save perturbed field to set the configuration back to reference
    // at the end of the test
    vector<Real > perturb(nLocalDoF, 0.0);
    srand( time(NULL) );
    for(int a = 0; a < nodeNum; a++) {
      for(int i = 0; i < _nodeDoF; i++) {
        Real randomNum =  perturbationFactor*(Real(rand())/RAND_MAX - 0.5);
        perturb[a*_nodeDoF + i] = randomNum;
        R->linearizedUpdate(a*_nodeDoF + i, randomNum);
      }
    }

    Real error = 0.0, norm = 0.0;

    set<MechanicsMaterial *> UNIQUEmaterials;
    for (int i = 0; i < _materials.size(); i++)
    UNIQUEmaterials.insert(_materials[i]);

    cout << "Number of unique materials = " <<  UNIQUEmaterials.size() << endl;
    R->setRequest(8); // First compute gradg and Hg numerically
    this->compute(R);

    // Test gradg //
    R->setRequest(2); // Reset result request so that only forces are computed
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

	Real RTRplus = 0.0;
	for (int i = 0; i < PbDoF; i++) {
	  RTRplus += square(R->getResidual(i));
	}
        // Perturb -2hM the material property alpha
        MatProp[m] -= 2.0*hM;
        // Reset matProp in the materials with MatProp[m]
        (*itMat)->setMaterialParameters(MatProp);
        // Compute R
        this->compute(R);
        Real RTRminus = 0;
	for (int i = 0; i < PbDoF; i++) {
	  RTRminus += square(R->getResidual(i));
	}
        // Bring back to original value of alpha
        MatProp[m] += hM;
        // Reset matProp in all the materials with MatProp[m]
        (*itMat)->setMaterialParameters(MatProp);

        error += pow( (RTRplus-RTRminus)/(2.0*hM) -
		      R->getGradg( MatID*MatProp.size() + m), 2.0 );
        norm  += pow( R->getGradg(  MatID*MatProp.size() + m), 2.0 );

        // cout << (*R->_Gradg)( MatID*MatProp.size() + m) << " " << (RTRplus-RTRminus)/(2.0*hM) << endl;
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
    R->setRequest(8);
    this->compute(R);
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
            Real GradPlus = R->getGradg( MatIDA*MatPropA.size() + mA);

            // Perturb -2hM the material property alpha
            MatPropB[mB] -= 2.0*hM;
            // Reset matProp in the materials with MatProp[m]
            (*itMatB)->setMaterialParameters(MatPropB);
            // Compute R
            this->compute(R);
            Real GradMinus = R->getGradg( MatIDA*MatPropA.size() + mA);

            // Bring back to original value of alpha
            MatPropB[mB] += hM;
            // Reset matProp in all the materials with MatPropB[m]
            (*itMatB)->setMaterialParameters(MatPropB);

            error += pow( (GradPlus - GradMinus)/(2.0*hM) -
            R->getStiffness( MatIDA*MatPropA.size() + mA, MatIDB*MatPropB.size() + mB ), 2.0);
            norm += pow( R->getStiffness( MatIDA*MatPropA.size() + mA, MatIDB*MatPropB.size() + mB ) , 2.0 );

            // cout << (*R->_Gradg)( MatID*MatProp.size() + m) << " " << (RTRplus-RTRminus)/(2.0*hM) << endl;

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
    for(int a = 0; a < nodeNum; a++) {
      for(int i = 0; i < _nodeDoF; i++) {
        R->linearizedUpdate(a*_nodeDoF + i, -perturb[a*_nodeDoF + i]);
      }
    }

  } // Check consistency of gradg and Hg - checkDmat

 





  // Writing output
  void MechanicsBody::writeOutputVTK(const string OutputFile, int step, Result* R)
  {
    /////
    // Todo: NEED TO BE REWRITTEN TAKING INTO ACCOUNT MULTIPLE QUADRATURE POINTS PER ELEMENT !!!
    // Todo: Figure out how to handle mixed meshes
    /////

    // Rewrite it with VTK Libraries
    // Create outputFile name
    string outputFileName = OutputFile + boost::lexical_cast<string>(step) + ".vtu";
    vtkSmartPointer<vtkUnstructuredGrid> newUnstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // Insert Points:
    int NumNodes = _myMesh->getNumberOfNodes();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    // points->SetNumberOfPoints(NumNodes);
    for (int i = 0; i < NumNodes; i++) {
      float x_point = 0.0; float y_point = 0.0; float z_point = 0.0;
      x_point = _myMesh->getX(i)(0);
      if (_myMesh->getDimension() > 1) y_point = _myMesh->getX(i)(1);
      if (_myMesh->getDimension() > 2) z_point = _myMesh->getX(i)(2);
      points->InsertNextPoint(x_point, y_point, z_point);
      // points->InsertPoint(i, x_point, y_point, z_point);
    }
    newUnstructuredGrid->SetPoints(points);
    
    // Element Connectivity:
    vector <GeomElement*> elements = _myMesh->getElements();
    int NumEl = elements.size();
    int NodePerEl = (elements[0])->getNodesPerElement();
    string ElType = _myMesh->getElementType();
    int dim = _myMesh->getDimension();

    VTKCellType cellType;
    // Set Cell Type: http://www.vtk.org/doc/nightly/html/vtkCellType_8h.html
    if (ElType == "C3D4") {
      // Full integration linear tetrahedral element
      cellType = VTK_TETRA; }
    else if (ElType == "C3D10") {
      // Full integration quadratic tetrahedral element
      cellType = VTK_QUADRATIC_TETRA; }
    else {
      cout << "3D Element type not implemented in MechanicsBody writeOutput." << endl;
      exit(EXIT_FAILURE);
    }

    for (int el_iter = 0; el_iter < NumEl; el_iter++) {
      vtkSmartPointer<vtkIdList> elConnectivity = vtkSmartPointer<vtkIdList>::New();

      const vector<int > & NodesID = (elements[el_iter])->getNodesID();
      for (int n = 0; n < NodePerEl; n++) {
        elConnectivity->InsertNextId(NodesID[n]);
      }
      newUnstructuredGrid->InsertNextCell(cellType, elConnectivity);
    }

    // ** BEGIN: POINT DATA ** //
    // ~~ BEGIN: DISPLACEMENTS ~~ //
    vtkSmartPointer<vtkDoubleArray> displacements = vtkSmartPointer<vtkDoubleArray>::New();
    displacements->SetNumberOfComponents(dim);
    displacements->SetName("displacement");
    displacements->SetComponentName(0, "X");
    if (dim > 1) displacements->SetComponentName(1, "Y");
    if (dim > 2) displacements->SetComponentName(2, "Z");

    for (int i = 0; i < NumNodes; i++ ) {
      double x[dim];
      VectorXd X = _myMesh->getX(i);
      for (int j = 0; j < dim; j++) {
        x[j] = R->getField(i*dim + j) - X(j);
      }
      displacements->InsertNextTuple(x);
    }
    newUnstructuredGrid->GetPointData()->AddArray(displacements);
    // ~~ END: DISPLACEMENTS ~~ //
    
    // ~~ BEGIN: RESIDUALS ~~ //
    vtkSmartPointer<vtkDoubleArray> residuals = vtkSmartPointer<vtkDoubleArray>::New();
    residuals->SetNumberOfComponents(dim);
    residuals->SetName("residual");
    residuals->SetComponentName(0, "X"); residuals->SetComponentName(1, "Y");
    if (dim > 2) residuals->SetComponentName(2, "Z");

    // Compute Residual
    int PbDoF = ( _myMesh->getNumberOfNodes())*this->getDoFperNode();
    R->setRequest(FORCE);
    this->compute(R);

    for (int i = 0; i < NumNodes; i++ ) {
      double res[dim];
      for (int j = 0; j < dim; j++) {
        res[j] = R->getResidual(i*dim + j);
      }
      residuals->InsertNextTuple(res);
    }
    newUnstructuredGrid->GetPointData()->AddArray(residuals);
    // ~~ END: RESIDUALS ~~ //
    // ** END: POINT DATA ** //
    
    // ** BEGIN: CELL DATA ** //
    // ~~ BEGIN: \alpha MATERIAL PROPERTY (MAT_PARAM_ID) ~~ //
    vtkSmartPointer<vtkDoubleArray> alpha = vtkSmartPointer<vtkDoubleArray>::New();
    alpha->SetNumberOfComponents(2);
    alpha->SetName("alpha");
    alpha->SetComponentName(0, "Alpha_1"); alpha->SetComponentName(1, "Alpha_2");
    for (int e = 0; e < NumEl; e++) {
      GeomElement* geomEl = elements[e];
      const int numQP = geomEl->getNumberOfQuadPoints();
      double alpha_arr[2];
      Real AvgMatProp_alpha1 = 0.0;
      Real AvgMatProp_alpha2 = 0.0;
      for (int q = 0; q < numQP; q++) {
        vector <Real> MatProp = _materials[e*numQP + q]->getMaterialParameters();
	if (!MatProp.empty()) {
	  if (MatProp.size() > 0) AvgMatProp_alpha1 += MatProp[0];
	  if (MatProp.size() > 1) AvgMatProp_alpha2 += MatProp[1];
        }
      }
      AvgMatProp_alpha1 /= double(numQP); AvgMatProp_alpha2 /= double(numQP);
      alpha_arr[0] = AvgMatProp_alpha1; alpha_arr[1] = AvgMatProp_alpha2;
      alpha->InsertNextTuple(alpha_arr);
    }
    newUnstructuredGrid->GetCellData()->AddArray(alpha);
    // ~~ END: \alpha MATERIAL PROPERTY (MAT_PARAM_ID) ~~ //
    
    // ~~ BEGIN: INTERNAL VARIABLES ~~ //
    // TODO: This method assumes the same material throughout the entire body
    int numInternalVariables = (_materials[0]->getInternalParameters()).size();
    if (numInternalVariables > 0) {
      vtkSmartPointer <vtkDoubleArray> internalVariables = vtkSmartPointer<vtkDoubleArray>::New();
      internalVariables->SetName("Material_Internal_Variables");
      internalVariables->SetNumberOfComponents(numInternalVariables);
      for (int i = 0; i < numInternalVariables; i++) {
        string tempName = "Internal_Variable_" +  boost::lexical_cast<string>(i);
        internalVariables->SetComponentName(i, tempName.c_str());
      }
      for (int e = 0; e < NumEl; e++) {
        GeomElement* geomEl = elements[e];
        const int numQP = geomEl->getNumberOfQuadPoints();
        vector<Real> IntProp = _materials[e*numQP]->getInternalParameters();
        if (!IntProp.size() == numInternalVariables) {
	  cout << "Internal Variables output for multi-materials not supported yet." << endl;
          // double* tempIntProp = new double[numInternalVariables]();
          double tempIntProp[numInternalVariables];
          for (int p = 0; p < numInternalVariables; p++) tempIntProp[p] = -123.4; // Some error value. NaN is better.
	  internalVariables->InsertNextTuple(tempIntProp);
	  // delete tempIntProp;
	  continue;
        }

        fill(IntProp.begin(), IntProp.end(), 0.0);

        // Get Internal Properties from each quad point.
        for (int q = 0; q < numQP; q++) {
          vector <Real> IntPropQuad = _materials[e*numQP + q]->getInternalParameters();
          for (int p = 0; p < numInternalVariables; p++) IntProp[p] = IntProp[p] + IntPropQuad[p];
        }
        // Average over quad points for cell data.
        // double* intPropArr = new double(numInternalVariables);
        double intPropArr[numInternalVariables];
        for (int i = 0; i < numInternalVariables; i++) 
  	  intPropArr[i] = IntProp[i]/numQP;
        internalVariables->InsertNextTuple(intPropArr);
	// delete [] intPropArr;
      }
      newUnstructuredGrid->GetCellData()->AddArray(internalVariables);
    }
    // ~~ END: INTERNAL VARIABLES ~~ //
  
    // ~~ BEGIN: 1ST PIOLA-KIRCHHOFF STRESS ~~ //
    // Loop through elements, also through material points array, which is unrolled
    // int index = 0;
    
    vtkSmartPointer <vtkDoubleArray> FirstPKStress = vtkSmartPointer<vtkDoubleArray>::New();
    FirstPKStress->SetName("First_PK_Stress");
    FirstPKStress->SetNumberOfComponents(dim*dim);
    for (int i = 0; i < dim; i++) { // row
      for (int j = 0; j < dim; j++) {
        string tempCompName = "P" + boost::lexical_cast<string>(i+1) + boost::lexical_cast<string>(j+1);
        FirstPKStress->SetComponentName(i*3 + j, tempCompName.c_str());
      }
    }

    MechanicsMaterial::FKresults FKres;
    FKres.request = 2;
    for(int e = 0; e < NumEl; e++)
    {
      // The dynamic allocation results in a memory leak + Seg fault. Not great.
      // double* eleFirstPKStress = new (nothrow) double(dim*dim);
      double eleFirstPKStress[dim * dim];
      for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
          eleFirstPKStress[i*dim + j] = 0.0;

      // double eleFirstPKStress[9] = {0.0};;
      GeomElement* geomEl = elements[e];
      const int numQP = geomEl->getNumberOfQuadPoints();
      
      // F at each quadrature point are computed at the same time in one element
      vector<Matrix3d > Flist(numQP, Matrix3d::Zero());
      // Compute deformation gradients for current element
      this->computeDeformationGradient(Flist, geomEl, R);

      // Loop over quadrature Points
      // for (int q = 0; q < numQP; q++) {
      //   _materials[e]->compute(FKres, Flist[q], &Fiber);
      /*
      WARNING
      THIS ONLY WORKS WITH 1QP PER ELEMENT - BAD - NEEDS TO BE CHANGED!!
      */
      // }
      // Read through this on how to visualize data at integration points:
      // http://www.vtk.org/Wiki/VTK/VTK_integration_point_support

      _materials[e*numQP + 0]->compute(FKres, Flist[0]);
      
      for (int i = 0; i < dim; i++)
	for (int j = 0; j < dim; j++) 
	  eleFirstPKStress[i*3 + j] = FKres.P(i,j);
      
      FirstPKStress->InsertNextTuple(eleFirstPKStress);
      // delete eleFirstPKStress;
    }
    newUnstructuredGrid->GetCellData()->AddArray(FirstPKStress);
    // ~~ END: 1ST PIOLA-KIRCHHOFF STRESS ~~ //
    // ** END: CELL DATA ** //

    // Write file
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(outputFileName.c_str());
    writer->SetInput(newUnstructuredGrid);
    writer->Write();
    
  } // writeOutput

  

  



} // namespace voom
