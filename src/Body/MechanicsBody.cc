#include "MechanicsBody.h"
#include <string>

namespace voom {

  // Constructor
  MechanicsBody::MechanicsBody(Mesh* myMesh, State* myState,
			       vector<MechanicsMaterial * > Materials):
    Body(myMesh, myState), _materials(Materials)
    {
      // Initialize _field vector for the part corresponding to this Body
      this->initializeField(1);
    } // Constructor





  void MechanicsBody::initializeField(Real fact) {
    vector<GeomElement* > Elements = _myMesh->getElements();
    for (int e = 0; e < _myMesh->getNumberOfElements(); e++) { // Loop Thhrough all elements - some nodes will be initiliazed multiple times - Otherwise we can first compile a set of unique nodal ids for this body
      vector<int > ElNodesID = Elements[e]->getNodesID();
      for (int n = 0; n < ElNodesID.size(); n++) {
	int nodeNum = ElNodesID[n];
	for (int i = 0; i < 3; i++) {
	  _myState->setPhi( nodeNum, i, _myState->getX(nodeNum)(i) );
	}
      }
    }
  } // initializeField





  // Compute deformation gradient
  void MechanicsBody::computeDeformationGradient(vector<Matrix3d > & Flist, GeomElement* geomEl)
  {
    // Compute F at all quadrature points
    const int numQP = geomEl->getNumberOfQuadPoints();
    const vector<int > & NodesID = geomEl->getNodesID();
    const int nodeNum = NodesID.size(), dim = 3; // We assume Mechanics Body is in 3D. Need to generalize this.

    for(int q = 0; q < numQP; q++) {
      // Initialize F to zero
      Flist[q] = Matrix3d::Zero();
      for(int i = 0; i < dim; i++) {
	for(int J = 0; J < dim; J++) {
	  for(int a = 0; a < nodeNum; a++) {
	    Flist[q](i,J) += _myState->getPhi(NodesID[a], i) * geomEl->getDN(q, a, J);
	  }
	}
      }
    } // loop over quadrature points
  };  // computeDeformationGradient





  void MechanicsBody::computeGreenLagrangianStrainTensor(vector<Matrix3d> & Elist, GeomElement* geomEl) 
  {
    const int numQP = geomEl->getNumberOfQuadPoints();
    vector<Matrix3d > Flist(numQP, Matrix3d::Zero());
    this->computeDeformationGradient(Flist, geomEl);

    for (int q = 0; q < numQP; q++) {
      Elist[q] = 0.5 * (Flist[q].transpose() * Flist[q] - Matrix3d::Identity());
    }
  };  // computeGreenLagrangianStrainTensor





  MechanicsBody::InvStruct MechanicsBody::computeInvariants(GeomElement* geomEl, int ElNum) 
  {
    InvStruct ElInvariants;
    // Compute deformation gradients
    const int numQP = geomEl->getNumberOfQuadPoints();
    vector<Matrix3d > Flist(numQP, Matrix3d::Zero());
    this->computeDeformationGradient(Flist, geomEl);
    
    // Setup names and storage for invariants
    ElInvariants.InvName.push_back("J");    vector<Real > J;
    int NumDir = vector<Vector3d> (_materials[0]-> getDirectionVectors()).size();
    for (int i = 0; i < NumDir; i++) {
      std::ostringstream oss;
      oss << "E_" << i << i;
      ElInvariants.InvName.push_back(oss.str()); 
    }
    vector<vector<Real > > TempInv(NumDir, vector<Real > (numQP, 0.0));

    // Loop through quadrature points and per each QP compute an invariant
    for (int q = 0; q < numQP; q++) {
      Matrix3d F = Flist[q], E = Matrix3d::Zero();
      E = 0.5 * (F.transpose() * F - Matrix3d::Identity());
      J.push_back( F.determinant() );

      vector<Vector3d > Dir = _materials[ElNum*numQP + q]->getDirectionVectors();
      for (int i = 0; i < NumDir; i++) {
	TempInv[i][q] = Dir[i].dot(E*Dir[i]);
      }
      
    }
    ElInvariants.InvValue.push_back(J);
    for (int i = 0; i < NumDir; i++) {
      ElInvariants.InvValue.push_back(TempInv[i]);
    }

    return ElInvariants;
     
  };  // computeInvariants





  // Compute Function - Compute Energy, Force, Stiffness
  void MechanicsBody::compute(Result* R)
  {
    const vector<GeomElement* > elements = _myMesh->getElements();
    const int AvgNodePerEl = ((elements[0])->getNodesID()).size(); // It is just to pre-allocate memory space
    const int NumEl = elements.size();
    const int dim = 3;
    
    int TotNumMatProp = R->getNumMatProp();
    int PbDoF = R->getPbDoF();
    vector<VectorXd > dRdalpha;
    VectorXd BodyResidual = VectorXd::Zero(PbDoF);
    if ( R->getRequest() & DMATPROP ) {
      dRdalpha.assign( TotNumMatProp, VectorXd::Zero(PbDoF) ); 
    }
    int NumPropPerMat = (_materials[0]->getMaterialParameters()).size(); // Assumes all materials in Body have the same number of material properties


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
      this->computeDeformationGradient(Flist, geomEl);

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
              R->addResidual(_myState->getGdof(NodesID[a])+i, tempResidual);
	      BodyResidual(_myState->getGdof(NodesID[a])+i) += tempResidual;
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
                      tempStiffness += FKres.K.get(i, M, j, N)*elements[e]->getDN(q, a, M)*elements[e]->getDN(q, b, N);
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
                (dRdalpha[ (_materials[e*numQP + q]->getMatID())*NumPropPerMat + alpha])( _myState->getGdof(NodesID[a]) + i ) += tempdRdalpha;
              } // i loop
            } // a loop
          } // alpha loop

        } // Compute DMATPROP

      } // QP loop

      if ( R->getRequest() & STIFFNESS ) {
        // Add Kele to global stiffness matrix
        for(int a = 0; a < numNodes; a++) {
          for(int i = 0; i < dim; i++) {
            for(int b = 0; b < numNodes; b++) {
              for(int j = 0; j < dim; j++) {
		R->addStiffness(_myState->getGdof(NodesID[a]) + i, _myState->getGdof(NodesID[b]) + j, Kele(a*dim + i, b*dim + j));
              }
            }
          }
        }
      }

    } // Element loop


 
    // Compute Gradg and Hg
    if ( R->getRequest() & DMATPROP ) {
      for (int alpha = 0; alpha < TotNumMatProp; alpha++) {
	R->addGradg(alpha, 2.0*dRdalpha[alpha].dot(BodyResidual) );
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
    vector<GeomElement* > Elements = _myMesh->getElements();
    set<int > UniqueBodyNodes;
    for (int e = 0; e < _myMesh->getNumberOfElements(); e++) {
      vector<int > ElNodesID = Elements[e]->getNodesID();
      for (int n = 0; n < ElNodesID.size(); n++) {
	UniqueBodyNodes.insert(n);
      }
    }

    // Perturb field randomly to change from reference configuration
    // Save perturbed field to set the configuration back to reference 
    // at the end of the test
    srand( time(NULL) );
    map<int, Real> perturb;
    for (set<int >::iterator it = UniqueBodyNodes.begin(); it != UniqueBodyNodes.end(); it++) {
      for (int i = 0; i < _myState->getNodeDof(*it); i++) {
	Real randomNum = perturbationFactor*(Real(rand())/RAND_MAX - 0.5);
	perturb.insert(pair<int, Real>(_myState->getGdof(*it) + i, randomNum));
	_myState->linearizedUpdate(_myState->getGdof(*it) + i, randomNum);
      }
    }

    Real error = 0.0, norm = 0.0;
    int PbDoF = R->getPbDoF();

    set<MechanicsMaterial *> UNIQUEmaterials;
    for (int i = 0; i < _materials.size(); i++) {
      UNIQUEmaterials.insert(_materials[i]); };

    cout << "Number of unique materials = " <<  UNIQUEmaterials.size() << endl;
    R->setRequest(DMATPROP); // First compute gradg and Hg numerically
    R->resetResults(DMATPROP);
    this->compute(R);

    // Test gradg //
    R->setRequest(FORCE); // Reset result request so that only forces are computed
    for ( set<MechanicsMaterial *>::iterator itMat =  UNIQUEmaterials.begin(); itMat != UNIQUEmaterials.end(); itMat++ )
    {
      vector<Real > MatProp = (*itMat)->getMaterialParameters();
      int MatID = (*itMat)->getMatID();
      for(int m = 0; m < MatProp.size(); m++) {
        // Perturb +hM the material property alpha
        MatProp[m] += hM;
        // Reset matProp in the materials with MatProp[m]
        (*itMat)->setMaterialParameters(MatProp);
        // Compute R
	R->resetResults(FORCE);
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
	R->resetResults(FORCE);
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
        norm  += pow( R->getGradg( MatID*MatProp.size() + m), 2.0 );

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
    R->setRequest(DMATPROP);
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
	    R->resetResults(DMATPROP);
            this->compute(R);
	    R->FinalizeHgAssembly();
            Real GradPlus = R->getGradg( MatIDA*MatPropA.size() + mA);

            // Perturb -2hM the material property alpha
            MatPropB[mB] -= 2.0*hM;
            // Reset matProp in the materials with MatProp[m]
            (*itMatB)->setMaterialParameters(MatPropB);
            // Compute R
	    R->resetResults(DMATPROP);
            this->compute(R);
	    R->FinalizeHgAssembly(); // Hg is not used at this point but it is computed since DMATPROP computes both Gradg and Hg
            Real GradMinus = R->getGradg( MatIDA*MatPropA.size() + mA);

            // Bring back to original value of alpha
            MatPropB[mB] += hM;
            // Reset matProp in all the materials with MatPropB[m]
            (*itMatB)->setMaterialParameters(MatPropB);
	    
	    // Need to compute Hg at the reference material properties
	    R->resetResults(DMATPROP);
            this->compute(R);
	    R->FinalizeHgAssembly(); // Hg is not used at this point but it is computed since DMATPROP computes both Gradg and Hg

            error += pow( (GradPlus - GradMinus)/(2.0*hM) -
            R->getHg( MatIDA*MatPropA.size() + mA, MatIDB*MatPropB.size() + mB ), 2.0);
            norm += pow( R->getHg( MatIDA*MatPropA.size() + mA, MatIDB*MatPropB.size() + mB ) , 2.0 );

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
    for (set<int >::iterator it = UniqueBodyNodes.begin(); it != UniqueBodyNodes.end(); it++) {
      for (int i = 0; i < _myState->getNodeDof(*it); i++) {
	int dof = _myState->getGdof(*it) + i;
	_myState->linearizedUpdate(dof, -perturb[dof]);
      }
    }

  } // Check consistency of gradg and Hg - checkDmat





  // Writing output
  void MechanicsBody::writeOutputVTK(const string OutputFile, int step)
  {
    /////
    // Todo: NEED TO BE REWRITTEN TAKING INTO ACCOUNT MULTIPLE QUADRATURE POINTS PER ELEMENT !!!
    /////

    // Rewrite it with VTK Libraries
    // Create outputFile name
    string outputFileName = OutputFile + boost::lexical_cast<string>(step) + ".vtu";
    vtkSmartPointer<vtkUnstructuredGrid> newUnstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // Go to local node numbering for VTK file
    vector<int > LtoG = _myMesh->getLocalToGlobal();
    map<int, int> GtoL;
    for (int i = 0; i < LtoG.size(); i++) {
      GtoL.insert(pair<int, int>(LtoG[i], i) );
    }



    // Insert Points:
    int NumNodes = LtoG.size();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    // points->SetNumberOfPoints(NumNodes);
    for (int n = 0; n < NumNodes; n++) 
      points->InsertNextPoint(_myState->getX(LtoG[n])(0), _myState->getX(LtoG[n])(1), _myState->getX(LtoG[n])(2) );
    newUnstructuredGrid->SetPoints(points);
    


    // Element Connectivity:
    // To-do: Figure out how to handle mixed meshes
    string ElType = _myMesh->getElementType();
    VTKCellType cellType;

    if (ElType == "C3D8") // Full integration hexahedral element
      cellType = VTK_HEXAHEDRON;
    else if (ElType == "C3D8R") // Reduced integration hexahedral element
      cellType = VTK_HEXAHEDRON;
    else if (ElType == "C3D4") // Full integration linear tetrahedral element
      cellType = VTK_TETRA;
    else if (ElType == "C3D10") // Full integration quadratic tetrahedral element
      cellType = VTK_QUADRATIC_TETRA;
    else if (ElType == "TD3")   // Full integration linear triangular element
      cellType = VTK_TRIANGLE;
    else if (ElType == "TD6")   // Full integration quadratic triangular element
      cellType = VTK_QUADRATIC_TRIANGLE;
    else if (ElType == "Q4")    // Full integration linear quadrilateral element
      cellType = VTK_QUAD;
    else {
      cout << "3D Element type not implemented in MechanicsModel writeOutput." << endl;
      exit(EXIT_FAILURE);
    }

    vector <GeomElement*> elements = _myMesh->getElements();
    for (int el_iter = 0; el_iter < elements.size(); el_iter++) {
      vtkSmartPointer<vtkIdList> elConnectivity = vtkSmartPointer<vtkIdList>::New();

      const vector<int > & NodesID = (elements[el_iter])->getNodesID();
      for (int n = 0; n < NodesID.size(); n++) {
        elConnectivity->InsertNextId(GtoL[NodesID[n]]);
      }
      newUnstructuredGrid->InsertNextCell(cellType, elConnectivity);
    }

 
    // ** BEGIN: POINT DATA ** //
    // ~~ BEGIN: DISPLACEMENTS ~~ //
    vtkSmartPointer<vtkDoubleArray> displacements = vtkSmartPointer<vtkDoubleArray>::New();
    displacements->SetNumberOfComponents(3); // Dimension is fixed to 3 for now
    displacements->SetName("displacement");
    displacements->SetComponentName(0, "X");
    displacements->SetComponentName(1, "Y");
    displacements->SetComponentName(2, "Z");

    for (int n = 0; n < NumNodes; n++ ) {
      double u[3];
      int node = LtoG[n];
      for (int i = 0; i < 3; i++) {
	u[i] = _myState->getPhi(node, i) - _myState->getX(node, i);
      }
      displacements->InsertNextTuple(u);
    }
    newUnstructuredGrid->GetPointData()->AddArray(displacements);
    // ~~ END: DISPLACEMENTS ~~ //


    // ** BEGIN CELL DATA ** //
    // ~~ INVARIANTS ~~ //
    vtkSmartPointer <vtkDoubleArray> InvariantsVTK = vtkSmartPointer<vtkDoubleArray>::New();
    InvariantsVTK->SetName("Invariants");
    InvStruct ElInvariants;
    ElInvariants = this->computeInvariants(elements[0], 0);
    int numInv = ElInvariants.InvName.size();

    InvariantsVTK->SetNumberOfComponents(numInv);
    for (int i = 0; i < numInv; i++) {
      InvariantsVTK->SetComponentName(i, (ElInvariants.InvName[i]).c_str());
    }

    for (int e = 0; e < elements.size(); e++) {
      // Setup an empty array to hold all Green strains
      double AvgElInvariants[numInv];

      ElInvariants = computeInvariants(elements[e], e);
      for (int i = 0; i < numInv; i++) {
	AvgElInvariants[i] = ElInvariants.AverageValue(i);
      }

      InvariantsVTK->InsertNextTuple(AvgElInvariants);
    }
    newUnstructuredGrid->GetCellData()->AddArray(InvariantsVTK);
    // ** END CELL DATA ** //
    

    // Write file
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(outputFileName.c_str());
    #if VTK_MAJOR_VERSION <= 5
        writer->SetInput(newUnstructuredGrid);
    #else
        writer->SetInputData(newUnstructuredGrid);
    #endif
    writer->Write();
    
  } // writeOutput
 


  void MechanicsBody::writeQPdataVTK(const string OutputFile, int step) 
  {
    string outputFileName = OutputFile + boost::lexical_cast<string>(step) + ".vtp";
    
    vtkSmartPointer<vtkPolyData> IntegrationPointGrid = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> IntegrationPoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkDoubleArray> IntegrationPointsDisplacements = vtkSmartPointer<vtkDoubleArray>::New();
    IntegrationPointsDisplacements->SetNumberOfComponents(3);
    IntegrationPointsDisplacements->SetName("Displacements");

    vector <GeomElement*> elements = _myMesh->getElements();
    for (int e = 0; e < elements.size(); e++) {
      GeomElement* geomEl = elements[e];
      const int numQP = geomEl->getNumberOfQuadPoints();
      const vector<int>& NodesID = geomEl->getNodesID();
      const uint numNodesOfEl = NodesID.size();

      for (int q = 0; q < numQP; q++) {
	Real Point[3] = {0.0};
	Real Displacement[3] = {0.0};
	for (int i = 0; i < 3; i++) {
	  for (int n = 0; n < numNodesOfEl; n++) {
	    Point[i]        += _myState->getX(NodesID[n], i) * geomEl->getN(q,n);
	    Displacement[i] += ( _myState->getPhi(NodesID[n], i) - _myState->getX(NodesID[n], i) ) * geomEl->getN(q, n);
	  }
	}
	IntegrationPoints->InsertNextPoint(Point);
	IntegrationPointsDisplacements->InsertNextTuple(Displacement);
      }
    }
    IntegrationPointGrid->SetPoints(IntegrationPoints);
    IntegrationPointGrid->GetPointData()->AddArray(IntegrationPointsDisplacements);

    
    // Plot material directions
    int NumDir = vector<Vector3d> (_materials[0]->getDirectionVectors()).size();
    vector<vtkSmartPointer<vtkDoubleArray> > directions; 
    vector<vtkSmartPointer<vtkDoubleArray> > alphas;

    for (int i = 0; i < NumDir; i++) {
      vtkSmartPointer<vtkDoubleArray> vtkP_A = vtkSmartPointer<vtkDoubleArray>::New();
      directions.push_back(vtkP_A);
      vtkSmartPointer<vtkDoubleArray> vtkP_B = vtkSmartPointer<vtkDoubleArray>::New();
      alphas.push_back(vtkP_B);
    
    }

    for (int i = 0; i < NumDir; i++) {
      std::ostringstream ossA, ossB;
      ossA << "v_" << i;
      directions[i]->SetNumberOfComponents(3);
      directions[i]->SetName(ossA.str().c_str());

      ossB << "alpha_" << i;
      alphas[i]->SetNumberOfComponents(1);
      alphas[i]->SetName(ossB.str().c_str());
    }
    
    // Jacobian
    vtkSmartPointer<vtkDoubleArray> Jacobian = vtkSmartPointer<vtkDoubleArray>::New();
    Jacobian->SetNumberOfComponents(1);
    Jacobian->SetName("J");
 
    // Other invariants
    vtkSmartPointer <vtkDoubleArray> InvariantsVTK = vtkSmartPointer<vtkDoubleArray>::New();
    InvariantsVTK->SetName("Invariants");
    InvStruct ElInvariants;
    ElInvariants = this->computeInvariants(elements[0], 0);
    int numInv = ElInvariants.InvName.size();

    InvariantsVTK->SetNumberOfComponents(numInv);
    for (int i = 0; i < numInv; i++) {
      InvariantsVTK->SetComponentName(i, (ElInvariants.InvName[i]).c_str());
    }
    


    for (int e = 0; e < elements.size(); e++) {
      GeomElement* geomEl = elements[e];
      const int numQP = geomEl->getNumberOfQuadPoints();
      vector<Matrix3d > Flist(numQP, Matrix3d::Zero());
      this->computeDeformationGradient(Flist, geomEl);

      InvStruct ElInvariants;
      ElInvariants = computeInvariants(elements[e], e);
      
      for (int q = 0; q < numQP; q++) {
	vector<Vector3d > Dir = _materials[e*numQP + q]->getDirectionVectors();
	double J[1] = {Flist[q].determinant()};
	Jacobian->InsertNextTuple(J);

	InvariantsVTK->InsertNextTuple((ElInvariants.InvQP(q)).data());

	for (int i = 0; i < NumDir; i++) {
	  Vector3d DirRotated;
	  DirRotated = Flist[q]*Dir[i];
	  DirRotated /= DirRotated.norm();  // Compute rotation of directions due to deformation
	  double v[3] = {DirRotated[0], DirRotated[1], DirRotated[2]};
	  directions[i]->InsertNextTuple(v);
	  Real angle = acos( Dir[i].dot(DirRotated) );
	  if (angle > 0.5*M_PI) { angle -= M_PI; };
	  alphas[i]->InsertNextTuple(&angle);
	}
      }
    }

    IntegrationPointGrid->GetPointData()->AddArray(Jacobian);
    IntegrationPointGrid->GetPointData()->AddArray(InvariantsVTK); 
    
    for (int i = 0; i < NumDir; i++) {
      IntegrationPointGrid->GetPointData()->AddArray(directions[i]);
      IntegrationPointGrid->GetPointData()->AddArray(alphas[i]);     // Plot direction angle change
    }
    
    

    // Write File
    vtkSmartPointer<vtkXMLPolyDataWriter> IntegrationPointWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    IntegrationPointWriter->SetFileName(outputFileName.c_str());
    #if VTK_MAJOR_VERSION <= 5
    	IntegrationPointWriter->SetInput(IntegrationPointGrid);
    #else	
        IntegrationPointWriter->SetInputData(IntegrationPointGrid);
    #endif    
    IntegrationPointWriter->Write();
  }




} // namespace voom
