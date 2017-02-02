#include "MechanicsModel.h"

namespace voom {

  // Constructor
  MechanicsModel::MechanicsModel(Mesh* aMesh, vector<MechanicsMaterial * > Materials,
    const uint NodeDoF,
    int PressureFlag, Mesh* SurfaceMesh,
    int NodalForcesFlag,
    int ResetFlag,
    int SpringBCflag):
    Model(aMesh, NodeDoF), _materials(Materials),
    _pressureFlag(PressureFlag), _pressure(0.0), _surfaceMesh(SurfaceMesh),
    _nodalForcesFlag(NodalForcesFlag), _forcesID(NULL), _forces(NULL), _resetFlag(ResetFlag), _springBCflag(SpringBCflag)
    {
      // THERE IS ONE MATERIAL PER ELEMENT - CAN BE CHANGED - DIFFERENT THAN BEFORE
      // Resize and initialize (default function) _field vector
      _field.resize(  (_myMesh->getNumberOfNodes() )*_nodeDoF );
      this->initializeField();

      // if (_pressureFlag == 1) {
      _prevField.resize( _field.size() );
      this->setPrevField();
      // }
    }



  // Compute deformation gradient
  void MechanicsModel::computeDeformationGradient(vector<Matrix3d > & Flist, GeomElement* geomEl)
  {
    // Compute F at all quadrature points
    const uint numQP = geomEl->getNumberOfQuadPoints();
    const uint dim   = _myMesh->getDimension();

    const vector<int > & NodesID = geomEl->getNodesID();
    const uint nodeNum = NodesID.size();

    // cout << "Nodes: " << _field[NodesID[0]*3 + 0] << "\t" << _field[NodesID[1]*3 + 0]  << "\t" << _field[NodesID[2]*3 + 0] << "\t" << _field[NodesID[3]*3 + 0]  << endl;
    // cout << "Nodes: " << _field[NodesID[0]*3 + 1] << "\t" << _field[NodesID[1]*3 + 1]  << "\t" << _field[NodesID[2]*3 + 1] << "\t" << _field[NodesID[3]*3 + 1]  << endl;
    // cout << "Nodes: " << _field[NodesID[0]*3 + 2] << "\t" << _field[NodesID[1]*3 + 2]  << "\t" << _field[NodesID[2]*3 + 2] << "\t" << _field[NodesID[3]*3 + 2]  << endl;
    // cout << "Nodes: " << NodesID[0] << "\t" << NodesID[1] << "\t" << NodesID[2] << "\t" << NodesID[3] << endl;

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


  void MechanicsModel::computeGreenLagrangianStrainTensor(vector<Matrix3d> & Elist, GeomElement* geomEl) {
    const int sizeofE = Elist.size();
    vector<Matrix3d> Flist(sizeofE, Matrix3d::Zero());
    this->computeDeformationGradient(Flist, geomEl);

    for (int i = 0; i < sizeofE; i++)
      Elist[i] = 0.5 * (Flist[i].transpose() * Flist[i] - Matrix3d::Identity());
  }





  // Compute Function - Compute Energy, Force, Stiffness
  void MechanicsModel::compute(Result * R)
  {
    const vector<GeomElement* > elements = _myMesh->getElements();
    const int AvgNodePerEl = ((elements[0])->getNodesID()).size();
    const int NumEl = elements.size();
    const int dim = _myMesh->getDimension();
    vector<Triplet<Real > > KtripletList, HgtripletList;

    int PbDoF = R->getPbDoF();
    int TotNumMatProp = R->getNumMatProp();
    vector<VectorXd > dRdalpha;
    int NumPropPerMat = (_materials[0]->getMaterialParameters()).size(); // Assume all materials have the same number of material properties


    // Reset values in result struct
    if ( R->getRequest() & ENERGY ) {
      R->setEnergy(0.0);
    }
    if ( R->getRequest() & FORCE || R->getRequest() & DMATPROP )  {
      R->resetResidualToZero();
    }
    if ( R->getRequest() & STIFFNESS ) {
      R->resetStiffnessToZero();
      KtripletList.reserve(dim*dim*NumEl*AvgNodePerEl*AvgNodePerEl);
    }

    if ( R->getRequest() & DMATPROP ) {
      dRdalpha.assign( TotNumMatProp, VectorXd::Zero(PbDoF) );
      if ( _resetFlag == 1 ) {
        R->resetGradgToZero();
        R->resetHgToZero();
        HgtripletList.reserve(TotNumMatProp*TotNumMatProp);
      }
    }



    // Loop through elements, also through material points array, which is unrolled
    // uint index = 0;
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
          for(uint a = 0; a < numNodes; a++) {
            for(uint i = 0; i < dim; i++) {
              Real tempResidual = 0.0;
              for (uint J = 0; J < dim; J++) {
                tempResidual += FKres.P(i,J) * geomEl->getDN(q, a, J);
              } // J loop
              tempResidual *= Vol;
              R->addResidual(NodesID[a]*dim+i, tempResidual);
            } // i loop
          } // a loop
        } // Internal force loop

        // Compute Stiffness
        if (R->getRequest() & STIFFNESS) {
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
                  // R->addStiffness(NodesID[a]*dim + i, NodesID[b]*dim + j, tempStiffness);
                  Kele(a*dim + i, b*dim + j) += tempStiffness;
                } // j loop
              } // b loop
            } // i loop
          } // a loop
        } // Compute stiffness matrix

        if ( R->getRequest() & DMATPROP ) {
          // Assemble dRdalpha
          for (uint alpha = 0; alpha < NumPropPerMat; alpha++) {
            for(uint a = 0; a < numNodes; a++) {
              for(uint i = 0; i < dim; i++) {
                Real tempdRdalpha = 0.0;
                for (uint J = 0; J < dim; J++) {
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

    // Insert BC terms from spring
    if (_springBCflag == 1) {
      vector<Triplet<Real > >  KtripletList_FromSpring = this->applySpringBC(*R);
      KtripletList.insert( KtripletList.end(), KtripletList_FromSpring.begin(), KtripletList_FromSpring.end() );
    }

    if (_torsionalSpringBCflag == 1) {
      vector<Triplet<Real> > KtripletList_FromTorsionalSpring = this->applyTorsionalSpringBC(*R);
      KtripletList.insert(KtripletList.end(), KtripletList_FromTorsionalSpring.begin(), KtripletList_FromTorsionalSpring.end());
    }

    // Sum up all stiffness entries with the same indices
    if ( R->getRequest() & STIFFNESS ) {
      R->setStiffnessFromTriplets(KtripletList);
      R->FinalizeGlobalStiffnessAssembly();
      // cout << "Stiffness assembled" << endl;
    }

    // Add pressure/external load if any
    if (_pressureFlag == 1) {
      this->applyPressure(R);
    }

    // Add external load if any
    if (_nodalForcesFlag == 1 && (R->getRequest() & FORCE || R->getRequest() & DMATPROP) ) {
      for (int i = 0; i < _forcesID->size(); i++) {
        R->addResidual( (*_forcesID)[i], (*_forces)[i] );
      }
    }

    // Compute Gradg and Hg
    if ( R->getRequest() & DMATPROP ) {

      // First extract residual
      VectorXd Residual = VectorXd::Zero(PbDoF);
      for (int i = 0; i < PbDoF; i++) {
        Residual(i) = R->getResidual(i); // local copy
      }

      if (_resetFlag == 1) {
        for (uint alpha = 0; alpha < TotNumMatProp; alpha++) {
          R->addGradg(alpha, 2.0*dRdalpha[alpha].dot(Residual) );
          for (uint beta = 0; beta < TotNumMatProp; beta++) {
            // !!!
            // WARNING : assume W is linear in alpha!!!!!
            // !!!
            HgtripletList.push_back( Triplet<Real >( alpha, beta,  2.0*dRdalpha[alpha].dot(dRdalpha[beta]) ) );
          } // alpha loop
        } // beta loop
        R->setHgFromTriplets(HgtripletList);
      }
      else {
        for (uint alpha = 0; alpha < TotNumMatProp; alpha++) {
          R->addGradg(alpha, 2.0*dRdalpha[alpha].dot(Residual) );
          for (uint beta = 0; beta < TotNumMatProp; beta++) {
            // !!!
            // WARNING : assume W is linear in alpha!!!!!
            // !!!
            R->addHg(alpha, beta,  2.0*dRdalpha[alpha].dot(dRdalpha[beta]) );
          } // alpha loop
        } // beta loop
      }




    } // Compute Gradg and Hg

  } // Compute Mechanics Model







  void MechanicsModel::applyPressure(Result * R) {

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
        if (R->getRequest() & ENERGY) {
          R->addEnergy( _pressure*Area*a3.dot(u) );
        }

        // Compute Residual
        if ( (R->getRequest() & FORCE) || (R->getRequest() & DMATPROP) ) {
          for(uint a = 0; a < NodesID.size(); a++) {
            for(uint i = 0; i < 3; i++) {
              R->addResidual(NodesID[a]*3+i, _pressure * Area * a3(i) * geomEl->getN(q, a));
              // cout << NodesID[a] << " " << geomEl->getN(q, a) << endl;
            } // i loop
          } // a loop
        } // Internal force loop

      } // loop over QP
    } // loop over elements

  } // apply pressure




  vector<Triplet<Real > > MechanicsModel::applySpringBC(Result & R) {

    vector<Triplet<Real > > KtripletList_FromSpring;

    // Set previous field for every solution step
    this->setPrevField();

    // Recompute normals - no change if _prevField has not changed.
    // this->computeNormals();

    // Loop through _spNodes
    for(int n = 0; n < _spNodes.size(); n++)
    {
      int NodeID = _spNodes[n];
      Vector3d xa_prev, xa_curr;
      xa_prev << _prevField[NodeID*3], _prevField[NodeID*3+1], _prevField[NodeID*3+2];
      xa_curr << _field[NodeID*3], _field[NodeID*3+1], _field[NodeID*3+2];

      // Compute energy
      if (R.getRequest() & ENERGY) {
        R.addEnergy( 0.5* _springK*pow( (xa_curr - xa_prev).dot(_spNormals[n]), 2.0) );
      }

      // Compute Residual
      if ( (R.getRequest() & FORCE) || (R.getRequest() & DMATPROP) ) {
        for(uint i = 0; i < 3; i++) {
          R.addResidual(NodeID*3+i,  _springK*_spNormals[n](i)*(xa_curr - xa_prev).dot(_spNormals[n]) );
        } // i loop
      } // Internal force loop

      // Compute stiffness matrix
      if ( R.getRequest() & STIFFNESS ) {
        for(uint i = 0; i < 3; i++) {
          for(uint j = 0; j < 3; j++) {
            KtripletList_FromSpring.push_back(Triplet<Real >( NodeID*3+i, NodeID*3+j, _springK*_spNormals[n](i)*_spNormals[n](j) ));
          } // j loop
        } // i loop
      } // Stiffness loop

    } // Spring nodes loop

    return  KtripletList_FromSpring;
  } // apply SpringBC



  void MechanicsModel::computeNormals() {

    // First compute normal of any element in _spMesh
    const vector<GeomElement* > elements = _spMesh->getElements();

    vector<Vector3d > ElNormals(elements.size(), Vector3d::Zero());

    // Loop over elements
    for(int e = 0; e < elements.size(); e++)
    {
      GeomElement* geomEl = elements[e];
      const vector<int  >& NodesID = geomEl->getNodesID();
      const int numQP    = geomEl->getNumberOfQuadPoints();
      const int numNodes = NodesID.size();

      // Loop over quadrature points
      for(int q = 0; q < numQP; q++) {

        // Compute normal based on _prevField
        Vector3d a1 = Vector3d::Zero(), a2 = Vector3d::Zero(), a3 = Vector3d::Zero();
        for (int a = 0; a < NodesID.size(); a++) {
          int nodeID = NodesID[a];
          Vector3d xa_prev;
          xa_prev << _prevField[nodeID*3], _prevField[nodeID*3+1], _prevField[nodeID*3+2];
          a1 += xa_prev*geomEl->getDN(q, a, 0);
          a2 += xa_prev*geomEl->getDN(q, a, 1);
        }
        ElNormals[e] += a1.cross(a2); // Not normalized with respect to area (elements with larger area count more)
      } // loop over QP

    } // loop over elements

    // loop over _spNodes
    for (int n = 0; n < _spNodes.size(); n++) {
      // Reset normal to zero
      _spNormals[n] = Vector3d::Zero();
      // Loop over all elements sharing that node
      for (int m = 0; m < _spNodesToEle[n].size(); m++) {
        _spNormals[n] += ElNormals[_spNodesToEle[n][m]];
      }
      Real normFactor = 1.0/_spNormals[n].norm();
      _spNormals[n] *= normFactor;

      // For testing only - to be commented out
      // cout << _spNormals[n](0) << " " << _spNormals[n](1) << " " << _spNormals[n](2) << endl;
    }

  } // compute Normals



  void MechanicsModel::initSpringBC(const string SpNodes, Mesh* SpMesh, Real SpringK) {

    // Activate flag
    _springBCflag = 1;

    // Store mesh information
    _spMesh = SpMesh;

    // Store spring constant
    _springK = SpringK;

    // Store node number on the outer surface
    ifstream inp(SpNodes.c_str());
    if (!inp.is_open()) {
        cout << "** Unable to open spring file " << SpNodes << ".\n ** Exiting..." << endl;
        exit(EXIT_FAILURE);
    }

    int numSpringNodes = 0;
    inp >> numSpringNodes;

    int nodeNum = 0;
    while (inp >> nodeNum)
      _spNodes.push_back(nodeNum);

    // Checking that the number at the top of the file corresponds to the number of nodes
    assert(numSpringNodes == _spNodes.size());
    cout << "** Applying the Spring BC to " << numSpringNodes << " nodes." << endl;

    // Collect elements that share a node in _spNodes
    const vector<GeomElement* > elements = _spMesh->getElements();

    for (int n = 0; n < _spNodes.size(); n++) {
      vector<int > connected;
      for(int e = 0; e < elements.size(); e++) {
        const vector<int >& NodesID = elements[e]->getNodesID();
        for (int m = 0; m < NodesID.size(); m++) {
          if (NodesID[m] == _spNodes[n]) {
            connected.push_back(e);
            break;
          } //  check if node belong to element
        } // loop over nodes of the element
      } // loop over elements
      _spNodesToEle.push_back(connected);

      // For testing only - to be commented out
      // cout << _spNodes[n] << " ";
      // for (int i = 0; i < connected.size(); i++) {
      // 	cout << connected[i] << " ";
      // }
      // cout << endl;

    } // loop over _spNodes

    _spNormals.resize(_spNodes.size(), Vector3d::Zero());
    // Compute initial node normals
    this->computeNormals();

  } // InitSpringBC



  // *** BEGIN: MATERIAL PARAMETER ID FUNCTIONS *** //
  void MechanicsModel::checkDmat(EigenResult * R, Real perturbationFactor, Real hM, Real tol)
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
        Real RTRplus = (R->_residual)->dot(*R->_residual);
        // Perturb -2hM the material property alpha
        MatProp[m] -= 2.0*hM;
        // Reset matProp in the materials with MatProp[m]
        (*itMat)->setMaterialParameters(MatProp);
        // Compute R
        this->compute(R);
        Real RTRminus = (R->_residual)->dot(*R->_residual);

        // Bring back to original value of alpha
        MatProp[m] += hM;
        // Reset matProp in all the materials with MatProp[m]
        (*itMat)->setMaterialParameters(MatProp);

        error += pow( (RTRplus-RTRminus)/(2.0*hM) -
        (*R->_Gradg)( MatID*MatProp.size() + m), 2.0 );
        norm += pow( (*R->_Gradg)(  MatID*MatProp.size() + m), 2.0 );

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
    SparseMatrix<Real > HgAn = *R->_Hg;
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
            HgAn.coeff( MatIDA*MatPropA.size() + mA, MatIDB*MatPropB.size() + mB ), 2.0);
            norm += pow( HgAn.coeff( MatIDA*MatPropA.size() + mA, MatIDB*MatPropB.size() + mB ) , 2.0 );

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
    for(uint a = 0; a < nodeNum; a++) {
      for(uint i = 0; i < _nodeDoF; i++) {
        this->linearizedUpdate(a, i, -perturb[a*_nodeDoF + i]);
      }
    }

  } // Check consistency of gradg and Hg - checkDmat

  // *** END: MATERIAL PARAMETER ID FUNCTIONS *** //



  // Compute volume functions - current and reference volumes
  Real MechanicsModel::computeRefVolume()
  {
    Real RefVol = 0.0;
    const vector<GeomElement* > elements = _myMesh->getElements();
    const int NumEl = elements.size();

    // Loop through elements, also through material points array, which is unrolled
    for(int e = 0; e < NumEl; e++)
    {
      GeomElement* geomEl = elements[e];
      const int numQP    = geomEl->getNumberOfQuadPoints();
      // Loop over quadrature points
      for(int q = 0; q < numQP; q++) {
        // Volume associated with QP q
        RefVol += geomEl->getQPweights(q);
      } // Loop over QP
    } // Loop over elements

    return RefVol;
  }

  Real MechanicsModel::computeCurrentVolume()
  {
    Real CurrVol = 0.0;
    const vector<GeomElement* > elements = _myMesh->getElements();
    const int NumEl = elements.size();

    // Loop through elements, also through material points array, which is unrolled
    for(int e = 0; e < NumEl; e++)
    {
      GeomElement* geomEl = elements[e];
      const int numQP    = geomEl->getNumberOfQuadPoints();

      // F at each quadrature point are computed at the same time in one element
      vector<Matrix3d > Flist(numQP, Matrix3d::Zero());
      // Compute deformation gradients for current element
      this->computeDeformationGradient(Flist, geomEl);

      // Loop over quadrature points
      for(int q = 0; q < numQP; q++) {
        // Volume associated with QP q
        if (isnan(Flist[q].determinant()))
        {
          // this->computeDeformationGradient(Flist, geomEl);
          // cout << Flist[q] << endl;
          cout << "Element with nan: " << e << endl;
        }
        CurrVol += ( geomEl->getQPweights(q) * Flist[q].determinant() );
      } // Loop over QP
    } // Loop over elements

    return CurrVol;
  }

  // Functions for applying Torsional Spring BC
  void MechanicsModel::initTorsionalSpringBC(const string torsionalSpringNodes, Real torsionalSpringK) {
    // Activate flag
     _torsionalSpringBCflag = 1;

    // Store node number on the outer surface
    ifstream inp(torsionalSpringNodes.c_str());
    if (!inp.is_open()) {
	cout << "** Unable to open Torsional spring file " << torsionalSpringNodes << ".\n ** Exiting..." << endl;
	exit(EXIT_FAILURE);
    }
    int numTorsionalSpringNodes = 0;
    inp >> numTorsionalSpringNodes;

    int nodeNum = 0;
    while (inp >> nodeNum)
      _torsionalSpringNodes.push_back(nodeNum);

    // Checking that the number at the top of the file corresponds to the number of nodes
    assert(numTorsionalSpringNodes == _torsionalSpringNodes.size());
    cout << "** Applying the Torsional Spring BC to " << numTorsionalSpringNodes << " nodes." << endl;

    _torsionalSpringK = torsionalSpringK;
    
    // Compute Centroid to compute the tangential vector
    computeCentroid();

    // Compute the tangential vectors
    computeTangents();
  }

  void MechanicsModel::computeCentroid() {
    double xavg = 0.0;
    double yavg = 0.0;
    double zavg = 0.0;
    for (int n = 0; n < _myMesh->getNumberOfNodes(); n++) {
      xavg += _myMesh->getX(n)(0);
      yavg += _myMesh->getX(n)(1);
      zavg += _myMesh->getX(n)(2);
    }
    xavg = xavg/_myMesh->getNumberOfNodes();
    yavg = yavg/_myMesh->getNumberOfNodes();
    zavg = zavg/_myMesh->getNumberOfNodes();
    _centroidLocation(0) = xavg; _centroidLocation(1) = yavg; _centroidLocation(2) = zavg;
  }

  void MechanicsModel::computeTangents() {
    for(int n = 0; n < _torsionalSpringNodes.size(); n++)
    {
      int NodeID = _torsionalSpringNodes[n];
      Vector3d xa_reference(_myMesh->getX(NodeID)(0), _myMesh->getX(NodeID)(1), _myMesh->getX(NodeID)(2));

      Vector2d normal = Vector2d::Zero();
      Vector3d tempNormal = xa_reference -_centroidLocation;
      normal(0) = tempNormal(0); normal(1) = tempNormal(1);
      normal = normal/normal.norm();

      Vector3d tangent(-1.0 * normal(1), normal(0), 0.0);
      _spTangents.push_back(tangent);
    }
  }

  vector<Triplet<Real> > MechanicsModel::applyTorsionalSpringBC(Result & R) {
    vector<Triplet<Real > > KtripletList_FromSpring;

    // Loop through _spNodes
    for(int n = 0; n < _torsionalSpringNodes.size(); n++)
    {
      int NodeID = _torsionalSpringNodes[n];
      Vector3d xa_reference, xa_curr;
      xa_reference << _myMesh->getX(NodeID)(0), _myMesh->getX(NodeID)(1), _myMesh->getX(NodeID)(2);
      xa_curr << _field[NodeID*3], _field[NodeID*3+1], _field[NodeID*3+2];

      // Compute energy
      if (R.getRequest() & ENERGY) {
        R.addEnergy( 0.5* _torsionalSpringK*pow( (xa_curr - xa_reference).dot(_spTangents[n]), 2.0) );
      }

      // Compute Residual
      if ( (R.getRequest() & FORCE) || (R.getRequest() & DMATPROP) ) {
        for(uint i = 0; i < 3; i++) {
          R.addResidual(NodeID*3+i,  _torsionalSpringK*_spTangents[n](i)*(xa_curr - xa_reference).dot(_spTangents[n]) );
        } // i loop
      } // Internal force loop

      // Compute stiffness matrix
      if ( R.getRequest() & STIFFNESS ) {
        for(uint i = 0; i < 3; i++) {
          for(uint j = 0; j < 3; j++) {
            KtripletList_FromSpring.push_back(Triplet<Real >( NodeID*3+i, NodeID*3+j, _torsionalSpringK*_spTangents[n](i)*_spTangents[n](j) ));
          } // j loop
        } // i loop
      } // Stiffness loop
    } // Spring nodes loop

    return  KtripletList_FromSpring;
  }


  // Writing output
  void MechanicsModel::writeOutputVTK(const string OutputFile, int step)
  {
    /////
    // Todo: NEED TO BE REWRITTEN TAKING INTO ACCOUNT MULTIPLE QUADRATURE POINTS PER ELEMENT !!!
    /////

    {
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
    // To-do: Figure out how to handle mixed meshes
    // To-do: It would be better to select based on Abaqus element names
    vector <GeomElement*> elements = _myMesh->getElements();
    int NumEl = elements.size();
    int NodePerEl = (elements[0])->getNodesPerElement();
    int dim = _myMesh->getDimension();

    VTKCellType cellType;
    // Set Cell Type: http://www.vtk.org/doc/nightly/html/vtkCellType_8h.html
    switch (dim) {
      case 3: // 3D
	switch (NodePerEl) {
	  case 4: // Linear Tetrahedron
	    cellType = VTK_TETRA;
	    break;
	  case 10: // Quadratic Tetrahedron
	    cellType = VTK_QUADRATIC_TETRA;
	  default:
	    cout << "3D Element type not implemented in MechanicsModel writeOutput." << endl;
	    exit(EXIT_FAILURE);
	}
	break;
      default:
        cout << "This element has not been implemented in MechanicsModel writeOutput." << endl;
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
    displacements->SetComponentName(1, "Y");
    if (dim > 2) displacements->SetComponentName(2, "Z");

    for (int i = 0; i < NumNodes; i++ ) {
      double x[dim];
      VectorXd X = _myMesh->getX(i);
      for (int j = 0; j < dim; j++) {
        x[j] = _field[i*dim + j] - X(j);
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
    uint PbDoF = ( _myMesh->getNumberOfNodes())*this->getDoFperNode();
    EigenResult myResults(PbDoF, 2);
    int myRequest = 2;
    myResults.setRequest(myRequest);
    this->compute(&myResults);
    VectorXd R = *(myResults._residual);

    for (int i = 0; i < NumNodes; i++ ) {
      double res[dim];
      for (int j = 0; j < dim; j++) {
        res[j] = R(i*dim + j);
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
    // uint index = 0;
    
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
      this->computeDeformationGradient(Flist, geomEl);

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
    
    // ** BEGIN: SUPPORT FOR INTEGRATION POINTS ** //
    
    // ~~ BEGIN: PLOT PRESSURE NORMALS ~~ //
    if (_pressureFlag) {
      string outputIntegrationPointDataFileName = OutputFile + "_IntPointData" + boost::lexical_cast<string>(step) + ".vtp";
      vtkSmartPointer<vtkPolyData> IntegrationPointGrid = vtkSmartPointer<vtkPolyData>::New();
      vtkSmartPointer<vtkPoints> IntegrationPoints = vtkSmartPointer<vtkPoints>::New();
      vtkSmartPointer<vtkDoubleArray> IntegrationPointsDisplacements = vtkSmartPointer<vtkDoubleArray>::New();
      IntegrationPointsDisplacements->SetNumberOfComponents(3);
      IntegrationPointsDisplacements->SetName("Displacements");

      vector <GeomElement*> elements2D = _surfaceMesh->getElements();
      for (int e = 0; e < elements2D.size(); e++) {
	GeomElement* geomEl = elements2D[e];
	const int numQP = geomEl->getNumberOfQuadPoints();
	const vector<int>& NodesID = geomEl->getNodesID();
	const uint numNodesOfEl = NodesID.size();

	for (int q = 0; q < numQP; q++) {
	  float tempDisplacement[3] = {0.0}; float tempPoint[3] = {0.0};
	  for (int d = 0; d < dim; d++) {
	    for (int n = 0; n < numNodesOfEl; n++) {
	      tempPoint[d] += _surfaceMesh->getX(n)(d) * geomEl->getN(q,n);
	      tempDisplacement[d] += (_field[NodesID[n]*dim + d] - _surfaceMesh->getX(n)(d)) * geomEl->getN(q, n);
	    }
	  }
	  IntegrationPoints->InsertNextPoint(tempPoint);
	  IntegrationPointsDisplacements->InsertNextTuple(tempDisplacement);
	}
      }
      IntegrationPointGrid->SetPoints(IntegrationPoints);
      IntegrationPointGrid->GetPointData()->AddArray(IntegrationPointsDisplacements);


      // Compute Normals:
      vtkSmartPointer<vtkDoubleArray> pressureNormal = vtkSmartPointer<vtkDoubleArray>::New();
      pressureNormal->SetNumberOfComponents(3);
      pressureNormal->SetName("Pressure_Normals");
      for(int e = 0; e < elements2D.size(); e++) {
	GeomElement* geomEl = elements2D[e];
	const vector<int  >& NodesID = geomEl->getNodesID();
	const int numQP    = geomEl->getNumberOfQuadPoints();
	const int numNodes = NodesID.size();
	
	// Loop over quadrature points
	for(int q = 0; q < numQP; q++) {
	  // Compute normal based on _prevField and displacement
	  Vector3d a1 = Vector3d::Zero(), a2 = Vector3d::Zero(), a3 = Vector3d::Zero();
	  for (int a = 0; a < NodesID.size(); a++) {
	    int nodeID = NodesID[a];
	    Vector3d xa_prev, xa_curr;
	    xa_prev << _prevField[nodeID*3], _prevField[nodeID*3+1], _prevField[nodeID*3+2];
	    xa_curr << _field[nodeID*3], _field[nodeID*3+1], _field[nodeID*3+2];
	    a1 += xa_prev*geomEl->getDN(q, a, 0);
	    a2 += xa_prev*geomEl->getDN(q, a, 1);
 	  }
	  a3 = a1.cross(a2);
	  
	  // Surface associated with QP q
	  Real Area = a3.norm();
	  a3 /= Area;
	  double tempNormal[3] = {a3[0], a3[1], a3[2]};
	  pressureNormal->InsertNextTuple(tempNormal);
	}
      }
      IntegrationPointGrid->GetPointData()->AddArray(pressureNormal);

      // Write File
      vtkSmartPointer<vtkXMLPolyDataWriter> IntegrationPointWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
      IntegrationPointWriter->SetFileName(outputIntegrationPointDataFileName.c_str());
      IntegrationPointWriter->SetInput(IntegrationPointGrid);
      IntegrationPointWriter->Write();
      // ** END: SUPPORT FOR INTEGRATION POINTS ** //
    }
    // ~~ END: PLOT PRESSURE NORMALS ~~ //

    
 
  }
    return;   

    

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
      out << _myMesh->getX(i)(0) << " " << _myMesh->getX(i)(1) << " " << _myMesh->getX(i)(2) << endl;
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
      x -= _myMesh->getX(i);
      out << x(0) << " " << x(1) << " " << x(2) << endl;
    }

    // Residuals
    out << "VECTORS residual double" << endl;

    // Compute Residual
    uint PbDoF = ( _myMesh->getNumberOfNodes())*this->getDoFperNode();
    EigenResult myResults(PbDoF, 2);
    int myRequest = 2;
    myResults.setRequest(myRequest);
    this->compute(&myResults);
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
      GeomElement* geomEl = elements[e];
      const int numQP = geomEl->getNumberOfQuadPoints();
      Real AvgMatProp = 0.0;
      for ( int q = 0; q < numQP; q++ ) {
        vector<Real > MatProp = _materials[e*numQP + q]->getMaterialParameters();
        if (!MatProp.empty())
        AvgMatProp += MatProp[0];
      }
      out << AvgMatProp/double(numQP) << endl;
    }

    // Alpha_2 material property
    out << "SCALARS Alpha2 double" << endl
    << "LOOKUP_TABLE default" << endl;
    for ( int e = 0; e < NumEl; e++ ) {
      GeomElement* geomEl = elements[e];
      const int numQP = geomEl->getNumberOfQuadPoints();
      Real AvgMatProp = 0.0;
      for ( int q = 0; q < numQP; q++ ) {
        vector<Real > MatProp = _materials[e*numQP + q]->getMaterialParameters();
        if (!MatProp.empty())
        AvgMatProp += MatProp[1];
      }
      out << AvgMatProp/double(numQP) << endl;

    }

    // Material internal variable
    out << "SCALARS InternalVariable double" << endl
    << "LOOKUP_TABLE default" << endl;
    for ( int e = 0; e < NumEl; e++ ) {
      vector<Real > IntProp = _materials[e]->getInternalParameters();
      // for ( int i = 0; i < IntProp.size(); i++ ) {
      /*
      WARNING
      THIS ONLY PRINTS THE FIRST INTERNAL VARIABLE
      BAD - NEED TO BE CHANGED!!
      */
      if ( IntProp.size() > 0 ) {
        out << IntProp[0] << endl;
      } else {
        out << 0.0 << endl;
      }
    }

    // Material stress tensor
    out << "TENSORS P double" << endl;

    // Loop through elements, also through material points array, which is unrolled
    // uint index = 0;
    MechanicsMaterial::FKresults FKres;
    FKres.request = 2;
    for(int e = 0; e < NumEl; e++)
    {
      GeomElement* geomEl = elements[e];
      const int numQP = geomEl->getNumberOfQuadPoints();

      // F at each quadrature point are computed at the same time in one element
      vector<Matrix3d > Flist(numQP, Matrix3d::Zero());
      // Compute deformation gradients for current element
      this->computeDeformationGradient(Flist, geomEl);

      // Loop over quadrature points
      // for(int q = 0; q < numQP; q++) {
      // _materials[e]->compute(FKres, Flist[q], &Fiber);
      /*
      WARNING
      THIS ONLY WORKS WITH 1QP PER ELEMENT
      BAD - NEED TO BE CHANGED!!
      */
      // }
      // Only print stress tensor at the first QP !!!
      _materials[e*numQP]->compute(FKres, Flist[0]);
      out << FKres.P(0,0) << " " <<  FKres.P(0,1) << " " << FKres.P(0,2) << endl <<
      FKres.P(1,0) << " " <<  FKres.P(1,1) << " " << FKres.P(1,2) << endl <<
      FKres.P(2,0) << " " <<  FKres.P(2,1) << " " << FKres.P(2,2) << endl;
    }

    // Close file
    out.close();
  } // writeOutput

} // namespace voom
