#include "FoundationModel.h"

namespace voom {

  // Constructor
  FoundationModel::FoundationModel(Mesh* aMesh, const uint NodeDoF, const string SpNodes, Real SpringK):
  Model(aMesh, NodeDoF), _springK(SpringK)
  {
    _field.resize(  (_myMesh->getNumberOfNodes() )*_nodeDoF );
    this->initializeField();
    this->initSpringBC(SpNodes);
  }

  // Compute Function - Compute Energy, Force, Stiffness
  void FoundationModel::compute(Result * R)
  {
    const vector<GeomElement* > elements = _myMesh->getElements();
    const int AvgNodePerEl = ((elements[0])->getNodesID()).size();
    const int NumEl = elements.size();
    const int dim = _myMesh->getDimension();
    vector<Triplet<Real > > KtripletList;

    int PbDoF = R->getPbDoF();

    // Reset values in result struct
    // Todo: Once we have the wrapper, we won't want to do this
    if ( R->getRequest() & ENERGY ) {
      R->setEnergy(0.0);
    }
    if ( R->getRequest() & FORCE)  {
      R->resetResidualToZero();
    }
    if ( R->getRequest() & STIFFNESS ) {
      // R->resetStiffnessToZero();
      KtripletList.reserve(dim*dim*NumEl*AvgNodePerEl*AvgNodePerEl);
    }

    // Loop through nodes
    MechanicsMaterial::FKresults FKres;
    FKres.request = R->getRequest();

    //** Insert stuff from the spring here **//
    vector <Triplet<Real > > KtripletList_FromSpring;

    // Loop through _spNodes
    for(int n = 0; n < _spNodes.size(); n++)
    {
      int NodeID = _spNodes[n];
      Vector3d xa_prev, xa_curr;
      xa_prev << _prevField[NodeID*3], _prevField[NodeID*3+1], _prevField[NodeID*3+2];
      xa_curr << _field[NodeID*3], _field[NodeID*3+1], _field[NodeID*3+2];

      // Compute energy
      if (R->getRequest() & ENERGY) {
        R->addEnergy( 0.5* _springK*pow( (xa_curr - xa_prev).dot(_spNormals[n]), 2.0) );
      }

      // Compute Residual
      if ((R->getRequest() & FORCE)) {
        for(uint i = 0; i < 3; i++) {
          R->addResidual(NodeID*3+i,  _springK*_spNormals[n](i)*(xa_curr - xa_prev).dot(_spNormals[n]) );
        } // i loop
      } // Internal force loop

      // Compute stiffness matrix
      if ( R->getRequest() & STIFFNESS ) {
        for(uint i = 0; i < 3; i++) {
          for(uint j = 0; j < 3; j++) {
            KtripletList_FromSpring.push_back(Triplet<Real >( NodeID*3+i, NodeID*3+j, _springK*_spNormals[n](i)*_spNormals[n](j) ));
          } // j loop
        } // i loop
        KtripletList.insert( KtripletList.end(), KtripletList_FromSpring.begin(), KtripletList_FromSpring.end() );
        R->setStiffnessFromTriplets(KtripletList);
        R->FinalizeGlobalStiffnessAssembly();
      } // Stiffness loop
    } // Spring nodes loop
  } // Compute Mechanics Model


  void FoundationModel::computeNormals() {

    // First compute normal of any element in _myMesh
    const vector<GeomElement* > elements = _myMesh->getElements();

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
    }

  } // compute Normals



  void FoundationModel::initSpringBC(const string SpNodes) {
    // Store node number on the outer surface
    ifstream inp(SpNodes.c_str());
    int nodeNum = 0;
    while (inp >> nodeNum)
    _spNodes.push_back(nodeNum);

    // Collect elements that share a node in _spNodes
    const vector<GeomElement* > elements = _myMesh->getElements();

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

    } // loop over _spNodes

    _spNormals.resize(_spNodes.size(), Vector3d::Zero());
    // Compute initial node normals
    this->computeNormals();

  } // InitSpringBC


  // Writing output
  void FoundationModel::writeOutputVTK(const string OutputFile, int step) {

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
    out << endl << "POINT_DATA " << NumNodes << endl << "VECTORS displacements double" << endl;

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
