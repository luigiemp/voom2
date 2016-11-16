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
    // Need to fill this in

  } // writeOutput

} // namespace voom
