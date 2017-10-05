#include "PotentialBody.h"

namespace voom {

  // Constructor
  PotentialBody::PotentialBody(Mesh* myMesh, State* myState, 
			       PairMaterial* Potential, 
			       const vector<int > & BodyNodes, Real SearchR):
    Body(myMesh, myState), _potential(Potential), _bodyNodes(BodyNodes), _searchR(SearchR)
  {
    this->getUniqueMeshNodes();
    this->setInteractions();
    this->setNodesToEle();
    this->computeNormals();
  };

  



  // Get unique set of nodes belonging to the mesh of this Body
  void PotentialBody::getUniqueMeshNodes() {
    vector<GeomElement* > Elements = _myMesh->getElements();
    for (int e = 0; e < _myMesh->getNumberOfElements(); e++) { // Loop Thhrough all elements
      vector<int > ElNodesID = Elements[e]->getNodesID();
      copy(ElNodesID.begin(), ElNodesID.end(), inserter(_uniqueNodes, _uniqueNodes.end())); 
    }
     
  } // getUniqueMeshNodes()





  void PotentialBody::setInteractions() { // From body nodes to mesh
    
    _interactions.clear();
    // For each node in BodyNodes, get the nodes in the boundary mesh within the search radius
    for (int i = 0; i < _bodyNodes.size(); i++) {
      Vector3d BodyNode;
      int nodeI = _bodyNodes[i];
      BodyNode << _myState->getPhi(nodeI, 0) , _myState->getPhi(nodeI, 1) , _myState->getPhi(nodeI, 2);
      vector<int > OneNodeInteractions; // Per each node in the body, find the nodes on the rigid boundary
      // Loop through all mesh nodes
      for (set<int >::iterator it = _uniqueNodes.begin(); it != _uniqueNodes.end(); it++) {
	Vector3d MeshNode = _myMesh->getX(*it);
	Vector3d Diff = BodyNode - MeshNode;
	if (Diff.norm() <= _searchR) {
	  OneNodeInteractions.push_back(*it);
	  // cout << *it << " ";
	}
      }
      _interactions.push_back(OneNodeInteractions);
      // cout << endl;
    }
     
  } // setInteractions





  void PotentialBody::setNodesToEle() { 

    vector<GeomElement* > Elements = _myMesh->getElements();

    // Collect elements that share a node
    for (set<int >::iterator it = _uniqueNodes.begin(); it != _uniqueNodes.end(); it++) {
      vector<int > connected;
      for(int e = 0; e < Elements.size(); e++) {
        const vector<int >& NodesID = Elements[e]->getNodesID();
        for (int m = 0; m < NodesID.size(); m++) {
          if (*it == NodesID[m]) {
            connected.push_back(e);
	    // cout << e << " ";
            break;
          } //  check if node belong to element
        } // loop over nodes of the element
      } // loop over elements
      // cout << *it << endl;
      _nodesToEle.insert(make_pair(*it, connected));
    } // Loop over all nodes in the mesh
				   
  } // setNodesToEle





  void PotentialBody::computeNormals() { // Per each BCnode

    // First compute normal of all elements in _myMesh
    const vector<GeomElement* > Elements = _myMesh->getElements();

    vector<Vector3d > ElNormals(Elements.size(), Vector3d::Zero());

    // Loop over elements
    for(int e = 0; e < Elements.size(); e++)
    {
      GeomElement* geomEl = Elements[e];
      const vector<int  >& NodesID = geomEl->getNodesID();
      const int numQP    = geomEl->getNumberOfQuadPoints();
      const int numNodes = NodesID.size();

      // Loop over quadrature points
      for(int q = 0; q < numQP; q++) {

        // Compute normal based on Mesh X
        Vector3d a1 = Vector3d::Zero(), a2 = Vector3d::Zero(), a3 = Vector3d::Zero();
        for (int a = 0; a < NodesID.size(); a++) {
          int nodeID = NodesID[a];
          Vector3d xa = _myMesh->getX(nodeID);
          a1 += xa*geomEl->getDN(q, a, 0);
          a2 += xa*geomEl->getDN(q, a, 1);
        }
        ElNormals[e] += a1.cross(a2); // Not normalized with respect to area (elements with larger area count more)
      } // loop over QP

    } // loop over elements

    // loop over mesh nodes and assign normal
    _normals.clear();
    for (set<int >::iterator it = _uniqueNodes.begin(); it != _uniqueNodes.end(); it++) {
      // Reset normal to zero
      Vector3d TempNormal = Vector3d::Zero();
      // Loop over all elements sharing that node
      for (int e = 0; e < _nodesToEle[*it].size(); e++) {
        TempNormal += ElNormals[_nodesToEle[*it][e]];
      }
      Real normFactor = 1.0/TempNormal.norm();
      _normals.insert(make_pair(*it, TempNormal *= normFactor));

      // For testing only
      // cout << *it << " " << _normals[*it](0) << " " << _normals[*it](1) << " " << _normals[*it](2) << endl;
    }

  } // computeNormals





  // Compute Function - Compute Energy, Force, Stiffness
  void PotentialBody::compute(Result* R)
  {
    // Initialize PairMaterial results
    PairMaterial::PairMresults PMR;
    PMR.request = (R->getRequest());

    // Loop through all interactions
    for(int i = 0; i < _bodyNodes.size(); i++) {
      Vector3d BodyNode;
      int nodeI = _bodyNodes[i];
      BodyNode << _myState->getPhi(nodeI, 0) , _myState->getPhi(nodeI, 1) , _myState->getPhi(nodeI, 2);
      for(int j = 0; j < _interactions[i].size(); j++) {
	// Compute poential gap
	int nodeJ = _interactions[i][j];
	Vector3d MeshNode; MeshNode = _myMesh->getX(nodeJ);
	Vector3d GapVec, Normal; GapVec = BodyNode - MeshNode; Normal = _normals[nodeJ];
	// cout << nodeJ << " " << Normal << endl;
	Real gap; gap = GapVec.dot(Normal);
	if (gap < 0) { 
	  gap = -gap;
	  Normal *= -1.0;
	}

	// Compute pair material
	_potential->compute(PMR, gap);
	// cout << gap << endl;

	// Energy
	if (R->getRequest() & ENERGY) {
	  R->addEnergy(PMR.W);
	}

	// Force
	if (R->getRequest() & FORCE) {
	  R->addResidual( _myState->getGdof(nodeI), PMR.F*Normal(0) );
	  R->addResidual( _myState->getGdof(nodeI) + 1, PMR.F*Normal(1) );
	  R->addResidual( _myState->getGdof(nodeI) + 2, PMR.F*Normal(2) );
	}

	// Stiffness
	if (R->getRequest() & STIFFNESS) {
	  for (int m = 0; m < 3; m++) {
	    for (int n = 0; n < 3; n++) {
	      R->addStiffness(_myState->getGdof(nodeI) + m, _myState->getGdof(nodeI) + n, PMR.K*Normal(m)*Normal(n));
	    }
	  }
	}
	  
      } // Loop over all interaction nodes for nodeI
    } // Loop over all BodyNodes at the boundaries

  } // PotentialBody::compute



 

  // Writing output
  void PotentialBody::writeOutputVTK(const string OutputFile, int step)
  {
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
      points->InsertNextPoint(_myMesh->getX(LtoG[n])(0), _myMesh->getX(LtoG[n])(1), _myMesh->getX(LtoG[n])(2) );
    newUnstructuredGrid->SetPoints(points);
    


    // Element Connectivity:
    // To-do: Figure out how to handle mixed meshes
    string ElType = _myMesh->getElementType();
    VTKCellType cellType;

     if (ElType == "TD3")   // Full integration linear triangular element
      cellType = VTK_TRIANGLE;
    else if (ElType == "TD6")   // Full integration quadratic triangular element
      cellType = VTK_QUADRATIC_TRIANGLE;
    else if (ElType == "Q4")    // Full integration linear quadrilateral element
      cellType = VTK_QUAD;
    else {
      cout << "Element type not implemented in MechanicsModel writeOutput." << endl;
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



    // Write file
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(outputFileName.c_str());
    #if VTK_MAJOR_VERSION <= 5
  	writer->SetInput(newUnstructuredGrid);
    #else
  	writer->SetInputData(newUnstructuredGrid);
    #endif    
    writer->Write();





    string neighborsFileName = OutputFile + "_link" + boost::lexical_cast<string>(step) + ".vtp";
    // Create a vtkPoints object and store the points in it
    vtkSmartPointer<vtkPoints> Rpoints = vtkSmartPointer<vtkPoints>::New();

    for(int i = 0; i < _myState->getXsize(); i++) {
      double BodyNode[3] = { _myState->getPhi(i, 0) , _myState->getPhi(i, 1) , _myState->getPhi(i, 2)};
      Rpoints->InsertNextPoint(BodyNode);
    }

    // Create a cell array to store the lines in and add the lines to it
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();

    for(int i = 0; i < _bodyNodes.size(); i++) {
      for(int j = 0; j < _interactions[i].size(); j++) {
	vtkSmartPointer<vtkPolyLine> polyLine = vtkSmartPointer<vtkPolyLine>::New();
	polyLine->GetPointIds()->SetNumberOfIds(2);
	polyLine->GetPointIds()->SetId(0, _bodyNodes[i]); 
	polyLine->GetPointIds()->SetId(1, _interactions[i][j]);
	cells->InsertNextCell(polyLine);
      }
    }
      
    // Create a polydata to store everything in
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
 
    // Add the points to the dataset
    polyData->SetPoints(Rpoints);
 
    // Add the lines to the dataset
    polyData->SetLines(cells);
      
    vtkSmartPointer<vtkXMLPolyDataWriter> polyDataWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    polyDataWriter->SetFileName(neighborsFileName.c_str());
    #if VTK_MAJOR_VERSION <= 5
  	polyDataWriter->SetInput(polyData);
    #else
  	polyDataWriter->SetInputData(polyData);
    #endif    
    polyDataWriter->Write();
    
  } // writeOutput




} // namespace voom
