#include "TorsionalBody.h"

namespace voom {

  // Constructor
  TorsionalBody::TorsionalBody(Mesh* myMesh, State* myState, 
			       const vector<int > & BodyNodes, 
			       Real TorK, Vector3d Center):
    Body(myMesh, myState), _bodyNodes(BodyNodes), _torK(TorK), _center(Center)
  {
    int BodyNodesSize = _bodyNodes.size();
    _Xprev.reserve(BodyNodesSize);
    for (int i = 0; i < BodyNodesSize; i++) {
      _Xprev.push_back( _myState->getX(_bodyNodes[i]) );
    }
  };





  void TorsionalBody::updateXprev() {
    for(int i = 0; i < _bodyNodes.size(); i++) {
      for (int j = 0; j < 3; j++) {
	_Xprev[i](j) = _myState->getPhi(_bodyNodes[i], j);
      }
    }
  }





  // Compute Function - Compute Energy, Force, Stiffness
  void TorsionalBody::compute(Result* R)
  {
    // Loop through all body nodes
    Real XC = _center(0), YC = _center(1);
    for(int i = 0; i < _bodyNodes.size(); i++) {
      Vector3d BodyNodeX, BodyNodePhi, Tangent, Displacement;
      int nodeI = _bodyNodes[i];
      BodyNodeX = _Xprev[i];
      BodyNodePhi << _myState->getPhi(nodeI, 0) , _myState->getPhi(nodeI, 1) , _myState->getPhi(nodeI, 2);
      Tangent <<  BodyNodeX(0) - XC, BodyNodeX(1) - YC, 0.0;
      Displacement = BodyNodePhi - BodyNodeX;
      Real disp = Displacement.dot(Tangent);

      // Energy
      if (R->getRequest() & ENERGY) {
	R->addEnergy( 0.5 * _torK * square(disp) );
      }

      // Force
      if (R->getRequest() & FORCE) {
	R->addResidual( _myState->getGdof(nodeI), _torK * disp * Tangent(0) );
	R->addResidual( _myState->getGdof(nodeI) + 1, _torK * disp * Tangent(1) );
      }

      // Stiffness
      if (R->getRequest() & STIFFNESS) {
	for (int m = 0; m < 2; m++) {
	  for (int n = 0; n < 2; n++) {
	    R->addStiffness(_myState->getGdof(nodeI) + m, _myState->getGdof(nodeI) + n, _torK*Tangent(m)*Tangent(n));
	  }
	}
      }
	  
    } // Loop over all BodyNodes

  } // TorsionalBody::compute



 

  // Writing output
  void TorsionalBody::writeOutputVTK(const string OutputFile, int step)
  {
    /*
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
    writer->SetInput(newUnstructuredGrid);
    writer->Write();
    */
    
  } // writeOutput




} // namespace voom
