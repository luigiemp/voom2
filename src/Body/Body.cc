//-*-C++-*-
#include "Body.h"

namespace voom {

  // Consistency Checks //
  void Body::checkConsistency(Result* R, Real perturbationFactor, int request,
			      Real h, Real tol)
  {
    // Find nodes in this body only
    vector<GeomElement* > Elements = _myMesh->getElements();
    set<int > UniqueBodyNodes;
    for (int e = 0; e < _myMesh->getNumberOfElements(); e++) {
      vector<int > ElNodesID = Elements[e]->getNodesID();
      for (int n = 0; n < ElNodesID.size(); n++) {
	UniqueBodyNodes.insert(ElNodesID[n]);
      }
    }

   

    // Perturb field randomly to change from reference configuration
    // Save perturbed field to set the configuration back to reference 
    // at the end of the test
    srand( time(NULL) );
    map<int, Real> perturb;
    for (set<int >::iterator it = UniqueBodyNodes.begin(); it != UniqueBodyNodes.end(); it++) {
      // cout << "Unique node = " << *it << endl;
      for (int i = 0; i < _myState->getNodeDof(*it); i++) {
	Real randomNum = perturbationFactor*(Real(rand())/RAND_MAX - 0.5);
	perturb.insert(pair<int, Real>(_myState->getGdof(*it) + i, randomNum));
	_myState->linearizedUpdate(_myState->getGdof(*it) + i, randomNum);
      }
    }
    


    // Force Check
    if ( request & FORCE ) {
      Real error = 0.0, norm = 0.0;

      R->setRequest(FORCE); // First compute forces numerically
      R->resetResults(FORCE);
      this->compute(R);

      R->setRequest(ENERGY); // Reset result request so that only energy is computed 
      R->resetResults(ENERGY);
      this->compute(R);

      cout << "Body energy at test start = " <<  R->getEnergy() << endl;

      for (set<int >::iterator it = UniqueBodyNodes.begin(); it != UniqueBodyNodes.end(); it++) {
	for (int i = 0; i < _myState->getNodeDof(*it); i++) {
	  int Gdof = _myState->getGdof(*it) + i;
	  // Perturb +h
	  _myState->linearizedUpdate(Gdof, h);
	  R->resetResults(ENERGY);
	  this->compute(R);
	  Real Wplus = R->getEnergy();
	  
	  // Perturb -2h
	  _myState->linearizedUpdate(Gdof, -2.0*h);
	  R->resetResults(ENERGY);
	  this->compute(R);
	  Real Wminus = R->getEnergy();
	  
	  // Bring back to original position
	  _myState->linearizedUpdate(Gdof, h);
	  
	  error += pow( (Wplus-Wminus)/(2.0*h) - 
			R->getResidual(Gdof), 2);
	  norm  += pow(R->getResidual(Gdof), 2);
	} // Loop over DoF
      }
      error = sqrt(error);
      norm  = sqrt(norm);

      if ( abs(error) < norm * tol) {
	cout << "** Elliptic Body Force consistency check PASSED" << endl;
	cout << "** Error: " << error << " Norm: " << norm  << " Norm*tol: " << norm*tol << endl;
      }
      else {
	cout << "** Elliptic Body Force consistency check FAILED" << endl;
	cout << "** Error: " << error << " Norm: " << norm << " Norm*tol: " << norm*tol << endl;
      }
    } // Check Forces loop



    // Stiffness check
    if ( request & STIFFNESS ) {
      Real error = 0.0, norm = 0.0;

      R->setRequest(STIFFNESS); // First compute stiffness numerically
      R->resetResults(STIFFNESS);
      this->compute(R);
      R->FinalizeGlobalStiffnessAssembly();

      R->setRequest(FORCE); // Reset result request so that only forces are computed 


      for (set<int >::iterator itA = UniqueBodyNodes.begin(); itA != UniqueBodyNodes.end(); itA++) {
	for (int i = 0; i < _myState->getNodeDof(*itA); i++) {
	  for (set<int >::iterator itB = UniqueBodyNodes.begin(); itB != UniqueBodyNodes.end(); itB++) {
	    for (int j = 0; j < _myState->getNodeDof(*itB); j++) {
	      // Perturb +h
	      int GdofI = _myState->getGdof(*itA) + i;
	      int GdofJ = _myState->getGdof(*itB) + j;
	      _myState->linearizedUpdate(GdofJ, h);
	      R->resetResults(FORCE);
	      this->compute(R);
	      Real Fplus = R->getResidual(GdofI);
	  
	      // Perturb -2h
	      _myState->linearizedUpdate(GdofJ, -2.0*h);
	      R->resetResults(FORCE);
	      this->compute(R);
	      Real Fminus = R->getResidual(GdofI);
	  
	      // Bring back to original position
	      _myState->linearizedUpdate(GdofJ, h);
	  
	      // Computing Error and Norm;
	      error += pow((Fplus - Fminus)/(2.*h) - R->getStiffness(GdofI, GdofJ), 2.0);
	      norm  += pow( R->getStiffness(GdofI, GdofJ), 2.0); 
	      // cout << R->getStiffness(GdofI, GdofJ) << "\t" << (Fplus - Fminus)/(2.*h) << endl;
	    } // loop over dof of node j
	  } // loop over nodes
	} // loop over dof of node i
      }// loop over nodes
      
      error = sqrt(error);
      norm  = sqrt(norm);
      if ( abs(error) < norm * tol) {
	cout << "** Elliptic Body Hessian consistency check PASSED" << endl;
	cout << "** Error: " << error << " Norm: " << norm << " Norm*tol: " << norm*tol << endl << endl;
      }
      else {
	cout << "** Elliptic Body Hessian consistency check FAILED" << endl;
	cout << "** Error: " << error << " Norm: " << norm << " Norm*tol: " << norm*tol << endl << endl;
      }
    } // Check Stiffness    

    // Reset field to initial values
    for (set<int >::iterator it = UniqueBodyNodes.begin(); it != UniqueBodyNodes.end(); it++) {
      for (int i = 0; i < _myState->getNodeDof(*it); i++) {
	int Gdof = _myState->getGdof(*it) + i;
	_myState->linearizedUpdate(Gdof, -perturb[Gdof]);
      }
    }
    
  } // Check consistency



  // Writing output
  void Body::writeOutputVTK(const string OutputFile, int step)
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

  } // writeOutput

  
} // Namespace voom
