#include "GelModel.h"
#include "Filament.h"
namespace voom {

  // Constructor: read the gelmesh to construct filament objects with material prop
  GelModel::GelModel(GelMesh* aGelMesh, vector<FilamentMaterial * > springs,
		     vector<FilamentMaterial * > angleSprings,
		     const uint NodeDoF,
		     int NodalForcesFlag,
		     int ResetFlag):
    Model( NodeDoF),_myGelMesh(aGelMesh) , _springs(springs),_angleSprings(angleSprings), 
    _nodalForcesFlag(NodalForcesFlag), _resetFlag(ResetFlag)
  {
    // THERE IS ONE MATERIAL PER ELEMENT - CAN BE CHANGED - DIFFERENT THAN BEFORE
    // Resize and initialize (default function) _field vector

    _field.resize(  (_myGelMesh->getNumberOfNodes() )*_nodeDoF );
    this->initializeField();

    const vector<GelElement* > filElement = _myGelMesh->getElements();
    uint NumFil = _myGelMesh->getNumberOfElements();

    // loop over the filament mesh elements to construct Filament object with materials:
    for (int nFil = 0 ; nFil < NumFil ; nFil++){
      // Construct new filament with material, add to _gel:
      Filament* Fil =  new Filament(filElement[nFil], _springs[nFil],_angleSprings[nFil]);
      _gel.push_back(Fil);
    }
    
  }
  

  void GelModel::getFilamentx(vector<Vector3d > & xlist,const vector<int> & NodesID)
  {
    // Get all filament node current position
    const uint dim   = _myGelMesh->getDimension();
    const uint nodeNum = NodesID.size();
    
    for(uint n = 0; n < nodeNum; n++)
      {
	for(uint i = 0; i < dim; i++)
	  {
	    xlist[n](i) = _field[NodesID[n]*dim+i];
	  }
      }
  }
  
  void GelModel::compute(Result & R)
  {
    //const vector<GelElement* > filElement = _myGelMesh->getFilaments();
    uint NumFil = _gel.size();

    if ( R.getRequest() & ENERGY ) {
      R.setEnergy(0.0);
    }
    if ( R.getRequest() & FORCE )  {
      R.resetResidualToZero();
    } 

    FilamentMaterial::Filresults Rf;
    Rf.request = R.getRequest();

    //Compute filament stretching and bending energy
    for (int nFil = 0 ; nFil < NumFil ; nFil++){
      // Get the filament nodes
      const vector<int> & FilNodesID = _gel[nFil]->getNodesID();
      int FilNodeNum = FilNodesID.size();
      // Get the current coords
      vector<Vector3d > Filxlist(FilNodeNum,Vector3d::Zero());
      this->getFilamentx(Filxlist,FilNodesID);
      
      //compute energy
      _gel[nFil]->compute(R, Filxlist);
      
    }

    //Compute cross link stretching energy
    // Crosslinks use the same classes as filament, the only difference
    // here is that they have only 2 nodes and their bending energy is zero

    uint NumCl = _crosslinks.size();
    
    for (int nCl = 0 ; nCl < NumCl ; nCl++){
      // Get the crosslink nodes
      const vector<int> & ClNodesID = _crosslinks[nCl]->getNodesID();
      int ClNodeNum =  _crosslinks[nCl]->getNumberOfNodes();
      vector<Vector3d > Clxlist(ClNodeNum,Vector3d::Zero());
      // Get the current coords
      this->getFilamentx(Clxlist,ClNodesID);
      //compute energy
      _crosslinks[nCl]->compute(R, Clxlist);
      
    }
    

  }

  void GelModel::writeOutput(const string OutputFile, int step) 
  {
    // Create outputFile name
    stringstream NodeFileNameStream;
    NodeFileNameStream << OutputFile << step << ".node";
    ofstream Nodeout;
    Nodeout.open( (NodeFileNameStream.str()).c_str() );
    int NumNodes = _myGelMesh->getNumberOfNodes();
    int dim = _myGelMesh->getDimension();

    for (int i = 0; i < NumNodes; i++ ) {
      for (int j = 0; j < dim; j++) {
	Nodeout << _myGelMesh->getX(i,j) << " ";
      }
      Nodeout << endl;
    }
    Nodeout.close();

    // stringstream ConFileNameStream;
    // ConFileNameStream << OutputFile << step << ".";
    // ofstream Conout;
    // Conout.open( (ConFileNameStream.str()).c_str() );
    // int NumNodes = _myGelMesh->getNumberOfNodes();
    // int dim = _myGelMesh->getDimension();

    // for (int i = 0; i < NumNodes; i++ ) {
    //   for (int j = 0; j < dim; j++) {
    // 	Conout << _myGelMesh->getX(i,j) << " ";
    //   }
    //   Conout << endl;
    // }
    // Conout.close();
  }
  
  
  // Writing output
  void GelModel::writeOutputVTK(const string OutputFile, int step) 
  {
   
  } // writeOutput
    






} // namespace voom
