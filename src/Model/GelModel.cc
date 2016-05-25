#include "GelModel.h"
#include "Filament.h"
namespace voom {

  // Constructor: read the gelmesh to construct filament objects with material prop
  GelModel::GelModel(GelMesh* aGelMesh,GelInput* input, vector<FilamentMaterial * > springs,
		     vector<FilamentMaterial * > angleSprings,
		     const uint NodeDoF,
		     PeriodicBox* box,
		     int NodalForcesFlag,
		     int ResetFlag):
    Model( NodeDoF),_myGelMesh(aGelMesh) , _springs(springs),_angleSprings(angleSprings), 
    _nodalForcesFlag(NodalForcesFlag),_box(box),_resetFlag(ResetFlag), _input(input)
  {
    // THERE IS ONE MATERIAL PER ELEMENT - CAN BE CHANGED - DIFFERENT THAN BEFORE
    // Resize and initialize (default function) _field vector

    _field.resize(  (_myGelMesh->getNumberOfNodes() )*_nodeDoF );
    this->initializeField();

    const vector<GelElement* > filElement = _myGelMesh->getElements();
    uint NumFil = _myGelMesh->getNumberOfElements();

    // loop over the filament mesh elements to construct Filament object with material prop:
    for (int nFil = 0 ; nFil < NumFil ; nFil++){
      // Construct new filament with material, add to _filaments:
      Filament* Fil =  new Filament(filElement[nFil], _springs[nFil],_angleSprings[nFil],_box);
      _filaments.push_back(Fil);
    }
    
  }
  
  //TO DO: transfer this function to filament and crosslink classes
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
    uint NumFil = _filaments.size();

    if ( R.getRequest() & ENERGY ) {
      R.setEnergy(0.0);
    }
    if ( R.getRequest() & FORCE )  {
      R.resetResidualToZero();
    } 

    //Compute filament stretching and bending energy
    for (int nFil = 0 ; nFil < NumFil ; nFil++){
      // Get the filament nodes //TO DO: transfer this to filament class
      const vector<int> & FilNodesID = _filaments[nFil]->getNodesID();
      int FilNodeNum = FilNodesID.size();
      // Get the current coords
      vector<Vector3d > Filxlist(FilNodeNum,Vector3d::Zero());
      this->getFilamentx(Filxlist,FilNodesID);
      
      //compute energy
      _filaments[nFil]->compute(R, Filxlist);
      
    }

    //Compute cross link stretching energy
    // Crosslinks use the same classes as filament, the only difference
    // here is that they have only 2 nodes and their bending energy is zero

    uint NumCl = _crosslinks.size();
    
    for (int nCl = 0 ; nCl < NumCl ; nCl++){
      // Get the crosslink nodes //TO DO: transfer this to crosslink class
      const vector<int> & ClNodesID = _crosslinks[nCl]->getNodesID();
      int ClNodeNum =  _crosslinks[nCl]->getNumberOfNodes();
      vector<Vector3d > Clxlist(ClNodeNum,Vector3d::Zero());
      
      // Get the current coords
      this->getFilamentx(Clxlist,ClNodesID);

      //compute energy
      _crosslinks[nCl]->compute(R, Clxlist);
      

      //Find crosslink that map accross the top/bottom boundary (for periodic BC)
      Vector3d dx = Clxlist[1]-Clxlist[0];

      if(_box->crossTop(dx)){
	
      }

    }
    

  }

  void GelModel::writeField(const string OutputFile, int step) {
      // Create outputFile name
      stringstream FileNameStream;
      FileNameStream <<_input->getValueOfKey<string>("savePath","results/") << OutputFile << step << ".node";
      ofstream out;
      out.open( (FileNameStream.str()).c_str() );
      int NumNode = _myGelMesh->getNumberOfNodes();
      int dim = _myGelMesh->getDimension();
      
      for (uint i = 0; i < NumNode; i++) {
	for(uint k = 0 ; k < dim; k++){
	  out << setprecision(15) << _field[i*dim+k] << " ";
	}
	out << endl;
      }
      out.close();
    }

  void GelModel::writeClConnectivity(const string OutputFile, int step){
     // Create outputFile name
    stringstream FileNameStream;
    FileNameStream <<_input->getValueOfKey<string>("savePath","results/") << OutputFile << step << ".ele";
    ofstream out;
    out.open( (FileNameStream.str()).c_str() );
    int NumCl = _crosslinks.size();
    for(int i = 0; i<NumCl; i++){
      out << _crosslinks[i]->getNodesID()[0] << " " <<  _crosslinks[i]->getNodesID()[1] << endl;
    }

  }

   void GelModel::writeFilConnectivity(const string OutputFile, int step){
     // Create outputFile name
    stringstream FileNameStream;
    FileNameStream <<_input->getValueOfKey<string>("savePath","results/") << OutputFile << step << ".ele";
    ofstream out;
    out.open( (FileNameStream.str()).c_str() );
    int NumCl = _filaments.size();
    for(int i = 0; i<NumCl; i++){
      int NumNodeFil = _filaments[i]->getNodesID().size();
      for(int k = 0 ; k< NumNodeFil ; k++){
	out << _filaments[i]->getNodesID()[k] << " ";
      }
      out << endl;
    }

  }

  void GelModel::writeX(const string OutputFile, int step) 
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
  }
  
  
  // Writing output
  void GelModel::writeOutputVTK(const string OutputFile, int step) 
  {
   
  } // writeOutput
    






} // namespace voom
