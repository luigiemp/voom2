#include "CrossLinker.h"
#include "Filament.h"
namespace voom {

  // Constructor: find initial crosslinks
  crosslinker::Crosslinker(GelModel* aGelModel, FilamentMaterial *  clMat):
    _myModel(aGelModel) , _clMat(clMat)
  {
    // THERE IS ONE MATERIAL PER ELEMENT - CAN BE CHANGED - DIFFERENT THAN BEFORE
    // Resize and initialize (default function) _field vector
    // _field.resize(  (_myGelMesh->getNumberOfNodes() )*_nodeDoF );
    // this->initializeField();

    // VectorXd _X0 = _myGelMesh->getX();

    // const vector<GeomFilament* > filElement = _myGelMesh->getFilaments();
    // uint NumFil = _myGelMesh->getNumberOfFilaments();
    
    // for (int nFil = 0 ; nFil < NumFil ; nFil++){

    //   const vector<int > & NodesID = filElement[nFil]->getNodesID();
    //   const uint nodeNum = NodesID.size();
    //   vector<Vector3d > xlist(nodeNum,Vector3d::Zero());

    //   getFilamentx(xlist,filElement[nFil]);
      
    //   Filament* Fil =  new Filament(filElement[nFil], _springs[nFil],_angleSprings[nFil]);

    //   _gel.push_back(Fil);
    // }
    
  }
  

  void crossLinker::updateCrosslinks()
  {
   
  }
 




} // namespace voom
