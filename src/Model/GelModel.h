//-*-C++-*-
#ifndef __GelModel_h__
#define __GelModel_h__
#include "Model.h"
#include "Result.h"
#include "FilamentMaterial.h"
#include "GelMesh.h"
#include "Filament.h"
#include "CrossLink.h"
#include "GelInput.h"
#include "PeriodicBox.h"

namespace voom{

  // Model Results
  class GelModel: public Model {

  public:

    //! Basic Constructor
    /*! Construct from basic data structures defining the mesh, materials, BCs. 
     */
    GelModel(GelMesh* aGelMesh,GelInput* input,vector<FilamentMaterial * > _springs,
	     vector<FilamentMaterial * > _angleSprings,
	     const uint NodeDoF,
	     PeriodicBox* box,
	     int NodalForcesFlag = 0,
	     int _resetFlag = 1);
    
		     
    //! Destructor
    
    ~GelModel() {
      set<FilamentMaterial *> UNIQUEmaterials;
      for (uint i = 0; i < _springs.size(); i++) 
	UNIQUEmaterials.insert(_springs[i]);
      for (uint i = 0; i < _angleSprings.size(); i++)
        UNIQUEmaterials.insert(_angleSprings[i]);

      for (set<FilamentMaterial *>::iterator it = UNIQUEmaterials.begin();
	   it != UNIQUEmaterials.end(); it++) 
		delete (*it);
      set<Filament *> UNIQUEfilaments;
      for (uint i = 0; i < _filaments.size(); i++) 
	UNIQUEfilaments.insert(_filaments[i]);

      for (set<Filament *>::iterator it = UNIQUEfilaments.begin();
	   it != UNIQUEfilaments.end(); it++) 
		delete (*it);

      set<CrossLink *> UNIQUEcrosslinks;
      for (uint i = 0; i < _crosslinks.size(); i++) 
	UNIQUEcrosslinks.insert(_crosslinks[i]);

      for (set<CrossLink *>::iterator it = UNIQUEcrosslinks.begin();
      	   it != UNIQUEcrosslinks.end(); it++) 
	delete (*it);
	
    };
    
    // Get number of Nodes:
    int getNumberOfNodes(){return _myGelMesh->getNumberOfNodes();}

    //! Initialize field
    // From constant value
    void initializeField(const Real value = 1.0) {
      const uint numNodes = _myGelMesh->getNumberOfNodes();
      const uint dim = _myGelMesh->getDimension();

      for (uint i = 0; i < numNodes; i++) 
	for (uint j = 0; j < dim; j++)
	  _field[i*dim+j] = _myGelMesh->getX(i,j)*value; // value = isotropic expansion/shrinking
    };

    //! From array
    void initializeField(const Real* value) {
      _field.assign(value, value+_field.size());
    };

    //! Linearized update
    void linearizedUpdate(const Real* localValues, Real fact) {
      const int nLocalDof = (_myGelMesh->getNumberOfNodes())*_nodeDoF;
      for(uint i = 0; i < nLocalDof; i++)
	_field[i] += fact*localValues[i];
    };

    // One value at the time (Node ID, dof index, value)
    void linearizedUpdate(const int id, const int dof, const Real value) {
      // const uint dim = _myGelMesh->getDimension();
      // assert( id < _field.size() && dof < dim );
      _field[id*_nodeDoF + dof] += value;
    }

    // One value at the time (Node ID, dof index, value)
    void linearizedUpdate(const int dof, const Real value) {
      _field[dof] += value;
    }
    
    void setField(uint dof, Real value) {
      _field[dof] = value;
    }
    void setField(const Real* value) {
      _field.assign(value, value+_field.size());
    };
    void getField(vector<double> & x) {
      assert(x.size() == _field.size());
      x = _field;
    }

    VectorXd getX()
    {
      VectorXd X =_myGelMesh->getX();
      return X;
    }
    
    Real getX(const uint nodeId, const uint dof){
      return _myGelMesh->getX(nodeId,dof);
    }
    
    void setPrevField() {
      _prevField = _field;
    };

    void printField() {
      int i = 0;
      while (i < _field.size()) {
	for (uint j = 0; j < _nodeDoF; j++) {
	  cout << _field[i] << " ";
	  i++;
	}
	cout << endl;
      }
    }

    void writeField(const string OutputFile, int step);

    uint getNumMat() {
      set<FilamentMaterial *> UNIQUEmaterials;
      for (uint i = 0; i < _springs.size(); i++) 
	UNIQUEmaterials.insert(_springs[i]);
      for (uint i = 0; i < _angleSprings.size(); i++)
        UNIQUEmaterials.insert(_angleSprings[i]);

      return UNIQUEmaterials.size();
    }

    uint getTotNumMatProp() {
      set<FilamentMaterial *> UNIQUEmaterials;
      for (uint i = 0; i < _springs.size(); i++)
        UNIQUEmaterials.insert(_springs[i]);
      for (uint i = 0; i < _angleSprings.size(); i++)
	UNIQUEmaterials.insert(_angleSprings[i]);

      return UNIQUEmaterials.size();

    }

    vector<FilamentMaterial * > getSprings() {
      return _springs;
    }

    void setResetFlag(int ResetFlag) {
      _resetFlag = ResetFlag;
    }

    
    void setNodalForcesFlag(int NodalForcesFlag) {
      _nodalForcesFlag = NodalForcesFlag;
    }
    
    //! Write output
    void writeOutputVTK(const string OutputFile, int step); 
    
    void writeX(const string OutputFile, int step); 

    void writeClConnectivity(const string OutputFile, int step);
    
    void writeFilConnectivity(const string OutputFile, int step);

    //! Solve the system
    void compute(Result & R);

    
    // Update nodal forces
    void updateNodalForces(vector<int > * ForcesID, vector<Real > * Forces) {
      _forcesID = ForcesID; _forces = Forces;
    }

    
    void getFilamentx(vector<Vector3d > & xlist, const vector<int> & NodesID);

    void addCrosslink(CrossLink * crosslink){
      _crosslinks.push_back(crosslink);
    }

    int getDimension(){return _myGelMesh->getDimension();}
    
    void deleteCrosslink(int i){
      _crosslinks.erase(_crosslinks.begin()+i);
    }

    //https://en.wikipedia.org/wiki/Erase%E2%80%93remove_idiom
    void deleteCrosslink(CrossLink * cl){
      delete(cl);
      _crosslinks.erase(remove(_crosslinks.begin(),_crosslinks.end(),cl),_crosslinks.end());
    }

    int getNumberOfCl(){ return _crosslinks.size();}

    int getNumberOfFil(){ return _filaments.size();}
    
    vector<CrossLink * > getCrossLinks(){return _crosslinks;}

  protected:
    
    //! List of Material data at each element in the model
    vector<FilamentMaterial * > _springs;
    vector<FilamentMaterial * > _angleSprings;
    
    vector<Real > _field;
    vector<VectorXd> _X0;

    GelInput* _input;
    
    vector<Filament * > _filaments;
    vector<CrossLink * > _crosslinks; 

    int _nodalForcesFlag;
    vector<int > * _forcesID;
    vector<Real > * _forces;
    GelMesh* _myGelMesh;
    vector<Real > _prevField;
    PeriodicBox* _box;
    int _resetFlag;
  };

} // namespace voom

#endif
