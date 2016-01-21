//-*-C++-*-
#ifndef __LBModel_h__
#define __LBModel_h__

#include "EllipticModel.h"
#include "LandauBrazovskii.h"
#include "ShellGeometry.h"
#include "EigenEllipticResult.h"

namespace voom{

  // Model Results
  class LBModel: public EllipticModel {

  public:

    //! Basic Constructor
    /*! Construct from basic data structures defining the mesh, materials, BCs. 
     */
    LBModel(Mesh* aMesh, vector<LandauBrazovskii * > _materials, 
		   const uint NodeDoF);
  
    //! Input-file-based Constructor
    // LBModel(Mesh* myMesh, const string inputFile, const uint NodeDoF);

    //! Destructor
    ~LBModel() {
      set<LandauBrazovskii *> UNIQUEmaterials;
      for (uint i = 0; i < _materials.size(); i++) 
	UNIQUEmaterials.insert(_materials[i]);
	
      for (set<LandauBrazovskii *>::iterator it = UNIQUEmaterials.begin();
	   it != UNIQUEmaterials.end(); it++) 
		delete (*it);
    };



    //! Initialize field
    // From constant value
    void initializeField(const Real value) {
      const uint numNodes = _myMesh->getNumberOfNodes();
      const uint dim = _myMesh->getDimension();
      // Case: Only scalar field DOF
      for (uint i = 0; i < numNodes; i++) 
	_field[i] =  value;       
    };

    // Random
    void initializeField() {
      const uint numNodes = _myMesh->getNumberOfNodes();
      const uint dim = _myMesh->getDimension();
      // Case: Only scalar field DOF
      for (uint i = 0; i < numNodes; i++) {
	_field[i] = 0.01*(0.5-double(rand())/double(RAND_MAX));   
	//_prevField[i] = 0.0; //*********************************** Implicit Time stepping *****************
      }
    };

    //! From array
    void initializeField(const Real* value) {
      _field.assign(value, value+_field.size());
    };



    //! Linearized update
    void linearizedUpdate(const Real* localValues, Real fact) {
      const int nLocalDof = (_myMesh->getNumberOfNodes())*_nodeDoF;
      for(uint i = 0; i < nLocalDof; i++)
	_field[i] += fact*localValues[i];
    };

    // One value at the time (Node ID, dof index, value)
    void linearizedUpdate(const int id, const int dof, const Real value) {
      const uint dim = _myMesh->getDimension();
      // assert( id < _field.size() && dof < dim );
      _field[id] += value;
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
    void setPrevField() {
      _prevField = _field;
    };

    void printField() {
      int i = 0;
      while (i < _field.size()) {
	cout << _field[i] << endl;
	i ++;
      }
    }

    void writeField(string FileName) {
      ofstream out;
      out.open( FileName.c_str() );
      out << _field.size() << endl;
      for (uint i = 0; i < _field.size(); i++) {
	out << setprecision(15) << _field[i] << endl;
      }
      out.close();
    }

    uint getNumMat() {
      set<LandauBrazovskii *> UNIQUEmaterials;
      for (uint i = 0; i < _materials.size(); i++) 
    	UNIQUEmaterials.insert(_materials[i]);
	
      return UNIQUEmaterials.size();
    }

    uint getTotNumMatProp() {
      // set<LandauBrazovskii *> UNIQUEmaterials;
      // for (uint i = 0; i < _materials.size(); i++) 
      // 	UNIQUEmaterials.insert(_materials[i]);
	
      // return ( UNIQUEmaterials.size() * (_materials[0]->getMaterialParameters()).size() );
    }

    vector<LandauBrazovskii * > getMaterials() {
      return _materials;
    }
    
    //! Write output
    void writeOutputVTK(const string OutputFile, int step); 

    //! Solve the system
    void compute(EllipticResult & R);
    
    void setImplicitDynamicsFlag(bool flag){
      
      _implicitDynamics == flag;
      if(_implicitDynamics == true)
      _prevField.resize( _myMesh->getNumberOfNodes()  ); //********************** Implicity time stepping **************
    
    }

  protected:
    //! Compute Gradients
    void computeGradients(vector<ShellGeometry> defGeomVec, GeomElement* geomEl);

    //! List of Material data at each element in the model
    // (need to be modified for history dependent materials, e.g. plasticity)
    vector<LandauBrazovskii * > _materials;

    //! Solution value at all nodes, local and ghost
    //! Displacement are stored unrolled, [phi_x, phi_y, phi_z]

    //vector<Vector3d> _displacements;

    vector<Real > _field;
    vector<Real > _prevField;

    bool _implicitDynamics;
  };

} // namespace voom

#endif
