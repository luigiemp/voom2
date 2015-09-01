//-*-C++-*-
#ifndef __LoopShellModel_h__
#define __LoopShellModel_h__

#include "EllipticModel.h"
#include "SCElastic.h"
#include "ShellGeometry.h"
#include "EigenEllipticResult.h"

namespace voom{

  // Model Results
  class LoopShellModel: public EllipticModel {

  public:

    //! Basic Constructor
    /*! Construct from basic data structures defining the mesh, materials, BCs. 
     */
    LoopShellModel(Mesh* aMesh, vector<SCElastic * > _materials, 
		   const uint NodeDoF);
  
    //! Input-file-based Constructor
    // LoopShellModel(Mesh* myMesh, const string inputFile, const uint NodeDoF);

    //! Destructor
    ~LoopShellModel() {
      set<SCElastic *> UNIQUEmaterials;
      for (uint i = 0; i < _materials.size(); i++) 
	UNIQUEmaterials.insert(_materials[i]);
	
      for (set<SCElastic *>::iterator it = UNIQUEmaterials.begin();
	   it != UNIQUEmaterials.end(); it++) 
		delete (*it);
    };



    //! Initialize field
    // From constant value
    void initializeField(const Real value = 1.0) {
      const uint numNodes = _myMesh->getNumberOfNodes();
      const uint dim = _myMesh->getDimension();
      if (_nodeDoF == dim ){
	//Case: Only displacement DOF
	for (uint i = 0; i < numNodes; i++) 
	  for (uint j = 0; j < dim; j++)
	    _field[i*dim+j] = _myMesh->getX(i,j)*value; // value = isotropic expansion/shrinking
      }
      else if (_nodeDoF == 1){
	// Case: Only scalar field DOF
	for (uint i = 0; i < numNodes; i++) 
	  _field[i] = 0; 
      }
      else if (_nodeDoF == 4){
	// Case: Displacement and scalar DOFs
	for (uint i = 0; i < numNodes; i++){ 
	  for (uint j = 0; j < dim; j++){
	    _field[i*(dim+1)+j] = _myMesh->getX(i,j)*value; // value = isotropic expansion/shrinking      
	  }
	  _field[i*(dim+1)+dim] = 0;
	}
      }
      else{
	cout << "Unknown initalization for _field"<<endl;
	exit(-1);
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
      _field[id*dim + dof] += value;
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
      // Assume dim = 3 - need to be changed
      int i = 0;
      while (i < _field.size()) {
	for (uint j = 0; j < 3; j++) {
	  cout << _field[i] << " ";
	  i++;
	}
	cout << endl;
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
      set<SCElastic *> UNIQUEmaterials;
      for (uint i = 0; i < _materials.size(); i++) 
	UNIQUEmaterials.insert(_materials[i]);
	
      return UNIQUEmaterials.size();
    }

    uint getTotNumMatProp() {
      set<SCElastic *> UNIQUEmaterials;
      for (uint i = 0; i < _materials.size(); i++) 
	UNIQUEmaterials.insert(_materials[i]);
	
      return ( UNIQUEmaterials.size() * (_materials[0]->getMaterialParameters()).size() );
    }

    vector<SCElastic * > getMaterials() {
      return _materials;
    }
    
    //! Write output
    void writeOutputVTK(const string OutputFile, int step); 

    //! Solve the system
    void compute(EllipticResult & R);


  protected:
    //! Compute Gradients
    //void computeGradients(vector<ShellGeometry> defGeomVec, GeomElement* geomEl);

    //! List of Material data at each element in the model
    // (need to be modified for history dependent materials, e.g. plasticity)
    vector<SCElastic * > _materials;

    //! Solution value at all nodes, local and ghost
    //! Displacement are stored unrolled, [phi_x, phi_y, phi_z]

    //vector<Vector3d> _displacements;

    vector<Real > _field;
    vector<Real > _prevField;
  };

} // namespace voom

#endif
