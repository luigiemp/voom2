//-*-C++-*-
#ifndef __PoissonModel_h__
#define __PoissonModel_h__

#include "Model.h"
#include "DiffusionMaterial.h"
#include "EigenResult.h"

namespace voom{

  class PoissonModel: public Model {
  public:

    //! Constructor
    PoissonModel(Mesh* myMesh, 
		 vector<DiffusionMaterial * > Materials, 
		 const uint NodeDoF);

    //! Destructor
    ~PoissonModel() {
      set<DiffusionMaterial *> UNIQUEmaterials;
      for (uint i = 0; i < _materials.size(); i++) 
	UNIQUEmaterials.insert(_materials[i]);
	
      for (set<DiffusionMaterial *>::iterator it = UNIQUEmaterials.begin();
	   it != UNIQUEmaterials.end(); it++) 
		delete (*it);
    };

    //! Initialize field
    // From constant value
    void initializeField(const Real value = 0.0) {
      _field.assign(_field.size(), value);
    };

    //! From array
    void initializeField(const Real* value) {
      _field.assign(value, value+_field.size());
    };

    //! Linearized update
    // From solution array
    void linearizedUpdate(const Real* localValues, Real fact=1.0) {
      const int nLocalDof = (_myMesh->getNumberOfNodes())*_nodeDoF;
      for(uint i = 0; i < nLocalDof; i++)
	_field[i] += localValues[i]*fact;
    };

    // One value at the time (Node ID, dof index, value)
    void linearizedUpdate(const int id, const int dof, const Real value) {
      // assert( id < _field.size() && dof == 0 );
      _field[id] += value;
    };

    void linearizedUpdate(const int dof, const Real value) {
      _field[dof] += value;
    };

    void setField(uint dof, Real value) {
      _field[dof] = value;
    };

    void setField(const Real* value) {
      _field.assign(value, value+_field.size());
    };

    void setPrevField() {
      _prevField = _field;
    };

    void getField(vector<double> & x) {
      assert(x.size() == _field.size());
      x = _field;
    }

    void printField() {
      for (uint i = 0; i < _field.size(); i++)
	cout << _field[i] << endl;
    };

    uint getNumMat() {
      set<DiffusionMaterial *> UNIQUEmaterials;
      for (uint i = 0; i < _materials.size(); i++) 
	UNIQUEmaterials.insert(_materials[i]);
	
      return UNIQUEmaterials.size();
    }

    //! Write output
    void writeOutputVTK(const string OutputFile, int step);  
    void writeOutput(const string OutputFile, const string format);

    //! Solve the system
    void compute(Result& result);

   

protected:
    // Default constructor cannot be used, Material need to be provided
    PoissonModel();
    //! List of material data at each quad point in the model
    vector<DiffusionMaterial* > _materials;

    //! Solution value at all nodes, local and ghost
    vector<Real > _field;
    vector<Real > _prevField;

  };

} // namespace voom

#endif
